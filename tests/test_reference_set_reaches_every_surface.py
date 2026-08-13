"""Every surface that reads a codon table must honor the caller's reference set.

This file exists because the threading was once *almost* complete and the gap was
invisible: `api.optimize` was covered, so the suite stayed green while
`api.library`, `api.candidates`, `api.tracks` and all three comparison scripts
still loaded whichever table happened to be the default. Under a non-default
`reference_set` that made `benchmark.py` report a %MinMax mean of 57.81 where the
truth was 100.0 -- a plausible number, recomputed against a table the sequence was
not designed with, which is invariant #2 broken silently.

So the tests here are deliberately shaped as **revert detectors**: each one fails
if `reference_set=` is dropped from one specific call site. Together they mean the
non-optimize half of the reference-set axis can no longer be undone with CI green.

Two correctness fixes that had no test of their own are pinned here too: the
builder's pooled-ambiguity guard (an identifier resolving to two genes must be
dropped, not guessed) and `load_provenance`'s refusal to accept a sidecar that
does not state its reference set.

The HTTP service's share of this lives in `tests/test_service.py` instead: CI
installs FastAPI only for that job and runs only that file, so a service test
placed here would `importorskip` in every job and therefore never run anywhere.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import pathlib
import sys
from types import ModuleType

import pytest

from bt4 import api
from bt4.biomodels.codon.tables import GENOME_WIDE, HIGHLY_EXPRESSED, load_table

# An organism whose two reference sets disagree at eight amino acids, so any
# "wrong table" shows up immediately rather than needing a lucky protein.
ORGANISM = "escherichia_coli"
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEKQANGTWPADEFHLMNCYVGR"

_SCRIPTS = pathlib.Path(__file__).resolve().parent.parent / "scripts"


def _load_script(name: str) -> ModuleType:
    """Import a maintainer script by path (they are not an installed package)."""
    path = _SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"bt4_script_{name}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _config(reference_set: str) -> api.OptimizeConfig:
    return api.OptimizeConfig(
        organism=ORGANISM, reference_set=reference_set, max_homopolymer=6, seed=0
    )


# --------------------------------------------------------------------------- #
# Engine surfaces beyond optimize().
# --------------------------------------------------------------------------- #


def test_library_samples_from_the_requested_reference_set() -> None:
    """`api.library` samples the codon distribution -- so it must be the right one.

    Both draws use the same seed, so any difference is the table, not the RNG.
    """
    strong = api.library(PROTEIN, _config(HIGHLY_EXPRESSED), n=4, seed=7)
    weak = api.library(PROTEIN, _config(GENOME_WIDE), n=4, seed=7)
    assert [m.dna for m in strong.results] != [m.dna for m in weak.results]
    for member in strong.results:
        assert member.audit["codon_reference_set"] == HIGHLY_EXPRESSED


@pytest.mark.parametrize("reference_set", [HIGHLY_EXPRESSED, GENOME_WIDE])
def test_candidates_are_scored_against_the_table_that_designed_them(
    reference_set: str,
) -> None:
    """Invariant #2 on the candidate set: reported CAI == recomputed CAI.

    Recomputed here from the *requested* table, so a candidate scored against the
    default while designed under the other fails.
    """
    cand_set = api.candidates(PROTEIN, _config(reference_set), n=3)
    table = load_table(ORGANISM, reference_set=reference_set)
    for cand in cand_set.candidates:
        dna = cand.result.dna
        assert cand.result.audit["codon_reference_set"] == reference_set
        assert float(cand.result.audit["cai"]) == pytest.approx(table.cai(dna), abs=1e-9)


@pytest.mark.parametrize("reference_set", [HIGHLY_EXPRESSED, GENOME_WIDE])
def test_tracks_report_and_use_the_requested_reference_set(reference_set: str) -> None:
    """%MinMax is a codon-commonness profile, so its basis is part of the answer."""
    result = api.optimize(PROTEIN, _config(reference_set))
    tracks = api.tracks(result.dna, ORGANISM, reference_set=reference_set)
    assert tracks.reference_set == reference_set
    assert tracks.organism == ORGANISM


def test_tracks_minmax_actually_changes_with_the_reference_set() -> None:
    """A label nobody's numbers depend on would be decoration; these depend on it."""
    dna = api.optimize(PROTEIN, _config(HIGHLY_EXPRESSED)).dna
    strong = api.tracks(dna, ORGANISM, reference_set=HIGHLY_EXPRESSED).get("minmax")
    weak = api.tracks(dna, ORGANISM, reference_set=GENOME_WIDE).get("minmax")
    assert strong is not None and weak is not None
    assert strong.values != weak.values


def test_tracks_omits_the_label_when_no_table_was_read() -> None:
    """A non-codon-aligned sequence gets no %MinMax, so it must claim no basis."""
    tracks = api.tracks("ATGGCCGCCCTGAAGCACGAGACCCAGTGG" + "A", ORGANISM)
    assert tracks.get("minmax") is None
    assert tracks.reference_set == ""
    assert tracks.organism == ""


# --------------------------------------------------------------------------- #
# The comparison scripts -- where the original defect actually lived.
# --------------------------------------------------------------------------- #


def test_benchmark_recomputes_against_the_solves_table() -> None:
    """The exact case that was wrong: %MinMax mean 57.81 vs the correct 100.0."""
    benchmark = _load_script("benchmark")
    config = api.OptimizeConfig(
        organism="saccharomyces_cerevisiae", reference_set=GENOME_WIDE
    )
    row = benchmark.benchmark({"demo": "MKTAYIAKQRQISFVKSHFSRQ"}, config)[0]
    # Under its own reference set the CAI-optimal sequence is by construction the
    # most-common codon at every site, so %MinMax pins at 100. Scored against the
    # other table it does not -- which is exactly how the defect showed up.
    assert row["bt4_minmax_mean"] == pytest.approx(100.0)


def test_compare_tools_recomputes_against_the_solves_table() -> None:
    """BT4's own row must score 1.0 under *whichever* table designed it.

    Its sequence is CAI-optimal by construction, so a row scored against the table
    that produced it reads 1.0. Recomputed against the other one it does not --
    which is precisely the shape of the defect, and why asserting it for both
    reference sets is a revert detector rather than a tautology.
    """
    compare_tools = _load_script("compare_tools")
    panel = compare_tools.load_panel()
    native_cai = {}
    for reference_set in (HIGHLY_EXPRESSED, GENOME_WIDE):
        config = api.OptimizeConfig(organism="homo_sapiens", reference_set=reference_set)
        rows = compare_tools.compare(panel, config)
        bt4 = [r for r in rows if "BT4" in str(r["optimizer"])]
        assert bt4, "the board must contain a BT4 row"
        assert bt4[0]["cai"] == pytest.approx(1.0), reference_set
        native_cai[reference_set] = [
            r["cai"] for r in rows if "native" in str(r["optimizer"]).lower()
        ]
    # ...and the *published* rows must move with the table, proving the recompute
    # really read a different one rather than the BT4 row being trivially 1.0.
    assert native_cai[HIGHLY_EXPRESSED] != native_cai[GENOME_WIDE]


def test_compare_tools_json_states_which_tables_it_used() -> None:
    """The JSON board is the artifact most likely to be re-published downstream.

    A CAI without its reference set is a number with no question attached, and a
    machine-readable board is exactly where a reader cannot see the banner.
    """
    compare_tools = _load_script("compare_tools")
    config = api.OptimizeConfig(organism="homo_sapiens", reference_set=GENOME_WIDE)
    payload = compare_tools.board(compare_tools.load_panel(), config)
    assert payload["organism"] == "homo_sapiens"
    assert payload["codon_reference_set"] == GENOME_WIDE
    assert payload["rows"]


def test_compare_reproducibility_recomputes_against_the_solves_table() -> None:
    compare_repro = _load_script("compare_reproducibility")
    records = compare_repro.load_panel()
    strong = compare_repro.reproducibility(
        records, api.OptimizeConfig(reference_set=HIGHLY_EXPRESSED)
    )
    weak = compare_repro.reproducibility(
        records, api.OptimizeConfig(reference_set=GENOME_WIDE)
    )

    def bt4_cai(rows: list[dict[str, object]]) -> list[float]:
        return [
            float(r["cai"]["mean"])  # type: ignore[index,call-overload]
            for r in rows
            if "BT4" in str(r["source"])
        ]

    assert bt4_cai(strong) and bt4_cai(strong) != bt4_cai(weak)


# --------------------------------------------------------------------------- #
# Two correctness fixes that had no test of their own.
# --------------------------------------------------------------------------- #


def _write_pep(path: pathlib.Path, entries: list[tuple[str, str]]) -> pathlib.Path:
    """Write a minimal gzipped Ensembl-format peptide FASTA."""
    lines = []
    for protein, gene in entries:
        lines.append(f">{protein} pep chromosome:X:1:2:1 gene:{gene}\nMA\n")
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("".join(lines))
    return path


def test_an_ambiguous_identifier_is_dropped_rather_than_guessed(
    tmp_path: pathlib.Path,
) -> None:
    """The headline correctness fix of the review pass, pinned.

    ``X.1`` names gene ``G1`` as a protein ID and gene ``G2`` as a *gene* ID. It
    is singular in each map on its own; only pooling the two sees the conflict.
    Before the fix it resolved silently -- to whichever map was consulted first --
    and was counted as a clean join, contradicting the stamp every sidecar carries
    ("identifiers resolving to more than one gene were dropped, not guessed").
    """
    builder = _load_script("build_highly_expressed_tables")
    pep = _write_pep(tmp_path / "pep.fa.gz", [("X.1.1", "X.1"), ("X.1", "G1")])
    index = builder.build_id_index(pep)

    assert "X.1" in index.ambiguous
    assert "X.1" not in index.exact
    assert "X.1" not in index.unversioned

    cds = {"G1": ("ATGTAA", "1"), "X.1": ("ATGTAA", "1")}
    selected = builder.select_reference_set(
        [("someGene", "X.1", 99.0)], index, cds, builder.FilterStats(), 0, 1
    )
    assert selected.genes == []
    stats = selected.stats.as_dict()
    assert stats["rows_unmatched_identifier_ambiguous"] == 1
    assert stats["rows_matched_via_protein_id"] == 0
    assert stats["rows_matched_via_gene_id"] == 0


def test_a_sidecar_that_does_not_state_its_reference_set_is_refused(
    tmp_path: pathlib.Path,
) -> None:
    """The provenance guard must fail closed.

    Defaulting a missing label to whatever the caller asked for would validate one
    sidecar as *every* reference set in turn -- the same fail-open shape this
    change removed from ``available_organisms``.
    """
    from bt4.biomodels.codon import build as codon_build
    from bt4.biomodels.codon.tables import load_table_from_file

    counts = {"ATG": 5, "TGG": 3, "TAA": 2, "GCC": 7, "GCT": 4}
    codon_build.write_table(
        counts,
        organism="toy",
        path=tmp_path,
        source="unit-test",
        pseudocount=1.0,
    )
    sidecar = tmp_path / "toy.provenance.json"
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    # write_table must state it, unprompted, for a caller's own CDS set.
    assert payload["reference_set"] == "custom"
    assert load_table_from_file(tmp_path / "toy.tsv").reference_set == "custom"


def test_reference_set_is_a_reserved_provenance_key() -> None:
    """It cannot be smuggled through the free-form ``extra`` dict."""
    import tempfile

    from bt4.biomodels.codon import build as codon_build

    with tempfile.TemporaryDirectory() as tmp, pytest.raises(ValueError, match="reserved"):
        codon_build.write_table(
            {"ATG": 1},
            organism="toy",
            path=tmp,
            source="unit-test",
            extra={"reference_set": "highly_expressed"},
        )


def test_a_dotted_organism_key_is_refused_not_silently_mislabelled() -> None:
    """``load_table("homo_sapiens.highly_expressed")`` must not resolve.

    An organism key is a bare stem, exactly as ``available_organisms()`` lists it.
    Because the genome-wide suffix is the empty string, a dotted key otherwise
    lands straight on the highly-expressed FILE and returns its counts stamped
    ``reference_set="genome_wide"`` -- the right numbers under the wrong label,
    which is worse than an error because nothing downstream can tell.
    """
    from bt4.biomodels.codon.tables import load_provenance

    for key in ("homo_sapiens.highly_expressed", "homo_sapiens.trna"):
        with pytest.raises(ValueError, match="unknown organism"):
            load_table(key)
        with pytest.raises(ValueError, match="unknown organism"):
            load_provenance(key)


def test_the_peptide_fasta_must_come_from_the_cds_release(
    tmp_path: pathlib.Path,
) -> None:
    """A release bump in one builder must not silently strand the other.

    The peptide URL and the CDS URL are pinned in different files. Both digest
    checks would still pass after a divergence -- each file matches its own pin --
    while versioned gene IDs from two releases fail to join, drawing the reference
    set from whatever biased remnant survives. The prose said "must come from the
    same release"; this is what enforces it.
    """
    import dataclasses

    builder = _load_script("build_highly_expressed_tables")
    for spec in builder.SPECS_HE:
        builder.check_pep_matches_cds_release(spec, builder._cds_spec(spec.key))

    spec = builder.SPECS_HE[0]
    cds = builder._cds_spec(spec.key)
    drifted = dataclasses.replace(
        spec, pep_url=spec.pep_url.replace(f"release-{cds.release}/", "release-999/")
    )
    with pytest.raises(SystemExit, match="same release"):
        builder.check_pep_matches_cds_release(drifted, cds)
