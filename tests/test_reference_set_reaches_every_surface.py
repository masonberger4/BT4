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
from bt4.biomodels.codon.tables import (
    GENOME_WIDE,
    HIGHLY_EXPRESSED,
    load_provenance,
    load_table,
)

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


# A repetitive protein whose CAI-optimal seed carries a repeat longer than the
# limit below, so candidate assembly must produce repeat-refined variants.
_REPEATY_PROTEIN = "M" + "EAAAK" * 6 + "GGGGSGGGGS" + "HHHHHH"


def test_repeat_refined_candidates_use_the_run_table_not_the_default() -> None:
    """The candidates revert-detector must reach the ``repeat_refined`` members.

    ``pipeline/candidates.py`` builds repeat-refined variants from a table loaded
    at one call site; the other candidates come from ``run_frontier``, which
    threads the reference set on its own. So a config with **no** GLOBAL rule
    produces zero repeat-refined variants and never exercises that call site --
    which is exactly why the earlier suite stayed green when it was reverted.
    This test forces the variants (a GLOBAL ``max_repeat_length`` the seed
    violates) under ``genome_wide`` and checks invariant #2 for every member.

    MUTATION THAT MUST FAIL THIS: drop ``reference_set=config.reference_set`` from
    the ``load_table`` call in ``assemble_and_rank_candidates``. Verified: the
    repeat-refined member then reports ``highly_expressed`` and its CAI recomputes
    against the wrong table.
    """
    config = api.OptimizeConfig(
        organism=ORGANISM, reference_set=GENOME_WIDE, max_repeat_length=10, seed=0
    )
    cand_set = api.candidates(_REPEATY_PROTEIN, config, n=6)

    sources = {c.source for c in cand_set.candidates}
    assert "repeat_refined" in sources, (
        "premise broken: no repeat_refined variant produced, so this test would "
        f"pass vacuously; sources={sources}"
    )

    table = load_table(ORGANISM, reference_set=GENOME_WIDE)
    for cand in cand_set.candidates:
        assert cand.result.audit["codon_reference_set"] == GENOME_WIDE, cand.source
        assert float(cand.result.audit["cai"]) == pytest.approx(
            table.cai(cand.result.dna), abs=1e-9
        ), cand.source


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


def test_write_table_stamps_a_custom_reference_set(tmp_path: pathlib.Path) -> None:
    """A caller's own CDS set is honestly labelled ``custom``, unprompted.

    (This does NOT test the load-time guard -- see the ``load_*`` tests below.
    An earlier version of this test was named for the guard but only checked
    ``write_table``/``load_table_from_file``, so a reverted guard passed it: the
    exact "test that would pass if the thing it names were broken" the project
    warns about. Renamed to what it actually checks.)
    """
    from bt4.biomodels.codon import build as codon_build
    from bt4.biomodels.codon.tables import load_table_from_file

    counts = {"ATG": 5, "TGG": 3, "TAA": 2, "GCC": 7, "GCT": 4}
    codon_build.write_table(
        counts, organism="toy", path=tmp_path, source="unit-test", pseudocount=1.0
    )
    payload = json.loads((tmp_path / "toy.provenance.json").read_text(encoding="utf-8"))
    assert payload["reference_set"] == "custom"
    assert load_table_from_file(tmp_path / "toy.tsv").reference_set == "custom"


# --------------------------------------------------------------------------- #
# The sidecar reference-set guard -- tested at the loader, not the writer.
# --------------------------------------------------------------------------- #
#
# These use a monkeypatched data directory so a crafted sidecar never touches the
# real package data. Each is paired with the exact one-line mutation that must
# make it fail; a test that stays green under its mutation is worthless (the
# project's own review caught one such vacuous test, which is why these exist).

_TSV_HEADER = "amino_acid\tcodon\tfrequency\n"


def _write_fake_table(
    directory: pathlib.Path, stem: str, sidecar: dict[str, object] | None
) -> None:
    """Write a valid <stem>.tsv (all 64 codons) and, optionally, a sidecar."""
    from bt4.domain.genetic_code import CODON_TABLE

    rows = [f"{aa}\t{codon}\t{i + 1}" for i, (codon, aa) in enumerate(sorted(CODON_TABLE.items()))]
    (directory / f"{stem}.tsv").write_text(_TSV_HEADER + "\n".join(rows) + "\n", encoding="utf-8")
    if sidecar is not None:
        (directory / f"{stem}.provenance.json").write_text(
            json.dumps(sidecar), encoding="utf-8"
        )


def _point_data_dir_at(monkeypatch: pytest.MonkeyPatch, directory: pathlib.Path) -> None:
    """Make the tables module read its bundled data from ``directory``."""
    from bt4.biomodels.codon import tables as tables_mod

    monkeypatch.setattr(tables_mod, "files", lambda _pkg: directory)


def _full_sidecar(reference_set: str | None) -> dict[str, object]:
    payload: dict[str, object] = {
        "source": "unit-test",
        "build": "unit-test",
        "cds_count": 1,
        "retrieved": "2026-01-01",
        "sha256": "0" * 64,
        "note": "unit-test",
    }
    if reference_set is not None:
        payload["reference_set"] = reference_set
    return payload


def test_load_provenance_refuses_a_keyless_sidecar(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A sidecar with no ``reference_set`` must be refused, not defaulted.

    MUTATION THAT MUST FAIL THIS: in ``_validate_stamped_reference_set``, replace
    the missing-key ``raise`` with ``stamped = str(data.get("reference_set",
    resolved))``. Verified: with that revert, this test fails and the suite
    otherwise stays green -- which is the whole point.
    """
    _write_fake_table(tmp_path, "keyless", _full_sidecar(None))
    _point_data_dir_at(monkeypatch, tmp_path)
    with pytest.raises(ValueError, match="does not state its reference set"):
        load_provenance("keyless")


def test_load_table_refuses_a_keyless_sidecar(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """load_table enforces the same guard as load_provenance (findings #5/#6)."""
    _write_fake_table(tmp_path, "keyless", _full_sidecar(None))
    _point_data_dir_at(monkeypatch, tmp_path)
    with pytest.raises(ValueError, match="does not state its reference set"):
        load_table("keyless")


def test_load_table_refuses_a_mismatched_sidecar(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A ``custom``-stamped sidecar at ``<organism>.tsv`` must not load as genome_wide.

    This is finding #6's own reproduction: a user's ``build-table`` output
    (honestly ``custom``) dropped into the data dir. Before the fix, load_table
    derived the label from the filename and returned ``genome_wide``; load_provenance
    already refused, so the two disagreed and ``tracks`` printed the false label.
    """
    _write_fake_table(tmp_path, "mismatch", _full_sidecar("custom"))
    _point_data_dir_at(monkeypatch, tmp_path)
    with pytest.raises(ValueError, match="claims reference set 'custom'"):
        load_table("mismatch")
    with pytest.raises(ValueError, match="claims reference set 'custom'"):
        load_provenance("mismatch")


def test_tracks_refuses_rather_than_printing_a_false_reference_set(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The tracks path must agree with optimize: refuse, never mislabel.

    ``api.tracks`` reads only ``load_table``; before the fix it would have
    published ``codon_reference_set="genome_wide"`` for a ``custom`` table while
    ``api.optimize`` (which reads the sidecar) refused. Now both refuse.
    """
    _write_fake_table(tmp_path, "mismatch", _full_sidecar("custom"))
    _point_data_dir_at(monkeypatch, tmp_path)
    dna = "ATG" + "GCC" * 8 + "TAA"
    with pytest.raises(ValueError, match="claims reference set 'custom'"):
        api.tracks(dna, "mismatch")


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
