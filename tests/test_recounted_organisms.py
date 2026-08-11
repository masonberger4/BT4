"""Tests for the six codon tables recounted from release-pinned Ensembl CDS sets.

Unlike the three older bundled tables (human / *E. coli* / yeast), which are
honestly labeled *representative published summaries*, these six are **genome-wide
recounts**: every number is a codon occurrence count from a public,
release-pinned CDS FASTA, produced by ``scripts/build_organism_tables.py``. So
the tests here assert two things the older tables cannot claim:

* **A complete re-derivation trail** -- source URL, the source file's own
  SHA-256, assembly, database release, and the filter tally -- so a third party
  can rebuild the exact shipped bytes (CLAUDE.md §8).
* **External ground truth** (§8, "not just self-consistency"): the tables must
  reproduce well-established, independently-published facts about each species'
  codon bias. A fabricated or mis-parsed table would fail these; internal
  consistency alone would not catch it.

These six are also why the bundled GtRNAdb tRNA tables became reachable at all:
tAI needs an organism you can actually select, which needs a codon table.
"""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from bt4 import api
from bt4.biomodels.codon.tables import (
    CodonUsageTable,
    available_organisms,
    load_provenance,
    load_table,
    sha256_hex,
)
from bt4.biomodels.codon.tai import available_tai_organisms
from bt4.domain.genetic_code import CODON_TABLE

RECOUNTED: tuple[str, ...] = (
    "arabidopsis_thaliana",
    "caenorhabditis_elegans",
    "danio_rerio",
    "drosophila_melanogaster",
    "mus_musculus",
    "rattus_norvegicus",
)

_STOPS = ("TAA", "TAG", "TGA")


def _raw_provenance(name: str) -> dict[str, object]:
    """The provenance sidecar as raw JSON (``TableProvenance`` drops extra keys)."""
    text = files("bt4.biomodels.codon.data").joinpath(f"{name}.provenance.json")
    return dict(json.loads(text.read_text(encoding="utf-8")))


def _top_codon(table: CodonUsageTable, amino_acid: str) -> str:
    """The highest-frequency codon among ``amino_acid``'s synonyms."""
    synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
    return max(synonyms, key=table.weight)


def _gc3(table: CodonUsageTable) -> float:
    """Fraction of sense codons ending in G or C -- the classic GC3 readout."""
    sense = {c: v for c, v in table.frequency.items() if c not in _STOPS}
    return sum(v for c, v in sense.items() if c[2] in "GC") / sum(sense.values())


# --------------------------------------------------------------------------- #
# The tables load, are advertised, and are complete.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", RECOUNTED)
def test_table_loads_and_is_advertised(name: str) -> None:
    table = load_table(name)
    assert table.organism == name
    assert name in available_organisms()
    # A genome-scale CDS set observes every codon; the builder refuses to write
    # a table that does not, rather than smoothing an invented value into it.
    assert set(table.frequency) == set(CODON_TABLE)


@pytest.mark.parametrize("name", RECOUNTED)
def test_counts_are_positive_integers(name: str) -> None:
    table = load_table(name)
    for codon, value in table.frequency.items():
        assert value == int(value), f"{name} {codon} is not a whole count"
        assert value >= 1


@pytest.mark.parametrize("name", RECOUNTED)
def test_group_max_codon_has_unit_weight(name: str) -> None:
    # Leucine is a six-box amino acid; its most-used synonym must normalize to 1.
    assert load_table(name).weight(_top_codon(load_table(name), "L")) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# Provenance: a complete, checkable re-derivation trail.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_sha256_matches_shipped_tsv(name: str) -> None:
    raw = files("bt4.biomodels.codon.data").joinpath(f"{name}.tsv").read_bytes()
    assert load_provenance(name).sha256 == sha256_hex(raw)


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_is_a_real_recount_not_a_summary(name: str) -> None:
    prov = load_provenance(name)
    # The load-bearing distinction from the older bundled tables: a real CDS
    # count stands behind these numbers, and the note says so plainly.
    assert prov.cds_count is not None
    assert prov.cds_count > 5_000, "a genome-wide CDS set, not a sample"
    assert "Real genome-wide codon counts" in prov.note
    assert "REPRESENTATIVE" not in prov.note


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_carries_the_rederivation_trail(name: str) -> None:
    prov = _raw_provenance(name)
    # Content-hashing our own output is not enough: without the SOURCE's hash and
    # URL nobody else can reproduce the table, which is what §8 asks for.
    for key in (
        "source_url",
        "source_sha256",
        "assembly",
        "database",
        "database_release",
        "total_codons_counted",
        "filters",
        "rebuild_command",
    ):
        assert key in prov, f"{name} provenance is missing {key!r}"
    assert str(prov["source_url"]).startswith("https://")
    assert len(str(prov["source_sha256"])) == 64
    assert f"release-{prov['database_release']}" in str(prov["source_url"])


@pytest.mark.parametrize("name", RECOUNTED)
def test_stamped_totals_match_the_shipped_counts(name: str) -> None:
    prov = _raw_provenance(name)
    table = load_table(name)
    assert sum(table.frequency.values()) == prov["total_codons_counted"]
    filters = prov["filters"]
    assert isinstance(filters, dict)
    assert filters["cds_counted"] == load_provenance(name).cds_count
    assert filters["records_in_source"] >= filters["cds_counted"]


# --------------------------------------------------------------------------- #
# External ground truth: independently-published facts about each species.
# --------------------------------------------------------------------------- #


def test_gc3_ordering_matches_known_genome_composition() -> None:
    """GC3 must order the species the way the literature does.

    Drosophila is strongly GC3-biased; mammals sit near 0.55-0.60; *C. elegans*
    and *Arabidopsis* are AT-rich coding genomes. Getting this ordering right is
    a real constraint on the data -- a mis-parsed or invented table would not.
    """
    gc3 = {name: _gc3(load_table(name)) for name in RECOUNTED}
    gc3["homo_sapiens"] = _gc3(load_table("homo_sapiens"))

    assert gc3["drosophila_melanogaster"] > gc3["danio_rerio"]
    assert gc3["danio_rerio"] > gc3["arabidopsis_thaliana"]
    assert gc3["arabidopsis_thaliana"] > gc3["caenorhabditis_elegans"]
    # AT-rich coding genomes sit well below half.
    for name in ("caenorhabditis_elegans", "arabidopsis_thaliana"):
        assert gc3[name] < 0.45
    # ...and the GC3-rich fly well above.
    assert gc3["drosophila_melanogaster"] > 0.60


def test_mammals_agree_with_each_other() -> None:
    """Mouse and rat were counted independently, yet must land near human.

    Three tables built from three separate CDS sets converging on the mammalian
    GC3 range is evidence the pipeline measured something real, not an artifact
    of one download.
    """
    human = _gc3(load_table("homo_sapiens"))
    for name in ("mus_musculus", "rattus_norvegicus"):
        assert abs(_gc3(load_table(name)) - human) < 0.05


@pytest.mark.parametrize(
    ("name", "amino_acid", "expected"),
    [
        # CTG is the preferred Leu codon across GC3-rich genomes...
        ("mus_musculus", "L", "CTG"),
        ("rattus_norvegicus", "L", "CTG"),
        ("danio_rerio", "L", "CTG"),
        ("drosophila_melanogaster", "L", "CTG"),
        # ...while AT-rich coding genomes prefer CTT.
        ("caenorhabditis_elegans", "L", "CTT"),
        ("arabidopsis_thaliana", "L", "CTT"),
    ],
)
def test_preferred_leucine_codon(name: str, amino_acid: str, expected: str) -> None:
    assert _top_codon(load_table(name), amino_acid) == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        # TGA dominates in GC-richer vertebrate/plant genomes, TAA in AT-rich ones.
        ("mus_musculus", "TGA"),
        ("rattus_norvegicus", "TGA"),
        ("danio_rerio", "TGA"),
        ("arabidopsis_thaliana", "TGA"),
        ("caenorhabditis_elegans", "TAA"),
        ("drosophila_melanogaster", "TAA"),
    ],
)
def test_preferred_stop_codon(name: str, expected: str) -> None:
    table = load_table(name)
    assert max(_STOPS, key=lambda c: table.frequency[c]) == expected


# --------------------------------------------------------------------------- #
# The point of the exercise: tAI is reachable for every organism that has it.
# --------------------------------------------------------------------------- #


def test_every_trna_table_has_a_selectable_organism() -> None:
    """No bundled tRNA table may be stranded without a codon table.

    tAI is only offered for an organism the user can actually select, and
    selection requires a codon-usage table. Before these six tables shipped, six
    of the eight bundled GtRNAdb tables were unreachable for exactly this reason
    -- this test is what keeps that from silently recurring.
    """
    stranded = sorted(set(available_tai_organisms()) - set(available_organisms()))
    assert stranded == [], f"tRNA data with no selectable organism: {stranded}"


@pytest.mark.parametrize("name", RECOUNTED)
def test_cai_in_unit_interval(name: str) -> None:
    table = load_table(name)
    dna = "ATG" + "GCA" + "CTT" + "CGG" + "AAA" + "TAA"
    assert 0.0 < table.cai(dna) <= 1.0


@pytest.mark.parametrize("name", RECOUNTED)
def test_table_content_hash_reaches_the_run_manifest(name: str) -> None:
    """Each new table's *content* hash must enter the provenance stamp (#9).

    Hashing config field names -- BT4's cautionary tale from BT3 (§10.10) -- would
    make two organisms stamp identically. Six new tables are six new chances to
    regress that, so assert the actual TSV digest appears in the manifest.
    """
    result = api.optimize(
        "MAALKHETQW", api.OptimizeConfig(organism=name, max_homopolymer=5)
    )
    manifest = json.loads(api.result_to_json(result))["audit"]["manifest"]
    assert manifest["inputs"]["codon_table_sha256"] == load_provenance(name).sha256


def test_each_organism_stamps_a_distinct_manifest() -> None:
    """Swapping the organism must change the stamp (invariant #9)."""
    stamps = set()
    for name in RECOUNTED:
        result = api.optimize(
            "MAALKHETQW", api.OptimizeConfig(organism=name, max_homopolymer=5)
        )
        manifest = json.loads(api.result_to_json(result))["audit"]["manifest"]
        stamps.add(json.dumps(manifest, sort_keys=True))
    assert len(stamps) == len(RECOUNTED)
