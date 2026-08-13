"""Tests for the BT4 codon-usage table biomodel."""

from __future__ import annotations

import math
from importlib.resources import files
from pathlib import Path

import pytest

from bt4.biomodels.codon.tables import (
    GENOME_WIDE,
    HIGHLY_EXPRESSED,
    REFERENCE_SET_SUFFIX,
    CodonUsageTable,
    available_organisms,
    available_reference_sets,
    load_provenance,
    load_table,
    load_table_from_file,
    sha256_hex,
)
from bt4.domain.genetic_code import CODON_TABLE

_HEADER = "amino_acid\tcodon\tfrequency\n"


def _full_tsv_rows() -> list[str]:
    """One row per codon with a distinct positive frequency."""
    rows: list[str] = []
    for i, (codon, aa) in enumerate(sorted(CODON_TABLE.items())):
        rows.append(f"{aa}\t{codon}\t{float(i + 1)}")
    return rows


def _write_tsv(path: Path, rows: list[str]) -> Path:
    path.write_text(_HEADER + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_homo_sapiens_loads_and_is_covered() -> None:
    table = load_table("homo_sapiens")
    assert table.organism == "homo_sapiens"
    assert "homo_sapiens" in available_organisms()


def test_available_organisms_excludes_trna_tables() -> None:
    # The tRNA tables (``<organism>.trna.tsv``) are tAI data, not codon-usage
    # tables. They must not leak into the codon-table organism list, where they
    # would appear as bogus organisms like "homo_sapiens.trna" (and then fail to
    # load as a codon table). Regression test for that filter.
    orgs = available_organisms()
    assert not any(name.endswith(".trna") for name in orgs)
    assert "homo_sapiens.trna" not in orgs
    assert "homo_sapiens" in orgs  # the real codon organism is still there


def test_alias_resolution() -> None:
    assert load_table("human").organism == "homo_sapiens"
    assert load_table("Human").organism == "homo_sapiens"


def test_weight_is_freq_over_group_max() -> None:
    """w(codon) is its frequency divided by its synonymous group's maximum.

    Derived from the table under test rather than pinned to literal frequencies:
    the numbers belong to the *data*, which is rebuilt from source, while this
    test is about the *formula*. Hard-coding one table's values made it fail on
    any recount for a reason that had nothing to do with the code.
    """
    table = load_table("homo_sapiens")
    leucines = [c for c, aa in CODON_TABLE.items() if aa == "L"]
    group_max = max(table.frequency[c] for c in leucines)
    for codon in leucines:
        assert table.weight(codon) == pytest.approx(
            table.frequency[codon] / group_max
        )
    # The most-used synonym normalizes to exactly 1...
    assert max(table.weight(c) for c in leucines) == pytest.approx(1.0)
    # ...as does a single-box amino acid (Met), trivially.
    assert table.weight("ATG") == pytest.approx(1.0)


def test_relative_adaptiveness_matches_weight() -> None:
    table = load_table("homo_sapiens")
    w = table.relative_adaptiveness()
    assert w["GCC"] == pytest.approx(table.weight("GCC"))


def test_cai_geometric_mean_known_sequence() -> None:
    table = load_table("homo_sapiens")
    # M A L K stop: only A, L, K are scored (M is single-box, stop excluded).
    # Use non-optimal codons so w < 1 and the geometric mean is discriminating.
    dna = "ATG" + "GCA" + "CTT" + "AAA" + "TAA"
    expected = math.exp(
        (
            math.log(table.weight("GCA"))
            + math.log(table.weight("CTT"))
            + math.log(table.weight("AAA"))
        )
        / 3
    )
    assert table.cai(dna) == pytest.approx(expected)
    assert 0.0 < table.cai(dna) < 1.0


def test_cai_all_nondegenerate_is_one() -> None:
    table = load_table("homo_sapiens")
    assert table.cai("ATGTGG") == pytest.approx(1.0)  # Met + Trp only


def test_cai_bad_length_raises() -> None:
    table = load_table("homo_sapiens")
    with pytest.raises(ValueError):
        table.cai("ATGG")


def test_cai_unknown_codon_raises() -> None:
    table = load_table("homo_sapiens")
    with pytest.raises(ValueError):
        table.cai("ATGXYZ")


def test_unknown_organism_raises() -> None:
    with pytest.raises(ValueError):
        load_table("martian")


def test_missing_amino_acid_raises(tmp_path: Path) -> None:
    rows = [r for r in _full_tsv_rows() if not r.startswith("W\t")]  # drop Trp
    p = _write_tsv(tmp_path / "broken.tsv", rows)
    with pytest.raises(ValueError, match="W"):
        load_table_from_file(p)


def test_load_table_from_file_roundtrip(tmp_path: Path) -> None:
    p = _write_tsv(tmp_path / "toy.tsv", _full_tsv_rows())
    table = load_table_from_file(p)
    assert table.organism == "toy"


@pytest.mark.parametrize("reference_set", [GENOME_WIDE, HIGHLY_EXPRESSED])
def test_provenance_sha256_matches_shipped_tsv(reference_set: str) -> None:
    prov = load_provenance("homo_sapiens", reference_set=reference_set)
    stem = f"homo_sapiens{REFERENCE_SET_SUFFIX[reference_set]}"
    raw = files("bt4.biomodels.codon.data").joinpath(f"{stem}.tsv").read_bytes()
    assert prov.sha256 == sha256_hex(raw)
    assert prov.reference_set == reference_set


def test_provenance_alias_respects_the_reference_set() -> None:
    for reference_set in available_reference_sets("human"):
        by_alias = load_provenance("human", reference_set=reference_set)
        by_key = load_provenance("homo_sapiens", reference_set=reference_set)
        assert by_alias == by_key


def test_provenance_alias() -> None:
    assert load_provenance("human").source == load_provenance("homo_sapiens").source


def test_construction_rejects_nonpositive_frequency() -> None:
    freq = {codon: 1.0 for codon in CODON_TABLE}
    freq["GCC"] = 0.0
    with pytest.raises(ValueError):
        CodonUsageTable(organism="x", frequency=freq)
