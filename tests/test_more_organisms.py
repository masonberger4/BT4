"""Tests for the additional bundled codon-usage tables (E. coli, S. cerevisiae).

These organisms are auto-discovered from packaged TSV data. Each test asserts
the table loads, is advertised by :func:`available_organisms`, covers every
amino acid plus stops (enforced by the table constructor), carries the expected
strong codon bias, and has a provenance SHA-256 that matches the shipped bytes.
"""

from __future__ import annotations

from importlib.resources import files

import pytest

from bt4.biomodels.codon.tables import (
    CodonUsageTable,
    available_organisms,
    load_provenance,
    load_table,
    sha256_hex,
)
from bt4.domain.genetic_code import CODON_TABLE

_NEW_ORGANISMS: tuple[str, ...] = ("escherichia_coli", "saccharomyces_cerevisiae")


def _group_max_codon(table: CodonUsageTable, amino_acid: str) -> str:
    """Return the highest-frequency codon among ``amino_acid``'s synonyms."""
    synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
    return max(synonyms, key=table.weight)


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_table_loads_with_canonical_organism(name: str) -> None:
    table = load_table(name)
    assert table.organism == name


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_appears_in_available_organisms(name: str) -> None:
    assert name in available_organisms()


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_covers_every_amino_acid_and_stop(name: str) -> None:
    # A successful load proves coverage: the constructor rejects any table that
    # is missing an amino acid or a stop codon. Confirm all 64 codons are present.
    table = load_table(name)
    assert set(table.frequency) == set(CODON_TABLE)


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_group_max_codon_has_unit_weight(name: str) -> None:
    table = load_table(name)
    # Leucine is a six-box amino acid; its most-used synonymous codon must
    # normalize to a relative adaptiveness of exactly 1.0.
    top_leu = _group_max_codon(table, "L")
    assert table.weight(top_leu) == pytest.approx(1.0)


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_provenance_sha256_matches_shipped_tsv(name: str) -> None:
    prov = load_provenance(name)
    raw = files("bt4.biomodels.codon.data").joinpath(f"{name}.tsv").read_bytes()
    assert prov.sha256 == sha256_hex(raw)


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_provenance_is_labeled_representative(name: str) -> None:
    prov = load_provenance(name)
    # Honesty invariant: representative summary, not a per-genome recount.
    assert prov.cds_count is None
    assert "REPRESENTATIVE" in prov.note
    assert "not an authoritative" in prov.note.lower()


@pytest.mark.parametrize("name", _NEW_ORGANISMS)
def test_cai_in_unit_interval(name: str) -> None:
    table = load_table(name)
    # A short mixed peptide (M A L R K stop) exercises degenerate residues.
    dna = "ATG" + "GCA" + "CTT" + "CGG" + "AAA" + "TAA"
    cai = table.cai(dna)
    assert 0.0 < cai <= 1.0


def test_expected_strong_biases() -> None:
    # E. coli strongly prefers CTG for Leu; S. cerevisiae strongly prefers AGA
    # for Arg. These well-known biases should surface as the unit-weight codon.
    ecoli = load_table("escherichia_coli")
    yeast = load_table("saccharomyces_cerevisiae")
    assert ecoli.weight("CTG") == pytest.approx(1.0)
    assert yeast.weight("AGA") == pytest.approx(1.0)
