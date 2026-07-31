"""Tests for the standard genetic code."""

from __future__ import annotations

import pytest

from bt4.domain.genetic_code import (
    AA_TO_CODONS,
    AMINO_ACIDS,
    CODON_TABLE,
    STOP,
    is_stop,
    synonymous_codons,
    translate,
)


def test_table_is_complete() -> None:
    # 61 sense codons + 3 stops = 64.
    assert len(CODON_TABLE) == 64
    assert sum(1 for aa in CODON_TABLE.values() if aa == STOP) == 3
    assert frozenset("ACDEFGHIKLMNPQRSTVWY") == AMINO_ACIDS


def test_every_amino_acid_has_codons_and_they_translate_back() -> None:
    for aa in AMINO_ACIDS | {STOP}:
        codons = synonymous_codons(aa)
        assert codons, f"{aa} has no codons"
        assert codons == tuple(sorted(codons)), "codons must be sorted (determinism)"
        for codon in codons:
            assert CODON_TABLE[codon] == aa


def test_met_and_trp_are_single_codon() -> None:
    assert synonymous_codons("M") == ("ATG",)
    assert synonymous_codons("W") == ("TGG",)


def test_translate_and_stop() -> None:
    assert translate("ATGGCCTGA") == "MA" + STOP
    assert is_stop("TAA") and is_stop("tag") and is_stop("TGA")
    assert not is_stop("ATG")


def test_translate_rejects_bad_length_and_unknown_codon() -> None:
    with pytest.raises(ValueError):
        translate("ATGG")
    with pytest.raises(ValueError):
        translate("ATGXYZ")


def test_case_insensitivity() -> None:
    assert translate("atggcc") == translate("ATGGCC")
    assert synonymous_codons("a") == AA_TO_CODONS["A"]
