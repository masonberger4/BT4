"""Tests for the concrete objective terms (CLAUDE.md invariant #4).

The load-bearing property here is ``delta == score``: the running sum of a
term's per-codon ``delta`` (as the DP would accumulate it, with real growing
prefixes) must equal the term's whole-sequence ``score``. We also pin the CAI
term's orientation (larger is better; the max-``w`` codon is optimal) and the
GC-target validation.
"""

from __future__ import annotations

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4.biomodels.codon.tables import load_table
from bt4.domain.genetic_code import (
    AMINO_ACIDS,
    CODON_TABLE,
    STOP,
    synonymous_codons,
)
from bt4.objectives import CaiTerm, GcProximityTerm, iter_codons

_TABLE = load_table("homo_sapiens")
_W = _TABLE.relative_adaptiveness()
_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=120)


def _backtranslate(protein: str) -> str:
    """Deterministically pick the first synonymous codon per residue."""
    return "".join(synonymous_codons(aa)[0] for aa in protein)


def _accumulated_delta(term: CaiTerm | GcProximityTerm, dna: str) -> float:
    """Sum the term's deltas over growing real prefixes (as the DP would)."""
    acc = 0.0
    prefix = ""
    for pos, codon in iter_codons(dna):
        acc += term.delta(prefix, codon, pos)
        prefix += codon
    return acc


def test_names_and_context() -> None:
    cai = CaiTerm(_W)
    gc = GcProximityTerm(0.5)
    assert cai.name == "cai_logw"
    assert gc.name == "gc_proximity"
    assert cai.context_len() == 0
    assert gc.context_len() == 0


@given(protein=_PROTEIN)
def test_cai_delta_equals_score(protein: str) -> None:
    term = CaiTerm(_W)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


@given(protein=_PROTEIN, target=st.floats(min_value=0.0, max_value=1.0))
def test_gc_delta_equals_score(protein: str, target: float) -> None:
    term = GcProximityTerm(target)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


def test_cai_delta_zero_for_non_degenerate_and_stop() -> None:
    term = CaiTerm(_W)
    # Met (ATG), Trp (TGG) and every stop codon carry no coding choice.
    assert term.delta("", "ATG", 0) == 0.0
    assert term.delta("", "TGG", 0) == 0.0
    for stop in synonymous_codons(STOP):
        assert term.delta("", stop, 0) == 0.0


def test_cai_delta_equals_log_weight_for_degenerate() -> None:
    term = CaiTerm(_W)
    for aa in AMINO_ACIDS:
        codons = synonymous_codons(aa)
        if len(codons) == 1 or aa in {"M", "W"}:
            continue
        for codon in codons:
            assert term.delta("", codon, 0) == pytest.approx(math.log(_TABLE.weight(codon)))


def test_cai_max_weight_codon_maximizes_delta() -> None:
    term = CaiTerm(_W)
    for aa in AMINO_ACIDS:
        codons = synonymous_codons(aa)
        if len(codons) == 1 or aa in {"M", "W"}:
            continue
        best_by_weight = max(codons, key=_TABLE.weight)
        best_by_delta = max(codons, key=lambda c: term.delta("", c, 0))
        assert term.delta("", best_by_weight, 0) == pytest.approx(
            term.delta("", best_by_delta, 0)
        )
        # The max-weight codon has w == 1, so its log-weight delta is 0 (the top).
        assert term.delta("", best_by_weight, 0) == pytest.approx(0.0)


def test_cai_delta_is_case_insensitive() -> None:
    term = CaiTerm(_W)
    assert term.delta("", "ctg", 0) == pytest.approx(term.delta("", "CTG", 0))


@pytest.mark.parametrize("codon,expected_gc", [("AAA", 0), ("ATG", 1), ("GCG", 3), ("GGC", 3)])
def test_gc_delta_matches_codon_gc(codon: str, expected_gc: int) -> None:
    target = 0.5
    term = GcProximityTerm(target)
    expected = -((expected_gc / 3.0 - target) ** 2)
    assert term.delta("", codon, 0) == pytest.approx(expected)


def test_gc_perfect_match_scores_zero() -> None:
    # target == 1/3 exactly matches a codon with a single G/C.
    term = GcProximityTerm(1.0 / 3.0)
    assert term.delta("", "ATG", 0) == pytest.approx(0.0)


@pytest.mark.parametrize("bad", [-0.1, 1.1, -1.0, 2.0])
def test_gc_target_out_of_range_rejected(bad: float) -> None:
    with pytest.raises(ValueError, match="gc target"):
        GcProximityTerm(bad)


@pytest.mark.parametrize("good", [0.0, 0.5, 1.0])
def test_gc_target_in_range_accepted(good: float) -> None:
    assert GcProximityTerm(good).target == good


def test_frozen_terms_are_immutable() -> None:
    cai = CaiTerm(_W)
    with pytest.raises((AttributeError, TypeError)):
        cai.name = "other"  # type: ignore[misc]
    gc = GcProximityTerm(0.5)
    with pytest.raises((AttributeError, TypeError)):
        gc.target = 0.1  # type: ignore[misc]


def test_score_over_scored_codons_matches_table_cai() -> None:
    # log-CAI * n_scored == score, so exp(score / n_scored) == table.cai(dna).
    term = CaiTerm(_W)
    dna = _backtranslate("MKLAADEFGHIKLPQRSTVWY")
    n_scored = sum(
        1
        for _, codon in iter_codons(dna)
        if CODON_TABLE[codon] != STOP and CODON_TABLE[codon] not in {"M", "W"}
    )
    assert n_scored > 0
    assert math.exp(term.score(dna) / n_scored) == pytest.approx(_TABLE.cai(dna))
