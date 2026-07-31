"""Tests for the codon-pair bias table and objective term.

Covers:

* :func:`build_codon_pair_table` on a hand-checkable synthetic CDS set: an
  over-represented pair scores ``CPS > 0`` and an under-represented pair scores
  ``CPS < 0`` (constructed so the sign is robust to the default pseudocount), and
  every score is finite;
* the length-multiple-of-three guard;
* :class:`CpbTerm` metadata and the ``delta == score`` invariant (#4) over random
  back-translated proteins;
* determinism of table construction.
"""

from __future__ import annotations

import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4.biomodels.codon.pairs import CodonPairTable, build_codon_pair_table
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons
from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons
from bt4.objectives.codon_pair import CpbTerm

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=120)

# Hand-checkable corpus over the amino-acid pair Phe(F)-Tyr(Y). Every CDS is
# "F Y stop". The codon marginals are balanced (each F and Y codon appears 4
# times), but the codon pairs are anti-correlated with the "matched" diagonal:
# TTT prefers TAC and TTC prefers TAT, while TTT-TAT is deliberately rare.
#
# For pair (A, B) with X = aa(A), Y = aa(B):
#     expected = (N(A) * N(B)) / (N(X) * N(Y)) * N(XY)
# Counts (no pseudocount): N(TTT)=N(TTC)=N(TAT)=N(TAC)=4, N(F)=N(Y)=8, N(F,Y)=8.
#   (TTT, TAC): observed 3, expected (4*4)/(8*8)*8 = 2.0  -> CPS > 0 (preferred)
#   (TTT, TAT): observed 1, expected 2.0                  -> CPS < 0 (avoided)
# Both signs survive the default pseudocount = 1.0 (see the module docstring).
_OVER_UNDER_CDS: list[str] = (
    ["TTTTACTAA"] * 3  # TTT TAC TAA  (F Y stop) -> over-represented (TTT, TAC)
    + ["TTCTATTAA"] * 3  # TTC TAT TAA  (F Y stop)
    + ["TTTTATTAA"] * 1  # TTT TAT TAA  (F Y stop) -> under-represented (TTT, TAT)
    + ["TTCTACTAA"] * 1  # TTC TAC TAA  (F Y stop)
)

# A broader fixed corpus for the objective-term tests (content is irrelevant to
# delta == score, which holds for any table, but a real table is used per spec).
_TERM_CDS: list[str] = [
    "ATGAAACTGGCAGCAGATGAATAA",  # M K L A A D E stop
    "ATGTTTTATTGCGGCGGACCGTAA",  # M F Y C G G P stop
    "ATGCTGCTGGAAGATAAACAGTAA",  # M L L E D K Q stop
]
_TERM_TABLE = build_codon_pair_table(_TERM_CDS)


def _backtranslate(protein: str) -> str:
    """Deterministically pick the first synonymous codon per residue."""
    return "".join(synonymous_codons(aa)[0] for aa in protein)


def _accumulated_delta(term: CpbTerm, dna: str) -> float:
    """Sum the term's deltas over growing real prefixes (as the DP would)."""
    acc = 0.0
    prefix = ""
    for pos, codon in iter_codons(dna):
        acc += term.delta(prefix, codon, pos)
        prefix += codon
    return acc


def test_build_returns_codon_pair_table() -> None:
    table = build_codon_pair_table(_OVER_UNDER_CDS)
    assert isinstance(table, CodonPairTable)
    assert table.name == "codon_pair"


def test_over_and_under_represented_signs() -> None:
    table = build_codon_pair_table(_OVER_UNDER_CDS)
    # Over-represented pair (observed 3 vs expected 2): positive CPS.
    assert table.score("TTT", "TAC") > 0.0
    # Under-represented pair (observed 1 vs expected 2): negative CPS.
    assert table.score("TTT", "TAT") < 0.0


def test_all_scores_are_finite() -> None:
    table = build_codon_pair_table(_OVER_UNDER_CDS)
    assert table.scores  # non-empty
    assert all(math.isfinite(cps) for cps in table.scores.values())


def test_score_is_case_insensitive_and_zero_for_unknown() -> None:
    table = build_codon_pair_table(_OVER_UNDER_CDS)
    assert table.score("ttt", "tac") == pytest.approx(table.score("TTT", "TAC"))
    # A pair that never occurs in the corpus is neutral.
    assert table.score("GGG", "CCC") == 0.0


def test_pseudocount_zero_still_finite_for_observed_pairs() -> None:
    table = build_codon_pair_table(_OVER_UNDER_CDS, pseudocount=0.0)
    assert all(math.isfinite(cps) for cps in table.scores.values())
    # Signs are even sharper without smoothing.
    assert table.score("TTT", "TAC") > 0.0
    assert table.score("TTT", "TAT") < 0.0


def test_length_not_multiple_of_three_raises() -> None:
    with pytest.raises(ValueError, match="multiple of three"):
        build_codon_pair_table(["ATGAAA", "ATGAA"])


def test_negative_pseudocount_raises() -> None:
    with pytest.raises(ValueError, match="pseudocount"):
        build_codon_pair_table(_OVER_UNDER_CDS, pseudocount=-1.0)


def test_invalid_dna_raises() -> None:
    with pytest.raises(ValueError, match="ACGT"):
        build_codon_pair_table(["ATGNNNTAA"])


def test_term_metadata() -> None:
    term = CpbTerm(_TERM_TABLE.scores)
    assert term.name == "codon_pair"
    assert term.scope() is Scope.PAIRWISE
    assert term.context_len() == 3


def test_term_first_codon_delta_is_zero() -> None:
    term = CpbTerm(_TERM_TABLE.scores)
    assert term.delta("", "ATG", 0) == 0.0


@given(protein=_PROTEIN)
def test_term_delta_equals_score(protein: str) -> None:
    term = CpbTerm(_TERM_TABLE.scores)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


def test_term_delta_reads_previous_codon() -> None:
    term = CpbTerm(_TERM_TABLE.scores)
    # delta at pos>0 uses only the trailing (previous) codon of the prefix.
    prev = "TTT"
    codon = "TAC"
    expected = _TERM_TABLE.score(prev, codon)
    assert term.delta("AAA" + prev, codon, 1) == pytest.approx(expected)
    assert term.delta("GGGCCC" + prev, codon, 5) == pytest.approx(expected)


def test_frozen_term_is_immutable() -> None:
    term = CpbTerm(_TERM_TABLE.scores)
    with pytest.raises((AttributeError, TypeError)):
        term.name = "other"  # type: ignore[misc]


def test_build_is_deterministic() -> None:
    first = build_codon_pair_table(_OVER_UNDER_CDS)
    second = build_codon_pair_table(_OVER_UNDER_CDS)
    assert dict(first.scores) == dict(second.scores)
