"""Tests for the repeat constraints (CLAUDE.md invariant #3).

The critical property, mirrored from :mod:`tests.test_constraints`: a sequence
built respecting ``ok_suffix`` must contain zero hard violations under
``validate`` (``ok_suffix <=> validate``), and the declared ``context_len`` must
actually suffice for the veto -- so a repeat that straddles a codon boundary is
still caught. BT3 documented this agreement and then broke it with a silent
context cap; here it is property-tested for tandem and inverted repeats.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.base import Constraint
from bt4.constraints.repeats import InvertedRepeatConstraint, TandemRepeatConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)


def _build_feasible(protein: str, constraints: list[Constraint]) -> str:
    """Greedily first-fit a codon per residue that passes every ``ok_suffix``.

    Returns the DNA built so far; stops early (returns the valid prefix) at the
    first residue where no synonymous codon is feasible.
    """
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if all(c.ok_suffix(dna, codon) for c in constraints):
                dna += codon
                break
        else:  # no synonymous codon worked - dead end, stop with the valid prefix.
            break
    return dna


def _hard(constraint: Constraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(
    protein=_PROTEIN,
    unit_len=st.integers(min_value=1, max_value=3),
    copies=st.integers(min_value=2, max_value=4),
)
def test_tandem_ok_suffix_implies_validate_clean(
    protein: str, unit_len: int, copies: int
) -> None:
    constraint = TandemRepeatConstraint(unit_len, copies)
    dna = _build_feasible(protein, [constraint])
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


@given(
    protein=_PROTEIN,
    stem=st.integers(min_value=2, max_value=5),
    loop_max=st.integers(min_value=0, max_value=3),
)
def test_inverted_ok_suffix_implies_validate_clean(
    protein: str, stem: int, loop_max: int
) -> None:
    constraint = InvertedRepeatConstraint(stem, loop_max)
    dna = _build_feasible(protein, [constraint])
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    assert translate(dna) == protein[: len(dna) // 3]


@given(
    protein=_PROTEIN,
    unit_len=st.integers(min_value=1, max_value=2),
    copies=st.integers(min_value=2, max_value=3),
    stem=st.integers(min_value=2, max_value=4),
    loop_max=st.integers(min_value=0, max_value=2),
)
def test_combined_constraints_ok_suffix_implies_validate_clean(
    protein: str, unit_len: int, copies: int, stem: int, loop_max: int
) -> None:
    constraints: list[Constraint] = [
        TandemRepeatConstraint(unit_len, copies),
        InvertedRepeatConstraint(stem, loop_max),
    ]
    dna = _build_feasible(protein, constraints)
    assume(len(dna) >= 3)
    for constraint in constraints:
        assert _hard(constraint, dna) == []


# --------------------------------------------------------------------------- #
# Positive detection: validate must catch real repeats.
# --------------------------------------------------------------------------- #


def test_tandem_validate_flags_mononucleotide_run() -> None:
    # unit_len=1, copies=3 bans any base repeated 3x. "AAAA" contains two such
    # length-3 windows (starts 0 and 1).
    constraint = TandemRepeatConstraint(1, 3)
    violations = list(constraint.validate("GGAAAAGG"))
    spans = sorted((v.start, v.end) for v in violations)
    assert spans == [(2, 5), (3, 6)]
    assert all(v.constraint == "tandem_repeat" for v in violations)
    assert all(v.severity is Severity.HARD for v in violations)


def test_tandem_validate_flags_dinucleotide_repeat() -> None:
    # unit_len=2, copies=2 bans a dinucleotide repeated twice, e.g. "ATAT".
    constraint = TandemRepeatConstraint(2, 2)
    violations = list(constraint.validate("GGATATGG"))
    assert [(v.start, v.end) for v in violations] == [(2, 6)]


def test_tandem_validate_clean_when_below_copies() -> None:
    # "ATAT" is only two copies; requiring three copies leaves it feasible.
    assert list(TandemRepeatConstraint(2, 3).validate("GGATATGG")) == []


def test_inverted_validate_flags_palindrome() -> None:
    # stem=2, loop_max=0: "GATC" is X="GA", revcomp("GA")="TC".
    constraint = InvertedRepeatConstraint(2, 0)
    violations = list(constraint.validate("GGGATCGG"))
    assert [(v.start, v.end) for v in violations] == [(2, 6)]
    assert violations[0].constraint == "inverted_repeat"


def test_inverted_validate_flags_hairpin_with_loop() -> None:
    # stem=2, loop_max=2: "GA" + "TT" (loop) + "TC" = "GATTTC" is a hairpin.
    constraint = InvertedRepeatConstraint(2, 2)
    starts = sorted(v.start for v in constraint.validate("CCGATTTCCC"))
    assert 2 in starts  # the GATTTC hairpin opens at index 2.


def test_inverted_validate_clean_without_stem() -> None:
    # A homopolymer has no arm whose reverse complement (all-T) appears
    # downstream, so no hairpin is detected.
    assert list(InvertedRepeatConstraint(3, 0).validate("AAAAAAAAA")) == []


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing repeats are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_tandem_ok_suffix_vetoes_boundary_crossing_run() -> None:
    constraint = TandemRepeatConstraint(1, 3)
    assert constraint.context_len() == 2
    # prefix ends in AA; the incoming codon's leading A completes AAA across the
    # seam and must be vetoed.
    assert constraint.ok_suffix("GGAA", "AGG") is False
    # A codon that does not extend the run is fine.
    assert constraint.ok_suffix("GGAA", "CGG") is True


def test_tandem_ok_suffix_vetoes_boundary_crossing_dinucleotide() -> None:
    constraint = TandemRepeatConstraint(2, 2)
    assert constraint.context_len() == 3
    # prefix ends in ATA; the codon 'TCC' completes ATAT across the seam.
    assert constraint.ok_suffix("GGATA", "TCC") is False
    assert constraint.ok_suffix("GGATA", "GCC") is True


def test_tandem_ok_suffix_allows_repeat_fully_in_prefix() -> None:
    # AAA lives entirely inside the already-feasible prefix; the incoming codon
    # neither extends nor introduces a repeat, so it is not vetoed.
    constraint = TandemRepeatConstraint(1, 3)
    assert constraint.ok_suffix("AAA", "CGT") is True


def test_inverted_ok_suffix_vetoes_boundary_crossing_hairpin() -> None:
    constraint = InvertedRepeatConstraint(3, 0)
    assert constraint.context_len() == 5
    # prefix ends in GGGCC; the codon 'CAA' supplies the final C, completing the
    # palindrome GGGCCC across the seam (revcomp("GGG") == "CCC").
    assert constraint.ok_suffix("TTGGGCC", "CAA") is False
    # A codon that does not close the stem is fine.
    assert constraint.ok_suffix("TTGGGCC", "AAA") is True


def test_inverted_ok_suffix_vetoes_boundary_crossing_hairpin_with_loop() -> None:
    constraint = InvertedRepeatConstraint(2, 2)
    # context_len = 2*2 + 2 - 1 = 5.
    assert constraint.context_len() == 5
    # prefix ends in GATT; codon 'TCC' supplies the final C completing GA+TT+TC.
    assert constraint.ok_suffix("CCGATT", "TCC") is False
    assert constraint.ok_suffix("CCGATT", "GCC") is True


# --------------------------------------------------------------------------- #
# Configuration / degenerate cases.
# --------------------------------------------------------------------------- #


def test_tandem_rejects_bad_params() -> None:
    with pytest.raises(ValueError, match="unit_len"):
        TandemRepeatConstraint(0, 3)
    with pytest.raises(ValueError, match="copies"):
        TandemRepeatConstraint(2, 1)


def test_inverted_rejects_bad_params() -> None:
    with pytest.raises(ValueError, match="stem"):
        InvertedRepeatConstraint(0, 2)
    with pytest.raises(ValueError, match="loop_max"):
        InvertedRepeatConstraint(3, -1)


def test_names_scope_and_penalty() -> None:
    tandem = TandemRepeatConstraint(2, 3)
    inverted = InvertedRepeatConstraint(4, 1)
    assert tandem.name == "tandem_repeat"
    assert inverted.name == "inverted_repeat"
    assert tandem.scope() is Scope.LOCAL
    assert inverted.scope() is Scope.LOCAL
    assert tandem.penalty("AAA", "CCC") == 0.0
    assert inverted.penalty("AAA", "CCC") == 0.0


def test_defaults() -> None:
    assert TandemRepeatConstraint(2).copies == 3
    assert InvertedRepeatConstraint(4).loop_max == 0
