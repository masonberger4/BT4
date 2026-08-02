"""Tests for :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint`.

The load-bearing property (CLAUDE.md invariant #3, "ok_suffix <=> validate"): a
sequence built respecting ``ok_suffix`` contains zero hard violations under
``validate``. This works precisely because ``MaxRepeatConstraint`` is honest
about being GLOBAL -- its ``ok_suffix`` reads the *whole* prefix rather than a
bounded window (CLAUDE.md section 10.1), so a repeat whose earlier copy lies far
back is still vetoed. The remaining tests pin the k = max_length + 1 detection of
direct, inverted, and palindromic repeats and the sufficiency of that single
k-mer length.
"""

from __future__ import annotations

import sys

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.base import Constraint
from bt4.constraints.max_repeat import MaxRepeatConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=60)


def _build_feasible(protein: str, constraint: Constraint) -> str:
    """Greedily first-fit a codon per residue that passes ``ok_suffix``.

    Returns the DNA built so far; stops early (returns the valid prefix) at the
    first residue where no synonymous codon is feasible.
    """
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if constraint.ok_suffix(dna, codon):
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


@given(protein=_PROTEIN, max_length=st.integers(min_value=4, max_value=9))
def test_ok_suffix_implies_validate_clean(protein: str, max_length: int) -> None:
    constraint = MaxRepeatConstraint(max_length)
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


# --------------------------------------------------------------------------- #
# Positive detection: validate must catch each repeat type.
# --------------------------------------------------------------------------- #


def test_validate_flags_direct_repeat() -> None:
    # max_length=3 => k=4. "AAAC" appears at 0 and 6; no other 4-mer repeats and
    # no reverse-complement (RC(AAAC)="GTTT") appears, so both occurrences are
    # flagged as a direct repeat.
    constraint = MaxRepeatConstraint(3)
    violations = list(constraint.validate("AAACCCAAAC"))
    assert [(v.start, v.end) for v in violations] == [(0, 4), (6, 10)]
    assert all(v.constraint == "max_repeat" for v in violations)
    assert all(v.severity is Severity.HARD for v in violations)
    assert all("direct" in v.detail for v in violations)


def test_validate_flags_inverted_repeat() -> None:
    # max_length=3 => k=4. "AAGG" at 0 and RC("AAGG")="CCTT" at 6 form an inverted
    # repeat; both arms are flagged and nothing else is.
    constraint = MaxRepeatConstraint(3)
    violations = list(constraint.validate("AAGGACCCTT"))
    assert [(v.start, v.end) for v in violations] == [(0, 4), (6, 10)]
    assert all("inverted" in v.detail for v in violations)


def test_validate_flags_palindrome() -> None:
    # "AATT" equals its own reverse complement: a palindromic 4-mer.
    constraint = MaxRepeatConstraint(3)
    violations = list(constraint.validate("AATT"))
    assert [(v.start, v.end) for v in violations] == [(0, 4)]
    assert "palindrome" in violations[0].detail
    assert violations[0].constraint == "max_repeat"


def test_validate_clean_sequence_passes() -> None:
    # No 4-mer repeats, no palindromic 4-mer, no reverse-complement coincidence.
    assert list(MaxRepeatConstraint(3).validate("ACGGATTCAG")) == []


# --------------------------------------------------------------------------- #
# k = max_length + 1 sufficiency: exactly max_length is allowed, +1 is not.
# --------------------------------------------------------------------------- #


def test_repeat_of_length_max_length_is_allowed() -> None:
    # "AAC" (length 3) repeats at 0 and 4, but no 4-mer repeats, so with
    # max_length=3 the sequence is clean.
    assert list(MaxRepeatConstraint(3).validate("AACGAAC")) == []


def test_repeat_of_length_max_length_plus_one_is_flagged() -> None:
    # "AACG" (length 4 = max_length + 1) repeats at 0 and 5 -> a direct repeat.
    violations = list(MaxRepeatConstraint(3).validate("AACGCAACG"))
    assert [(v.start, v.end) for v in violations] == [(0, 4), (5, 9)]
    assert all("direct" in v.detail for v in violations)


# --------------------------------------------------------------------------- #
# ok_suffix reads the whole prefix (global): a far-apart repeat is vetoed.
# --------------------------------------------------------------------------- #


def test_ok_suffix_vetoes_far_apart_direct_repeat() -> None:
    constraint = MaxRepeatConstraint(3)
    # prefix carries "AAGG" at its very start and ends in "A"; the codon "AGG"
    # completes a second "AAGG" ten bases later -- caught only by reading the
    # whole prefix.
    assert constraint.ok_suffix("AAGGTTCACGA", "AGG") is False
    # A codon that completes no repeat is accepted.
    assert constraint.ok_suffix("AAGGTTCACGA", "TTG") is True


def test_ok_suffix_vetoes_boundary_crossing_palindrome() -> None:
    constraint = MaxRepeatConstraint(3)
    # prefix ends in "AAT"; the codon supplying the final "T" completes the
    # palindrome "AATT" across the seam.
    assert constraint.ok_suffix("CGAAT", "TCG") is False
    assert constraint.ok_suffix("CGAAT", "GCG") is True


# --------------------------------------------------------------------------- #
# Contract surface: scope, context, name, penalty, configuration.
# --------------------------------------------------------------------------- #


def test_scope_is_global() -> None:
    assert MaxRepeatConstraint(6).scope() is Scope.GLOBAL


def test_context_len_is_unbounded_sentinel() -> None:
    # ok_suffix reads the whole prefix; context_len is the unbounded sentinel.
    assert MaxRepeatConstraint(6).context_len() == sys.maxsize


def test_name_and_penalty() -> None:
    constraint = MaxRepeatConstraint(6)
    assert constraint.name == "max_repeat"
    assert constraint.penalty("AAA", "CCC") == 0.0


def test_rejects_bad_params() -> None:
    with pytest.raises(ValueError, match="max_length"):
        MaxRepeatConstraint(0)
