"""Tests for :class:`~bt4.constraints.gc_run.GcRunConstraint` (invariant #3).

The critical property mirrors the other local constraints: a sequence built
respecting ``ok_suffix`` must contain zero hard violations under ``validate``
(``ok_suffix <=> validate``), and the declared ``context_len`` must actually
suffice for the veto - so a GC run that straddles a codon boundary is still
caught. BT3 documented such an agreement and then broke it with a silent 12-nt
context cap; here it is property-tested.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.gc_run import GcRunConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)


def _build_feasible(protein: str, constraint: GcRunConstraint) -> str:
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


def _hard(constraint: GcRunConstraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN, max_run=st.integers(min_value=3, max_value=6))
def test_gc_run_ok_suffix_implies_validate_clean(protein: str, max_run: int) -> None:
    constraint = GcRunConstraint(max_run)
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


# --------------------------------------------------------------------------- #
# Positive detection: validate must catch real over-long GC runs.
# --------------------------------------------------------------------------- #


def test_validate_flags_over_long_run_with_correct_span() -> None:
    constraint = GcRunConstraint(3)
    # ATG CGC AAA -> positions 2..5 are G,C,G,C, a GC run of length 4.
    violations = list(constraint.validate("ATGCGCAAA"))
    assert len(violations) == 1
    v = violations[0]
    assert v.constraint == "gc_run"
    assert v.severity is Severity.HARD
    assert (v.start, v.end) == (2, 6)


def test_validate_flags_each_maximal_run_once() -> None:
    constraint = GcRunConstraint(2)
    # GGG (0-3) and CCCC (5-9) are two maximal over-long GC runs, split by the
    # non-GC "AT" in the middle.
    spans = sorted((v.start, v.end) for v in constraint.validate("GGGATCCCC"))
    assert spans == [(0, 3), (5, 9)]


def test_validate_leaves_at_limit_run_alone() -> None:
    # A GC run exactly equal to max_run is allowed.
    assert list(GcRunConstraint(3).validate("ATGCGA")) == []


def test_run_of_exactly_max_run_allowed_but_one_more_flagged() -> None:
    # AGCGCA: positions 1..4 are G,C,G,C -> a GC run of exactly 4.
    assert list(GcRunConstraint(4).validate("AGCGCA")) == []
    # AGCGCGA: positions 1..5 are G,C,G,C,G -> a GC run of 5 > 4.
    over = list(GcRunConstraint(4).validate("AGCGCGA"))
    assert [(v.start, v.end) for v in over] == [(1, 6)]


def test_mixed_gc_counts_as_one_run() -> None:
    # A fully mixed G/C stretch is a single run, not many short ones.
    six = list(GcRunConstraint(5).validate("GCGCGC"))
    assert [(v.start, v.end) for v in six] == [(0, 6)]
    # The same run is fine when the limit admits its full length.
    assert list(GcRunConstraint(6).validate("GCGCGC")) == []


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing runs are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_context_len_equals_max_run() -> None:
    assert GcRunConstraint(5).context_len() == 5
    assert GcRunConstraint(1).context_len() == 1


def test_ok_suffix_vetoes_boundary_crossing_run() -> None:
    constraint = GcRunConstraint(3)
    # prefix ends in GCG (a run of 3); a codon starting with C makes a run of 4.
    assert constraint.ok_suffix("ATGCG", "CAA") is False
    # A codon that does not extend the GC run is fine.
    assert constraint.ok_suffix("ATGCG", "AAA") is True


def test_ok_suffix_vetoes_run_formed_across_seam() -> None:
    constraint = GcRunConstraint(3)
    # prefix ends in GC (a run of 2); GCA on the seam completes GCGC -> 4.
    assert constraint.ok_suffix("ATGC", "GCA") is False


def test_ok_suffix_vetoes_run_inside_incoming_codon() -> None:
    constraint = GcRunConstraint(2)
    # The incoming codon GCG is itself a mixed GC run of length 3 > 2.
    assert constraint.ok_suffix("ATA", "GCG") is False


# --------------------------------------------------------------------------- #
# Names, scope, and configuration guards.
# --------------------------------------------------------------------------- #


def test_name_scope_and_penalty() -> None:
    constraint = GcRunConstraint(4)
    assert constraint.name == "gc_run"
    assert constraint.scope() is Scope.LOCAL
    assert constraint.penalty("GCG", "CGC") == 0.0


def test_rejects_non_positive_max_run() -> None:
    with pytest.raises(ValueError, match="max_run"):
        GcRunConstraint(0)
    with pytest.raises(ValueError, match="max_run"):
        GcRunConstraint(-3)
