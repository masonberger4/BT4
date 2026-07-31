"""Tests for the concrete constraints (CLAUDE.md invariant #3).

The critical property: a sequence built respecting ``ok_suffix`` must contain
zero hard violations under ``validate`` (``ok_suffix <=> validate``), and the
declared ``context_len`` must actually suffice for the veto - so a violation
that straddles a codon boundary is still caught. BT3 documented this agreement
and then broke it with a silent 12-nt context cap; here it is property-tested.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.base import Constraint
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)
# A small pool of forbidden motifs (kept short so feasible builds usually exist).
_MOTIFS = ("GAATTC", "GGATCC", "AAG", "TTT", "CG")


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


@given(protein=_PROTEIN, max_run=st.integers(min_value=3, max_value=5))
def test_homopolymer_ok_suffix_implies_validate_clean(protein: str, max_run: int) -> None:
    constraint = HomopolymerConstraint(max_run)
    dna = _build_feasible(protein, [constraint])
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


@given(protein=_PROTEIN, rc=st.booleans())
def test_forbidden_motif_ok_suffix_implies_validate_clean(protein: str, rc: bool) -> None:
    constraint = ForbiddenMotifConstraint(_MOTIFS, reverse_complement=rc)
    dna = _build_feasible(protein, [constraint])
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []


@given(
    protein=_PROTEIN,
    max_run=st.integers(min_value=3, max_value=5),
    rc=st.booleans(),
)
def test_combined_constraints_ok_suffix_implies_validate_clean(
    protein: str, max_run: int, rc: bool
) -> None:
    constraints: list[Constraint] = [
        HomopolymerConstraint(max_run),
        ForbiddenMotifConstraint(_MOTIFS, reverse_complement=rc),
    ]
    dna = _build_feasible(protein, constraints)
    assume(len(dna) >= 3)
    for constraint in constraints:
        assert _hard(constraint, dna) == []


# --------------------------------------------------------------------------- #
# Positive detection: validate must catch real violations.
# --------------------------------------------------------------------------- #


def test_homopolymer_validate_flags_long_run() -> None:
    constraint = HomopolymerConstraint(3)
    violations = list(constraint.validate("ATGAAAACCC"))  # AAAA is a run of 4.
    assert len(violations) == 1
    v = violations[0]
    assert v.constraint == "homopolymer"
    assert v.severity is Severity.HARD
    assert (v.start, v.end) == (3, 7)


def test_homopolymer_validate_flags_each_maximal_run_once() -> None:
    constraint = HomopolymerConstraint(2)
    # AAA (0-3) and TTTT (6-10) are two maximal over-long runs.
    violations = list(constraint.validate("AAAGCTTTT"))
    spans = sorted((v.start, v.end) for v in violations)
    assert spans == [(0, 3), (5, 9)]


def test_homopolymer_validate_clean_at_boundary_length() -> None:
    # A run exactly equal to max_run is allowed.
    assert list(HomopolymerConstraint(3).validate("AAAGGG")) == []


def test_forbidden_motif_validate_flags_motif_and_rc() -> None:
    # AAG (fwd) and its reverse complement CTT should both be caught with rc=True.
    constraint = ForbiddenMotifConstraint(("AAG",), reverse_complement=True)
    fwd = list(constraint.validate("GGGAAGGGG"))
    assert [v.constraint for v in fwd] == ["forbidden_motif"]
    rc = list(constraint.validate("GGGCTTGGG"))
    assert [v.constraint for v in rc] == ["forbidden_motif"]
    # Without rc, the reverse complement is NOT banned.
    no_rc = ForbiddenMotifConstraint(("AAG",), reverse_complement=False)
    assert list(no_rc.validate("GGGCTTGGG")) == []


def test_forbidden_motif_validate_reports_all_occurrences() -> None:
    constraint = ForbiddenMotifConstraint(("CG",))
    starts = sorted(v.start for v in constraint.validate("CGATCGATCG"))
    assert starts == [0, 4, 8]


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing violations are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_homopolymer_ok_suffix_vetoes_boundary_crossing_run() -> None:
    constraint = HomopolymerConstraint(3)
    # prefix ends in AAA; adding a codon starting with A would make a run of 4.
    assert constraint.ok_suffix("GGGAAA", "AGG") is False
    # A codon that does not extend the run is fine.
    assert constraint.ok_suffix("GGGAAA", "CGG") is True


def test_forbidden_motif_ok_suffix_vetoes_boundary_crossing_motif() -> None:
    constraint = ForbiddenMotifConstraint(("GAATTC",))
    # prefix ends in GAATT; the codon 'CAA' completes GAATTC across the seam.
    assert constraint.ok_suffix("GGGGAATT", "CAA") is False
    # context_len is exactly maxlen - 1 == 5, which is what the veto inspects.
    assert constraint.context_len() == 5
    # A codon that does not complete the motif is allowed.
    assert constraint.ok_suffix("GGGGAATT", "GAA") is True


def test_forbidden_motif_ok_suffix_allows_motif_fully_in_prefix() -> None:
    # A motif entirely inside the (already-feasible) prefix, not touching the new
    # codon, is not this extension's fault - ok_suffix only judges next_codon.
    constraint = ForbiddenMotifConstraint(("CG",))
    assert constraint.ok_suffix("ACG", "AAA") is True
    # But a motif overlapping the incoming codon is vetoed.
    assert constraint.ok_suffix("AAC", "GAA") is False


# --------------------------------------------------------------------------- #
# Degenerate / configuration cases.
# --------------------------------------------------------------------------- #


def test_empty_motifs_is_inert() -> None:
    constraint = ForbiddenMotifConstraint(())
    assert constraint.context_len() == 0
    assert constraint.ok_suffix("", "ATG") is True
    assert constraint.ok_suffix("GAATTC", "GAA") is True
    assert list(constraint.validate("GAATTCGGATCC")) == []


def test_homopolymer_rejects_non_positive_max_run() -> None:
    with pytest.raises(ValueError, match="max_run"):
        HomopolymerConstraint(0)
    with pytest.raises(ValueError, match="max_run"):
        HomopolymerConstraint(-2)


def test_forbidden_motif_rejects_invalid_motif() -> None:
    with pytest.raises(ValueError):
        ForbiddenMotifConstraint(("GAXTTC",))


def test_names_and_scope() -> None:
    from bt4.domain.scope import Scope

    homo = HomopolymerConstraint(4)
    motif = ForbiddenMotifConstraint(("GAATTC",))
    assert homo.name == "homopolymer"
    assert motif.name == "forbidden_motif"
    assert homo.scope() is Scope.LOCAL
    assert motif.scope() is Scope.LOCAL
    assert homo.penalty("AAA", "CCC") == 0.0
    assert motif.penalty("AAA", "CCC") == 0.0


def test_reverse_complement_adds_distinct_motif() -> None:
    constraint = ForbiddenMotifConstraint(("AAG",), reverse_complement=True)
    # AAG -> CTT; both must be banned by ok_suffix on the seam.
    assert constraint.ok_suffix("GGGA", "AGG") is False  # completes AAG
    assert constraint.ok_suffix("GGGC", "TTG") is False  # completes CTT
