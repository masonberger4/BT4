"""Tests for graceful constraint relaxation (CLAUDE.md §4.2, defect A.4).

:mod:`bt4.domain.relax` turns a hard constraint into a visible-but-non-binding
soft one so an otherwise-infeasible instance degrades gracefully instead of
dead-ending. These tests pin the contract: the soft form never vetoes, keeps its
occurrences visible (downgraded to SOFT), and relaxation is opt-in per constraint.
"""

from __future__ import annotations

from bt4.constraints.kozak import InternalStartConstraint
from bt4.constraints.rules import HomopolymerConstraint
from bt4.domain import Severity, SoftConstraint, is_relaxable, relax_constraint


def test_soft_constraint_never_vetoes() -> None:
    hard = InternalStartConstraint()
    soft = SoftConstraint(hard)
    # A strong internal ATG the hard rule vetoes must NOT be vetoed by the soft form.
    prefix, codon = "GCCACC", "ATG"  # A at -3, ATG whose +4 (next base) can be G
    # The hard rule vetoes some extension; the soft rule never does.
    assert soft.ok_suffix(prefix, codon) is True
    assert soft.name == hard.name
    assert soft.scope() == hard.scope()
    assert soft.context_len() == hard.context_len()


def test_soft_constraint_downgrades_violations_to_soft() -> None:
    hard = InternalStartConstraint()
    soft = SoftConstraint(hard)
    dna = "GCCACCATGGCC"  # a strong internal ATG at index 6
    hard_viols = list(hard.validate(dna))
    soft_viols = list(soft.validate(dna))
    assert hard_viols, "sanity: the hard rule flags this sequence"
    assert all(v.severity is Severity.HARD for v in hard_viols)
    # Same occurrences, but visible as SOFT rather than dropped -- nothing hidden.
    assert len(soft_viols) == len(hard_viols)
    assert all(v.severity is Severity.SOFT for v in soft_viols)
    assert {v.start for v in soft_viols} == {v.start for v in hard_viols}


def test_relaxation_is_opt_in() -> None:
    # A constraint that defines relax() opts in; one that does not is never relaxed.
    assert is_relaxable(InternalStartConstraint()) is True
    assert is_relaxable(HomopolymerConstraint(6)) is False


def test_relax_constraint_uses_declared_relax_then_falls_back() -> None:
    # Opt-in constraint: relax_constraint returns its own soft form.
    kozak = InternalStartConstraint()
    relaxed = relax_constraint(kozak)
    assert isinstance(relaxed, SoftConstraint)
    assert relaxed.name == "internal_start"
    # Non-opt-in constraint: the generic wrapper still relaxes it if asked directly.
    generic = relax_constraint(HomopolymerConstraint(6))
    assert isinstance(generic, SoftConstraint)
    assert generic.ok_suffix("AAAAAA", "AAA") is True  # would-be homopolymer veto lifted
