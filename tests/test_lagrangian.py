"""Tests for the Lagrangian-relaxation backend (one global budget + local rules).

The backend's promise is the union the exact DP and CP-SAT each miss: a single
*whole-sequence* count budget (here a GC-count budget) enforced together with the
*local* sequence constraints CP-SAT drops. The tests below pin the load-bearing
behaviour -- budget respected, local constraints still honoured, an honest
certificate, determinism, and a proven-infeasible budget raising rather than
returning a bogus answer -- and cross-check a no-budget call against
:func:`~bt4.optimize.exact_dp.solve_exact` directly.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from bt4._accel import gc_count, max_homopolymer_run
from bt4.constraints.rules import HomopolymerConstraint
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.genetic_code import STOP, synonymous_codons, translate
from bt4.optimize.exact_dp import InfeasibleError, solve_exact
from bt4.optimize.lagrangian import solve_lagrangian

ScalarDelta = Callable[[str, str, int], float]


def gc_reward() -> ScalarDelta:
    """A simple scalar objective that rewards GC-rich codons.

    Making the unconstrained optimum GC-maximal guarantees an upper GC budget
    genuinely binds, so the relaxation has real work to do.
    """

    def delta(prefix: str, codon: str, pos: int) -> float:
        return float(gc_count(codon))

    return delta


def min_achievable_gc(residues: Sequence[str]) -> int:
    """Lowest total GC count attainable over all synonymous assignments."""
    return sum(min(gc_count(c) for c in synonymous_codons(r)) for r in residues)


def max_achievable_gc(residues: Sequence[str]) -> int:
    """Highest total GC count attainable over all synonymous assignments."""
    return sum(max(gc_count(c) for c in synonymous_codons(r)) for r in residues)


def test_gc_upper_budget_is_respected_with_honest_certificate() -> None:
    residues = [*"LRALRSLRA", STOP]
    delta = gc_reward()

    lo = min_achievable_gc(residues)
    hi = max_achievable_gc(residues)
    # A budget strictly between the min and the (GC-maximal) unconstrained optimum
    # so it genuinely binds but stays feasible.
    budget_max = (lo + hi) // 2
    assert lo <= budget_max < hi

    result = solve_lagrangian(
        residues,
        scalar_delta=delta,
        constraints=[],
        amount=gc_count,
        budget_max=budget_max,
    )

    # Budget enforced on the delivered sequence.
    assert gc_count(result.dna) <= budget_max
    # Still a valid back-translation.
    assert translate(result.dna) == "LRALRSLRA*"
    # Honest certificate: named solver, budget flagged as relaxed, real gap.
    cert = result.certificate
    assert cert.solver == "lagrangian"
    assert cert.status is OptimalityStatus.GAP_BOUNDED
    assert cert.relaxed_terms == ("budget",)
    assert cert.gap is not None
    assert 0.0 <= cert.gap <= 1.0
    # objective_scalar is the true objective recomputed from the codons.
    true_obj = sum(
        delta("", result.dna[i : i + 3], i // 3)
        for i in range(0, len(result.dna), 3)
    )
    assert result.objective_scalar == true_obj


def test_no_budget_delegates_to_exact_dp_unchanged() -> None:
    residues = [*"MARKLE", STOP]
    delta = gc_reward()

    relaxed = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_count
    )
    exact = solve_exact(residues, scalar_delta=delta, constraints=[])

    # With no budget nothing is relaxed: identical result and certificate.
    assert relaxed.dna == exact.dna
    assert relaxed.objective_scalar == exact.objective_scalar
    assert relaxed.certificate.is_proven_optimal
    assert relaxed.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_local_constraint_is_honored_alongside_budget() -> None:
    # Poly-glycine is homopolymer-prone (GGG|GGG spells GGGGGG); this is exactly
    # the local constraint CP-SAT's budget backend would silently ignore.
    residues = [*"GGGGGG", STOP]
    delta = gc_reward()
    homopolymer = HomopolymerConstraint(3)

    budget_max = max_achievable_gc(residues) - 2  # bind the GC budget too

    result = solve_lagrangian(
        residues,
        scalar_delta=delta,
        constraints=[homopolymer],
        amount=gc_count,
        budget_max=budget_max,
    )

    # The whole point: the local homopolymer rule is enforced *and* the budget.
    assert max_homopolymer_run(result.dna) <= 3
    assert gc_count(result.dna) <= budget_max
    assert translate(result.dna) == "GGGGGG*"
    # And validate() (the whole-sequence audit) agrees there are no violations.
    assert list(homopolymer.validate(result.dna)) == []


def test_deterministic_identical_calls_match() -> None:
    residues = [*"LRALRSLRA", STOP]
    delta = gc_reward()
    budget_max = (min_achievable_gc(residues) + max_achievable_gc(residues)) // 2

    first = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_count, budget_max=budget_max
    )
    second = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_count, budget_max=budget_max
    )

    assert first.dna == second.dna
    assert first.objective_scalar == second.objective_scalar
    assert first.certificate.gap == second.certificate.gap


def test_infeasible_budget_raises() -> None:
    residues = [*"GGGGGG", STOP]
    delta = gc_reward()

    # A GC ceiling strictly below the minimum achievable GC: no assignment can
    # satisfy it, and the exact extremal probe proves that immediately.
    impossible = min_achievable_gc(residues) - 1

    raised = False
    try:
        solve_lagrangian(
            residues,
            scalar_delta=delta,
            constraints=[],
            amount=gc_count,
            budget_max=impossible,
            budget_name="gc_max",
        )
    except InfeasibleError as exc:
        raised = True
        assert "gc_max" in exc.constraints
    assert raised
