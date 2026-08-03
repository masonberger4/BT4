"""Tests for the budgeted trellis backend (one global budget + local rules).

The backend's promise is the union the exact DP and CP-SAT each miss: a single
*whole-sequence* count budget (here a GC-count budget) enforced together with the
*local* sequence constraints CP-SAT drops. It does this with an amount-bucketed
exact DP, so a budget is enforced *exactly* -- there is no Lagrangian duality gap
and no interior budget window can be skipped (the confirmed audit bug).

The tests below pin the load-bearing behaviour: an interior two-sided budget is
feasible (the reproduction), one-sided budgets bind, local constraints are still
honoured, the certificate is honestly ``PROVEN_OPTIMAL`` when the solve is exact
and never proven under a beam, a genuinely infeasible budget raises rather than
returning a bogus answer, a beam never falsely reports infeasibility, and the
solve is deterministic. A no-budget call is cross-checked against
:func:`~bt4.optimize.exact_dp.solve_exact` directly.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from bt4 import api
from bt4._accel import gc_count, max_homopolymer_run
from bt4.constraints.rules import HomopolymerConstraint
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.genetic_code import STOP, synonymous_codons, translate
from bt4.optimize.exact_dp import InfeasibleError, solve_exact
from bt4.optimize.lagrangian import solve_lagrangian

ScalarDelta = Callable[[str, str, int], float]


def gc_amount(prefix: str, codon: str) -> int:
    """GC count of ``codon`` in the context-aware ``amount(prefix, codon)`` shape.

    GC is a per-base count that never straddles a codon boundary, so the prefix is
    ignored and ``budget_context`` is 0 -- this simply adapts ``gc_count`` to the
    signature :func:`~bt4.optimize.lagrangian.solve_lagrangian` now expects.
    """
    return gc_count(codon)


# A protein whose achievable total GC count (including the trailing stop) spans
# [17, 37]. The interior band [26, 28] excludes *both* extremes -- exactly the
# case the old linear Lagrangian penalty jumped over and falsely called
# infeasible (the confirmed audit reproduction).
_REPRO = "MKAILVDEQTRSFYNWHGP"


def gc_reward() -> ScalarDelta:
    """A simple scalar objective that rewards GC-rich codons.

    Making the unconstrained optimum GC-maximal guarantees an upper GC budget
    genuinely binds, so the budget has real work to do.
    """

    def delta(prefix: str, codon: str, pos: int) -> float:
        return float(gc_count(codon))

    return delta


def zero_delta() -> ScalarDelta:
    """A flat objective: every codon scores the same.

    With a flat objective the solver has no gradient pulling it toward the budget
    window, so landing inside the band is purely the budget DP's doing -- exactly
    what the interior-window bug broke.
    """

    def delta(prefix: str, codon: str, pos: int) -> float:
        return 0.0

    return delta


def min_achievable_gc(residues: Sequence[str]) -> int:
    """Lowest total GC count attainable over all synonymous assignments."""
    return sum(min(gc_count(c) for c in synonymous_codons(r)) for r in residues)


def max_achievable_gc(residues: Sequence[str]) -> int:
    """Highest total GC count attainable over all synonymous assignments."""
    return sum(max(gc_count(c) for c in synonymous_codons(r)) for r in residues)


def test_interior_two_sided_window_is_feasible_reproduction() -> None:
    # The confirmed audit reproduction: a two-sided band that excludes both
    # extremes of the achievable GC range. The old backend raised InfeasibleError
    # here; the bucketed DP must return a feasible sequence with GC in the band.
    residues = [*_REPRO, STOP]
    lo = min_achievable_gc(residues)
    hi = max_achievable_gc(residues)
    assert (lo, hi) == (17, 37)
    assert lo < 26 and hi > 28  # band is strictly interior

    result = solve_lagrangian(
        residues,
        scalar_delta=zero_delta(),
        constraints=[],
        amount=gc_amount,
        budget_min=26,
        budget_max=28,
        budget_name="gc_budget",
    )

    assert 26 <= gc_count(result.dna) <= 28
    assert translate(result.dna) == _REPRO + "*"
    # Enforced exactly => honestly proven optimal subject to the budget.
    cert = result.certificate
    assert cert.solver == "lagrangian"
    assert cert.status is OptimalityStatus.PROVEN_OPTIMAL
    assert cert.relaxed_terms == ()


def test_interior_two_sided_window_via_api_end_to_end() -> None:
    # The exact end-to-end reproduction command from the audit. The default
    # config carries max_homopolymer=6, so the pipeline routes this through the
    # Lagrangian/bucketed backend.
    result = api.optimize(_REPRO, api.OptimizeConfig(gc_min=26, gc_max=28))

    assert 26 <= gc_count(result.dna) <= 28
    assert translate(result.dna) == _REPRO + "*"
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert result.metrics.hard_violations == 0


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
        amount=gc_amount,
        budget_max=budget_max,
    )

    # Budget enforced on the delivered sequence.
    assert gc_count(result.dna) <= budget_max
    # Still a valid back-translation.
    assert translate(result.dna) == "LRALRSLRA*"
    # Honest certificate: the budget is enforced exactly, so nothing is relaxed
    # and the solve is proven optimal (no duality-gap hand-waving).
    cert = result.certificate
    assert cert.solver == "lagrangian"
    assert cert.status is OptimalityStatus.PROVEN_OPTIMAL
    assert cert.gap in (None, 0.0)
    assert cert.relaxed_terms == ()
    # objective_scalar is the true objective accumulated over the codons.
    true_obj = sum(
        delta("", result.dna[i : i + 3], i // 3)
        for i in range(0, len(result.dna), 3)
    )
    assert result.objective_scalar == true_obj


def test_lower_budget_only_binds_and_is_proven() -> None:
    residues = [*"LRALRSLRA", STOP]
    lo = min_achievable_gc(residues)
    hi = max_achievable_gc(residues)
    # A lower bound strictly above the minimum but below the maximum: it binds and
    # is feasible. A flat objective proves the budget DP lands in-band on its own.
    budget_min = (lo + hi) // 2 + 1
    assert lo < budget_min <= hi

    result = solve_lagrangian(
        residues,
        scalar_delta=zero_delta(),
        constraints=[],
        amount=gc_amount,
        budget_min=budget_min,
        budget_name="gc_min",
    )

    assert gc_count(result.dna) >= budget_min
    assert translate(result.dna) == "LRALRSLRA*"
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_both_bounds_interior_band_excludes_extremes() -> None:
    residues = [*"LRALRSLRA", STOP]
    lo = min_achievable_gc(residues)
    hi = max_achievable_gc(residues)
    # An interior band that excludes both the min and the max achievable GC.
    band_lo = lo + 1
    band_hi = hi - 1
    assert lo < band_lo <= band_hi < hi

    result = solve_lagrangian(
        residues,
        scalar_delta=zero_delta(),
        constraints=[],
        amount=gc_amount,
        budget_min=band_lo,
        budget_max=band_hi,
        budget_name="gc_budget",
    )

    assert band_lo <= gc_count(result.dna) <= band_hi
    assert translate(result.dna) == "LRALRSLRA*"
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_no_budget_delegates_to_exact_dp_unchanged() -> None:
    residues = [*"MARKLE", STOP]
    delta = gc_reward()

    relaxed = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_amount
    )
    exact = solve_exact(residues, scalar_delta=delta, constraints=[])

    # With no budget nothing is bucketed: identical result and certificate.
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
        amount=gc_amount,
        budget_max=budget_max,
    )

    # The whole point: the local homopolymer rule is enforced *and* the budget.
    assert max_homopolymer_run(result.dna) <= 3
    assert gc_count(result.dna) <= budget_max
    assert translate(result.dna) == "GGGGGG*"
    # And validate() (the whole-sequence audit) agrees there are no violations.
    assert list(homopolymer.validate(result.dna)) == []
    # Both budget and local constraint enforced exactly => proven optimal.
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_deterministic_identical_calls_match() -> None:
    residues = [*"LRALRSLRA", STOP]
    delta = gc_reward()
    budget_max = (min_achievable_gc(residues) + max_achievable_gc(residues)) // 2

    first = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_amount, budget_max=budget_max
    )
    second = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_amount, budget_max=budget_max
    )

    assert first.dna == second.dna
    assert first.objective_scalar == second.objective_scalar
    assert first.certificate.status == second.certificate.status
    assert first.certificate.gap == second.certificate.gap


def test_certificate_proven_when_exact_but_never_proven_under_beam() -> None:
    residues = [*"LRALRSLRA", STOP]
    delta = gc_reward()
    budget_max = (min_achievable_gc(residues) + max_achievable_gc(residues)) // 2

    exact = solve_lagrangian(
        residues, scalar_delta=delta, constraints=[], amount=gc_amount, budget_max=budget_max
    )
    assert exact.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert gc_count(exact.dna) <= budget_max

    # A beam narrow enough to genuinely truncate: the certificate must not claim
    # optimality, yet the budget stays enforced exactly and the answer is valid.
    beamed = solve_lagrangian(
        residues,
        scalar_delta=delta,
        constraints=[],
        amount=gc_amount,
        budget_max=budget_max,
        beam=1,
    )
    assert beamed.certificate.status is OptimalityStatus.BEAM_TRUNCATED
    assert beamed.certificate.status is not OptimalityStatus.PROVEN_OPTIMAL
    assert gc_count(beamed.dna) <= budget_max
    assert translate(beamed.dna) == "LRALRSLRA*"


def test_beam_never_falsely_infeasible_on_interior_window() -> None:
    # The audit's second confirmed issue: under a beam, the old code could hit
    # ``if best is None: raise InfeasibleError`` without having proven anything.
    # The bucketed DP's sound amount prune runs before the beam cap, so a beam can
    # only trade optimality, never feasibility -- an interior, genuinely feasible
    # band must still yield an in-band solution under a beam.
    residues = [*_REPRO, STOP]

    result = solve_lagrangian(
        residues,
        scalar_delta=zero_delta(),
        constraints=[],
        amount=gc_amount,
        budget_min=26,
        budget_max=28,
        beam=1,
        budget_name="gc_budget",
    )

    assert 26 <= gc_count(result.dna) <= 28
    assert translate(result.dna) == _REPRO + "*"
    assert result.certificate.status is not OptimalityStatus.PROVEN_OPTIMAL


def test_infeasible_budget_raises() -> None:
    residues = [*"GGGGGG", STOP]
    delta = gc_reward()

    # A GC ceiling strictly below the minimum achievable GC: no assignment can
    # satisfy it, and the sound amount prune proves that immediately.
    impossible = min_achievable_gc(residues) - 1

    raised = False
    try:
        solve_lagrangian(
            residues,
            scalar_delta=delta,
            constraints=[],
            amount=gc_amount,
            budget_max=impossible,
            budget_name="gc_max",
        )
    except InfeasibleError as exc:
        raised = True
        assert "gc_max" in exc.constraints
    assert raised


def test_infeasible_budget_still_raises_under_beam() -> None:
    # Budget infeasibility is proven by the *sound* amount prune, which runs
    # before any beam cap, so it must still raise (correctly) even with a beam --
    # this is genuine infeasibility, not a beam artefact.
    residues = [*"GGGGGG", STOP]
    impossible = min_achievable_gc(residues) - 1

    raised = False
    try:
        solve_lagrangian(
            residues,
            scalar_delta=gc_reward(),
            constraints=[],
            amount=gc_amount,
            budget_max=impossible,
            budget_name="gc_max",
            beam=1,
        )
    except InfeasibleError as exc:
        raised = True
        assert "gc_max" in exc.constraints
    assert raised
