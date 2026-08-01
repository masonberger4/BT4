"""Lagrangian-relaxation backend: local constraints AND one global count budget.

The exact DP (:mod:`bt4.optimize.exact_dp`) is the honest workhorse for local,
context-bounded objectives and constraints, but it cannot carry a
*whole-sequence* count budget -- a total GC count between ``gc_min``/``gc_max``,
or a total CpG count under a cap. A running global total is not a
bounded-trailing-context state, so folding it into the trellis would blow the
state space up (one state per attainable count per prefix). The CP-SAT backend
(:mod:`bt4.optimize.cpsat`) can model such a budget as a single linear
constraint, but it drops *all* local sequence constraints (homopolymer runs,
forbidden or restriction motifs, repeats). This backend fills exactly that gap:
it enforces one global count budget **and** the local constraints together, by
dualizing the budget into the per-codon objective and reusing the exact DP.

How it stays honest (CLAUDE.md invariant #6):

- The budget is moved into the objective with a multiplier ``lambda >= 0`` and
  the relaxed problem is solved by :func:`~bt4.optimize.exact_dp.solve_exact`
  with the *same* constraints and objective context. Every candidate the DP
  returns therefore already satisfies all local constraints -- that is the whole
  point of routing through the DP rather than an ILP.
- Each exact relaxed solve also yields a *valid dual bound*: for the maximization
  ``max f(seq) s.t. amount(seq) <= budget_max``, ``f(seq) - lambda*(amount(seq) -
  budget_max)`` is an upper bound on the true constrained optimum for any
  ``lambda >= 0`` (and symmetrically for a lower budget). We keep the tightest
  (minimum) such bound seen and pair it with the best budget-feasible candidate's
  *true* objective (recomputed from the codons, never trusted from the relaxed
  accumulator) to report a real, non-negative optimality gap.
- Integer problems have a Lagrangian duality gap, so we never claim
  ``PROVEN_OPTIMAL``: the certificate is ``GAP_BOUNDED`` with the honest computed
  gap (possibly ``0.0``), ``solver="lagrangian"``, and ``relaxed_terms`` naming
  the dualized budget. If a ``beam`` actually truncated any relaxed solve the
  dual bound is no longer trustworthy, so the certificate degrades to
  ``BEAM_TRUNCATED`` and carries no gap. With no budget at all nothing is
  relaxed and we return the exact DP's own result and certificate unchanged.

Determinism (invariant #7): the multiplier follows a fixed diminishing
subgradient schedule (``step = initial_step / (k + 1)``), the underlying DP is
deterministic, and ties are broken toward the lexicographically smaller DNA, so
identical inputs yield byte-identical output. No randomness, no wall-clock.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.optimize.exact_dp import InfeasibleError, SolveResult, solve_exact

__all__ = ["solve_lagrangian"]

_INITIAL_STEP: float = 4.0
"""Initial subgradient step size for the multiplier update.

The subgradient is normalized by the number of residues so it is O(1) regardless
of sequence length, and the step diminishes as ``_INITIAL_STEP / (k + 1)``. A
value of four lets ``lambda`` sweep the price range that flips realistic codon
choices within the default iteration budget; because the best budget-feasible
candidate is retained across all iterations, an over-large step only widens
exploration and never corrupts the returned solution.
"""

ScalarDelta = Callable[[str, str, int], float]


def solve_lagrangian(
    residues: Sequence[str],
    *,
    scalar_delta: ScalarDelta,
    constraints: Sequence[Constraint],
    amount: Callable[[str], int],
    budget_min: int | None = None,
    budget_max: int | None = None,
    beam: int | None = None,
    objective_context: int = 0,
    budget_name: str = "budget",
    max_iters: int = 50,
) -> SolveResult:
    """Solve the additive core plus one global count budget via Lagrangian dual.

    The budgeted quantity is ``amount(seq) = sum(amount(codon) for codon in seq)``
    and the budget is ``budget_min <= amount(seq) <= budget_max`` for whichever
    bounds are given. A single multiplier ``lambda >= 0`` dualizes the currently
    binding side: for an upper bound the relaxed per-codon objective is
    ``scalar_delta(...) - lambda*amount(codon)``, for a lower bound it is
    ``scalar_delta(...) + lambda*amount(codon)``. Each iteration solves the
    relaxed problem exactly with :func:`~bt4.optimize.exact_dp.solve_exact`
    (so all local constraints stay enforced), measures the delivered
    ``amount``, records the candidate if it is budget-feasible, and takes a
    projected subgradient step on ``lambda`` based on the budget violation. When
    both bounds are given, whichever side the delivered candidate violates is the
    side dualized next; a feasible candidate leaves the current side's multiplier
    to decay.

    Args:
        residues: Amino-acid letters to back-translate, **including** a trailing
            ``"*"`` as the final residue (the stop). Residue ``pos`` may use any
            codon in ``synonymous_codons(residues[pos])``.
        scalar_delta: Incremental *true* objective of placing a codon, called as
            ``scalar_delta(prefix, codon, pos)`` and oriented so larger is
            better. The returned ``objective_scalar`` is this term summed over the
            delivered codons -- never the relaxed accumulator.
        constraints: Hard local-feasibility rules, honored on every candidate
            because each relaxed solve runs through the exact DP with them in
            force. This is the advantage over CP-SAT, which drops them.
        amount: Per-codon contribution to the budgeted whole-sequence quantity,
            e.g. ``gc_count`` for a GC budget. Must depend only on the codon.
        budget_min: If given, require ``amount(seq) >= budget_min``.
        budget_max: If given, require ``amount(seq) <= budget_max``.
        beam: Passed straight through to the exact DP as a speed knob. If a beam
            actually truncates any relaxed solve the dual bound is unreliable and
            the certificate degrades to ``BEAM_TRUNCATED`` (no gap).
        objective_context: Trailing DNA context ``scalar_delta`` depends on, in
            characters (e.g. ``3`` for a pairwise objective). Passed through to
            the DP; the dualized ``amount`` term adds no context of its own.
        budget_name: Human-readable name of the dualized budget, used in the
            certificate's ``relaxed_terms`` and in ``InfeasibleError``.
        max_iters: Maximum number of subgradient iterations.

    Returns:
        A :class:`~bt4.optimize.exact_dp.SolveResult` whose ``dna`` translates
        back to ``residues`` and satisfies both the local constraints and the
        budget, whose ``objective_scalar`` is the true objective recomputed from
        the codons, and whose certificate is ``GAP_BOUNDED`` with the honest
        computed gap (or ``BEAM_TRUNCATED`` if a beam truncated a relaxed solve,
        or the DP's own certificate when no budget is relaxed).

    Raises:
        InfeasibleError: If a relaxed solve finds the local constraints
            infeasible (propagated from the DP), or if the budget itself admits
            no assignment. With ``beam`` unset the budget infeasibility is proven
            by an exact extremal solve; the error names ``budget_name``.
    """
    # No budget => nothing to relax; the exact DP is already the honest answer.
    if budget_min is None and budget_max is None:
        return solve_exact(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )

    scale = max(1.0, float(len(residues)))
    exact = beam is None

    def amount_total(dna: str) -> int:
        return sum(amount(dna[i : i + 3]) for i in range(0, len(dna), 3))

    def true_objective(dna: str) -> float:
        total = 0.0
        prefix = ""
        for pos, i in enumerate(range(0, len(dna), 3)):
            codon = dna[i : i + 3]
            total += scalar_delta(prefix, codon, pos)
            prefix += codon
        return total

    def is_feasible(amt: int) -> bool:
        if budget_min is not None and amt < budget_min:
            return False
        return not (budget_max is not None and amt > budget_max)

    # The best budget-feasible candidate seen: (true objective, DNA). A single
    # optional tuple keeps the two coupled fields narrowable together.
    best: tuple[float, str] | None = None

    def consider(dna: str) -> None:
        nonlocal best
        amt = amount_total(dna)
        if not is_feasible(amt):
            return
        value = true_objective(dna)
        # Prefer higher true objective; break ties toward the lexicographically
        # smaller DNA so the result is deterministic (invariant #7).
        if best is None or value > best[0] or (value == best[0] and dna < best[1]):
            best = (value, dna)

    beam_truncated = False

    # Extremal feasibility probes: the min-amount (and max-amount) locally-feasible
    # solutions. When the DP is exact these *prove* budget infeasibility, and
    # whenever an extreme happens to be budget-feasible it seeds a guaranteed
    # incumbent so a subgradient that never lands feasible cannot force a false
    # InfeasibleError.
    if budget_max is not None:
        low = solve_exact(
            residues,
            scalar_delta=_neg_amount(amount),
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )
        beam_truncated = beam_truncated or _truncated(low.certificate)
        if exact and amount_total(low.dna) > budget_max:
            raise InfeasibleError([budget_name])
        consider(low.dna)

    if budget_min is not None:
        high = solve_exact(
            residues,
            scalar_delta=_pos_amount(amount),
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )
        beam_truncated = beam_truncated or _truncated(high.certificate)
        if exact and amount_total(high.dna) < budget_min:
            raise InfeasibleError([budget_name])
        consider(high.dna)

    lam = 0.0
    side = "max" if budget_max is not None else "min"
    best_dual: float | None = None

    for k in range(max_iters):
        relaxed = _relaxed_delta(scalar_delta, amount, side, lam)
        result = solve_exact(
            residues,
            scalar_delta=relaxed,
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )
        proven = result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
        beam_truncated = beam_truncated or not proven
        amt = amount_total(result.dna)

        # A valid upper bound on the true constrained optimum -- but only when the
        # relaxed solve was exact (a truncated solve may undershoot the relaxed
        # max and hence is not a bound).
        if proven:
            if side == "max":
                assert budget_max is not None
                dual = result.objective_scalar + lam * float(budget_max)
            else:
                assert budget_min is not None
                dual = result.objective_scalar - lam * float(budget_min)
            if best_dual is None or dual < best_dual:
                best_dual = dual

        consider(result.dna)

        # Dualize whichever side the delivered candidate violates; a feasible
        # candidate keeps the current side so its multiplier decays.
        over = budget_max is not None and amt > budget_max
        under = budget_min is not None and amt < budget_min
        if over:
            nxt = "max"
        elif under:
            nxt = "min"
        else:
            nxt = side
        if nxt == "max":
            assert budget_max is not None
            subgrad = amt - budget_max
        else:
            assert budget_min is not None
            subgrad = budget_min - amt
        step = _INITIAL_STEP / (k + 1)
        lam = max(0.0, lam + step * (subgrad / scale))
        side = nxt

    if best is None:
        raise InfeasibleError([budget_name])
    best_true, best_dna = best

    if beam_truncated:
        certificate = OptimalityCertificate(
            status=OptimalityStatus.BEAM_TRUNCATED,
            solver="lagrangian",
            relaxed_terms=(budget_name,),
            detail=(
                f"dualized budget {budget_name!r} over {len(residues)} residues; "
                f"a beam (width {beam}) truncated a relaxed solve, so no dual "
                "bound is claimed"
            ),
        )
    else:
        if best_dual is None:
            gap = 1.0
        else:
            gap = (best_dual - best_true) / max(abs(best_dual), 1.0)
            gap = min(1.0, max(0.0, gap))
        certificate = OptimalityCertificate(
            status=OptimalityStatus.GAP_BOUNDED,
            solver="lagrangian",
            gap=gap,
            relaxed_terms=(budget_name,),
            detail=(
                f"dualized budget {budget_name!r} over {len(residues)} residues; "
                f"subgradient dual bound gives relative gap {gap:.3g}"
            ),
        )

    return SolveResult(dna=best_dna, objective_scalar=best_true, certificate=certificate)


def _relaxed_delta(
    scalar_delta: ScalarDelta, amount: Callable[[str], int], side: str, lam: float
) -> ScalarDelta:
    """Build the relaxed per-codon objective for one dualized budget side.

    Args:
        scalar_delta: The true incremental objective (larger is better).
        amount: Per-codon contribution to the budgeted quantity.
        side: ``"max"`` to dualize an upper bound (subtract ``lambda*amount``),
            ``"min"`` to dualize a lower bound (add ``lambda*amount``).
        lam: The current non-negative multiplier.

    Returns:
        A ``scalar_delta``-shaped callable folding the dualized budget into the
        per-codon score. The ``amount`` term reads only the codon, so it adds no
        trailing context to the trellis state.
    """
    sign = -1.0 if side == "max" else 1.0

    def relaxed(prefix: str, codon: str, pos: int) -> float:
        return scalar_delta(prefix, codon, pos) + sign * lam * float(amount(codon))

    return relaxed


def _neg_amount(amount: Callable[[str], int]) -> ScalarDelta:
    """Return a ``scalar_delta`` that maximizes ``-amount`` (minimizes amount)."""

    def delta(prefix: str, codon: str, pos: int) -> float:
        return -float(amount(codon))

    return delta


def _pos_amount(amount: Callable[[str], int]) -> ScalarDelta:
    """Return a ``scalar_delta`` that maximizes ``+amount`` (maximizes amount)."""

    def delta(prefix: str, codon: str, pos: int) -> float:
        return float(amount(codon))

    return delta


def _truncated(certificate: OptimalityCertificate) -> bool:
    """Return ``True`` iff a beam actually dropped states in the given solve."""
    return certificate.status is OptimalityStatus.BEAM_TRUNCATED
