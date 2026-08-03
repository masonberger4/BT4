"""Budgeted trellis backend: local constraints AND one global count budget.

The exact DP (:mod:`bt4.optimize.exact_dp`) is the honest workhorse for local,
context-bounded objectives and constraints, but it cannot on its own carry a
*whole-sequence* count budget -- a total GC count between ``gc_min``/``gc_max``,
or a total CpG/UpA dinucleotide count in a target range. The CP-SAT backend
(:mod:`bt4.optimize.cpsat`) models such a budget as a single linear constraint
but drops *all* local sequence constraints (homopolymer runs, forbidden or
restriction motifs, repeats) and any pairwise objective term. This backend fills
exactly that gap: it enforces one global count budget **and** the local
constraints (and any pairwise objective) together.

The budgeted quantity's per-codon contribution may be *context-aware*: a
dinucleotide (e.g. CpG) can straddle a codon boundary -- the ``C`` ending one
codon and the ``G`` starting the next -- so the contribution of a codon depends
on the last ``budget_context`` bases of the prefix, exactly as
:meth:`bt4.objectives.dinucleotide.DinucleotideTerm.delta` attributes each
occurrence to the codon holding its END base. ``budget_context`` is folded into
the state's context length so prefixes that differ only in their trailing bases
stay unmerged and a straddling count is attributed exactly once (a per-base GC
count is context-free, so it uses ``budget_context == 0`` and this reduces to the
original behaviour).

How it works -- an amount-bucketed exact DP:

- The trellis carries one layer per residue (including the trailing stop). Its
  state is the pair *(trailing DNA context, cumulative budgeted amount)* where
  the context length is the union of every constraint's declared ``context_len``,
  the objective's ``context_len``, and the budget's ``budget_context``. Because
  the running amount is part of the state, the budget is enforced *exactly*, not
  relaxed: at the final layer only states whose total amount lands in
  ``[budget_min, budget_max]`` survive, and the highest-objective survivor is
  reconstructed.
- Two prefixes that share a state key are interchangeable for every future
  feasibility decision (same trailing context => same ``ok_suffix`` behaviour and
  same objective ``delta``) *and* for the budget (same cumulative amount), so
  keeping only the best of each is exact -- the same merge argument the exact DP
  makes, extended by the amount dimension.
- A sound bound prune keeps the amount dimension small: at each layer a state is
  dropped only when, given the min/max amount still attainable over the remaining
  residues, it *provably cannot* reach the budget window. This never drops a
  state that could still become budget-feasible, so it preserves optimality while
  bounding the bucket count to (roughly) the width of the budget window. When the
  amount is context-aware, the per-residue min/max are taken over *all* codons of
  the residue **and all possible trailing contexts** (see :func:`_budget_contexts`),
  so the bound stays valid under any preceding context -- possibly loose, never
  too tight, which is exactly what soundness requires.

How it stays honest (CLAUDE.md invariants #6 and #7):

- With ``beam is None`` the DP explores the full (sound-pruned) state space, so
  the returned sequence is the *proven* optimum subject to the budget and every
  local constraint. The certificate is ``PROVEN_OPTIMAL`` -- nothing is relaxed.
- A ``beam`` caps the states kept per layer (a speed knob); if it actually drops
  states the certificate degrades to ``BEAM_TRUNCATED``. Crucially, the sound
  amount prune runs *before* the beam cap and guarantees every surviving state
  can still reach the budget window, so a beam can never leave a *feasible*
  instance with no in-window solution -- it only trades away optimality, never
  feasibility. A beam that never bites still reports ``PROVEN_OPTIMAL`` (matching
  the exact DP).
- ``InfeasibleError`` is raised only when infeasibility is genuinely proven: when
  a layer empties under the constraints' ``ok_suffix`` vetoes (local
  infeasibility, always exact) or when the sound amount prune -- which never
  drops a reachable-feasible state -- empties a layer (budget infeasibility, sound
  even under a beam). It is never raised merely because a beam truncated the
  search.
- Determinism (#7): states are iterated in sorted order, codons in
  :func:`~bt4.domain.genetic_code.synonymous_codons` order, and every tie is
  broken toward the lexicographically smaller DNA. No randomness, no wall-clock.

The certificate keeps ``solver="lagrangian"`` as a stable backend label (the
pipeline and its tests key on it); the certificate ``detail`` states plainly that
the mechanism is an amount-bucketed exact budget DP.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons
from bt4.optimize.exact_dp import InfeasibleError, SolveResult, solve_exact

__all__ = ["solve_lagrangian"]

ScalarDelta = Callable[[str, str, int], float]
# Per-codon contribution to the budgeted quantity, given the running prefix (only
# its last ``budget_context`` chars are read) and the codon being placed.
AmountFn = Callable[[str, str], int]

# One trellis layer: (trailing context, cumulative amount) -> (best scalar, DNA).
_StateKey = tuple[str, int]
_Layer = dict[_StateKey, tuple[float, str]]


def _budget_contexts(budget_context: int) -> list[str]:
    """Every trailing-context string an amount function might read.

    ``amount(prefix, codon)`` reads only the last ``budget_context`` chars of
    ``prefix``. Across a whole run those trailing views range over *all* ACGT
    strings of length ``0..budget_context`` -- the shorter ones occur at the
    start, before ``budget_context`` bases have been placed (Python's
    ``prefix[-k:]`` on a short prefix is the whole prefix). Bounding the
    per-residue amounts over this full set makes the reachability prune sound
    under any preceding context (possibly loose, never too tight).

    Args:
        budget_context: Number of trailing prefix chars the amount reads. ``0``
            means the amount is context-free (a per-base count like GC), for
            which the only relevant context is the empty string.

    Returns:
        The list of candidate trailing contexts, always including ``""``.
    """
    if budget_context <= 0:
        return [""]
    contexts: list[str] = [""]
    current: list[str] = [""]
    for _ in range(budget_context):
        current = [c + base for c in current for base in "ACGT"]
        contexts.extend(current)
    return contexts


def solve_lagrangian(
    residues: Sequence[str],
    *,
    scalar_delta: ScalarDelta,
    constraints: Sequence[Constraint],
    amount: AmountFn,
    budget_min: int | None = None,
    budget_max: int | None = None,
    beam: int | None = None,
    objective_context: int = 0,
    budget_context: int = 0,
    budget_name: str = "budget",
    max_iters: int = 50,
) -> SolveResult:
    """Solve the additive core plus one global count budget, exactly.

    The budgeted quantity is ``amount(seq) = sum(amount(prefix_i, codon_i))`` over
    the codons, where ``prefix_i`` is the sequence placed before codon ``i``, and
    the budget is ``budget_min <= amount(seq) <= budget_max`` for whichever bounds
    are given. When neither bound is set nothing is budgeted and the call
    delegates unchanged to :func:`~bt4.optimize.exact_dp.solve_exact`. Otherwise
    an amount-bucketed exact DP (see the module docstring) enforces the budget
    *exactly* alongside every local constraint and the (possibly pairwise)
    additive objective, so unlike a linear Lagrangian relaxation there is no
    duality gap and no interior budget window can be skipped over.

    Args:
        residues: Amino-acid letters to back-translate, **including** a trailing
            ``"*"`` as the final residue (the stop). Residue ``pos`` may use any
            codon in ``synonymous_codons(residues[pos])``.
        scalar_delta: Incremental *true* objective of placing a codon, called as
            ``scalar_delta(prefix, codon, pos)`` and oriented so larger is
            better. It must depend only on the last ``objective_context``
            characters of ``prefix``. The returned ``objective_scalar`` is this
            term accumulated over the delivered codons.
        constraints: Hard local-feasibility rules, honored on every candidate via
            their ``ok_suffix`` veto. This is the advantage over CP-SAT, which
            drops them.
        amount: Per-codon contribution to the budgeted whole-sequence quantity,
            called as ``amount(prefix, codon)`` and depending only on the last
            ``budget_context`` chars of ``prefix``. A context-free per-base count
            (e.g. GC) ignores ``prefix`` and uses ``budget_context == 0``; a
            dinucleotide count reads one trailing base to catch a
            boundary-straddling occurrence (``budget_context == 1``).
        budget_min: If given, require ``amount(seq) >= budget_min``.
        budget_max: If given, require ``amount(seq) <= budget_max``.
        beam: If ``None``, run the full exact bucketed DP (proven optimal). An int
            caps the states kept per layer for speed; if it drops any state the
            certificate reports ``BEAM_TRUNCATED``. A beam never turns a feasible
            instance infeasible (see the module docstring).
        objective_context: Trailing DNA context ``scalar_delta`` depends on, in
            characters (e.g. ``3`` for a pairwise objective).
        budget_context: Trailing DNA context ``amount`` depends on, in characters
            (``0`` for a context-free count, ``1`` for a dinucleotide). Folded into
            the state's context length so straddling contributions stay exact.
        budget_name: Human-readable name of the budget, used in the certificate
            and in ``InfeasibleError`` when the budget itself is infeasible.
        max_iters: Retained for backward compatibility with the previous
            subgradient implementation; the exact bucketed DP takes no iteration
            budget and ignores it.

    Returns:
        A :class:`~bt4.optimize.exact_dp.SolveResult` whose ``dna`` translates
        back to ``residues`` and satisfies both the local constraints and the
        budget, whose ``objective_scalar`` is the true objective accumulated over
        the delivered codons, and whose certificate is ``PROVEN_OPTIMAL`` (exact
        solve, budget in force) or ``BEAM_TRUNCATED`` (a beam dropped states).

    Raises:
        InfeasibleError: If the local constraints admit no assignment (named by
            the offending constraints) or the budget admits no assignment (named
            ``budget_name``). Both are genuinely proven -- never raised merely
            because a beam truncated the search.
    """
    # No budget => nothing to bucket; the exact DP is already the honest answer.
    if budget_min is None and budget_max is None:
        return solve_exact(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )
    return _solve_budgeted(
        residues,
        scalar_delta=scalar_delta,
        constraints=constraints,
        amount=amount,
        budget_min=budget_min,
        budget_max=budget_max,
        beam=beam,
        objective_context=objective_context,
        budget_context=budget_context,
        budget_name=budget_name,
    )


def _solve_budgeted(
    residues: Sequence[str],
    *,
    scalar_delta: ScalarDelta,
    constraints: Sequence[Constraint],
    amount: AmountFn,
    budget_min: int | None,
    budget_max: int | None,
    beam: int | None,
    objective_context: int,
    budget_context: int,
    budget_name: str,
) -> SolveResult:
    """Amount-bucketed exact DP enforcing one global count budget exactly.

    See the module docstring for the state, the soundness of the amount prune,
    and the certificate semantics. This is the honest core the public
    :func:`solve_lagrangian` dispatches to whenever a budget is present.

    Args:
        residues: Amino-acid letters (with a trailing stop) to back-translate.
        scalar_delta: The true incremental objective (larger is better).
        constraints: Hard local-feasibility rules honored via ``ok_suffix``.
        amount: Per-codon contribution to the budgeted quantity, reading the last
            ``budget_context`` chars of the prefix.
        budget_min: Lower bound on the total amount, or ``None``.
        budget_max: Upper bound on the total amount, or ``None``.
        beam: Per-layer state cap (a speed knob), or ``None`` for exact.
        objective_context: Trailing context ``scalar_delta`` depends on.
        budget_context: Trailing context ``amount`` depends on.
        budget_name: Human-readable budget name for the certificate/error.

    Returns:
        The budget-feasible :class:`~bt4.optimize.exact_dp.SolveResult`.

    Raises:
        InfeasibleError: On genuinely proven local or budget infeasibility.
    """
    context_len = max(
        [objective_context, budget_context, *(c.context_len() for c in constraints)]
    )
    n = len(residues)

    # Per-residue codon choices and their min/max amount, plus suffix sums giving
    # the min/max amount still attainable over residues [i:]. These drive the
    # sound reachability prune. The amount can be context-aware (a dinucleotide
    # straddling a codon boundary), so each residue's min/max is taken over every
    # codon AND every possible trailing context (see :func:`_budget_contexts`):
    # bounding over all contexts keeps the reachability bound valid regardless of
    # what precedes the residue -- sound (possibly loose), never too tight.
    budget_ctxs = _budget_contexts(budget_context)
    codon_choices = [synonymous_codons(r) for r in residues]
    per_min: list[int] = []
    per_max: list[int] = []
    for cods in codon_choices:
        amts = [amount(ctx, c) for c in cods for ctx in budget_ctxs]
        per_min.append(min(amts))
        per_max.append(max(amts))
    remaining_min = [0] * (n + 1)
    remaining_max = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        remaining_min[i] = remaining_min[i + 1] + per_min[i]
        remaining_max[i] = remaining_max[i + 1] + per_max[i]

    def can_reach(amt: int, next_pos: int) -> bool:
        """Whether ``amt`` can still land in the budget over residues [next_pos:].

        Sound: returns ``False`` only when *no* assignment of the remaining
        residues can bring the total into the budget window, so pruning on it
        never drops a state that could still become budget-feasible.
        """
        lo = amt + remaining_min[next_pos]
        hi = amt + remaining_max[next_pos]
        if budget_max is not None and lo > budget_max:
            return False
        return not (budget_min is not None and hi < budget_min)

    layer: _Layer = {("", 0): (0.0, "")}
    truncated = False

    for pos in range(n):
        cods = codon_choices[pos]
        generated: _Layer = {}
        # Sorted iteration keeps the build order deterministic; the merge rule is
        # order-independent regardless.
        for (_ctx, amt), (score, dna) in sorted(layer.items()):
            for codon in cods:
                if not all(c.ok_suffix(dna, codon) for c in constraints):
                    continue
                new_dna = dna + codon
                new_amt = amt + amount(dna, codon)
                new_score = score + scalar_delta(dna, codon, pos)
                new_ctx = new_dna[-context_len:] if context_len > 0 else ""
                key = (new_ctx, new_amt)
                cur = generated.get(key)
                if cur is None or _wins(new_score, new_dna, cur):
                    generated[key] = (new_score, new_dna)

        # Empty under ok_suffix => genuine local infeasibility (always exact).
        if not generated:
            raise InfeasibleError([c.name for c in constraints])

        # Sound amount prune: drop states that provably cannot reach the budget.
        next_pos = pos + 1
        pruned: _Layer = {
            key: value
            for key, value in generated.items()
            if can_reach(key[1], next_pos)
        }
        # Empty after a sound prune => genuine budget infeasibility (sound even
        # under a beam, because the prune never drops a reachable-feasible state).
        if not pruned:
            raise InfeasibleError([budget_name])

        # Optional speed knob: keep only the top-``beam`` states per layer. The
        # sound prune above already guarantees every kept state can still reach
        # the window, so this only trades optimality, never feasibility.
        if beam is not None and len(pruned) > beam:
            kept = sorted(pruned.items(), key=lambda kv: (-kv[1][0], kv[1][1]))[:beam]
            pruned = dict(kept)
            truncated = True

        layer = pruned

    # At the final layer every surviving state's amount lies in the budget window
    # (``can_reach`` with no residues left forces ``budget_min <= amt <=
    # budget_max``), so the best-objective survivor is the answer. Ties break
    # toward the lexicographically smaller DNA (invariant #7).
    best_score, best_dna = min(layer.values(), key=lambda sv: (-sv[0], sv[1]))

    if truncated:
        certificate = OptimalityCertificate(
            status=OptimalityStatus.BEAM_TRUNCATED,
            solver="lagrangian",
            detail=(
                f"amount-bucketed budget DP over {n} residues under budget "
                f"{budget_name!r}; a beam (width {beam}) truncated a layer, so "
                "optimality is not proven (the budget is still enforced exactly)"
            ),
        )
    else:
        certificate = OptimalityCertificate.proven(
            "lagrangian",
            detail=(
                f"amount-bucketed budget DP over {n} residues: proven optimal "
                f"subject to the {budget_name!r} budget and all local constraints "
                f"(context K={context_len})"
            ),
        )

    return SolveResult(dna=best_dna, objective_scalar=best_score, certificate=certificate)


def _wins(new_score: float, new_dna: str, current: tuple[float, str]) -> bool:
    """Return True iff ``(new_score, new_dna)`` should replace ``current``.

    A candidate wins on a strictly higher scalar, or on an equal scalar with a
    lexicographically smaller DNA (deterministic tie-break, invariant #7).

    Args:
        new_score: The candidate's accumulated objective.
        new_dna: The candidate's coding-sequence prefix.
        current: The incumbent ``(score, dna)`` for the same state key.

    Returns:
        ``True`` if the candidate should replace the incumbent.
    """
    cur_score, cur_dna = current
    if new_score != cur_score:
        return new_score > cur_score
    return new_dna < cur_dna
