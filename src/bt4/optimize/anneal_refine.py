"""Incremental simulated-annealing refinement over synonymous codon swaps.

The exact DP (:mod:`bt4.optimize.exact_dp`), CP-SAT, and the Lagrangian dual all
handle objectives that decompose over a *bounded trailing context*. The genuinely
non-local terms BT4 cares about -- 5' folding free energy, whole-sequence splice
risk, a learned expression head -- do not. They are only ever *scored*, never
per-codon-decomposable, so BT4 optimizes them in a **refinement layer**: start
from a feasible sequence and take synonymous-only single-codon swaps, accepting
or rejecting each by Metropolis on a caller-supplied whole-sequence ``score``.

This module is that engine. It is honest and cheap by construction:

* **Round-trip preserved (invariant #1).** Every move replaces one codon with a
  *synonymous* codon of the same residue, so ``translate(dna)`` is invariant
  under refinement. The seed is checked to translate to ``residues`` at the door,
  so the returned DNA does too - by construction, never by luck.
* **No new hard violations (invariant #5).** A proposed swap is accepted only if
  it introduces no hard violation. The check is *incremental*: changing the codon
  at position ``pos`` can only affect a constraint's ``ok_suffix`` verdict for
  that codon and the few following codons whose declared context still reaches
  back into it (``reach = ceil(context/3)`` codons). We re-run ``ok_suffix`` over
  exactly that bounded window on the mutated sequence; if any constraint vetoes,
  the move is rejected. Because everything outside the window is unchanged and was
  already feasible, an accepted state is feasible everywhere, so the hard-
  violation count can only fall or stay flat - it never rises.
* **Reaching farther without weakening #5: block moves + parallel tempering.**
  A *single-codon* move with a strict global-constraint gate (never raise the
  whole-sequence hard count) cannot cross a barrier that needs a **coordinated
  multi-codon** change - it can leave a dispersed max-repeat / uORF in place
  (CLAUDE.md 9, Phase 3). Two opt-in mechanisms widen the reach:
  - **Block moves** (``block_size`` / ``block_prob``) propose synonymous swaps at
    several positions *at once*, so a repeat that only clears when two copies move
    together becomes reachable in one step. A block candidate is checked for local
    feasibility over the *union* of affected windows and then through the same
    global gate, so invariant #5 still holds move-by-move. **Block moves always
    full-``score`` re-score** (never ``delta_score``): summing per-position deltas
    is valid only for additive, disjoint-context terms, and the non-local terms
    this engine exists to serve (folding / splice / expression) are neither.
  - **Parallel tempering** (``replicas`` / ``temps`` / ``swap_every``) runs several
    replicas at different temperatures; hot replicas accept uphill moves and swap
    their (already-feasible) configuration into the cold chain, so a barrier is
    crossed without any single chain ever accepting a hard-count increase. Every
    replica passes **both** gates against **its own** current global count, and
    every configuration ever visited has a global count ``<=`` the seed's, so the
    delivered result (chosen best over all replicas, lower global count first)
    still satisfies invariant #5. **Feasibility floor (honest):** a repeat pinned
    to synonymously-immovable bases (Met ``ATG`` / Trp ``TGG``, or a base-locked
    degenerate position) is unremovable by *any* synonymous scheme - neither a
    block move nor a hot replica can clear it, and it is reported as a residual,
    never silently claimed clean.
  Both default off; with ``block_size=1``, ``block_prob=0.0``, ``replicas=1`` and
  ``swap_every=0`` the engine reproduces the single-chain trajectory below
  byte-for-byte (invariant #7).
* **Incremental scoring, O(context) per move (CLAUDE.md 7 / 10.8).** When the
  caller supplies ``delta_score`` - the score change of a single swap - each move
  costs only what that callable costs (O(context) for a local term, a bounded
  window for folding), never a full O(L) re-score. BT3's SA recomputed the whole
  objective *and* full validation per move; that quadratic blow-up is exactly what
  this design refuses. Without ``delta_score`` the engine falls back to calling
  the whole-sequence ``score`` once per proposal, which is O(L) per iteration -
  correct but quadratic overall; pass ``delta_score`` for real work.
* **Deterministic (invariant #7).** All randomness flows through a single
  ``random.Random(seed)``; no global RNG, no wall-clock. Identical inputs and seed
  yield byte-identical output.
* **Honest certificate.** Simulated annealing is a heuristic: it can escape local
  optima but never proves global optimality. The result therefore always carries
  an ``OptimalityCertificate`` of status ``HEURISTIC`` - it never claims
  ``PROVEN_OPTIMAL``.

The returned ``objective_scalar`` is the whole-sequence ``score`` recomputed on
the delivered DNA (reported == computed, invariant #2), not a drifting move
accumulator.
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable, Sequence

from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.optimize.exact_dp import SolveResult

__all__ = ["anneal_refine"]

Score = Callable[[str], float]
"""Whole-sequence objective, oriented so larger is better."""

DeltaScore = Callable[[str, int, str, str], float]
"""Incremental score of one synonymous swap.

Called as ``delta_score(dna, pos, old_codon, new_codon)`` where ``dna`` is the
*current* sequence (before the swap) and ``pos`` is the codon index. It must
return ``score(dna_with_new_codon) - score(dna)`` and should read only a bounded
window of ``dna`` around ``pos`` (that is what keeps a move O(context)).
"""


def anneal_refine(
    seed_dna: str,
    residues: Sequence[str],
    score: Score,
    constraints: Sequence[Constraint],
    *,
    global_constraints: Sequence[Constraint] = (),
    iterations: int = 2000,
    seed: int = 0,
    delta_score: DeltaScore | None = None,
    temp0: float = 1.0,
    cooling: float = 0.995,
    block_size: int = 1,
    block_prob: float = 0.0,
    replicas: int = 1,
    temps: Sequence[float] | None = None,
    swap_every: int = 0,
) -> SolveResult:
    """Refine a feasible coding sequence by synonymous-swap simulated annealing.

    Args:
        seed_dna: A feasible starting coding sequence, one codon per residue. It
            must translate exactly to ``residues`` (checked); this anchors the
            round-trip invariant so every accepted state translates back too.
        residues: Amino-acid letters, **including** the trailing ``"*"`` stop,
            used to enumerate synonymous options per position via
            :func:`~bt4.domain.genetic_code.synonymous_codons`.
        score: Whole-sequence objective, oriented so **larger is better**. This is
            where non-local terms (folding, splice, learned expression) plug in.
            It is called once on the seed and once on the delivered result; during
            the loop it is used only when ``delta_score`` is not supplied.
        constraints: Hard-feasibility rules (the ``Constraint`` contract). Each
            move is rejected if it would introduce a new hard violation, checked
            incrementally via ``ok_suffix`` over the affected window. These are the
            *local* (bounded-context) constraints; the seed must already satisfy
            them (checked at the door).
        global_constraints: Non-local (``Scope.GLOBAL``) constraints -- e.g. the
            dispersed :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint` --
            whose veto reads the whole sequence and so cannot be checked by a
            bounded ``ok_suffix`` window. A move is rejected if it *raises* the
            total number of hard violations these report (invariant #5), computed
            by a whole-sequence ``validate`` on the candidate. Unlike ``constraints``
            the seed **may** violate these (that is what refinement drives down);
            the count can only fall or stay flat, never rise. Empty by default, so
            callers that only refine local/folding objectives are unaffected.
        iterations: Number of proposed swaps (``>= 0``). ``0`` returns the seed
            unchanged (still a valid, honest ``HEURISTIC`` result).
        seed: Seed for the local ``random.Random`` driving proposals and
            acceptance; fixes the whole trajectory (invariant #7).
        delta_score: Optional incremental score of a single swap (see
            :data:`DeltaScore`). When given, each move is O(context); when
            ``None``, the engine re-scores the whole candidate with ``score``
            (O(L) per iteration - correct but quadratic overall).
        temp0: Initial annealing temperature (``>= 0``). ``0`` degenerates to
            greedy hill-climbing (only improving moves accepted).
        cooling: Geometric cooling factor applied to the temperature every
            iteration; must lie in ``(0, 1]``. Applied per replica, per step.
        block_size: Maximum number of positions a single **block move** may swap
            at once (``>= 1``). ``1`` (the default) disables block moves entirely
            (they need ``>= 2`` movable positions). When ``> 1`` and a block move
            is drawn, between 2 and ``min(block_size, #movable)`` distinct movable
            positions are proposed together and the candidate is **full-``score``
            re-scored** (never ``delta_score``), because summing per-position
            deltas is only valid for additive disjoint-context terms.
        block_prob: Probability in ``[0, 1]`` that any given step is a block move
            rather than a single-codon move (``0.0`` default = never; then no RNG
            draw is spent on the decision, so the single-chain trajectory is
            byte-identical to before this argument existed). Only has effect when
            ``block_size >= 2`` and at least two positions are movable.
        replicas: Number of parallel-tempering replicas (``>= 1``). ``1`` (default)
            is a single chain identical to before. Each replica keeps its own
            ``(dna, score, global-count, temperature)`` and passes both feasibility
            gates against **its own** current global count, so invariant #5 holds
            per replica; the delivered result is the best over all replicas (lower
            global count first, then higher score).
        temps: Optional explicit per-replica starting temperatures (length must
            equal ``replicas``, each ``>= 0``). When ``None``, every replica starts
            at ``temp0`` (so ``replicas=1`` reproduces the single-chain schedule);
            supply an ascending ladder (e.g. ``(0.1, 0.5, 2.0)``) for a real
            temperature spread. Each replica's temperature is still cooled by
            ``cooling`` on every step it takes.
        swap_every: Attempt adjacent replica-exchange swaps every ``swap_every``
            iterations (``0`` default = never; also a no-op when ``replicas == 1``).
            Exchanges use the standard replica-exchange Metropolis criterion and
            only relabel already-feasible configurations, so no swap can introduce
            a hard violation.

    Returns:
        A :class:`~bt4.optimize.exact_dp.SolveResult` whose ``dna`` is the best
        (highest-``score``) feasible sequence visited, whose ``objective_scalar``
        is ``score`` recomputed on that DNA, and whose certificate has status
        :attr:`~bt4.domain.certificate.OptimalityStatus.HEURISTIC` and solver
        ``"anneal_refine"``.

    Raises:
        ValueError: If ``seed_dna`` does not translate to ``residues``, if its
            length is not ``3 * len(residues)``, or if ``iterations``, ``temp0``,
            ``cooling``, ``block_size``, ``block_prob``, ``replicas``, or ``temps``
            are out of range.
    """
    n = len(residues)
    if len(seed_dna) != 3 * n:
        raise ValueError(
            f"seed_dna length {len(seed_dna)} != 3 * {n} residues"
        )
    if translate(seed_dna) != "".join(residues):
        raise ValueError("seed_dna does not translate to residues (round-trip broken)")
    # Invariant #5 is stated relative to a feasible seed; enforce that precondition
    # unconditionally (one O(L) pass) so "never raises the hard-violation count"
    # holds absolutely, not just when the caller happens to pass a feasible seed.
    if any(v.severity is Severity.HARD for c in constraints for v in c.validate(seed_dna)):
        raise ValueError("seed_dna is not constraint-feasible (has hard violations)")
    if iterations < 0:
        raise ValueError(f"iterations must be >= 0, got {iterations}")
    if temp0 < 0.0:
        raise ValueError(f"temp0 must be >= 0, got {temp0}")
    if not 0.0 < cooling <= 1.0:
        raise ValueError(f"cooling must be in (0, 1], got {cooling}")
    if block_size < 1:
        raise ValueError(f"block_size must be >= 1, got {block_size}")
    if not 0.0 <= block_prob <= 1.0:
        raise ValueError(f"block_prob must be in [0, 1], got {block_prob}")
    if replicas < 1:
        raise ValueError(f"replicas must be >= 1, got {replicas}")
    if temps is not None:
        temps = tuple(temps)
        if len(temps) != replicas:
            raise ValueError(
                f"temps must have length replicas={replicas}, got {len(temps)}"
            )
        if any(t < 0.0 for t in temps):
            raise ValueError(f"temps must all be >= 0, got {temps}")

    context_len = max((c.context_len() for c in constraints), default=0)
    # Codons after ``pos`` whose ok_suffix window still reaches the changed codon.
    reach = (context_len + 2) // 3

    # Only positions with a genuine synonymous choice are worth proposing.
    movable = [i for i, res in enumerate(residues) if len(synonymous_codons(res)) > 1]

    rng = random.Random(seed)

    seed_score = score(seed_dna)
    # Whole-sequence hard-violation count for the non-local constraints. Refinement
    # may only drive this down, never up (invariant #5, global edition).
    seed_global = _global_hard(seed_dna, global_constraints)

    # Best visited across *all* replicas, ranked lower-global-count-first then
    # higher-score (see :func:`_better`) -- the mechanical guarantee that the
    # delivered result satisfies invariant #5 even with tempering swaps in play.
    best_dna = seed_dna
    best_score = seed_score
    best_global = seed_global

    if movable and iterations > 0:
        # Each replica carries its own (dna, score, global-count, temperature).
        # temps=None => every replica starts at temp0, so replicas=1 reproduces the
        # single-chain schedule exactly (invariant #7).
        start_temps = temps if temps is not None else (temp0,) * replicas
        rep_dna = [seed_dna] * replicas
        rep_score = [seed_score] * replicas
        rep_global = [seed_global] * replicas
        rep_temp = [float(t) for t in start_temps]

        # Block moves need >= 2 movable positions to be a "block"; with block_size
        # left at its default of 1 (or too few movable positions) they are off, and
        # -- because ``block_prob > 0.0`` short-circuits first -- no RNG draw is
        # spent deciding, so the single-codon trajectory is byte-identical.
        allow_block = block_size >= 2 and len(movable) >= 2

        for it in range(iterations):
            for r in range(replicas):
                is_block = allow_block and block_prob > 0.0 and rng.random() < block_prob
                if is_block:
                    cand_dna, cand_score, cand_global, accepted = _block_step(
                        rep_dna[r], rep_score[r], rep_global[r], rep_temp[r],
                        residues, movable, block_size, reach, context_len,
                        constraints, global_constraints, score, rng, n,
                    )
                else:
                    cand_dna, cand_score, cand_global, accepted = _single_step(
                        rep_dna[r], rep_score[r], rep_global[r], rep_temp[r],
                        residues, movable, reach, context_len,
                        constraints, global_constraints, score, delta_score, rng, n,
                    )
                if accepted:
                    rep_dna[r] = cand_dna
                    rep_score[r] = cand_score
                    rep_global[r] = cand_global
                    if _better(cand_global, cand_score, best_global, best_score):
                        best_global = cand_global
                        best_score = cand_score
                        best_dna = cand_dna
                rep_temp[r] *= cooling

            # Parallel-tempering exchange: swap adjacent replicas' (already-feasible)
            # configurations by the standard replica-exchange Metropolis rule. Only
            # relabels states, so no swap can introduce a hard violation; every
            # config keeps global-count <= seed's, preserving invariant #5.
            if swap_every > 0 and replicas > 1 and (it + 1) % swap_every == 0:
                for r in range(replicas - 1):
                    ti, tj = rep_temp[r], rep_temp[r + 1]
                    if ti <= 0.0 or tj <= 0.0:
                        continue  # replica-exchange criterion needs positive temps
                    delta = (1.0 / ti - 1.0 / tj) * (rep_score[r + 1] - rep_score[r])
                    if delta >= 0.0 or rng.random() < math.exp(delta):
                        rep_dna[r], rep_dna[r + 1] = rep_dna[r + 1], rep_dna[r]
                        rep_score[r], rep_score[r + 1] = rep_score[r + 1], rep_score[r]
                        rep_global[r], rep_global[r + 1] = rep_global[r + 1], rep_global[r]

    certificate = OptimalityCertificate(
        status=OptimalityStatus.HEURISTIC,
        solver="anneal_refine",
        detail=(
            f"simulated-annealing refinement over {iterations} iteration(s) from a "
            "feasible seed; heuristic - not proven optimal"
        ),
    )
    # Report the score recomputed on the delivered DNA, never the accumulator
    # (invariant #2, reported == computed). Guard against delta_score drift: if the
    # best visited (chosen by the accumulator) recomputes *worse than the seed by our
    # ranking*, we deliver the seed. Because every visited config has global-count
    # <= the seed's, the fallback only ever fires on an equal-count, lower-score
    # drift -- it can never *raise* the delivered hard-violation count (invariant #5).
    delivered_dna = best_dna
    delivered_score = score(best_dna)
    delivered_global = _global_hard(best_dna, global_constraints)
    if _better(seed_global, seed_score, delivered_global, delivered_score):
        delivered_dna = seed_dna
        delivered_score = seed_score
    return SolveResult(
        dna=delivered_dna, objective_scalar=delivered_score, certificate=certificate
    )


def _global_hard(dna: str, global_constraints: Sequence[Constraint]) -> int:
    """Return the total hard-violation count of ``dna`` under ``global_constraints``.

    A whole-sequence ``validate`` per constraint -- the only honest way to count a
    non-local (``Scope.GLOBAL``) constraint such as the dispersed max-repeat rule,
    whose two copies can lie any distance apart. ``0`` when there are none.

    Args:
        dna: The candidate coding sequence.
        global_constraints: The non-local constraints to audit.
    """
    return sum(
        1
        for c in global_constraints
        for v in c.validate(dna)
        if v.severity is Severity.HARD
    )


def _better(glob_a: int, score_a: float, glob_b: int, score_b: float) -> bool:
    """Return whether state ``a`` outranks state ``b`` for delivery.

    The ranking is **lower global hard-violation count first, then higher score** --
    the ordering that makes invariant #5 hold on the delivered result: a config that
    clears a dispersed max-repeat / uORF always beats one that does not, regardless
    of objective score, and ties on the count break toward the better objective.
    """
    if glob_a != glob_b:
        return glob_a < glob_b
    return score_a > score_b


def _single_step(
    current_dna: str,
    current_score: float,
    current_global: int,
    temp: float,
    residues: Sequence[str],
    movable: Sequence[int],
    reach: int,
    context_len: int,
    constraints: Sequence[Constraint],
    global_constraints: Sequence[Constraint],
    score: Score,
    delta_score: DeltaScore | None,
    rng: random.Random,
    n: int,
) -> tuple[str, float, int, bool]:
    """Propose and evaluate one single-codon synonymous swap for one replica.

    Draws (position, replacement) from ``rng`` in the same order the original
    single-chain loop did -- so with ``replicas=1`` and block moves off, the whole
    trajectory is byte-identical (invariant #7). Returns
    ``(dna, score, global_count, accepted)``: the candidate and its bookkeeping on
    acceptance, or the unchanged state with ``accepted=False`` on veto/rejection.
    Feasibility is checked incrementally (local) then by whole-sequence recount
    (global); a move that would raise the global hard count is rejected (#5).
    """
    pos = rng.choice(movable)
    base = pos * 3
    old_codon = current_dna[base : base + 3]
    options = [c for c in synonymous_codons(residues[pos]) if c != old_codon]
    new_codon = rng.choice(options)

    if not _move_feasible(current_dna, pos, new_codon, reach, context_len, constraints, n):
        return current_dna, current_score, current_global, False

    candidate: str | None
    if global_constraints:
        built = current_dna[:base] + new_codon + current_dna[base + 3 :]
        cand_global = _global_hard(built, global_constraints)
        if cand_global > current_global:
            return current_dna, current_score, current_global, False
        candidate = built
        change = score(built) - current_score
    elif delta_score is not None:
        cand_global = current_global
        change = delta_score(current_dna, pos, old_codon, new_codon)
        candidate = None  # built lazily only on acceptance
    else:
        cand_global = current_global
        candidate = current_dna[:base] + new_codon + current_dna[base + 3 :]
        change = score(candidate) - current_score

    if _accept(change, temp, rng):
        if candidate is None:
            candidate = current_dna[:base] + new_codon + current_dna[base + 3 :]
        return candidate, current_score + change, cand_global, True
    return current_dna, current_score, current_global, False


def _block_step(
    current_dna: str,
    current_score: float,
    current_global: int,
    temp: float,
    residues: Sequence[str],
    movable: Sequence[int],
    block_size: int,
    reach: int,
    context_len: int,
    constraints: Sequence[Constraint],
    global_constraints: Sequence[Constraint],
    score: Score,
    rng: random.Random,
    n: int,
) -> tuple[str, float, int, bool]:
    """Propose and evaluate one **block move** (2..block_size coordinated swaps).

    A block move swaps several movable positions at once, so a barrier that only
    clears when two codons move *together* (a dispersed repeat, a uORF) becomes
    reachable in a single step. Unlike a single move it **always full-``score``
    re-scores** the candidate -- summing per-position deltas is valid only for
    additive disjoint-context terms, which the non-local terms this engine serves
    are not. Local feasibility is checked over the *union* of affected windows and
    the global gate rejects any count increase, so invariant #5 holds move-by-move.
    """
    k = rng.randint(2, min(block_size, len(movable)))
    positions = sorted(rng.sample(list(movable), k))
    candidate = current_dna
    for pos in positions:
        base = pos * 3
        old_codon = candidate[base : base + 3]
        options = [c for c in synonymous_codons(residues[pos]) if c != old_codon]
        new_codon = rng.choice(options)
        candidate = candidate[:base] + new_codon + candidate[base + 3 :]

    if not _block_feasible(candidate, positions, reach, context_len, constraints, n):
        return current_dna, current_score, current_global, False

    if global_constraints:
        cand_global = _global_hard(candidate, global_constraints)
        if cand_global > current_global:
            return current_dna, current_score, current_global, False
    else:
        cand_global = current_global

    change = score(candidate) - current_score
    if _accept(change, temp, rng):
        return candidate, current_score + change, cand_global, True
    return current_dna, current_score, current_global, False


def _block_feasible(
    candidate: str,
    positions: Sequence[int],
    reach: int,
    context_len: int,
    constraints: Sequence[Constraint],
    n_codons: int,
) -> bool:
    """Return whether a multi-position block ``candidate`` adds no hard violation.

    The block candidate already has *all* its swaps applied, so we re-run
    ``ok_suffix`` over the **union** of the affected codon windows (each changed
    position plus the ``reach`` codons whose trailing context reaches it) against
    the candidate's own bases -- exact because ``ok_suffix`` reads only the last
    ``context_len`` characters. Any codon outside every affected window sees the
    same context it saw in the (feasible) pre-move sequence, so it cannot newly
    fail and need not be checked.
    """
    if not constraints:
        return True
    affected: set[int] = set()
    for pos in positions:
        for q in range(pos, min(n_codons, pos + reach + 1)):
            affected.add(q)
    for q in sorted(affected):
        q_off = q * 3
        tail = candidate[max(0, q_off - context_len) : q_off]
        codon_q = candidate[q_off : q_off + 3]
        for constraint in constraints:
            if not constraint.ok_suffix(tail, codon_q):
                return False
    return True


def _accept(change: float, temp: float, rng: random.Random) -> bool:
    """Return whether a move with score delta ``change`` is accepted (Metropolis).

    An improving-or-flat move (``change >= 0``) is always accepted. A worsening
    move is accepted with probability ``exp(change / temp)``; at ``temp <= 0`` it
    is always rejected (greedy). ``change`` is negative here, so the exponent is
    negative and ``exp`` stays in ``(0, 1)`` (underflowing to ``0.0`` harmlessly
    at tiny temperatures).

    Args:
        change: The score delta of the proposed move (larger is better).
        temp: The current annealing temperature.
        rng: The seeded random source (its ``random()`` is the only draw here).
    """
    if change >= 0.0:
        return True
    if temp <= 0.0:
        return False
    return rng.random() < math.exp(change / temp)


def _move_feasible(
    dna: str,
    pos: int,
    new_codon: str,
    reach: int,
    context_len: int,
    constraints: Sequence[Constraint],
    n_codons: int,
) -> bool:
    """Return whether swapping in ``new_codon`` at ``pos`` adds no hard violation.

    The check is incremental (CLAUDE.md 7). Replacing the codon at DNA offset
    ``pos * 3`` can only change a constraint's ``ok_suffix`` verdict for that codon
    and the following ``reach`` codons - those whose declared trailing context
    still overlaps the mutated bases. We build the bounded local region spanning
    the largest constraint context before ``pos`` through the last reachable codon,
    substitute ``new_codon`` into it, and re-run ``ok_suffix`` for every constraint
    at each codon in ``[pos, pos + reach]``. Passing only the ``context_len``-char
    trailing tail as the prefix is exact: the contract guarantees ``ok_suffix``
    reads only the last ``context_len`` characters. Any codon whose context does
    not reach the change sees unchanged, already-feasible bases, so re-checking it
    (which some constraints with shorter context trigger) can only return ``True``.

    Args:
        dna: The current (pre-swap) coding sequence.
        pos: Codon index being changed.
        new_codon: The synonymous codon proposed at ``pos``.
        reach: Number of trailing codons after ``pos`` whose ok_suffix window can
            reach the change (``ceil(context_len / 3)``).
        context_len: Largest constraint ``context_len`` (the trailing-tail width).
        constraints: Hard-feasibility rules to check.
        n_codons: Total number of codons (so we do not run past the end).

    Returns:
        ``True`` if no constraint vetoes any affected codon, else ``False``.
    """
    if not constraints:
        return True

    region_start = max(0, pos * 3 - context_len)
    region_end = min(len(dna), (pos + reach + 1) * 3)
    off = pos * 3 - region_start
    region = dna[region_start:region_end]
    region = region[:off] + new_codon + region[off + 3 :]

    last_q = min(n_codons, pos + reach + 1)
    for q in range(pos, last_q):
        q_off = q * 3 - region_start
        tail = region[max(0, q_off - context_len) : q_off]
        codon_q = region[q_off : q_off + 3]
        for constraint in constraints:
            if not constraint.ok_suffix(tail, codon_q):
                return False
    return True
