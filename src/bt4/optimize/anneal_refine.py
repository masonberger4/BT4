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
            iteration; must lie in ``(0, 1]``.

    Returns:
        A :class:`~bt4.optimize.exact_dp.SolveResult` whose ``dna`` is the best
        (highest-``score``) feasible sequence visited, whose ``objective_scalar``
        is ``score`` recomputed on that DNA, and whose certificate has status
        :attr:`~bt4.domain.certificate.OptimalityStatus.HEURISTIC` and solver
        ``"anneal_refine"``.

    Raises:
        ValueError: If ``seed_dna`` does not translate to ``residues``, if its
            length is not ``3 * len(residues)``, or if ``iterations``, ``temp0``,
            or ``cooling`` are out of range.
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

    context_len = max((c.context_len() for c in constraints), default=0)
    # Codons after ``pos`` whose ok_suffix window still reaches the changed codon.
    reach = (context_len + 2) // 3

    # Only positions with a genuine synonymous choice are worth proposing.
    movable = [i for i, res in enumerate(residues) if len(synonymous_codons(res)) > 1]

    rng = random.Random(seed)

    current_dna = seed_dna
    current_score = score(seed_dna)
    seed_score = current_score
    best_dna = seed_dna
    best_score = current_score
    # Whole-sequence hard-violation count for the non-local constraints. Refinement
    # may only drive this down, never up (invariant #5, global edition).
    current_global = _global_hard(seed_dna, global_constraints)

    if movable and iterations > 0:
        temp = temp0
        for _ in range(iterations):
            pos = rng.choice(movable)
            base = pos * 3
            old_codon = current_dna[base : base + 3]
            options = [c for c in synonymous_codons(residues[pos]) if c != old_codon]
            new_codon = rng.choice(options)

            if not _move_feasible(current_dna, pos, new_codon, reach, context_len, constraints, n):
                temp *= cooling
                continue

            # A global (whole-sequence) constraint's veto cannot be checked by a
            # bounded window, so we build the candidate and re-count its hard
            # violations; a move that would raise the count is rejected outright.
            if global_constraints:
                candidate = current_dna[:base] + new_codon + current_dna[base + 3 :]
                cand_global = _global_hard(candidate, global_constraints)
                if cand_global > current_global:
                    temp *= cooling
                    continue
                change = score(candidate) - current_score
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
                current_dna = candidate
                current_score += change
                current_global = cand_global
                if current_score > best_score:
                    best_score = current_score
                    best_dna = current_dna

            temp *= cooling

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
    # best visited (chosen by the accumulator) recomputes worse than the seed, we
    # deliver the seed, so refinement is never a regression even if a caller's
    # delta_score violated its contract.
    delivered_dna = best_dna
    delivered_score = score(best_dna)
    if delivered_score < seed_score:
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
