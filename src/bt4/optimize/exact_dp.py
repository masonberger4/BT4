"""Exact dynamic-programming solver over the codon trellis (CLAUDE.md §4.4).

The single load-bearing idea BT4 inherits from BT3 is that back-translation is
*constrained combinatorial optimization*, not greedy per-codon substitution.
This module realises the honest core of that idea: an exact trellis DP that,
when the state space is small enough, provably returns the maximum-scalar coding
sequence subject to every constraint's declared context.

The trellis keeps one layer per residue. A layer maps a *state key* - the last
``K`` characters of the chosen prefix, where ``K`` is the largest ``context_len``
any constraint declares - to the best (highest-scalar) prefix ending in that
context. Because every constraint inspects at most ``K`` trailing characters,
two prefixes sharing a state key are interchangeable for all future feasibility
decisions, so keeping only the best of each is exact. When ``beam`` is set the
solver may drop low-scoring states; it then reports ``BEAM_TRUNCATED`` rather
than pretending the answer is optimal (invariant #6, certificate honesty).

Ties are always broken toward the lexicographically smaller DNA so identical
inputs yield byte-identical output (invariant #7, determinism).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass

from bt4 import _accel
from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons

__all__ = ["InfeasibleError", "SolveResult", "frontier_solver", "solve_exact"]

# A reusable solver: given a scalar objective delta and a beam width, return the
# solve for that weighting. Used by the Pareto frontier, which solves the same
# trellis (same residues + constraints) once per objective scalarization.
ReusableSolver = Callable[[Callable[[str, str, int], float], "int | None"], "SolveResult"]

# Minimum residue count before the native (precomputed-table) trellis path is
# worth its precompute cost. Below this the pure-Python layer loop is already
# fast and the Python-side precompute would dominate (CLAUDE.md §7); the native
# path is a pure performance option, so the choice never changes the output.
_NATIVE_MIN_RESIDUES = 24

# Cap on the reachable-context count before the native precompute is abandoned in
# favor of the pure-Python DP. A large context (e.g. a wide GC-run window) can
# make the transition graph huge; past this bound the precompute would cost more
# than it saves, so callers fall back rather than pessimize (CLAUDE.md §7).
_NATIVE_MAX_CONTEXTS = 100_000


class InfeasibleError(ValueError):
    """Raised when no coding sequence satisfies every hard constraint.

    Attributes:
        constraints: Names of the constraints in force when feasibility was
            lost. The offending residue cannot be extended under their combined
            ``ok_suffix`` vetoes.
    """

    def __init__(self, constraints: Iterable[str]) -> None:
        """Record the constraint names and build a human-readable message.

        Args:
            constraints: The names of the constraints that were active.
        """
        self.constraints: tuple[str, ...] = tuple(constraints)
        names = ", ".join(self.constraints) if self.constraints else "(none)"
        super().__init__(f"no feasible codon under constraints: {names}")


@dataclass(frozen=True, slots=True)
class SolveResult:
    """The outcome of an exact (or beam-truncated) trellis solve.

    Attributes:
        dna: The chosen coding sequence, one codon per input residue.
        objective_scalar: The accumulated scalar objective of ``dna`` (the sum
            of ``scalar_delta`` over its codons; larger is better).
        certificate: How optimal the solve is - ``PROVEN_OPTIMAL`` for a full
            exact DP, ``BEAM_TRUNCATED`` when the beam dropped states.
    """

    dna: str
    objective_scalar: float
    certificate: OptimalityCertificate


# One trellis layer: state key (trailing K chars) -> (best scalar, best DNA).
_Layer = dict[str, tuple[float, str]]


def solve_exact(
    residues: Sequence[str],
    *,
    scalar_delta: Callable[[str, str, int], float],
    constraints: Sequence[Constraint],
    beam: int | None = None,
    objective_context: int = 0,
    position_independent: bool = False,
) -> SolveResult:
    """Solve for the maximum-scalar coding sequence over the codon trellis.

    Args:
        residues: Amino-acid letters to back-translate, **including** a trailing
            ``"*"`` as the final residue (the stop). Residue ``pos`` may use any
            codon in ``synonymous_codons(residues[pos])``.
        scalar_delta: Incremental objective of placing a codon: called as
            ``scalar_delta(prefix, codon, pos)`` and oriented so larger is
            better. Callers typically pass a weighted sum of objective-term
            deltas.
        constraints: Hard-feasibility rules. Each contributes a trailing-context
            veto via ``ok_suffix``; the DP state carries the union of their
            declared contexts (no global cap - CLAUDE.md §10.1).
        beam: If ``None``, run the full exact DP. Otherwise keep at most ``beam``
            states per layer (highest scalar first); dropping any state marks the
            result ``BEAM_TRUNCATED``.
        objective_context: Trailing DNA context (in characters) that
            ``scalar_delta`` itself depends on -- e.g. ``3`` for a pairwise
            objective that reads the previous codon. The trellis state carries
            ``max`` of this and the constraint contexts, so pairwise/positional
            objectives stay exact (states are only merged when both their future
            constraint *and* objective behavior are identical).
        position_independent: Caller's guarantee that ``scalar_delta`` is a pure
            function of ``(prefix[-K:], codon)`` and does **not** depend on the
            codon index ``pos`` (i.e. no ``POSITIONAL``-scope objective term is
            active). When ``True``, the Rust ``bt4_native.trellis_solve`` inner
            loop is used for large instances -- transitions are precomputed once
            per ``(context, codon)`` in Python (so the float summation order, and
            therefore the tie-break, stay bit-for-bit identical) and the layer DP
            runs in Rust. This is a **pure performance** switch: the output is
            byte-identical to the default pure-Python DP, and the pure-Python DP
            is used unchanged when the guarantee does not hold, the native
            extension is absent, or the instance is small (:data:`_NATIVE_MIN_RESIDUES`).

    Returns:
        A :class:`SolveResult` whose ``dna`` translates back to
        ``residues`` and whose certificate states the solve's optimality.

    Raises:
        InfeasibleError: If some residue admits no codon that keeps every
            constraint's ``ok_suffix`` satisfied.
    """
    context_len = max([objective_context, *(c.context_len() for c in constraints)])

    if (
        position_independent
        and _accel.ACCELERATED
        and len(residues) >= _NATIVE_MIN_RESIDUES
    ):
        # Native fast path: precompute the transition tables in Python (keeping
        # the float/tie-break behavior identical) and run the layer DP in Rust.
        return _solve_native(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,
            beam=beam,
            context_len=context_len,
        )

    layer: _Layer = {"": (0.0, "")}
    pruned = False

    for pos, residue in enumerate(residues):
        codons = synonymous_codons(residue)
        next_layer: _Layer = {}
        # Sorted iteration keeps the build order deterministic; the merge rule
        # below is order-independent regardless, but sorting makes it obvious.
        for _key, (score, dna) in sorted(layer.items()):
            for codon in codons:
                if not all(c.ok_suffix(dna, codon) for c in constraints):
                    continue
                new_dna = dna + codon
                new_score = score + scalar_delta(dna, codon, pos)
                new_key = new_dna[-context_len:] if context_len > 0 else ""
                current = next_layer.get(new_key)
                # Keep the higher scalar; break exact ties toward the
                # lexicographically smaller DNA so output is reproducible.
                if current is None or _wins(new_score, new_dna, current):
                    next_layer[new_key] = (new_score, new_dna)

        if not next_layer:
            raise InfeasibleError([c.name for c in constraints])

        if beam is not None and len(next_layer) > beam:
            kept = sorted(
                next_layer.items(), key=lambda kv: (-kv[1][0], kv[1][1])
            )[:beam]
            next_layer = dict(kept)
            pruned = True

        layer = next_layer

    best_score, best_dna = min(
        layer.values(), key=lambda sv: (-sv[0], sv[1])
    )

    if beam is not None and pruned:
        certificate = OptimalityCertificate(
            status=OptimalityStatus.BEAM_TRUNCATED,
            solver="beam_dp",
            detail=f"beam width {beam}",
        )
    else:
        certificate = OptimalityCertificate.proven(
            "exact_dp",
            detail=f"exact DP over {len(residues)} residues, context K={context_len}",
        )

    return SolveResult(dna=best_dna, objective_scalar=best_score, certificate=certificate)


# Sentinel for "key absent" in a cache that legitimately stores ``None`` values.
_MISSING = object()


@dataclass(frozen=True, slots=True)
class _TrellisStructure:
    """The position-independent transition graph of a codon trellis.

    The graph depends only on the residues, the constraints (via ``ok_suffix``),
    and the context length ``K`` -- **not** on the objective weights -- so it can
    be built once and reused across many objective scalarizations (the frontier
    sweep, CLAUDE.md §7 requirement to cache tables across grid points). Only the
    per-transition scalar deltas change with the weights, and those are cheap to
    recompute (see :func:`_delta_tables`).

    Context id ``0`` is the empty start context. ``layer_from/to/codon[l]`` are the
    parallel lists of the *allowed* transitions of layer ``l``: a transition places
    ``codons[layer_codon[l][t]]`` from context ``layer_from[l][t]``, arriving at the
    merged trailing-context ``layer_to[l][t]``. ``id_to_str`` maps a context id back
    to its string (needed only to recompute deltas). ``feasible`` is ``False`` when
    a layer became unreachable (the problem is infeasible).
    """

    context_len: int
    codons: tuple[str, ...]
    id_to_str: tuple[str, ...]
    layer_from: tuple[tuple[int, ...], ...]
    layer_to: tuple[tuple[int, ...], ...]
    layer_codon: tuple[tuple[int, ...], ...]
    feasible: bool


def _precompute_structure(
    residues: Sequence[str],
    constraints: Sequence[Constraint],
    context_len: int,
    *,
    max_contexts: int | None = None,
) -> _TrellisStructure | None:
    """Enumerate the reachable trailing-context transition graph of the trellis.

    Walks the reachable context states layer by layer -- exactly the states the
    pure-Python DP would reach -- assigning each distinct context string an integer
    id (``""`` is id ``0``). Feasibility of a ``(context, codon)`` transition is
    decided by ``ok_suffix``, which reads at most each constraint's ``context_len``
    (<= ``context_len``) trailing characters; the state key holds those characters,
    so evaluating on the truncated context equals evaluating on the full prefix
    (CLAUDE.md §5 #3). The objective is *not* consulted here (that is why the
    structure is reusable across weightings).

    Returns the :class:`_TrellisStructure`, or ``None`` when the reachable-context
    count would exceed ``max_contexts`` (the caller then falls back to the
    pure-Python DP rather than pay an unbounded precompute -- CLAUDE.md §7).
    """
    context_ids: dict[str, int] = {"": 0}
    id_to_str: list[str] = [""]
    codon_ids: dict[str, int] = {}
    codons: list[str] = []
    new_ctx_cache: dict[tuple[str, str], str | None] = {}

    def transition_ctx(ctx: str, codon: str) -> str | None:
        key = (ctx, codon)
        cached = new_ctx_cache.get(key, _MISSING)
        if cached is not _MISSING:
            return cached  # type: ignore[return-value]
        if all(c.ok_suffix(ctx, codon) for c in constraints):
            result: str | None = (ctx + codon)[-context_len:] if context_len > 0 else ""
        else:
            result = None
        new_ctx_cache[key] = result
        return result

    def codon_id(codon: str) -> int:
        cid = codon_ids.get(codon)
        if cid is None:
            cid = len(codons)
            codon_ids[codon] = cid
            codons.append(codon)
        return cid

    layer_from: list[tuple[int, ...]] = []
    layer_to: list[tuple[int, ...]] = []
    layer_codon: list[tuple[int, ...]] = []

    reachable: set[int] = {0}
    feasible = True
    for residue in residues:
        syn = synonymous_codons(residue)
        lf: list[int] = []
        lt: list[int] = []
        lc: list[int] = []
        nxt: set[int] = set()
        # Sorted iteration keeps the emitted order deterministic; the DP merge is
        # order-independent regardless.
        for c_id in sorted(reachable):
            ctx = id_to_str[c_id]
            for codon in syn:
                new_ctx = transition_ctx(ctx, codon)
                if new_ctx is None:
                    continue
                nid = context_ids.get(new_ctx)
                if nid is None:
                    nid = len(id_to_str)
                    context_ids[new_ctx] = nid
                    id_to_str.append(new_ctx)
                    if max_contexts is not None and nid >= max_contexts:
                        return None
                lf.append(c_id)
                lt.append(nid)
                lc.append(codon_id(codon))
                nxt.add(nid)
        layer_from.append(tuple(lf))
        layer_to.append(tuple(lt))
        layer_codon.append(tuple(lc))
        reachable = nxt
        if not reachable:
            # No feasible extension: this and every later layer are empty.
            feasible = False
            break

    return _TrellisStructure(
        context_len=context_len,
        codons=tuple(codons),
        id_to_str=tuple(id_to_str),
        layer_from=tuple(layer_from),
        layer_to=tuple(layer_to),
        layer_codon=tuple(layer_codon),
        feasible=feasible,
    )


def _delta_tables(
    structure: _TrellisStructure,
    scalar_delta: Callable[[str, str, int], float],
) -> list[list[float]]:
    """Per-layer scalar deltas for ``structure`` under one objective weighting.

    ``delta[l][t] = scalar_delta(context_string, codon, 0)`` for the ``t``-th
    transition of layer ``l``. Computed in Python so the float summation order
    (and hence the tie-break) is byte-for-byte identical to the pure-Python DP;
    memoized per distinct ``(from_context, codon)`` so a residue that recurs across
    layers evaluates the objective once. Under the caller's position-independence
    guarantee this equals the pure DP's ``scalar_delta(full_prefix, codon, pos)``
    (the objective reads at most ``context_len`` trailing characters and never
    ``pos``).
    """
    cache: dict[tuple[int, int], float] = {}
    id_to_str = structure.id_to_str
    codons = structure.codons
    out: list[list[float]] = []
    for froms, cods in zip(structure.layer_from, structure.layer_codon, strict=True):
        row: list[float] = []
        for from_id, codon_gid in zip(froms, cods, strict=True):
            key = (from_id, codon_gid)
            d = cache.get(key)
            if d is None:
                d = scalar_delta(id_to_str[from_id], codons[codon_gid], 0)
                cache[key] = d
            row.append(d)
        out.append(row)
    return out


def _solve_from_structure(
    structure: _TrellisStructure,
    scalar_delta: Callable[[str, str, int], float],
    constraints: Sequence[Constraint],
    beam: int | None,
    n_residues: int,
) -> SolveResult:
    """Run the native trellis DP over a prebuilt structure and one weighting.

    Byte-for-byte equivalent to the pure-Python loop in :func:`solve_exact` under
    the position-independence guarantee; only the inner loop runs in native code.
    """
    if not structure.feasible:
        raise InfeasibleError([c.name for c in constraints])
    layer_delta = _delta_tables(structure, scalar_delta)
    out = _accel.trellis_solve(
        list(structure.codons),
        [list(x) for x in structure.layer_from],
        [list(x) for x in structure.layer_to],
        [list(x) for x in structure.layer_codon],
        layer_delta,
        beam,
    )
    if out is None:
        raise InfeasibleError([c.name for c in constraints])
    best_dna, best_score, pruned = out
    if beam is not None and pruned:
        certificate = OptimalityCertificate(
            status=OptimalityStatus.BEAM_TRUNCATED,
            solver="beam_dp",
            detail=f"beam width {beam}",
        )
    else:
        certificate = OptimalityCertificate.proven(
            "exact_dp",
            detail=f"exact DP over {n_residues} residues, "
            f"context K={structure.context_len}",
        )
    return SolveResult(dna=best_dna, objective_scalar=best_score, certificate=certificate)


def _solve_native(
    residues: Sequence[str],
    *,
    scalar_delta: Callable[[str, str, int], float],
    constraints: Sequence[Constraint],
    beam: int | None,
    context_len: int,
) -> SolveResult:
    """Solve one instance via the precomputed tables + native (or twin) DP.

    Builds the structure and deltas and runs :func:`_accel.trellis_solve`. Byte-for
    -byte equivalent to :func:`solve_exact`'s pure-Python loop under the position
    -independence guarantee. Note that for a *single* solve the Python precompute
    dominates, so this is used by the frontier (where the structure is reused, see
    :func:`_solve_from_structure`) and by tests, not by single-shot optimization.
    """
    structure = _precompute_structure(residues, constraints, context_len)
    assert structure is not None  # no max_contexts cap here
    return _solve_from_structure(structure, scalar_delta, constraints, beam, len(residues))


def frontier_solver(
    residues: Sequence[str],
    constraints: Sequence[Constraint],
    *,
    objective_context: int = 0,
    position_independent: bool = False,
) -> ReusableSolver:
    """Build a reusable solver that shares one trellis structure across weightings.

    The Pareto frontier solves the *same* trellis -- identical residues and
    constraints -- once per objective scalarization, changing only the objective
    weights. This is the public seam that lets that sweep amortize the work: when
    the objective is ``position_independent`` (no ``POSITIONAL`` term) and the Rust
    accelerator is present on a large-enough instance, the reachable-context
    transition graph is precomputed **once** here and reused for every returned
    solve (only the cheap per-transition deltas are recomputed, and the DP inner
    loop runs in native code); otherwise every call falls back to the pure-Python
    :func:`solve_exact`. Either way each returned solve is byte-identical to the
    pure-Python exact DP (a proven-optimal, or explicitly beam-truncated, result).

    The native-path gating (:data:`_NATIVE_MIN_RESIDUES`, :data:`_NATIVE_MAX_CONTEXTS`)
    and the trellis-structure internals stay encapsulated here in the ``optimize``
    layer -- callers get a plain public callable and never touch a private symbol
    (CLAUDE.md §3/§10.9).

    Args:
        residues: Amino-acid letters including the trailing stop.
        constraints: Hard-feasibility rules (the shared trellis structure).
        objective_context: Extra trailing context the objective needs (as in
            :func:`solve_exact`).
        position_independent: Caller's guarantee that every weighting passed to the
            returned solver is a pure function of ``(prefix[-K:], codon)`` with no
            dependence on the codon index (see :func:`solve_exact`). Only then is
            the native shared-structure path eligible.

    Returns:
        A :data:`ReusableSolver`: call it as ``solve(scalar_delta, beam)`` to get
        the :class:`SolveResult` for that objective weighting.
    """
    context_len = max([objective_context, *(c.context_len() for c in constraints)])
    n_residues = len(residues)
    structure: _TrellisStructure | None = None
    if position_independent and _accel.ACCELERATED and n_residues >= _NATIVE_MIN_RESIDUES:
        structure = _precompute_structure(
            residues, constraints, context_len, max_contexts=_NATIVE_MAX_CONTEXTS
        )

    def solve(
        scalar_delta: Callable[[str, str, int], float], beam: int | None
    ) -> SolveResult:
        if structure is not None:
            return _solve_from_structure(
                structure, scalar_delta, constraints, beam, n_residues
            )
        return solve_exact(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,
            beam=beam,
            objective_context=objective_context,
        )

    return solve


def _wins(new_score: float, new_dna: str, current: tuple[float, str]) -> bool:
    """Return True iff ``(new_score, new_dna)`` should replace ``current``.

    A candidate wins on a strictly higher scalar, or on an equal scalar with a
    lexicographically smaller DNA (deterministic tie-break, invariant #7).
    """
    cur_score, cur_dna = current
    if new_score != cur_score:
        return new_score > cur_score
    return new_dna < cur_dna
