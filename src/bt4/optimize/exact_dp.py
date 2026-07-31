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

from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons

__all__ = ["InfeasibleError", "SolveResult", "solve_exact"]


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

    Returns:
        A :class:`SolveResult` whose ``dna`` translates back to
        ``residues`` and whose certificate states the solve's optimality.

    Raises:
        InfeasibleError: If some residue admits no codon that keeps every
            constraint's ``ok_suffix`` satisfied.
    """
    context_len = max([objective_context, *(c.context_len() for c in constraints)])

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


def _wins(new_score: float, new_dna: str, current: tuple[float, str]) -> bool:
    """Return True iff ``(new_score, new_dna)`` should replace ``current``.

    A candidate wins on a strictly higher scalar, or on an equal scalar with a
    lexicographically smaller DNA (deterministic tie-break, invariant #7).
    """
    cur_score, cur_dna = current
    if new_score != cur_score:
        return new_score > cur_score
    return new_dna < cur_dna
