"""CP-SAT (OR-Tools) backend for the additive core with a global GC budget.

The exact DP (:mod:`bt4.optimize.exact_dp`) is the honest workhorse for local,
context-bounded objectives and constraints. It cannot, however, carry a
*whole-sequence* GC-count budget: a global running total is not a
bounded-trailing-context state, so encoding it in the trellis would blow up the
state space (one state per attainable GC count per prefix). This module hands
that exact problem to OR-Tools CP-SAT, which models the global budget as a
single linear constraint and reports honest optimality (proven-optimal, or
gap-bounded when it hits the time limit).

Scope - what this backend does and, deliberately, does not do:

- **Does:** solve an *additive, context-free* per-codon objective (one scalar
  coefficient per ``(codon, pos)`` placement) subject to (a) exactly one codon
  per residue and (b) an optional global GC-count budget over the whole
  sequence. This is exactly the class of problems the DP cannot bound.
- **Does not (yet):** encode any *local sequence* constraint - homopolymer runs,
  forbidden or restriction motifs, tandem/inverted repeats, Kozak/uORF context,
  or any pairwise/positional objective term. Those read trailing sequence
  context and remain the exact DP's job. Do not route a problem carrying such
  constraints through this backend expecting them to be honored; they are simply
  absent from the model here.

Optimality is stated with respect to the *integer-scaled* objective. Float
coefficients are converted to integers via a fixed multiplier ``_SCALE`` (a
million) and rounded, because CP-SAT is an integer program; ``_SCALE`` is large
enough to preserve realistic codon-score resolution. The returned
``objective_scalar`` is always the *true* float objective, recomputed from the
chosen codons - never the scaled integer the solver optimized.

Determinism (invariant #7): the solver runs single-threaded with a fixed seed,
and codons are enumerated in the sorted order :func:`synonymous_codons`
guarantees, so identical inputs yield byte-identical output.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from ortools.sat.python import cp_model

from bt4._accel import gc_count
from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus
from bt4.domain.genetic_code import synonymous_codons
from bt4.optimize.exact_dp import InfeasibleError, SolveResult

__all__ = ["solve_cpsat"]

_SCALE: int = 1_000_000
"""Fixed float-to-integer multiplier for the CP-SAT objective.

CP-SAT optimizes integers, so each float coefficient ``codon_score(codon, pos)``
is stored as ``round(codon_score(codon, pos) * _SCALE)``. A million preserves
about six decimal digits of resolution, ample for codon scores in ``[-1, 1]``.
"""


def solve_cpsat(
    residues: Sequence[str],
    *,
    codon_score: Callable[[str, int], float],
    gc_min: int | None = None,
    gc_max: int | None = None,
    max_time_s: float = 10.0,
) -> SolveResult:
    """Solve the additive core plus an optional global GC budget with CP-SAT.

    Builds a CP-SAT model with one Boolean placement variable ``x[pos, codon]``
    per candidate codon, an ``AddExactlyOne`` clause per residue, and a maximized
    linear objective over the integer-scaled per-codon scores. When ``gc_min``
    and/or ``gc_max`` are given, a single linear constraint bounds the
    whole-sequence GC nucleotide count.

    Args:
        residues: Amino-acid letters to back-translate, **including** a trailing
            ``"*"`` as the final residue (the stop). Residue ``pos`` may use any
            codon in ``synonymous_codons(residues[pos])``.
        codon_score: The additive, context-free objective coefficient for placing
            ``codon`` at index ``pos``; called as ``codon_score(codon, pos)`` and
            oriented so larger is better. Only context-free scores are meaningful
            here - this backend does not read trailing sequence context.
        gc_min: If given, require the whole sequence to contain at least this many
            G/C nucleotides.
        gc_max: If given, require the whole sequence to contain at most this many
            G/C nucleotides.
        max_time_s: Wall-clock ceiling handed to the solver. On expiry the solver
            returns the best sequence found so far with a gap-bounded certificate
            rather than a proof of optimality.

    Returns:
        A :class:`SolveResult` whose ``dna`` translates back to ``residues``
        (exactly one codon per residue, in order), whose ``objective_scalar`` is
        the true float objective recomputed from the chosen codons, and whose
        certificate is ``PROVEN_OPTIMAL`` (solver proved optimality of the
        integer-scaled objective) or ``GAP_BOUNDED`` (time limit hit, with the
        solver's proven relative gap).

    Raises:
        InfeasibleError: If the GC budget (or the model generally) admits no
            assignment, or the solver returns a non-solution status. The named
            constraints identify the binding GC bound(s) when applicable.
    """
    model = cp_model.CpModel()

    # Per residue: the ordered list of (codon, BoolVar) candidates. Ortools
    # objects are untyped, so their containers are annotated with ``Any``.
    grid: list[list[tuple[str, Any]]] = []
    objective_terms: list[Any] = []
    gc_terms: list[Any] = []
    budget: bool = gc_min is not None or gc_max is not None

    for pos, residue in enumerate(residues):
        row: list[tuple[str, Any]] = []
        for codon in synonymous_codons(residue):
            var = model.NewBoolVar(f"x_{pos}_{codon}")
            row.append((codon, var))
            objective_terms.append(round(codon_score(codon, pos) * _SCALE) * var)
            if budget:
                gc_terms.append(gc_count(codon) * var)
        model.AddExactlyOne(var for _codon, var in row)
        grid.append(row)

    model.Maximize(sum(objective_terms))

    budget_names: list[str] = []
    if budget:
        gc_expr = sum(gc_terms)
        if gc_min is not None:
            model.Add(gc_expr >= gc_min)
            budget_names.append("gc_min")
        if gc_max is not None:
            model.Add(gc_expr <= gc_max)
            budget_names.append("gc_max")

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = max_time_s
    solver.parameters.random_seed = 0
    solver.parameters.num_search_workers = 1
    status = solver.Solve(model)

    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        if status == cp_model.INFEASIBLE:
            raise InfeasibleError(budget_names or ["gc_budget"])
        raise InfeasibleError([f"cpsat:{solver.StatusName(status)}"])

    chosen = [
        next(codon for codon, var in row if solver.Value(var)) for row in grid
    ]
    dna = "".join(chosen)
    objective_scalar = sum(
        codon_score(codon, pos) for pos, codon in enumerate(chosen)
    )

    if status == cp_model.OPTIMAL:
        certificate = OptimalityCertificate.proven(
            "cpsat",
            detail=(
                f"CP-SAT proven optimal over {len(residues)} residues "
                f"w.r.t. the integer-scaled objective (SCALE={_SCALE})"
            ),
        )
    else:
        bound = solver.BestObjectiveBound()
        obj = solver.ObjectiveValue()
        gap = (bound - obj) / max(abs(bound), 1.0)
        gap = min(1.0, max(0.0, gap))
        certificate = OptimalityCertificate(
            status=OptimalityStatus.GAP_BOUNDED,
            solver="cpsat",
            gap=gap,
            detail=(
                f"CP-SAT hit the {max_time_s}s time limit without proof; "
                f"relative gap on the integer-scaled objective is {gap:.3g}"
            ),
        )

    return SolveResult(
        dna=dna, objective_scalar=objective_scalar, certificate=certificate
    )
