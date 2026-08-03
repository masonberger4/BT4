"""The optimality certificate — BT4's honesty about how good a solve is.

BT3 printed "CAI-optimal" in a docstring while silently beam-pruning and
context-capping. BT4 attaches an ``OptimalityCertificate`` to every result that
states, in machine-readable form, exactly how optimal the solve is and what (if
anything) was relaxed. A run can never silently lose optimality.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

__all__ = ["OptimalityCertificate", "OptimalityStatus"]


class OptimalityStatus(Enum):
    """How optimal a returned solution is known to be.

    Ordered from strongest to weakest guarantee. ``PROVEN_OPTIMAL`` means an
    exact solver explored the full (constraint-declared) state space;
    ``GAP_BOUNDED`` means a relaxation/ILP returned a solution with a proven
    bound on its distance from the optimum; the remaining statuses are explicit
    admissions that a guarantee was traded for speed or tractability.
    ``SAMPLED`` is the odd one out: it is not a weaker optimum but a
    non-optimizing draw, and it never claims to be optimal at all.
    """

    PROVEN_OPTIMAL = "proven_optimal"
    GAP_BOUNDED = "gap_bounded"
    BEAM_TRUNCATED = "beam_truncated"
    CONTEXT_CAPPED = "context_capped"
    RELAXED = "relaxed"
    HEURISTIC = "heuristic"
    SAMPLED = "sampled"
    """A stochastic draw from the codon distribution (library/degenerate-design
    mode) -- **no optimality is claimed**. The sequence was sampled, not
    optimized: it round-trips and carries recomputed metrics, but it is neither
    a proven optimum nor an expression prediction. See :mod:`bt4.optimize.sample`
    and :mod:`bt4.pipeline.library`."""


@dataclass(frozen=True, slots=True)
class OptimalityCertificate:
    """A machine-readable statement of solve quality.

    Attributes:
        status: The strongest guarantee that holds for the returned solution.
        solver: Name of the solver backend that produced the solution
            (e.g. ``"exact_dp"``, ``"beam_dp"``, ``"cpsat_ilp"``, ``"lagrangian"``).
        gap: Proven relative optimality gap in ``[0, 1]`` when known (0.0 for
            ``PROVEN_OPTIMAL``), else ``None``.
        relaxed_terms: Names of objective/constraint terms that were relaxed
            (e.g. dualized global budgets), empty when nothing was relaxed.
        detail: Free-form human-readable note (e.g. beam width used, which
            constraint's context exceeded a cap).
    """

    status: OptimalityStatus
    solver: str
    gap: float | None = None
    relaxed_terms: tuple[str, ...] = ()
    detail: str = ""

    def __post_init__(self) -> None:
        if self.gap is not None and not (0.0 <= self.gap <= 1.0):
            raise ValueError(f"optimality gap must be in [0, 1], got {self.gap}")
        if self.status is OptimalityStatus.PROVEN_OPTIMAL:
            if self.gap not in (None, 0.0):
                raise ValueError("PROVEN_OPTIMAL cannot carry a non-zero gap")
            if self.relaxed_terms:
                raise ValueError("PROVEN_OPTIMAL cannot have relaxed terms")

    @property
    def is_proven_optimal(self) -> bool:
        """True iff the solution is certified globally optimal."""
        return self.status is OptimalityStatus.PROVEN_OPTIMAL

    @classmethod
    def proven(cls, solver: str, detail: str = "") -> OptimalityCertificate:
        """Construct a certificate for a provably optimal solve."""
        return cls(OptimalityStatus.PROVEN_OPTIMAL, solver, gap=0.0, detail=detail)


# Sentinel used before any solver has run (e.g. a validation-only result).
NOT_OPTIMIZED: OptimalityCertificate = OptimalityCertificate(
    status=OptimalityStatus.HEURISTIC,
    solver="none",
    detail="no optimization performed",
)
