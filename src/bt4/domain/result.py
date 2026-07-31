"""Immutable result types produced by the BT4 pipeline.

These are pure data carriers (frozen dataclasses, stdlib only). Every metric
reported here is defined to be *recomputed from the sequence* by the owning
model — never trusted from a solver accumulator (the ``reported == computed``
invariant).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum

from .certificate import OptimalityCertificate
from .objective import ObjectiveVector

__all__ = ["Metrics", "Result", "Severity", "Violation"]


class Severity(Enum):
    """Whether a constraint vetoes (HARD) or merely scores (SOFT)."""

    HARD = "hard"
    SOFT = "soft"


@dataclass(frozen=True, slots=True)
class Violation:
    """A single constraint violation found by a whole-sequence audit.

    Attributes:
        constraint: Registered name of the constraint that fired.
        severity: HARD (feasibility) or SOFT (quality).
        start: 0-based inclusive start position on the DNA.
        end: 0-based exclusive end position on the DNA.
        detail: Human-readable description of the violation.
    """

    constraint: str
    severity: Severity
    start: int
    end: int
    detail: str = ""


@dataclass(frozen=True, slots=True)
class Metrics:
    """Recomputed sequence metrics attached to a result.

    Every field is computed from ``Result.dna`` by the owning model, so it can
    be independently re-derived and checked. ``objective`` carries the full
    multi-objective vector; the scalar convenience fields mirror common terms.
    """

    objective: ObjectiveVector
    gc: float
    length_nt: int
    hard_violations: int = 0
    soft_violations: int = 0


@dataclass(frozen=True, slots=True)
class Result:
    """The outcome of one back-translation optimization.

    Attributes:
        protein: The input protein (stop-free, upper-cased).
        dna: The optimized coding DNA (includes the stop codon when appended).
        metrics: Recomputed metrics for ``dna``.
        certificate: How optimal the solve was, and what (if anything) was relaxed.
        violations: Whole-sequence audit findings (empty when fully feasible).
        audit: Free-form provenance/audit payload (manifest, seed, timings).
    """

    protein: str
    dna: str
    metrics: Metrics
    certificate: OptimalityCertificate
    violations: tuple[Violation, ...] = ()
    audit: Mapping[str, object] = field(default_factory=dict)

    @property
    def hard_violations(self) -> tuple[Violation, ...]:
        """The subset of ``violations`` with HARD severity."""
        return tuple(v for v in self.violations if v.severity is Severity.HARD)

    @property
    def is_feasible(self) -> bool:
        """True when the sequence has no hard violations."""
        return not self.hard_violations
