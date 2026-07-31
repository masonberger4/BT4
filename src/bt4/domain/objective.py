"""Multi-objective primitives: the objective vector and the Pareto frontier.

BT4's objective is a *vector*, not a scalar. These pure types let the pipeline
carry, compare, and expose trade-offs honestly — no single magic-weighted
number, and every delivered point knows where it sits on the frontier.

Convention: **every objective is expressed so that larger is better.** Terms
that are naturally costs (splice risk, CpG deviation, folding stability near the
start) are negated by their owning ``ObjectiveTerm`` before landing here.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType

__all__ = ["Frontier", "ObjectiveVector", "dominates", "pareto_front"]


@dataclass(frozen=True, slots=True)
class ObjectiveVector:
    """A named bundle of objective values, all oriented so larger is better.

    Stored as an immutable mapping ``term name -> value``. Two vectors are only
    comparable when they carry the same set of term names.
    """

    values: Mapping[str, float]

    def __post_init__(self) -> None:
        # Freeze into a read-only mapping so the vector is truly immutable.
        object.__setattr__(self, "values", MappingProxyType(dict(self.values)))

    def terms(self) -> frozenset[str]:
        """The set of objective term names present in this vector."""
        return frozenset(self.values)

    def get(self, term: str) -> float:
        """Return the value of ``term`` (raises ``KeyError`` if absent)."""
        return self.values[term]

    def scalarize(self, weights: Mapping[str, float]) -> float:
        """Weighted-sum scalarization over the given weights.

        Only terms present in ``weights`` contribute. This is a *reporting*
        convenience for ranking a frontier — never the primary objective.
        """
        return sum(self.values[t] * w for t, w in weights.items())


def _comparable(a: ObjectiveVector, b: ObjectiveVector) -> None:
    if a.terms() != b.terms():
        raise ValueError(
            f"objective vectors are not comparable: {sorted(a.terms())} "
            f"!= {sorted(b.terms())}"
        )


def dominates(a: ObjectiveVector, b: ObjectiveVector) -> bool:
    """Return True if ``a`` Pareto-dominates ``b`` (larger-is-better on all terms).

    ``a`` dominates ``b`` iff it is no worse on every term and strictly better on
    at least one.

    Raises:
        ValueError: If the two vectors do not share the same term set.
    """
    _comparable(a, b)
    no_worse = all(a.values[t] >= b.values[t] for t in a.values)
    strictly_better = any(a.values[t] > b.values[t] for t in a.values)
    return no_worse and strictly_better


def pareto_front(points: Iterable[ObjectiveVector]) -> tuple[ObjectiveVector, ...]:
    """Return the non-dominated subset of ``points`` (order preserved).

    A point is kept iff no other point dominates it. Duplicate-valued points are
    all kept (none strictly dominates another).
    """
    pts = list(points)
    keep: list[ObjectiveVector] = []
    for i, p in enumerate(pts):
        if not any(dominates(q, p) for j, q in enumerate(pts) if j != i):
            keep.append(p)
    return tuple(keep)


@dataclass(frozen=True, slots=True)
class Frontier:
    """A computed Pareto frontier plus the index of the delivered point.

    Attributes:
        points: The non-dominated objective vectors.
        chosen: Index into ``points`` of the point actually returned to the
            caller (``-1`` when no single point was selected).
    """

    points: tuple[ObjectiveVector, ...]
    chosen: int = -1

    def chosen_point(self) -> ObjectiveVector | None:
        """The delivered objective vector, or ``None`` if none was selected."""
        if 0 <= self.chosen < len(self.points):
            return self.points[self.chosen]
        return None
