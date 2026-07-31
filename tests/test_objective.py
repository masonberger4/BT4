"""Tests for the multi-objective primitives."""

from __future__ import annotations

import pytest

from bt4.domain.objective import (
    Frontier,
    ObjectiveVector,
    dominates,
    pareto_front,
)


def _v(**kw: float) -> ObjectiveVector:
    return ObjectiveVector(kw)


def test_vector_is_immutable_and_comparable() -> None:
    v = _v(cai=0.8, gc=0.5)
    assert v.terms() == frozenset({"cai", "gc"})
    assert v.get("cai") == 0.8
    with pytest.raises(TypeError):
        v.values["cai"] = 0.1  # type: ignore[index]


def test_dominance_larger_is_better() -> None:
    a = _v(cai=0.9, gc=0.6)
    b = _v(cai=0.8, gc=0.6)
    assert dominates(a, b)
    assert not dominates(b, a)
    # Equal vectors do not dominate each other.
    assert not dominates(a, _v(cai=0.9, gc=0.6))


def test_dominance_requires_matching_terms() -> None:
    with pytest.raises(ValueError):
        dominates(_v(cai=0.9), _v(gc=0.6))


def test_pareto_front_keeps_only_nondominated() -> None:
    pts = [
        _v(cai=0.9, gc=0.4),  # non-dominated
        _v(cai=0.5, gc=0.9),  # non-dominated
        _v(cai=0.4, gc=0.3),  # dominated by the first
    ]
    front = pareto_front(pts)
    assert pts[0] in front and pts[1] in front
    assert pts[2] not in front


def test_scalarize_reports_weighted_sum() -> None:
    v = _v(cai=1.0, gc=0.5)
    assert v.scalarize({"cai": 2.0, "gc": 4.0}) == pytest.approx(4.0)


def test_frontier_tracks_chosen_point() -> None:
    pts = (_v(cai=0.9, gc=0.4), _v(cai=0.5, gc=0.9))
    f = Frontier(points=pts, chosen=1)
    assert f.chosen_point() is pts[1]
    assert Frontier(points=pts).chosen_point() is None
