"""Tests for the ExpressionPredictor scaffold and the frontier-rerank hook.

The scaffold ships no calibrated model, so the load-bearing behaviours are
*honesty* behaviours: the default predictor reports ``calibrated is False`` and an
information-free score, and the rerank hook annotates a frontier with expression
scores but only re-picks the delivered point when the predictor is calibrated (an
uncalibrated score must never steer delivery -- CLAUDE.md §10.5).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from bt4 import api
from bt4.biomodels.expression import (
    ExpressionPredictor,
    ExpressionResult,
    NullExpressionModel,
    default,
)
from bt4.pipeline.optimize import OptimizeConfig
from bt4.pipeline.rerank import rerank_by_expression

_PROTEIN = "MKAILVDEQTRSFYNWHGP"


# --------------------------------------------------------------------------- #
# The default backend is the honest neutral placeholder.
# --------------------------------------------------------------------------- #


def test_default_is_uncalibrated_placeholder() -> None:
    model = default()
    assert isinstance(model, NullExpressionModel)
    assert isinstance(model, ExpressionPredictor)  # satisfies the protocol
    assert model.calibrated is False
    assert model.name == "null_expression"


def test_placeholder_scores_zero_and_labels_itself() -> None:
    model = default()
    r = model.score_sequence("ATGGCCTAA")
    assert isinstance(r, ExpressionResult)
    assert r.score == 0.0
    assert r.calibrated is False
    assert r.model_name == "null_expression"
    assert "placeholder" in r.units


# --------------------------------------------------------------------------- #
# The rerank hook: annotate always, steer delivery only when calibrated.
# --------------------------------------------------------------------------- #


def test_rerank_with_placeholder_annotates_but_never_steers() -> None:
    frontier = api.frontier(_PROTEIN, OptimizeConfig(gc_weight=1.0), steps=7)
    before = frontier.frontier.chosen
    reranked = rerank_by_expression(frontier)  # default placeholder
    # Delivery is unchanged: an uncalibrated score must not steer it.
    assert reranked.frontier.chosen == before
    # Every result is annotated for reporting, honestly flagged uncalibrated.
    assert reranked.results
    for r in reranked.results:
        assert r.audit["expression_calibrated"] is False
        assert r.audit["expression_score"] == 0.0
        assert r.audit["expression_model"] == "null_expression"


@dataclass(frozen=True, slots=True)
class _FakeCalibrated:
    """A stand-in calibrated head that prefers higher-GC sequences (for the test)."""

    name: str = field(default="fake_calibrated", init=False)

    @property
    def calibrated(self) -> bool:
        return True

    def score_sequence(self, dna: str) -> ExpressionResult:
        gc = sum(1 for b in dna.upper() if b in "GC")
        return ExpressionResult(float(gc), self.name, calibrated=True, units="test")


def test_rerank_with_calibrated_backend_may_steer_delivery() -> None:
    frontier = api.frontier(_PROTEIN, OptimizeConfig(gc_weight=1.0), steps=9)
    reranked = rerank_by_expression(frontier, _FakeCalibrated())
    scores = [r.audit["expression_score"] for r in reranked.results]
    # The chosen point is now the highest-scoring one under the calibrated head.
    assert reranked.frontier.chosen == max(range(len(scores)), key=lambda i: scores[i])
    for r in reranked.results:
        assert r.audit["expression_calibrated"] is True


def test_rerank_preserves_results_and_manifest() -> None:
    frontier = api.frontier(_PROTEIN, steps=5)
    reranked = rerank_by_expression(frontier)
    assert len(reranked.results) == len(frontier.results)
    assert reranked.manifest == frontier.manifest
    # Sequences are untouched (annotation only).
    assert [r.dna for r in reranked.results] == [r.dna for r in frontier.results]
