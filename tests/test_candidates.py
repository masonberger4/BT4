"""Tests for candidate-set assembly + expression rerank (design-flow step 3).

These pin the honesty/correctness rules hardened in design review: an uncalibrated
head only annotates (discovery order, solver-delivered chosen); a calibrated head
reorders and re-picks; the solver-delivered sequence is invariant to ``n`` (the cap
never drops it, and — calibrated — is applied after scoring); de-dup/cap counts are
reported; repeat-refined variants fire only when the seed violates a GLOBAL rule;
scoring uses the batch path when available; and everything is deterministic (#7).
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from bt4.api import candidates, frontier
from bt4.biomodels.expression import (
    BatchExpressionPredictor,
    ExpressionResult,
    NullExpressionModel,
)
from bt4.domain.genetic_code import translate
from bt4.pipeline.optimize import OptimizeConfig

_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAV"


@dataclass(frozen=True)
class _GcHead:
    """A deterministic *calibrated* fake head: score = G/C fraction (larger better)."""

    name: str = "gc_fake"
    calibrated: bool = True

    def score_sequence(self, dna: str) -> ExpressionResult:
        gc = sum(1 for c in dna if c in "GC") / len(dna)
        return ExpressionResult(score=gc, model_name=self.name, calibrated=True, units="gc")


@dataclass(frozen=True)
class _BatchGcHead(_GcHead):
    """The same head, but also batched -- structurally a BatchExpressionPredictor."""

    name: str = "gc_fake_batch"

    def score_many(self, dnas: list[str]) -> list[ExpressionResult]:
        return [self.score_sequence(d) for d in dnas]


# --------------------------------------------------------------------------- #
# Uncalibrated (default): annotate only, discovery order, delivered pinned.
# --------------------------------------------------------------------------- #


def test_default_is_discovery_order_and_no_op() -> None:
    cs = candidates(_PROTEIN, OptimizeConfig(), n=10)
    assert cs.calibrated is False
    assert cs.order_basis == "discovery"
    # Every score is the Null placeholder's 0.0 and is honestly labelled uncalibrated.
    assert all(c.expression_score == 0.0 for c in cs.candidates)
    assert all(c.expression_calibrated is False for c in cs.candidates)
    # The delivered candidate is the solver's own delivered frontier sequence.
    fr = frontier(_PROTEIN, OptimizeConfig())
    delivered = cs.delivered()
    assert delivered is not None
    assert delivered.result.dna == fr.delivered().dna
    assert cs.chosen == 0
    # And it round-trips to the protein (invariant #1).
    assert translate(delivered.result.dna).rstrip("*") == _PROTEIN


def test_default_predictor_is_null_placeholder() -> None:
    cs = candidates(_PROTEIN, OptimizeConfig(), n=5)
    assert cs.candidates[0].expression_model == NullExpressionModel().name


# --------------------------------------------------------------------------- #
# Calibrated: reorder by score, re-pick the top; solver-delivered still retained.
# --------------------------------------------------------------------------- #


def test_calibrated_reorders_and_repicks() -> None:
    cs = candidates(_PROTEIN, OptimizeConfig(), n=10, predictor=_GcHead())
    assert cs.calibrated is True
    assert cs.order_basis == "expression_rank"
    scores = [c.expression_score for c in cs.candidates]
    assert scores == sorted(scores, reverse=True)  # ranked best-first
    assert cs.chosen == 0  # delivered = top predicted expression
    # The solver-delivered sequence is retained for transparency even when a
    # calibrated head steers delivery elsewhere.
    fr = frontier(_PROTEIN, OptimizeConfig())
    assert any(c.result.dna == fr.delivered().dna for c in cs.candidates)


def test_delivered_retained_under_tight_cap() -> None:
    # n=1 must never drop the sequence that would ship.
    fr = frontier(_PROTEIN, OptimizeConfig())
    for predictor in (None, _GcHead()):
        cs = candidates(_PROTEIN, OptimizeConfig(), n=1, predictor=predictor)
        assert len(cs.candidates) == 1
        if predictor is None:
            # Uncalibrated: the retained one IS the solver-delivered sequence.
            assert cs.candidates[0].result.dna == fr.delivered().dna
        # Either way the delivered index is valid.
        assert cs.delivered() is not None


# --------------------------------------------------------------------------- #
# Batched scoring path.
# --------------------------------------------------------------------------- #


def test_batch_path_used_when_available_and_agrees() -> None:
    batch = _BatchGcHead()
    assert isinstance(batch, BatchExpressionPredictor)
    cs_batch = candidates(_PROTEIN, OptimizeConfig(), n=10, predictor=batch)
    cs_seq = candidates(_PROTEIN, OptimizeConfig(), n=10, predictor=_GcHead())
    assert cs_batch.scored_batched is True
    assert cs_seq.scored_batched is False
    # score_many and score_sequence must yield identical scores (batch contract).
    by_dna_batch = {c.result.dna: c.expression_score for c in cs_batch.candidates}
    by_dna_seq = {c.result.dna: c.expression_score for c in cs_seq.candidates}
    assert by_dna_batch == by_dna_seq


# --------------------------------------------------------------------------- #
# Repeat-refined variants: only when the seed violates a GLOBAL rule.
# --------------------------------------------------------------------------- #


def test_repeat_refined_variants_when_seed_violates_global_rule() -> None:
    # A repetitive protein under a tight max-repeat: the exact-DP seed violates it,
    # so refinement variants are drawn from the delivered seed.
    cfg = OptimizeConfig(max_repeat_length=8)
    cs = candidates("AAAAAAAAAAAAKKKKKKKKKK", cfg, n=20, repeat_variants=4)
    assert cs.n_repeat_refined > 0
    assert "repeat_refined" in {c.source for c in cs.candidates}
    assert "refined from the delivered seed" in cs.repeat_note


def test_no_repeat_variants_without_global_rule() -> None:
    cs = candidates(_PROTEIN, OptimizeConfig(), n=10, repeat_variants=4)
    assert cs.n_repeat_refined == 0
    assert "no GLOBAL rule active" in cs.repeat_note
    assert all(c.source == "frontier" for c in cs.candidates)


def test_repeat_variants_zero_is_noted() -> None:
    cfg = OptimizeConfig(max_repeat_length=8)
    cs = candidates("AAAAAAAAAAAAKKKKKKKKKK", cfg, n=20, repeat_variants=0)
    assert cs.n_repeat_refined == 0
    assert cs.repeat_note == "repeat_variants=0"


# --------------------------------------------------------------------------- #
# Bookkeeping, determinism, provenance, validation.
# --------------------------------------------------------------------------- #


def test_counts_are_consistent() -> None:
    cfg = OptimizeConfig(max_repeat_length=8)
    cs = candidates("AAAAAAAAAAAAKKKKKKKKKK", cfg, n=3, repeat_variants=4)
    scored = cs.n_frontier + cs.n_repeat_refined - cs.n_dedup_dropped
    # Everything scored is either kept or dropped by the cap; no silent loss.
    assert scored == len(cs.candidates) + cs.n_dropped_cap
    assert cs.n_dropped_cap >= 0
    assert len(cs.candidates) <= 3


def test_determinism() -> None:
    a = candidates(_PROTEIN, OptimizeConfig(), n=10, predictor=_GcHead())
    b = candidates(_PROTEIN, OptimizeConfig(), n=10, predictor=_GcHead())
    assert [c.result.dna for c in a.candidates] == [c.result.dna for c in b.candidates]
    assert a.manifest.config_hash == b.manifest.config_hash


def test_manifest_differs_by_predictor_n_and_steps() -> None:
    base = candidates(_PROTEIN, OptimizeConfig(), n=10, steps=11)
    diff_n = candidates(_PROTEIN, OptimizeConfig(), n=9, steps=11)
    diff_pred = candidates(_PROTEIN, OptimizeConfig(), n=10, steps=11, predictor=_GcHead())
    diff_steps = candidates(_PROTEIN, OptimizeConfig(), n=10, steps=7)
    assert base.manifest.config_hash != diff_n.manifest.config_hash
    assert base.manifest.config_hash != diff_pred.manifest.config_hash
    assert base.manifest.config_hash != diff_steps.manifest.config_hash


def test_calibrated_tight_cap_delivers_top_score() -> None:
    # n=1 with a calibrated head must deliver the head's TOP pick (not the
    # solver-delivered seed) and label the ordering honestly.
    cs = candidates(_PROTEIN, OptimizeConfig(), n=1, predictor=_GcHead())
    assert len(cs.candidates) == 1
    assert cs.order_basis == "expression_rank"
    full = candidates(_PROTEIN, OptimizeConfig(), n=50, predictor=_GcHead())
    top_score = max(c.expression_score for c in full.candidates)
    assert cs.candidates[0].expression_score == top_score  # the reranker's best, not dropped


def test_input_validation() -> None:
    with pytest.raises(ValueError, match="n must be"):
        candidates(_PROTEIN, OptimizeConfig(), n=0)
    with pytest.raises(ValueError, match="repeat_variants"):
        candidates(_PROTEIN, OptimizeConfig(), repeat_variants=-1)
    with pytest.raises(ValueError, match="steps"):
        candidates(_PROTEIN, OptimizeConfig(), steps=0)
    with pytest.raises(ValueError):
        candidates("not a protein!!", OptimizeConfig())
