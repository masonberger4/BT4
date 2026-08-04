"""Tests for the model-agnostic expression acceptance gate.

These pin the honesty behavior of :mod:`bt4.biomodels.expression.gate`: it splits
by group with no leakage, reports Spearman/Pearson/R^2 + conformal coverage, and
``passed`` requires both point accuracy and honest uncertainty. Nothing here needs
torch or a real model -- the gate is model-agnostic and runs on plain numbers or a
tiny fake predictor.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from bt4.biomodels.expression import (
    ExpressionEvalCase,
    ExpressionResult,
    NullExpressionModel,
    run_expression_gate,
    verify_expression_gate,
)


def _cases_from(
    triples: list[tuple[float, float, str]],
) -> list[ExpressionEvalCase]:
    return [ExpressionEvalCase(predicted=p, measured=m, group=g) for p, m, g in triples]


def test_perfect_predictor_passes() -> None:
    # predicted == measured across many groups: Spearman 1.0, residuals 0, coverage 1.0.
    cases = _cases_from([(float(i), float(i), f"g{i}") for i in range(12)])
    report = verify_expression_gate(cases, min_spearman=0.5, coverage_tolerance=0.2)
    assert report.passed is True
    assert report.spearman == pytest.approx(1.0)
    assert report.r2 == pytest.approx(1.0)
    assert report.n_groups == 12
    # Calibration residuals are all 0, so the interval is a point and still covers.
    assert report.empirical_coverage == pytest.approx(1.0)


def test_uncorrelated_predictor_fails_spearman() -> None:
    # Predictions unrelated to truth: Spearman near 0 < threshold => fail.
    triples = [
        (0.0, 5.0, "a"), (5.0, 1.0, "a"), (1.0, 4.0, "b"), (4.0, 2.0, "b"),
        (2.0, 3.0, "c"), (3.0, 0.0, "c"), (0.5, 4.5, "d"), (4.5, 0.5, "d"),
    ]
    report = verify_expression_gate(_cases_from(triples), min_spearman=0.5)
    assert report.passed is False


def test_grouped_split_has_no_leakage() -> None:
    # Two members per group; calibration and test groups must be disjoint, so the
    # per-fold counts reflect whole groups moving together, never a split group.
    triples = [(float(i), float(i), f"g{i // 2}") for i in range(12)]  # 6 groups x2
    report = verify_expression_gate(_cases_from(triples), calibration_fraction=0.5)
    assert report.n_calibration + report.n_test == 12
    assert report.n_calibration % 2 == 0  # whole groups (size 2) only
    assert report.n_test % 2 == 0


def test_requires_at_least_two_groups() -> None:
    cases = _cases_from([(1.0, 1.0, "only"), (2.0, 2.0, "only"), (3.0, 3.0, "only")])
    with pytest.raises(ValueError, match="two distinct groups"):
        verify_expression_gate(cases)


def test_requires_two_test_points() -> None:
    # Two groups but the test fold ends up with a single point => undefined corr.
    cases = _cases_from([(1.0, 1.0, "a"), (2.0, 2.0, "a"), (3.0, 3.0, "b")])
    with pytest.raises(ValueError, match="need >= 2"):
        verify_expression_gate(cases, calibration_fraction=0.5)


def test_argument_range_guards() -> None:
    cases = _cases_from([(float(i), float(i), f"g{i}") for i in range(6)])
    with pytest.raises(ValueError, match="target_coverage"):
        verify_expression_gate(cases, target_coverage=1.0)
    with pytest.raises(ValueError, match="calibration_fraction"):
        verify_expression_gate(cases, calibration_fraction=0.0)


def test_coverage_gap_can_fail_even_with_good_ranking() -> None:
    # Monotone predictions (Spearman 1.0) but a systematic offset on the TEST groups
    # that the calibration residuals underestimate => empirical coverage below target,
    # so the honesty (coverage) half of the gate fails despite perfect ranking. Uses
    # 10 calibration groups (>= 9, so the 90% conformal quantile is finite: with all
    # calibration residuals 0, q_hat = 0) and 4 offset test groups (residual 10).
    cal = [(float(i), float(i), f"cal{i}") for i in range(10)]  # zero residual
    test = [(float(i), float(i) + 10.0, f"test{i}") for i in range(4)]  # big offset
    report = verify_expression_gate(
        _cases_from(cal + test),
        calibration_fraction=0.7,  # first 10 (all cal*) groups => calibration
        min_spearman=0.5,
        target_coverage=0.9,
        coverage_tolerance=0.05,
    )
    assert report.n_calibration == 10
    assert report.n_test == 4
    assert report.conformal_half_width == pytest.approx(0.0)
    assert report.spearman == pytest.approx(1.0)  # ranking is perfect
    assert report.empirical_coverage == pytest.approx(0.0)  # but nothing is covered
    assert report.passed is False  # honesty half of the gate fails


@dataclass(frozen=True)
class _FakePredictor:
    """A deterministic calibrated-looking head for the wrapper test (measured+noise)."""

    name: str = "fake"
    calibrated: bool = True

    def score_sequence(self, dna: str) -> ExpressionResult:
        # Score by GC fraction so ranking is deterministic and sequence-driven.
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        return ExpressionResult(score=gc, model_name=self.name, calibrated=True, units="log-TE")


def test_run_expression_gate_wraps_a_predictor() -> None:
    # measured == the predictor's own GC score, so ranking is perfect.
    seqs = ["ATG", "GCG", "AAT", "CCG", "ATA", "GGC"]
    samples = []
    fake = _FakePredictor()
    for i, s in enumerate(seqs):
        measured = fake.score_sequence(s).score
        samples.append((s, measured, f"g{i}"))
    report = run_expression_gate(fake, samples, min_spearman=0.5, coverage_tolerance=0.2)
    assert report.spearman == pytest.approx(1.0)
    assert report.passed is True


def test_null_model_does_not_pass_the_gate() -> None:
    # The placeholder scores every sequence 0.0: zero-variance predictions => Spearman
    # 0.0, so it can never pass -- exactly the honesty guarantee (an uncalibrated
    # placeholder never self-certifies).
    seqs = ["ATG", "GCG", "AAT", "CCG", "ATA", "GGC"]
    samples = [(s, float(i), f"g{i}") for i, s in enumerate(seqs)]
    report = run_expression_gate(NullExpressionModel(), samples, min_spearman=0.3)
    assert report.spearman == pytest.approx(0.0)
    assert report.passed is False
