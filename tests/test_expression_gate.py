"""Tests for the model-agnostic expression acceptance gate.

These pin the honesty behavior of :mod:`bt4.biomodels.expression.gate`: it splits
by group with no leakage, reports Spearman/Pearson/R^2 + conformal coverage, and
``passed`` requires both point accuracy and honest uncertainty. Nothing here needs
torch or a real model -- the gate is model-agnostic and runs on plain numbers or a
tiny fake predictor.
"""

from __future__ import annotations

import math
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


# --- Defect A: pooled scoring credits between-group skill ----------------------
#
# THE load-bearing test of this whole exercise. A head that knows only "which protein
# is this" -- which is what training across natural genes teaches, and what BT4 cannot
# use -- must pass the pooled gate and FAIL the within-group one. If the second
# assertion ever stops holding, the strict bar is not being enforced and the gate can
# hand out a false pass.


def _gene_identity_panel() -> list[ExpressionEvalCase]:
    """Four proteins x four variants. Prediction = the protein's baseline, exactly.

    Within any one protein the prediction is constant, so there is no ordering
    information at all; across proteins it is perfect.
    """
    baselines = {"g1": 1.0, "g2": 5.0, "g3": 9.0, "g4": 13.0}
    cases: list[ExpressionEvalCase] = []
    for group, baseline in baselines.items():
        for offset in (-0.3, -0.1, 0.1, 0.3):
            cases.append(
                ExpressionEvalCase(
                    predicted=baseline,  # knows the gene, nothing about the variant
                    measured=baseline + offset,
                    group=group,
                )
            )
    return cases


def test_gene_identity_head_passes_pooled_but_fails_within_group() -> None:
    cases = _gene_identity_panel()

    pooled = verify_expression_gate(cases, min_spearman=0.3, bootstrap_resamples=0)
    strict = verify_expression_gate(
        cases, min_spearman=0.3, within_group=True, bootstrap_resamples=0
    )

    # Pooled: high rank agreement, purely from between-protein differences. (Not 1.0
    # only because the head ties every variant within a protein, and tie-averaged
    # ranks cost a little.)
    assert pooled.spearman > 0.8
    # Within-group: no ordering information survives centring, so it scores ~0 and the
    # Spearman half of the gate fails.
    assert strict.spearman == pytest.approx(0.0)
    assert strict.passed is False
    assert strict.within_group is True


def test_within_group_credits_a_head_that_orders_variants() -> None:
    # The converse: a head with NO between-protein skill but perfect within-protein
    # ordering is invisible to the pooled metric and correctly credited by the strict
    # one. This is exactly BT4's regime.
    cases: list[ExpressionEvalCase] = []
    for group, baseline in (("g1", 1.0), ("g2", 5.0), ("g3", 9.0), ("g4", 13.0)):
        for offset in (-0.3, -0.1, 0.1, 0.3):
            cases.append(
                ExpressionEvalCase(
                    predicted=offset,  # knows the variant, nothing about the gene
                    measured=baseline + offset,
                    group=group,
                )
            )

    strict = verify_expression_gate(
        cases, min_spearman=0.3, within_group=True, bootstrap_resamples=0
    )
    assert strict.spearman == pytest.approx(1.0)
    assert strict.n_groups_ranked == 2  # the two test-fold groups
    assert dict(strict.per_group_spearman) == {"g3": pytest.approx(1.0), "g4": pytest.approx(1.0)}


def test_within_group_skips_singleton_groups_rather_than_scoring_them_zero() -> None:
    # A singleton has no ordering to get right; counting it as 0.0 would silently drag
    # the aggregate down and make a good head look mediocre.
    cases = [
        ExpressionEvalCase(predicted=1.0, measured=1.0, group="a1"),
        ExpressionEvalCase(predicted=2.0, measured=2.0, group="a1"),
        ExpressionEvalCase(predicted=3.0, measured=3.0, group="a2"),
        ExpressionEvalCase(predicted=1.0, measured=1.0, group="b1"),
        ExpressionEvalCase(predicted=2.0, measured=2.0, group="b1"),
        ExpressionEvalCase(predicted=9.0, measured=9.0, group="b2"),
    ]
    report = verify_expression_gate(
        cases, min_spearman=0.3, within_group=True, bootstrap_resamples=0
    )
    assert report.n_groups_test == 2  # b1 and b2
    assert report.n_groups_ranked == 1  # only b1 has >= 2 members
    assert report.spearman == pytest.approx(1.0)  # not dragged toward 0 by b2


def test_within_group_raises_when_every_test_group_is_a_singleton() -> None:
    cases = [
        ExpressionEvalCase(predicted=float(i), measured=float(i), group=f"g{i}")
        for i in range(6)
    ]
    with pytest.raises(ValueError, match="no ordering to score"):
        verify_expression_gate(cases, within_group=True, bootstrap_resamples=0)


# --- Defect B: raw residuals across different units are vacuous ----------------


def _scaled_panel(slope: float, intercept: float) -> list[ExpressionEvalCase]:
    """A perfectly rank-correct head reporting on a completely different scale."""
    cases: list[ExpressionEvalCase] = []
    for group in range(8):
        for step in range(4):
            predicted = float(group * 4 + step) / 10.0
            cases.append(
                ExpressionEvalCase(
                    predicted=predicted,
                    measured=slope * predicted + intercept,
                    group=f"g{group:02d}",
                )
            )
    return cases


def test_recalibration_rescues_a_rank_perfect_head_on_the_wrong_scale() -> None:
    cases = _scaled_panel(slope=7.0, intercept=100.0)

    raw = verify_expression_gate(cases, min_spearman=0.3, bootstrap_resamples=0)
    linked = verify_expression_gate(
        cases, min_spearman=0.3, recalibrate=True, bootstrap_resamples=0
    )

    # Ranking is scale-free, so both see a perfect Spearman ...
    assert raw.spearman == pytest.approx(1.0)
    assert linked.spearman == pytest.approx(1.0)
    # ... but only the linked run produces an interval worth anything. Unrecalibrated,
    # the residuals are dominated by the ruler mismatch.
    assert raw.width_over_iqr > 1.0
    assert linked.width_over_iqr < 1e-6
    assert linked.slope == pytest.approx(7.0)
    assert linked.intercept == pytest.approx(100.0)


def test_a_negative_link_slope_is_reported_not_hidden() -> None:
    # A head that ranks BACKWARDS on this panel. The fitted link makes it usable, but
    # the sign must be visible in the report -- the sign convention of a CLR residual
    # against an assay readout is not guaranteed and must never be assumed.
    cases = _scaled_panel(slope=-3.0, intercept=2.0)
    report = verify_expression_gate(
        cases, min_spearman=0.3, recalibrate=True, bootstrap_resamples=0
    )
    assert report.slope < 0.0
    assert report.spearman == pytest.approx(-1.0)  # raw ranking really is inverted
    assert report.passed is False  # and the gate does not credit it


def test_link_slope_spread_flags_a_link_that_does_not_transfer() -> None:
    # Same head, but each protein needs a different slope. Coverage on this panel may
    # look fine; the spread is what says the link does not generalise.
    cases: list[ExpressionEvalCase] = []
    for group in range(8):
        group_slope = 1.0 + group  # a different ruler per protein
        for step in range(4):
            predicted = float(step)
            cases.append(
                ExpressionEvalCase(
                    predicted=predicted,
                    measured=group_slope * predicted,
                    group=f"g{group:02d}",
                )
            )
    report = verify_expression_gate(
        cases, min_spearman=0.3, recalibrate=True, bootstrap_resamples=0
    )
    assert report.link_slope_spread > 0.5


# --- vacuity: a constant predictor must not be able to game the report ---------


def test_a_constant_predictor_gets_valid_coverage_and_is_caught_anyway() -> None:
    # Split conformal is valid for ANY score function, so a constant predictor achieves
    # correct coverage with a uselessly wide interval. Coverage alone therefore cannot
    # be the gate; the report must catch this on both the rank and the width axes.
    cases: list[ExpressionEvalCase] = []
    for group in range(10):
        for step in range(4):
            cases.append(
                ExpressionEvalCase(
                    predicted=42.0,  # says the same thing about everything
                    measured=float(group) + step * 0.5,
                    group=f"g{group:02d}",
                )
            )
    report = verify_expression_gate(cases, min_spearman=0.3, bootstrap_resamples=0)

    assert report.spearman == pytest.approx(0.0)  # zero-variance predictions
    assert report.passed is False  # caught on the rank axis
    assert report.width_over_iqr > 1.0  # and visibly vacuous on the width axis


# --- cluster bootstrap ---------------------------------------------------------


def test_bootstrap_ci_is_deterministic_and_brackets_the_estimate() -> None:
    cases = _scaled_panel(slope=1.0, intercept=0.0)
    first = verify_expression_gate(cases, bootstrap_resamples=200, bootstrap_seed=7)
    again = verify_expression_gate(cases, bootstrap_resamples=200, bootstrap_seed=7)

    assert (first.spearman_ci_low, first.spearman_ci_high) == (
        again.spearman_ci_low,
        again.spearman_ci_high,
    )  # invariant #7
    assert first.bootstrap_resamples == 200
    assert first.spearman_ci_low <= first.spearman <= first.spearman_ci_high


def test_bootstrap_resamples_whole_groups_so_the_ci_is_not_falsely_narrow() -> None:
    # A noisy head: the CI must be wide enough to admit real uncertainty. Resampling
    # individual cases would treat one protein's variants as independent observations
    # and shrink this interval dramatically.
    measured = [0.0, 3.0, 1.0, 4.0, 2.0, 5.0, 7.0, 6.0, 9.0, 8.0, 11.0, 10.0]
    cases = [
        ExpressionEvalCase(
            predicted=float(i), measured=measured[i], group=f"g{i // 2:02d}"
        )
        for i in range(12)
    ]
    report = verify_expression_gate(cases, bootstrap_resamples=500, bootstrap_seed=1)
    assert report.bootstrap_resamples > 0
    assert report.spearman_ci_high - report.spearman_ci_low > 0.0


def test_bootstrap_is_disabled_honestly_rather_than_faked() -> None:
    cases = _scaled_panel(slope=1.0, intercept=0.0)
    report = verify_expression_gate(cases, bootstrap_resamples=0)
    assert report.bootstrap_resamples == 0
    assert math.isnan(report.spearman_ci_low)
    assert math.isnan(report.spearman_ci_high)


# --- the within-group scope caveat --------------------------------------------


def test_within_group_marks_its_coverage_as_anchor_conditional() -> None:
    # In within-group mode the target is an offset from the protein's own baseline, so
    # the interval is only achievable at design time if a member of that protein has
    # been measured. That is a narrower claim and must be stamped as one.
    pooled = verify_expression_gate(_scaled_panel(1.0, 0.0), bootstrap_resamples=0)
    strict = verify_expression_gate(
        _scaled_panel(1.0, 0.0), within_group=True, bootstrap_resamples=0
    )
    assert pooled.coverage_conditional_on_group_anchor is False
    assert strict.coverage_conditional_on_group_anchor is True


# --- batching -----------------------------------------------------------------


def test_run_expression_gate_uses_the_batch_surface_when_available() -> None:
    # Scoring a panel row-by-row through RiboNN would multiply the wall clock by the
    # row count and re-hash 90 weight files each time.
    calls: list[int] = []

    class _BatchHead:
        name = "batch"
        calibrated = False

        def score_sequence(self, dna: str) -> ExpressionResult:
            raise AssertionError("score_sequence must not be used when batching")

        def score_many(self, dnas: list[str]) -> list[ExpressionResult]:
            calls.append(len(dnas))
            return [
                ExpressionResult(
                    score=float(len(dna)), model_name="batch", calibrated=False, units="u"
                )
                for dna in dnas
            ]

    samples = [
        ("ATG" * (i + 1) + "TAA", float(i), f"g{i // 2:02d}") for i in range(8)
    ]
    report = run_expression_gate(_BatchHead(), samples, bootstrap_resamples=0)

    assert calls == [8]  # exactly one invocation for the whole panel
    assert report.n_calibration + report.n_test == 8
