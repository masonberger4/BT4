"""Tests for the binary-classification estimators in :mod:`bt4.biomodels._stats`.

These back the splice acceptance gate, where the task is per-position "is this a
splice site?" -- severely imbalanced, a handful of sites among thousands of
positions. CLAUDE.md section 6 names the metric set (**PR-AUC / MCC / ECE /
Brier, never bare accuracy**) precisely because that imbalance makes accuracy
meaningless: calling "not a site" everywhere scores >99%.

The tests are organised around the ways a metric can *lie* on this task:

* **PR-AUC** must be average precision, not trapezoidal -- linear interpolation in
  PR space is unachievable and optimistically biased (Davis & Goadrich 2006);
* **ROC-AUC** must stay tie-aware, and is pinned here *showing* that it flatters a
  model with hopeless precision, which is why it is never the headline;
* **Brier** looks tiny for a useless model on an imbalanced set, so
  **brier_skill_score** must score the base-rate predictor exactly 0;
* **ECE** cannot distinguish informative from vacuous calibration, so a
  base-rate-everywhere model is pinned as scoring *well* on it -- the trap the gate
  must report around;
* **top-k accuracy** must resolve ties pessimistically, or a model scoring
  everything equally would claim a perfect score.

Values are checked against hand-computed references (and, where noted, against
what scikit-learn's equivalents return) rather than against the implementation.
"""

from __future__ import annotations

import pytest

from bt4.biomodels._stats import (
    brier_score,
    brier_skill_score,
    expected_calibration_error,
    mcc,
    pr_auc,
    reliability_bins,
    roc_auc,
    top_k_accuracy,
)

# --------------------------------------------------------------------------
# PR-AUC


def test_pr_auc_perfect_and_reference_value() -> None:
    """Perfect ranking gives 1.0; a mixed case matches average precision exactly."""
    assert pr_auc([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]) == pytest.approx(1.0)
    # sklearn.metrics.average_precision_score([0,0,1,1],[0.1,0.4,0.35,0.8]) == 0.8333...
    assert pr_auc([0, 0, 1, 1], [0.1, 0.4, 0.35, 0.8]) == pytest.approx(5 / 6)


def test_pr_auc_is_average_precision_not_trapezoidal() -> None:
    """Pin the step-wise definition, which is *lower* than the trapezoidal one.

    Trapezoidal interpolation between PR operating points is not achievable by any
    classifier and inflates the score; on the imbalanced splice task that bias
    lands exactly where the honest number matters.
    """
    labels = [1, 0, 0, 0, 1]
    scores = [0.9, 0.8, 0.7, 0.6, 0.5]
    # Operating points: (R=0.5,P=1.0) then (R=1.0,P=0.4).
    # Average precision  = 0.5*1.0 + 0.5*0.4 = 0.70
    # Trapezoidal would be 0.5*1.0 + 0.5*(1.0+0.4)/2 = 0.85 -- inflated.
    assert pr_auc(labels, scores) == pytest.approx(0.70)


def test_pr_auc_groups_ties_into_one_operating_point() -> None:
    """Equal scores cannot be flattered by a lucky sort order."""
    tied = pr_auc([1, 0, 1, 0], [0.5, 0.5, 0.5, 0.5])
    assert tied == pytest.approx(0.5)  # the base rate, not 1.0


def test_pr_auc_without_positives_is_zero() -> None:
    """No positives means no skill is demonstrable -- report it, do not raise."""
    assert pr_auc([0, 0, 0], [0.9, 0.5, 0.1]) == 0.0


# --------------------------------------------------------------------------
# ROC-AUC


@pytest.mark.parametrize(
    ("labels", "scores", "expected"),
    [
        ([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], 1.0),
        ([0, 0, 1, 1], [0.9, 0.8, 0.2, 0.1], 0.0),
        ([0, 0, 1, 1], [0.5, 0.5, 0.5, 0.5], 0.5),
        ([0, 0, 1, 1], [0.1, 0.4, 0.35, 0.8], 0.75),
        ([0, 1], [0.5, 0.5], 0.5),
    ],
)
def test_roc_auc_reference_values(
    labels: list[int], scores: list[float], expected: float
) -> None:
    """Tie-aware ROC-AUC against hand-computed references."""
    assert roc_auc(labels, scores) == pytest.approx(expected)


def test_roc_auc_flatters_a_model_with_hopeless_precision() -> None:
    """The reason ROC-AUC is reported but never the headline.

    One true site among a thousand positions, ranked above every negative but with
    a hundred false positives above the rest: ROC-AUC is near-perfect while
    precision at that operating point is 1%.
    """
    labels = [1] + [0] * 999
    scores = [1.0] + [0.9] * 100 + [0.1] * 899
    assert roc_auc(labels, scores) > 0.94
    assert pr_auc(labels, scores) == pytest.approx(1.0)  # the single positive ranks first
    # ...but push it below the false positives and PR-AUC collapses while ROC stays fair:
    scores_worse = [0.5] + [0.9] * 100 + [0.1] * 899
    assert roc_auc(labels, scores_worse) > 0.85
    assert pr_auc(labels, scores_worse) < 0.02


def test_roc_auc_single_class_is_zero() -> None:
    """With one class absent, no discrimination is measurable."""
    assert roc_auc([1, 1], [0.1, 0.2]) == 0.0
    assert roc_auc([0, 0], [0.1, 0.2]) == 0.0


# --------------------------------------------------------------------------
# MCC


@pytest.mark.parametrize(
    ("labels", "preds", "expected"),
    [
        ([1, 1, 0, 0], [1, 1, 0, 0], 1.0),
        ([1, 1, 0, 0], [0, 0, 1, 1], -1.0),
        ([1, 1, 0, 0], [1, 0, 1, 0], 0.0),
    ],
)
def test_mcc_reference_values(
    labels: list[int], preds: list[int], expected: float
) -> None:
    """MCC over the standard confusion-matrix corner cases."""
    assert mcc(labels, preds) == pytest.approx(expected)


def test_mcc_punishes_the_all_negative_degenerate() -> None:
    """The metric that makes 'call nothing a site' score zero, not 99%."""
    labels = [1] + [0] * 99
    all_negative = [0] * 100
    accuracy = sum(1 for t, p in zip(labels, all_negative, strict=True) if t == p) / 100
    assert accuracy == pytest.approx(0.99)  # what bare accuracy would report
    assert mcc(labels, all_negative) == 0.0  # what MCC reports


def test_mcc_rejects_non_binary_predictions() -> None:
    """A thresholded prediction must be 0/1, not a score."""
    with pytest.raises(ValueError, match="predictions must be 0/1"):
        mcc([1, 0], [1, 2])


# --------------------------------------------------------------------------
# Brier and skill


def test_brier_bounds() -> None:
    """Perfect is 0, maximally wrong is 1."""
    assert brier_score([1, 0], [1.0, 0.0]) == pytest.approx(0.0)
    assert brier_score([1, 0], [0.0, 1.0]) == pytest.approx(1.0)


def test_brier_looks_tiny_for_a_useless_imbalanced_model() -> None:
    """Why the raw Brier score must not be read alone."""
    labels = [1] + [0] * 999
    base_rate = [0.001] * 1000
    assert brier_score(labels, base_rate) < 0.002  # flatteringly small
    assert brier_skill_score(labels, base_rate) == pytest.approx(0.0, abs=1e-9)


def test_brier_skill_score_is_the_vacuous_pass_detector() -> None:
    """0 for the base-rate predictor, 1 for perfect, negative for worse."""
    labels = [1] + [0] * 999
    assert brier_skill_score(labels, [1.0] + [0.0] * 999) == pytest.approx(1.0)
    assert brier_skill_score(labels, [0.5] * 1000) < 0.0


def test_brier_skill_score_single_class_is_zero() -> None:
    """A perfect baseline leaves no skill to demonstrate."""
    assert brier_skill_score([0, 0, 0], [0.0, 0.0, 0.0]) == 0.0


# --------------------------------------------------------------------------
# ECE and reliability


def test_ece_zero_when_calibrated() -> None:
    """A bin whose mean confidence equals its observed rate contributes nothing."""
    assert expected_calibration_error([1, 0], [1.0, 0.0]) == pytest.approx(0.0)
    assert expected_calibration_error([1] * 5 + [0] * 5, [0.5] * 10) == pytest.approx(0.0)


def test_ece_measures_overconfidence() -> None:
    """Predicting 0.9 where the rate is 0.1 is an ECE of 0.8."""
    assert expected_calibration_error([1] + [0] * 9, [0.9] * 10) == pytest.approx(0.8)


def test_ece_cannot_detect_a_vacuous_model() -> None:
    """The documented blind spot, pinned so it is never mistaken for a pass.

    A predictor emitting the base rate everywhere is perfectly calibrated and
    completely uninformative. ECE says it is excellent; PR-AUC and the Brier skill
    score are what expose it.
    """
    labels = [1] * 10 + [0] * 990
    vacuous = [0.01] * 1000
    assert expected_calibration_error(labels, vacuous) < 0.01  # looks great
    assert brier_skill_score(labels, vacuous) == pytest.approx(0.0, abs=1e-9)
    assert pr_auc(labels, vacuous) == pytest.approx(0.01)  # the base rate


def test_reliability_bins_omit_empty_bins() -> None:
    """No point is reported where there is no evidence."""
    bins = reliability_bins([1, 0], [0.95, 0.05], n_bins=10)
    assert len(bins) == 2
    assert all(count > 0 for _, _, count in bins)


def test_reliability_bins_counts_sum_to_n() -> None:
    """Every case lands in exactly one bin, including p at the boundaries."""
    labels = [1, 0, 1, 0, 1]
    probs = [0.0, 1.0, 0.5, 0.999, 0.001]
    assert sum(c for _, _, c in reliability_bins(labels, probs, n_bins=4)) == len(labels)


def test_ece_rejects_non_positive_bins() -> None:
    """A bin count must be meaningful."""
    with pytest.raises(ValueError, match="n_bins must be positive"):
        expected_calibration_error([1, 0], [0.9, 0.1], n_bins=0)


# --------------------------------------------------------------------------
# top-k accuracy


def test_top_k_accuracy_reference_values() -> None:
    """k equals the positive count, so precision and recall coincide."""
    assert top_k_accuracy([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]) == pytest.approx(1.0)
    assert top_k_accuracy([0, 0, 1, 1], [0.9, 0.2, 0.8, 0.1]) == pytest.approx(0.5)


def test_top_k_accuracy_resolves_ties_pessimistically() -> None:
    """A model scoring everything equally gets the base rate, never 1.0."""
    assert top_k_accuracy([0, 0, 1, 1], [0.5] * 4) == pytest.approx(0.5)
    assert top_k_accuracy([1] + [0] * 9, [0.5] * 10) == pytest.approx(0.1)


def test_top_k_accuracy_without_positives_is_zero() -> None:
    """No sites to find means no score to claim."""
    assert top_k_accuracy([0, 0], [0.9, 0.1]) == 0.0


# --------------------------------------------------------------------------
# Shared validation


@pytest.mark.parametrize(
    "fn", [pr_auc, roc_auc, brier_score, brier_skill_score, top_k_accuracy]
)
def test_length_mismatch_raises(fn) -> None:  # type: ignore[no-untyped-def]
    """Mismatched series are a caller bug, not something to silently zip-truncate."""
    with pytest.raises(ValueError, match="differ in length"):
        fn([1, 0, 1], [0.5, 0.5])


@pytest.mark.parametrize(
    "fn", [pr_auc, roc_auc, brier_score, brier_skill_score, top_k_accuracy]
)
def test_empty_input_raises(fn) -> None:  # type: ignore[no-untyped-def]
    """An empty evaluation is a mistake worth surfacing."""
    with pytest.raises(ValueError, match="at least one case"):
        fn([], [])


@pytest.mark.parametrize(
    "fn", [pr_auc, roc_auc, brier_score, brier_skill_score, top_k_accuracy]
)
def test_non_binary_labels_raise(fn) -> None:  # type: ignore[no-untyped-def]
    """Ground truth must be 0/1; a probability passed as a label is a real bug."""
    with pytest.raises(ValueError, match="labels must be 0/1"):
        fn([0, 2], [0.5, 0.5])
