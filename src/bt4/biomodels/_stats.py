"""Shared, dependency-free statistics for BT4's biomodel evaluation gates.

Small, pure-standard-library estimators used by more than one biomodel: rank and
linear correlation (the splice agreement report, :mod:`bt4.biomodels.splice.agreement`),
the coefficient of determination, and split-conformal prediction helpers (the
expression acceptance gate, :mod:`bt4.biomodels.expression.gate`). Kept here so the
two do not each carry their own copy (and so a single well-tested implementation
backs both). No numpy; :mod:`math` only.

The correlation functions return an honest ``0.0`` -- "no detectable relationship"
-- when a series has zero variance, rather than raising or returning ``nan``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

__all__ = [
    "brier_score",
    "brier_skill_score",
    "conformal_quantile",
    "empirical_coverage",
    "expected_calibration_error",
    "iqr",
    "linear_fit",
    "mcc",
    "pearson",
    "pr_auc",
    "pr_auc_skill",
    "r2_score",
    "reliability_bins",
    "roc_auc",
    "spearman",
    "top_k_accuracy",
]


def _ranks(values: Sequence[float]) -> list[float]:
    """Return fractional (tie-averaged) ranks of ``values``.

    Ties share the average of the ranks they span, so a rank correlation on the
    result is well defined in the presence of equal values.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0  # 1-based, averaged over the tie block
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    """Return the Pearson (linear) correlation of ``x`` and ``y``.

    Returns:
        The correlation in ``[-1, 1]``, or ``0.0`` when either series has zero
        variance (correlation is undefined; ``0.0`` is the honest "no detectable
        linear relationship" default).

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    if len(x) != len(y):
        raise ValueError(f"series differ in length: {len(x)} vs {len(y)}")
    n = len(x)
    if n < 2:
        raise ValueError("correlation needs at least two points")
    mean_x = math.fsum(x) / n
    mean_y = math.fsum(y) / n
    dx = [xi - mean_x for xi in x]
    dy = [yi - mean_y for yi in y]
    cov = math.fsum(a * b for a, b in zip(dx, dy, strict=True))
    var_x = math.fsum(a * a for a in dx)
    var_y = math.fsum(b * b for b in dy)
    if var_x == 0.0 or var_y == 0.0:
        return 0.0
    return cov / math.sqrt(var_x * var_y)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Return the Spearman rank correlation of ``x`` and ``y``.

    Computed as the Pearson correlation of the tie-averaged ranks, so it measures
    monotonic (rank) agreement rather than linear agreement.

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    return pearson(_ranks(x), _ranks(y))


def r2_score(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """Return the coefficient of determination R^2 of a prediction.

    ``R^2 = 1 - SS_res / SS_tot`` where ``SS_res = sum (true - pred)^2`` and
    ``SS_tot = sum (true - mean(true))^2``. It can be negative when the prediction
    is worse than predicting the mean -- reported honestly, not clamped.

    Returns:
        The R^2 value, or ``0.0`` when the truth has zero variance (``SS_tot == 0``,
        so R^2 is undefined).

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    if len(y_true) != len(y_pred):
        raise ValueError(f"series differ in length: {len(y_true)} vs {len(y_pred)}")
    n = len(y_true)
    if n < 2:
        raise ValueError("r2_score needs at least two points")
    mean_true = math.fsum(y_true) / n
    ss_tot = math.fsum((t - mean_true) ** 2 for t in y_true)
    ss_res = math.fsum((t - p) ** 2 for t, p in zip(y_true, y_pred, strict=True))
    if ss_tot == 0.0:
        return 0.0
    return 1.0 - ss_res / ss_tot


def linear_fit(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    """Return the least-squares ``(slope, intercept)`` of ``y`` on ``x``.

    Used to fit the **link** between a model's arbitrary-unit output and an assay's
    units, on a held-out fold, before residuals mean anything: a rank-correct
    predictor on the wrong scale has huge raw residuals that say nothing about its
    quality (:mod:`bt4.biomodels.expression.gate`).

    Returns:
        ``(slope, intercept)``. When ``x`` has zero variance no slope is identifiable,
        and the honest answer is the intercept-only fit ``(0.0, mean(y))`` -- the same
        convention as the correlation functions returning ``0.0`` for "nothing
        detectable" rather than raising or returning ``nan``.

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    if len(x) != len(y):
        raise ValueError(f"series differ in length: {len(x)} vs {len(y)}")
    n = len(x)
    if n < 2:
        raise ValueError("linear_fit needs at least two points")
    mean_x = math.fsum(x) / n
    mean_y = math.fsum(y) / n
    dx = [xi - mean_x for xi in x]
    cov = math.fsum(a * (yi - mean_y) for a, yi in zip(dx, y, strict=True))
    var_x = math.fsum(a * a for a in dx)
    if var_x == 0.0:
        return 0.0, mean_y
    slope = cov / var_x
    return slope, mean_y - slope * mean_x


def iqr(values: Sequence[float]) -> float:
    """Return the interquartile range of ``values`` (linear-interpolated quartiles).

    The scale a conformal interval must be judged against: an interval is only useful
    if it is narrow *relative to the spread of the labels*. Coverage alone cannot say
    that -- a constant predictor achieves exactly valid coverage with a uselessly wide
    interval -- so the gate reports median interval width divided by this.

    Returns:
        ``Q3 - Q1``, or ``0.0`` for fewer than two values.

    Raises:
        ValueError: If ``values`` is empty.
    """
    n = len(values)
    if n < 1:
        raise ValueError("iqr needs at least one value")
    if n < 2:
        return 0.0
    ordered = sorted(values)

    def _quantile(q: float) -> float:
        position = q * (n - 1)
        low = math.floor(position)
        high = min(low + 1, n - 1)
        weight = position - low
        return ordered[low] * (1.0 - weight) + ordered[high] * weight

    return _quantile(0.75) - _quantile(0.25)


def conformal_quantile(scores: Sequence[float], coverage: float) -> float:
    """Return the finite-sample split-conformal quantile of nonconformity ``scores``.

    For a target ``coverage`` (1 - alpha) and ``n`` calibration nonconformity
    scores, the conformal quantile is the ``k``-th smallest score with
    ``k = ceil((n + 1) * coverage)``. When ``k > n`` (too few calibration points
    to guarantee the coverage) the quantile is ``+inf`` -- the honest answer that
    the interval must cover everything to make the finite-sample guarantee.

    Args:
        scores: Calibration nonconformity scores (e.g. absolute residuals).
        coverage: Target coverage in the open interval ``(0, 1)``.

    Returns:
        The conformal quantile (``math.inf`` when ``n`` is too small for the level).

    Raises:
        ValueError: If ``scores`` is empty or ``coverage`` is not in ``(0, 1)``.
    """
    if not 0.0 < coverage < 1.0:
        raise ValueError(f"coverage must be in (0, 1), got {coverage}")
    n = len(scores)
    if n < 1:
        raise ValueError("conformal_quantile needs at least one calibration score")
    k = math.ceil((n + 1) * coverage)
    if k > n:
        return math.inf
    return sorted(scores)[k - 1]


def empirical_coverage(residuals: Sequence[float], half_width: float) -> float:
    """Return the fraction of ``residuals`` within ``half_width`` (empirical coverage).

    With a symmetric conformal interval of half-width ``q`` around each prediction,
    a test point is covered iff its absolute residual is ``<= q``. This reports the
    realized coverage on a held-out set, to compare against the target.

    Args:
        residuals: Absolute test residuals ``|true - pred|``.
        half_width: The conformal interval half-width (:func:`conformal_quantile`).

    Returns:
        The covered fraction in ``[0, 1]``.

    Raises:
        ValueError: If ``residuals`` is empty.
    """
    n = len(residuals)
    if n < 1:
        raise ValueError("empirical_coverage needs at least one residual")
    return sum(1 for r in residuals if r <= half_width) / n


# ---------------------------------------------------------------------------
# Binary-classification estimators.
#
# Added for the splice acceptance gate, where the task is per-position "is this a
# splice site?" -- severely imbalanced (a handful of sites among thousands of
# positions). CLAUDE.md section 6 names the metric set: **PR-AUC / MCC / ECE /
# Brier, never bare accuracy**, and the imbalance is exactly why: a model calling
# "not a site" everywhere scores >99% accuracy while finding nothing.


def _check_binary(labels: Sequence[int], scores: Sequence[float]) -> None:
    """Validate a labels/scores pair for the classification estimators."""
    if len(labels) != len(scores):
        raise ValueError(f"series differ in length: {len(labels)} vs {len(scores)}")
    if not labels:
        raise ValueError("need at least one case")
    bad = {int(v) for v in labels} - {0, 1}
    if bad:
        raise ValueError(f"labels must be 0/1, found {sorted(bad)}")


def pr_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return the area under the precision-recall curve (**average precision**).

    Computed as average precision -- the step-wise sum ``sum (R_k - R_{k-1}) * P_k``
    -- deliberately **not** trapezoidal interpolation between operating points.
    Linear interpolation in PR space is not achievable by any classifier and is
    optimistically biased (Davis & Goadrich, ICML 2006); on the heavily imbalanced
    per-position splice task that bias would flatter the model precisely where the
    honest number matters.

    Ties are handled by grouping equal scores into a single operating point, so a
    model that assigns many positions the same score cannot gain from an arbitrary
    ordering among them.

    Args:
        labels: Ground truth, 0 or 1.
        scores: Predicted scores; higher means more likely positive. Need not be
            probabilities -- PR-AUC is rank-based.

    Returns:
        Average precision in ``[0, 1]``. Returns ``0.0`` when there are no
        positives (nothing to retrieve, so no skill is demonstrable) -- the honest
        "nothing detectable" convention this module uses elsewhere.

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, scores)
    n_pos = sum(1 for v in labels if v == 1)
    if n_pos == 0:
        return 0.0

    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    total = 0.0
    tp = 0
    fp = 0
    prev_recall = 0.0
    i = 0
    while i < len(order):
        # Consume every case sharing this score as one operating point.
        score = scores[order[i]]
        while i < len(order) and scores[order[i]] == score:
            if labels[order[i]] == 1:
                tp += 1
            else:
                fp += 1
            i += 1
        recall = tp / n_pos
        precision = tp / (tp + fp)
        total += (recall - prev_recall) * precision
        prev_recall = recall
    return total


def pr_auc_skill(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return PR-AUC rescaled against its no-skill floor: ``(AP - p) / (1 - p)``.

    Average precision has a floor at the panel's **prevalence** ``p`` and a ceiling
    at 1, so a raw PR-AUC cannot be compared across panels with different positive
    rates -- and prevalence here is a *construction choice* (how many negatives the
    panel samples), not a constant of nature. Measured on this module's own
    :func:`pr_auc` with model quality held fixed, thinning negatives moved AP from
    0.88 to 0.98 while :func:`roc_auc` barely moved.

    A **ratio** (``AP / p``) does not fix this: its ceiling is ``1 / p``, so it drifts
    with prevalence too and systematically rewards the sparser panel. The skill form
    is 0.0 at no-skill and 1.0 at perfect for **every** prevalence, so it reads the
    same way as :func:`brier_skill_score` and the expression gate's
    ``width_over_iqr``.

    Returns:
        The skill score; negative means worse than predicting at random. Returns
        ``0.0`` when prevalence is 0 or 1 (no skill is demonstrable).

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, scores)
    prevalence = math.fsum(float(v) for v in labels) / len(labels)
    if prevalence in (0.0, 1.0):
        return 0.0
    return (pr_auc(labels, scores) - prevalence) / (1.0 - prevalence)


def roc_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return the area under the ROC curve, tie-aware.

    Computed from tie-averaged ranks (the Mann-Whitney U identity), so equal scores
    contribute 0.5 rather than depending on sort order.

    **Read this alongside, never instead of, :func:`pr_auc`.** On a task where
    positives are ~0.1% of positions, ROC-AUC stays high for a model with hopeless
    precision, because the false-positive rate is divided by an enormous negative
    count. It is reported for comparability with published numbers, not as the
    headline.

    Returns:
        ROC-AUC in ``[0, 1]``, or ``0.0`` when either class is absent (no
        discrimination is measurable).

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, scores)
    n_pos = sum(1 for v in labels if v == 1)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.0
    ranks = _ranks(scores)  # tie-averaged, ascending, already 1-based
    pos_rank_sum = math.fsum(r for r, v in zip(ranks, labels, strict=True) if v == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def mcc(labels: Sequence[int], predictions: Sequence[int]) -> float:
    """Return the Matthews correlation coefficient of a thresholded prediction.

    MCC is the honest single-number summary on an imbalanced task: it is high only
    when all four confusion-matrix cells are good, so the "call everything negative"
    degenerate scores 0 rather than 99%.

    Returns:
        MCC in ``[-1, 1]``, or ``0.0`` when a denominator vanishes (one predicted or
        one true class is empty, so no correlation is defined).

    Raises:
        ValueError: If the series differ in length, are empty, or contain non-0/1.
    """
    _check_binary(labels, [float(p) for p in predictions])
    bad = {int(v) for v in predictions} - {0, 1}
    if bad:
        raise ValueError(f"predictions must be 0/1, found {sorted(bad)}")

    tp = sum(1 for t, p in zip(labels, predictions, strict=True) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(labels, predictions, strict=True) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(labels, predictions, strict=True) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(labels, predictions, strict=True) if t == 1 and p == 0)
    denom = math.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    if denom == 0.0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def brier_score(labels: Sequence[int], probs: Sequence[float]) -> float:
    """Return the Brier score: mean squared error of probabilistic predictions.

    Unlike PR-AUC and ROC-AUC, this is **not** rank-based -- it penalises a model
    whose ordering is perfect but whose probabilities are systematically too
    confident. Lower is better; 0 is perfect.

    Because its scale depends on the base rate, a raw Brier score on an imbalanced
    task looks impressively small for a useless model. Read it through
    :func:`brier_skill_score`.

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, probs)
    return math.fsum((p - t) ** 2 for t, p in zip(labels, probs, strict=True)) / len(labels)


def brier_skill_score(labels: Sequence[int], probs: Sequence[float]) -> float:
    """Return the Brier skill score against the base-rate predictor.

    ``1 - BS_model / BS_baseline``, where the baseline predicts the observed
    prevalence at every position. This is the splice analogue of the expression
    gate's ``width_over_iqr``: the number that exposes a **vacuous pass**. On a task
    where 0.1% of positions are sites, a model emitting 0.001 everywhere achieves a
    tiny Brier score and a flattering ECE while being useless -- and scores exactly
    ``0.0`` here.

    Returns:
        The skill score. ``1.0`` is perfect, ``0.0`` is no better than predicting the
        base rate, and negative means **worse** than that baseline. Returns ``0.0``
        when the baseline itself is perfect (a single-class set), since no skill is
        demonstrable.

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, probs)
    prevalence = math.fsum(float(v) for v in labels) / len(labels)
    baseline = math.fsum((prevalence - t) ** 2 for t in labels) / len(labels)
    if baseline == 0.0:
        return 0.0
    return 1.0 - brier_score(labels, probs) / baseline


def reliability_bins(
    labels: Sequence[int], probs: Sequence[float], *, n_bins: int = 10
) -> list[tuple[float, float, int]]:
    """Return equal-width reliability bins as ``(mean_prob, observed_rate, count)``.

    The raw material of a reliability diagram, and of :func:`expected_calibration_error`.
    Empty bins are omitted rather than reported as ``(0, 0, 0)``, so a caller
    plotting the result does not draw a point where there is no evidence.

    Args:
        labels: Ground truth, 0 or 1.
        probs: Predicted probabilities in ``[0, 1]``.
        n_bins: Number of equal-width bins across ``[0, 1]``.

    Raises:
        ValueError: If the series differ in length, are empty, labels are not 0/1,
            or ``n_bins`` is not positive.
    """
    _check_binary(labels, probs)
    if n_bins <= 0:
        raise ValueError(f"n_bins must be positive, got {n_bins}")

    sums: list[float] = [0.0] * n_bins
    hits: list[int] = [0] * n_bins
    counts: list[int] = [0] * n_bins
    for t, p in zip(labels, probs, strict=True):
        clamped = min(max(p, 0.0), 1.0)
        idx = min(int(clamped * n_bins), n_bins - 1)
        sums[idx] += clamped
        hits[idx] += int(t)
        counts[idx] += 1
    return [
        (sums[i] / counts[i], hits[i] / counts[i], counts[i])
        for i in range(n_bins)
        if counts[i] > 0
    ]


def expected_calibration_error(
    labels: Sequence[int], probs: Sequence[float], *, n_bins: int = 10
) -> float:
    """Return the expected calibration error: count-weighted |confidence - accuracy|.

    **ECE has two well-known failure modes, and neither is fixed here -- they are
    reported around.** It is sensitive to ``n_bins`` (a coarse binning can hide
    miscalibration inside a wide bin), and it is *unable* to distinguish a
    well-calibrated informative model from a well-calibrated useless one: a
    predictor emitting the base rate everywhere has near-zero ECE. That is why the
    splice gate reports ECE **together with** :func:`brier_skill_score` and
    :func:`pr_auc`, and why a low ECE alone is never a pass.

    Returns:
        ECE in ``[0, 1]``; lower is better.

    Raises:
        ValueError: If the series differ in length, are empty, labels are not 0/1,
            or ``n_bins`` is not positive.
    """
    bins = reliability_bins(labels, probs, n_bins=n_bins)
    total = len(labels)
    return math.fsum(count * abs(mean_p - rate) for mean_p, rate, count in bins) / total


def top_k_accuracy(labels: Sequence[int], scores: Sequence[float]) -> float:
    """Return top-k accuracy, where ``k`` is the number of true positives.

    Implements the **pooled** construction used by Zeng & Li (2022) and the shared
    SpliceAI/Pangolin evaluation code: rank all candidate positions, take the top
    ``k`` where ``k`` is the number of labelled positives *in the set as passed*, and
    report the fraction that are real. Because ``k`` equals the positive count,
    precision and recall coincide, which is what makes it a single scalar.

    **There is no single "top-k accuracy".** OpenSpliceAI (*eLife* 2025) computes
    ``k`` per gene *and* per class then averages, which disagrees with this on the
    same panel. Cite whichever construction is used; never "the" definition. This
    function is **class-agnostic and does no pooling of its own** -- the caller owns
    what goes into one call, and pooling across genes versus per-gene averaging is
    that caller's decision, not this function's.

    Comparison anchors that use *this* construction: Pangolin top-1 79%, SpliceAI
    75% on their shared benchmark (Zeng & Li 2022); SpliceAI's own headline ~0.95 is
    on its GENCODE test set.

    Ties at the cutoff are resolved **pessimistically**: every case sharing the
    boundary score is considered, and the credit is the expected fraction under a
    random tie-break rather than the best case. A model scoring everything equally
    therefore gets the base rate, not 1.0. This is a deliberate deviation from the
    published implementations' ``np.argsort``, which breaks ties by array order and
    can be flattered by it -- so BT4's number may sit slightly *below* a published
    one wherever ties exist.

    Returns:
        Top-k accuracy in ``[0, 1]``, or ``0.0`` when there are no positives.

    Raises:
        ValueError: If the series differ in length, are empty, or labels are not 0/1.
    """
    _check_binary(labels, scores)
    k = sum(1 for v in labels if v == 1)
    if k == 0:
        return 0.0

    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    cutoff = scores[order[k - 1]]
    strictly_above = [i for i in order if scores[i] > cutoff]
    at_cutoff = [i for i in order if scores[i] == cutoff]

    hits = float(sum(1 for i in strictly_above if labels[i] == 1))
    remaining = k - len(strictly_above)
    if remaining > 0 and at_cutoff:
        tied_pos = sum(1 for i in at_cutoff if labels[i] == 1)
        hits += remaining * (tied_pos / len(at_cutoff))
    return hits / k
