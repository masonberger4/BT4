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
    "conformal_quantile",
    "empirical_coverage",
    "iqr",
    "linear_fit",
    "pearson",
    "r2_score",
    "spearman",
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
