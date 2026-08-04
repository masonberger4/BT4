"""Tests for the shared dependency-free statistics (:mod:`bt4.biomodels._stats`)."""

from __future__ import annotations

import math

import pytest

from bt4.biomodels._stats import (
    conformal_quantile,
    empirical_coverage,
    pearson,
    r2_score,
    spearman,
)


def test_pearson_perfect_and_anti() -> None:
    assert pearson([1, 2, 3, 4], [2, 4, 6, 8]) == pytest.approx(1.0)
    assert pearson([1, 2, 3, 4], [8, 6, 4, 2]) == pytest.approx(-1.0)


def test_pearson_zero_variance_is_zero_not_nan() -> None:
    assert pearson([1, 1, 1, 1], [1, 2, 3, 4]) == 0.0


def test_spearman_is_rank_monotone() -> None:
    # A monotone but non-linear relation: Spearman is 1.0, Pearson is not.
    x = [1, 2, 3, 4, 5]
    y = [1, 4, 9, 16, 25]
    assert spearman(x, y) == pytest.approx(1.0)
    assert pearson(x, y) < 1.0


def test_r2_perfect_and_mean_predictor_and_negative() -> None:
    y = [1.0, 2.0, 3.0, 4.0]
    assert r2_score(y, y) == pytest.approx(1.0)
    mean = sum(y) / len(y)
    assert r2_score(y, [mean] * len(y)) == pytest.approx(0.0)  # predicting the mean => 0
    # A prediction worse than the mean yields a negative R^2 (reported, not clamped).
    assert r2_score(y, [4.0, 3.0, 2.0, 1.0]) < 0.0


def test_r2_zero_variance_truth_is_zero() -> None:
    assert r2_score([2.0, 2.0, 2.0], [1.0, 2.0, 3.0]) == 0.0


def test_correlation_length_and_size_guards() -> None:
    with pytest.raises(ValueError):
        pearson([1.0], [1.0])  # fewer than 2
    with pytest.raises(ValueError):
        pearson([1.0, 2.0], [1.0])  # length mismatch
    with pytest.raises(ValueError):
        spearman([1.0], [1.0])  # spearman inherits the guard via pearson


def test_conformal_quantile_selects_kth_smallest() -> None:
    # n=9 scores, coverage 0.9 => k = ceil(10*0.9) = 9 => the 9th smallest (max).
    scores = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    assert conformal_quantile(scores, 0.9) == pytest.approx(0.9)
    # coverage 0.5 => k = ceil(10*0.5) = 5 => the 5th smallest.
    assert conformal_quantile(scores, 0.5) == pytest.approx(0.5)


def test_conformal_quantile_infinite_when_too_few_points() -> None:
    # n=5, coverage 0.9 => k = ceil(6*0.9) = 6 > 5 => +inf (must cover everything).
    assert conformal_quantile([0.1, 0.2, 0.3, 0.4, 0.5], 0.9) == math.inf


def test_conformal_quantile_argument_guards() -> None:
    with pytest.raises(ValueError):
        conformal_quantile([0.1], 0.0)
    with pytest.raises(ValueError):
        conformal_quantile([0.1], 1.0)
    with pytest.raises(ValueError):
        conformal_quantile([], 0.9)


def test_empirical_coverage() -> None:
    residuals = [0.1, 0.2, 0.9, 1.5]
    assert empirical_coverage(residuals, 0.5) == pytest.approx(0.5)  # 2 of 4 within 0.5
    assert empirical_coverage(residuals, math.inf) == 1.0
    with pytest.raises(ValueError):
        empirical_coverage([], 0.5)
