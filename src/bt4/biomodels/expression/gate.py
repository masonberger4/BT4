"""The acceptance gate a learned expression head must pass to claim calibration.

CLAUDE.md §6/§8/§10.6 make ``calibrated`` an *earned* property: an
:class:`~bt4.biomodels.expression.base.ExpressionPredictor` may report
``calibrated is True`` only after a held-out, grouped acceptance gate on data
**from the regime it serves**. This module is that gate -- model-agnostic and
dependency-free, so it can be written and tested now while the placeholder is
still the default, exactly as the splice fidelity/attestation machinery shipped
before any calibrated splice backend.

**What it measures (Q5 spec).** For a coding-sequence expression head predicting
**log-TE** (log translation efficiency), a regression target, the gate reports:

* **Spearman rank correlation** (primary) -- does the head rank sequences by
  expression the way the measured data does?
* **Pearson** and **R^2** (secondary) -- linear agreement and variance explained.
* **Split-conformal coverage** at a target level (default 90%) -- does a conformal
  interval calibrated on held-out data actually contain the truth at the promised
  rate? This is the honest uncertainty check (§8): a point estimate is never
  presented as ground truth.

**Grouped split (no leakage).** Calibration and test are split by **group**
(a homology cluster for designed variants, a chromosome for natural genes), so no
group's members straddle the split. This is the distribution-shift-aware
evaluation the maintainer's framing demands: a head validated only on natural-gene
TE has *not* earned calibration for the CDS-variant regime BT4 optimizes, and a
grouped split is what surfaces that gap rather than hiding it behind a random
split's optimistic leakage.

**This gate does not flip anything.** It computes an honest
:class:`ExpressionGateReport`; a maintainer promotes a head to ``calibrated=True``
only after it passes on real in-regime data (and never by assignment). The
thresholds (``min_spearman``, ``coverage_tolerance``) are inputs, set by the
maintainer at gate time -- this module does not bless a default that would let a
weak head self-certify.

Pure standard library; depends only on :mod:`bt4.biomodels._stats` and the
:class:`ExpressionPredictor` contract.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels._stats import (
    conformal_quantile,
    empirical_coverage,
    pearson,
    r2_score,
    spearman,
)
from bt4.biomodels.expression.base import ExpressionPredictor

__all__ = [
    "ExpressionEvalCase",
    "ExpressionGateReport",
    "run_expression_gate",
    "verify_expression_gate",
]


@dataclass(frozen=True, slots=True)
class ExpressionEvalCase:
    """One held-out evaluation point for the expression acceptance gate.

    Attributes:
        predicted: The head's predicted score for this sequence (log-TE).
        measured: The measured log-TE from the reference dataset.
        group: The leakage-control group id (homology cluster / chromosome). Cases
            sharing a group are never split across calibration and test.
    """

    predicted: float
    measured: float
    group: str


@dataclass(frozen=True, slots=True)
class ExpressionGateReport:
    """The honest result of an expression acceptance gate (nothing is flipped).

    Attributes:
        passed: ``True`` iff ``spearman >= min_spearman`` **and** the empirical
            conformal coverage is within ``coverage_tolerance`` of the target.
        n_calibration: Number of calibration cases (used for the conformal quantile).
        n_test: Number of test cases (metrics are computed on these).
        n_groups: Number of distinct leakage-control groups across all cases.
        spearman: Spearman rank correlation on the test fold (primary metric).
        pearson: Pearson correlation on the test fold.
        r2: Coefficient of determination on the test fold (may be negative).
        target_coverage: The requested conformal coverage level (1 - alpha).
        empirical_coverage: The realized coverage on the test fold.
        coverage_tolerance: Allowed absolute gap between target and empirical
            coverage for a pass.
        conformal_half_width: The split-conformal interval half-width from the
            calibration residuals (``math.inf`` if calibration was too small).
        min_spearman: The Spearman threshold required to pass.
    """

    passed: bool
    n_calibration: int
    n_test: int
    n_groups: int
    spearman: float
    pearson: float
    r2: float
    target_coverage: float
    empirical_coverage: float
    coverage_tolerance: float
    conformal_half_width: float
    min_spearman: float


def _grouped_split(
    cases: Sequence[ExpressionEvalCase], calibration_fraction: float
) -> tuple[list[ExpressionEvalCase], list[ExpressionEvalCase]]:
    """Split ``cases`` into (calibration, test) with **disjoint groups**.

    Groups are ordered deterministically (sorted by id); the first
    ``ceil(calibration_fraction * n_groups)`` groups form the calibration fold and
    the rest form the test fold, so no group leaks across the split and the split
    is reproducible (invariant #7 in spirit -- no RNG).
    """
    groups = sorted({c.group for c in cases})
    n_cal_groups = math.ceil(calibration_fraction * len(groups))
    # Keep at least one group on each side when there are >= 2 groups.
    n_cal_groups = max(1, min(n_cal_groups, len(groups) - 1))
    cal_groups = set(groups[:n_cal_groups])
    calibration = [c for c in cases if c.group in cal_groups]
    test = [c for c in cases if c.group not in cal_groups]
    return calibration, test


def verify_expression_gate(
    cases: Sequence[ExpressionEvalCase],
    *,
    target_coverage: float = 0.90,
    coverage_tolerance: float = 0.05,
    min_spearman: float = 0.3,
    calibration_fraction: float = 0.5,
) -> ExpressionGateReport:
    """Run the expression acceptance gate over held-out ``cases``.

    Splits the cases into calibration and test folds with **disjoint groups**,
    fits a split-conformal interval on the calibration residuals, and reports the
    Spearman/Pearson/R^2 point metrics and the realized conformal coverage on the
    test fold. ``passed`` requires the primary Spearman threshold *and* conformal
    coverage close to the target -- both point accuracy and honest uncertainty.

    Args:
        cases: Held-out ``(predicted, measured, group)`` evaluation points.
        target_coverage: Conformal coverage level (1 - alpha) in ``(0, 1)``.
        coverage_tolerance: Max absolute gap between target and empirical coverage
            for a pass.
        min_spearman: Minimum test Spearman rank correlation for a pass.
        calibration_fraction: Fraction of *groups* assigned to the calibration
            fold (the rest are the test fold).

    Returns:
        An :class:`ExpressionGateReport` (this function never flips ``calibrated``).

    Raises:
        ValueError: If there are fewer than two groups (a grouped split is
            impossible), if either fold ends up empty, if the test fold has fewer
            than two points (correlation is undefined), or if the level/fraction
            arguments are out of range.
    """
    if not 0.0 < target_coverage < 1.0:
        raise ValueError(f"target_coverage must be in (0, 1), got {target_coverage}")
    if not 0.0 < calibration_fraction < 1.0:
        raise ValueError(
            f"calibration_fraction must be in (0, 1), got {calibration_fraction}"
        )
    n_groups = len({c.group for c in cases})
    if n_groups < 2:
        raise ValueError(
            "expression gate needs at least two distinct groups for a leakage-free "
            f"split, got {n_groups}"
        )

    calibration, test = _grouped_split(cases, calibration_fraction)
    if not calibration or not test:
        raise ValueError("grouped split produced an empty calibration or test fold")
    if len(test) < 2:
        raise ValueError(
            f"test fold has {len(test)} case(s); need >= 2 for correlation metrics"
        )

    cal_residuals = [abs(c.predicted - c.measured) for c in calibration]
    q_hat = conformal_quantile(cal_residuals, target_coverage)

    test_pred = [c.predicted for c in test]
    test_meas = [c.measured for c in test]
    test_residuals = [abs(c.predicted - c.measured) for c in test]

    rho = spearman(test_meas, test_pred)
    r = pearson(test_meas, test_pred)
    r2 = r2_score(test_meas, test_pred)
    coverage = empirical_coverage(test_residuals, q_hat)

    passed = (
        rho >= min_spearman
        and abs(coverage - target_coverage) <= coverage_tolerance
    )
    return ExpressionGateReport(
        passed=passed,
        n_calibration=len(calibration),
        n_test=len(test),
        n_groups=n_groups,
        spearman=rho,
        pearson=r,
        r2=r2,
        target_coverage=target_coverage,
        empirical_coverage=coverage,
        coverage_tolerance=coverage_tolerance,
        conformal_half_width=q_hat,
        min_spearman=min_spearman,
    )


def run_expression_gate(
    predictor: ExpressionPredictor,
    samples: Sequence[tuple[str, float, str]],
    **kwargs: float,
) -> ExpressionGateReport:
    """Score ``samples`` with ``predictor`` and run :func:`verify_expression_gate`.

    A convenience wrapper: each sample is ``(dna, measured_log_te, group)``. The
    predictor's :meth:`~ExpressionPredictor.score_sequence` supplies the
    prediction, then the model-agnostic gate does the rest.

    Args:
        predictor: The expression head to evaluate (calibrated or not -- the gate
            reports honestly either way).
        samples: ``(dna, measured, group)`` triples from the deployment regime.
        **kwargs: Passed through to :func:`verify_expression_gate`
            (``target_coverage`` / ``coverage_tolerance`` / ``min_spearman`` /
            ``calibration_fraction``).

    Returns:
        An :class:`ExpressionGateReport`.
    """
    cases = [
        ExpressionEvalCase(
            predicted=predictor.score_sequence(dna).score,
            measured=measured,
            group=group,
        )
        for dna, measured, group in samples
    ]
    return verify_expression_gate(cases, **kwargs)
