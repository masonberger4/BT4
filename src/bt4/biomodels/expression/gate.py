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

**Two modes, because the default one can hand out a false pass.**

*Pooled* (the default, ``within_group=False``) computes one Spearman over the whole
test fold. When the groups are **proteins**, that fold mixes several proteins whose
baseline expression differs wildly -- so a head that knows nothing about codons but
recognises "this is a highly-expressed gene" scores well. That is exactly what a head
trained across natural genes learns, and exactly the skill BT4 cannot use, since BT4
chooses between synonymous variants *of one protein*.

*Within-group* (``within_group=True``) centres predictions and measurements inside
each group before scoring, and aggregates a per-group Spearman across groups. Only
"did you order this protein's own variants correctly?" can contribute. This is the
regime BT4 actually deploys in, and it is the strict bar (§10.6).

**The link (``recalibrate=True``), because raw residuals across different units are
meaningless.** A head whose output is in arbitrary units (RiboNN reports a CLR
compositional residual) cannot be compared to an assay's units by subtraction. Split
conformal stays *valid* under any units -- its guarantee holds for any score function
-- but the interval becomes uselessly wide, and a **constant predictor achieves exactly
valid coverage**. So the gate (a) optionally fits the affine link
``measured ~= slope * predicted + intercept`` on the *calibration* fold only, never on
the fold it then conformalizes, and (b) always reports **median interval width divided
by the label IQR**, which is the number that exposes a vacuous pass. Reporting a
"calibration slope of 1" after fitting the slope would be circular, so what is reported
instead is the *spread of the fitted link across calibration groups*.

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
import random
from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels._stats import (
    conformal_quantile,
    empirical_coverage,
    iqr,
    linear_fit,
    pearson,
    r2_score,
    spearman,
)
from bt4.biomodels.expression.base import (
    BatchExpressionPredictor,
    ExpressionPredictor,
)

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
        spearman: Primary rank metric on the test fold, computed on the head's
            **raw** (un-linked) predictions. Rank correlation needs no link, and this
            must describe what BT4 would really do -- it ranks candidates by the raw
            score and never applies a fitted link at design time. So a head that ranks
            backwards shows a **negative** value here and cannot be rescued by a
            negative fitted slope.
        pearson: Pearson correlation on the test fold, after the link (if fitted).
        r2: Coefficient of determination on the test fold, after the link (if fitted);
            may be negative. A bare R^2 on arbitrary-unit output is meaningless, which
            is why the link exists.
        target_coverage: The requested conformal coverage level (1 - alpha).
        empirical_coverage: The realized coverage on the test fold.
        coverage_tolerance: Allowed absolute gap between target and empirical
            coverage for a pass.
        conformal_half_width: The split-conformal interval half-width from the
            calibration residuals (``math.inf`` if calibration was too small).
        min_spearman: The Spearman threshold required to pass.
        within_group: Whether metrics were computed inside groups (the strict bar)
            rather than pooled across the test fold.
        recalibrate: Whether the affine link was fitted on the calibration fold.
        slope: Fitted link slope (``1.0`` when ``recalibrate`` is off). A **negative**
            slope means the head ranks backwards on this panel -- visible here rather
            than silently absorbed.
        intercept: Fitted link intercept (``0.0`` when ``recalibrate`` is off).
        link_slope_spread: Standard deviation of the link slope refitted within each
            calibration group. Reported *instead of* a calibration slope, which is
            ``1.0`` by construction once the link is fitted; a large spread means the
            link does not transfer between groups and the interval should not be
            trusted.
        per_group_spearman: ``(group, spearman)`` for every test group with at least
            two members. Empty in pooled mode.
        n_groups_test: Number of distinct groups in the test fold.
        n_groups_ranked: Test groups that actually contributed a rank correlation
            (those with >= 2 members). The effective sample size for a cross-group
            claim is this, **not** ``n_test``.
        width_over_iqr: Median conformal interval width divided by the interquartile
            range of the test labels. The vacuity check: valid coverage with a ratio
            >= 1 means the interval is as wide as the data and says nothing.
        spearman_ci_low: Lower bound of the cluster-bootstrap CI on ``spearman``
            (whole groups resampled, so the dependence between one protein's variants
            is respected).
        spearman_ci_high: Upper bound of the same CI.
        bootstrap_resamples: Resamples that yielded a defined statistic. ``0`` means no
            CI was computed and the bounds are ``nan``.
        coverage_conditional_on_group_anchor: ``True`` in within-group mode, where the
            target is each variant's offset from its group's mean. Such an interval is
            only achievable at design time if a member of that protein has already been
            measured to anchor it -- which is BT4's ``delta_logte`` framing, but is a
            **narrower claim** than an unconditional interval and must be reported as
            one.
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
    within_group: bool = False
    recalibrate: bool = False
    slope: float = 1.0
    intercept: float = 0.0
    link_slope_spread: float = 0.0
    per_group_spearman: tuple[tuple[str, float], ...] = ()
    n_groups_test: int = 0
    n_groups_ranked: int = 0
    width_over_iqr: float = math.inf
    spearman_ci_low: float = math.nan
    spearman_ci_high: float = math.nan
    bootstrap_resamples: int = 0
    coverage_conditional_on_group_anchor: bool = False


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


def _centre_within_groups(
    cases: Sequence[ExpressionEvalCase],
) -> list[ExpressionEvalCase]:
    """Subtract each group's own mean from both predictions and measurements.

    This is what turns the gate into a test of BT4's actual question. After centring, a
    case's values are *offsets from its own protein's baseline*, so between-protein
    differences -- the thing a natural-gene-trained head is good at and BT4 cannot use
    -- contribute exactly nothing.

    A singleton group centres to ``(0.0, 0.0)``; such a case carries no within-group
    information and is dropped from the rank aggregation by
    :func:`_per_group_spearman`, though it still counts in the fold sizes.
    """
    sums: dict[str, tuple[float, float, int]] = {}
    for case in cases:
        total_p, total_m, count = sums.get(case.group, (0.0, 0.0, 0))
        sums[case.group] = (total_p + case.predicted, total_m + case.measured, count + 1)
    means = {
        group: (total_p / count, total_m / count)
        for group, (total_p, total_m, count) in sums.items()
    }
    return [
        ExpressionEvalCase(
            predicted=case.predicted - means[case.group][0],
            measured=case.measured - means[case.group][1],
            group=case.group,
        )
        for case in cases
    ]


def _by_group(
    cases: Sequence[ExpressionEvalCase],
) -> dict[str, list[ExpressionEvalCase]]:
    """Bucket cases by group id, preserving input order within each bucket."""
    buckets: dict[str, list[ExpressionEvalCase]] = {}
    for case in cases:
        buckets.setdefault(case.group, []).append(case)
    return buckets


def _per_group_spearman(
    cases: Sequence[ExpressionEvalCase],
) -> list[tuple[str, float]]:
    """Return ``(group, spearman)`` for every group with at least two members.

    Groups are visited in sorted order so the result is deterministic (invariant #7).
    A group with one member has no ordering to get right and is omitted rather than
    contributing a fabricated ``0.0``.
    """
    out: list[tuple[str, float]] = []
    for group, members in sorted(_by_group(cases).items()):
        if len(members) < 2:
            continue
        measured = [c.measured for c in members]
        predicted = [c.predicted for c in members]
        out.append((group, spearman(measured, predicted)))
    return out


def _primary_spearman(cases: Sequence[ExpressionEvalCase], within_group: bool) -> float:
    """Return the primary rank metric for the active mode.

    Pooled mode: one Spearman over the fold. Within-group mode: the **unweighted mean**
    of the per-group Spearmans, so a protein with 30 variants does not outvote one with
    4 -- each protein is one independent observation of "can it rank variants?"
    (the aggregation ProteinGym uses).

    Raises:
        ValueError: In within-group mode when no group has two or more members, so no
            within-group ordering exists to score.
    """
    if not within_group:
        return spearman([c.measured for c in cases], [c.predicted for c in cases])
    per_group = _per_group_spearman(cases)
    if not per_group:
        raise ValueError(
            "within-group mode needs at least one test group with >= 2 members; "
            "every test group is a singleton, so there is no ordering to score"
        )
    return math.fsum(rho for _group, rho in per_group) / len(per_group)


def _cluster_bootstrap_ci(
    cases: Sequence[ExpressionEvalCase],
    within_group: bool,
    *,
    resamples: int,
    seed: int,
    level: float = 0.95,
) -> tuple[float, float, int]:
    """Percentile CI on the primary rank metric, resampling **whole groups**.

    Variants of one protein are a dependent cluster, so resampling individual cases
    would treat 30 variants of one protein as 30 independent observations and produce a
    CI far too narrow. Resampling groups with replacement respects the dependence, and
    the effective sample size is the number of groups.

    Deterministic from ``seed`` (invariant #7). Resamples whose statistic is undefined
    (e.g. every drawn group a singleton in within-group mode) are skipped rather than
    counted as zero.

    Returns:
        ``(low, high, n_valid)``. With fewer than 20 valid resamples the bounds are
        ``nan`` and ``n_valid`` is ``0`` -- an honest "not estimated" rather than a
        CI computed from noise.
    """
    if resamples <= 0:
        return math.nan, math.nan, 0
    buckets = _by_group(cases)
    group_ids = sorted(buckets)
    rng = random.Random(seed)
    statistics_: list[float] = []
    for _ in range(resamples):
        drawn = [rng.choice(group_ids) for _ in group_ids]
        resampled = [case for group in drawn for case in buckets[group]]
        try:
            statistics_.append(_primary_spearman(resampled, within_group))
        except ValueError:
            continue
    if len(statistics_) < 20:
        return math.nan, math.nan, 0
    statistics_.sort()
    alpha = (1.0 - level) / 2.0
    low_index = max(0, math.floor(alpha * (len(statistics_) - 1)))
    high_index = min(len(statistics_) - 1, math.ceil((1.0 - alpha) * (len(statistics_) - 1)))
    return statistics_[low_index], statistics_[high_index], len(statistics_)


def _link_slope_spread(cases: Sequence[ExpressionEvalCase]) -> float:
    """Standard deviation of the affine link refitted inside each calibration group.

    Reported instead of a "calibration slope", which is 1.0 by construction once the
    link has been fitted. A link that differs sharply between proteins does not
    transfer, so the interval built from it should not be trusted even if its coverage
    on this panel looks fine.
    """
    slopes = [
        linear_fit([c.predicted for c in members], [c.measured for c in members])[0]
        for _group, members in sorted(_by_group(cases).items())
        if len(members) >= 2
    ]
    if len(slopes) < 2:
        return 0.0
    mean = math.fsum(slopes) / len(slopes)
    variance = math.fsum((s - mean) ** 2 for s in slopes) / (len(slopes) - 1)
    return math.sqrt(variance)


def verify_expression_gate(
    cases: Sequence[ExpressionEvalCase],
    *,
    target_coverage: float = 0.90,
    coverage_tolerance: float = 0.05,
    min_spearman: float = 0.3,
    calibration_fraction: float = 0.5,
    within_group: bool = False,
    recalibrate: bool = False,
    bootstrap_resamples: int = 1000,
    bootstrap_seed: int = 0,
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
        within_group: Score inside each group instead of pooling across the test fold
            -- the strict bar. Between-group differences then cannot contribute, so a
            head that only recognises which gene it is looking at scores ~0.
        recalibrate: Fit ``measured ~= slope * predicted + intercept`` on the
            calibration fold and apply it before residuals. Required whenever the
            head's units differ from the assay's, which is the normal case.
        bootstrap_resamples: Cluster-bootstrap resamples (whole groups) for the CI on
            the primary metric. ``0`` disables it.
        bootstrap_seed: Seed for that bootstrap; identical inputs give identical
            bounds (invariant #7).

    Returns:
        An :class:`ExpressionGateReport` (this function never flips ``calibrated``).

    Raises:
        ValueError: If there are fewer than two groups (a grouped split is
            impossible), if either fold ends up empty, if the test fold has fewer
            than two points (correlation is undefined), if ``within_group`` is set but
            no test group has two or more members, or if the level/fraction arguments
            are out of range.
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

    # Centre first, then fit the link. In within-group mode the quantity of interest is
    # a variant's offset from its own protein's baseline, so the link that matters
    # relates *offsets* to *offsets*, not absolute levels.
    if within_group:
        calibration = _centre_within_groups(calibration)
        test = _centre_within_groups(test)

    # The link is fitted on the calibration fold ONLY. Fitting it on the fold that is
    # then conformalized would break the independence split conformal requires.
    slope, intercept = (1.0, 0.0)
    if recalibrate:
        slope, intercept = linear_fit(
            [c.predicted for c in calibration], [c.measured for c in calibration]
        )
    slope_spread = _link_slope_spread(calibration) if recalibrate else 0.0

    def _linked(case: ExpressionEvalCase) -> ExpressionEvalCase:
        return ExpressionEvalCase(
            predicted=slope * case.predicted + intercept,
            measured=case.measured,
            group=case.group,
        )

    # The RANK metric is computed on the *un-linked* predictions, deliberately.
    # Rank correlation needs no link (it is invariant to any strictly increasing map),
    # and more importantly it must describe what BT4 would actually do: BT4 ranks
    # candidates by the head's raw score and never applies a fitted link at design
    # time. Scoring the linked predictions instead would let a head that ranks
    # *backwards* be rescued by a negative fitted slope and reported as passing, when a
    # deployed BT4 would hand the user the worst candidate. A negative ``slope`` in the
    # report is the visible signal for that case.
    #
    # The link applies to everything that lives on the measurement scale -- Pearson,
    # R^2, the conformal residuals and the interval width -- because those are
    # meaningless when the head's units are arbitrary.
    cal_linked = [_linked(c) for c in calibration]
    test_linked = [_linked(c) for c in test]

    cal_residuals = [abs(c.predicted - c.measured) for c in cal_linked]
    q_hat = conformal_quantile(cal_residuals, target_coverage)

    test_meas = [c.measured for c in test_linked]
    test_linked_pred = [c.predicted for c in test_linked]
    test_residuals = [abs(c.predicted - c.measured) for c in test_linked]

    rho = _primary_spearman(test, within_group)
    per_group = tuple(_per_group_spearman(test)) if within_group else ()
    r = pearson(test_meas, test_linked_pred)
    r2 = r2_score(test_meas, test_linked_pred)
    coverage = empirical_coverage(test_residuals, q_hat)

    # Coverage alone is passable by a constant predictor, so always report the interval
    # width against the spread of the labels it is supposed to be informative about.
    label_iqr = iqr(test_meas)
    width_over_iqr = (2.0 * q_hat / label_iqr) if label_iqr > 0.0 else math.inf

    ci_low, ci_high, n_resamples = _cluster_bootstrap_ci(
        test, within_group, resamples=bootstrap_resamples, seed=bootstrap_seed
    )

    passed = (
        rho >= min_spearman
        and abs(coverage - target_coverage) <= coverage_tolerance
    )
    return ExpressionGateReport(
        passed=passed,
        n_calibration=len(cal_linked),
        n_test=len(test_linked),
        n_groups=n_groups,
        spearman=rho,
        pearson=r,
        r2=r2,
        target_coverage=target_coverage,
        empirical_coverage=coverage,
        coverage_tolerance=coverage_tolerance,
        conformal_half_width=q_hat,
        min_spearman=min_spearman,
        within_group=within_group,
        recalibrate=recalibrate,
        slope=slope,
        intercept=intercept,
        link_slope_spread=slope_spread,
        per_group_spearman=per_group,
        n_groups_test=len({c.group for c in test}),
        n_groups_ranked=len(_per_group_spearman(test)),
        width_over_iqr=width_over_iqr,
        spearman_ci_low=ci_low,
        spearman_ci_high=ci_high,
        bootstrap_resamples=n_resamples,
        coverage_conditional_on_group_anchor=within_group,
    )


def run_expression_gate(
    predictor: ExpressionPredictor,
    samples: Sequence[tuple[str, float, str]],
    *,
    target_coverage: float = 0.90,
    coverage_tolerance: float = 0.05,
    min_spearman: float = 0.3,
    calibration_fraction: float = 0.5,
    within_group: bool = False,
    recalibrate: bool = False,
    bootstrap_resamples: int = 1000,
    bootstrap_seed: int = 0,
) -> ExpressionGateReport:
    """Score ``samples`` with ``predictor`` and run :func:`verify_expression_gate`.

    A convenience wrapper: each sample is ``(dna, measured_log_te, group)``.

    **Scored in one batched call where the backend allows it.** A backend implementing
    :class:`~bt4.biomodels.expression.base.BatchExpressionPredictor` (RiboNN does) has
    a large fixed *per-invocation* cost -- for RiboNN, hashing 90 weight files, loading
    50 models, and a full forward pass of the whole set per model -- so scoring a panel
    one row at a time would multiply the wall clock by the number of rows and re-verify
    the weights that many times. Backends without the batch surface fall back to
    per-sequence scoring.

    Args:
        predictor: The expression head to evaluate (calibrated or not -- the gate
            reports honestly either way).
        samples: ``(dna, measured, group)`` triples from the deployment regime.
        target_coverage: See :func:`verify_expression_gate`.
        coverage_tolerance: See :func:`verify_expression_gate`.
        min_spearman: See :func:`verify_expression_gate`.
        calibration_fraction: See :func:`verify_expression_gate`.
        within_group: See :func:`verify_expression_gate` -- the strict bar.
        recalibrate: See :func:`verify_expression_gate`.
        bootstrap_resamples: See :func:`verify_expression_gate`.
        bootstrap_seed: See :func:`verify_expression_gate`.

    Returns:
        An :class:`ExpressionGateReport`.
    """
    dnas = [dna for dna, _measured, _group in samples]
    if isinstance(predictor, BatchExpressionPredictor):
        scores = [result.score for result in predictor.score_many(dnas)]
    else:
        scores = [predictor.score_sequence(dna).score for dna in dnas]
    cases = [
        ExpressionEvalCase(predicted=score, measured=measured, group=group)
        for score, (_dna, measured, group) in zip(scores, samples, strict=True)
    ]
    return verify_expression_gate(
        cases,
        target_coverage=target_coverage,
        coverage_tolerance=coverage_tolerance,
        min_spearman=min_spearman,
        calibration_fraction=calibration_fraction,
        within_group=within_group,
        recalibrate=recalibrate,
        bootstrap_resamples=bootstrap_resamples,
        bootstrap_seed=bootstrap_seed,
    )
