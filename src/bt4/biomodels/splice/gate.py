"""The acceptance gate a splice backend must pass to claim *statistical* calibration.

A wrapped CNN earns ``calibrated=True`` today by passing an **integration-fidelity**
gate: proof that BT4's adapter reproduces the published model bit-for-bit
(:mod:`bt4.biomodels.splice.attestation`). Pangolin has passed that. It says nothing
about whether the model's numbers *mean* what they claim -- whether a score of 0.5
corresponds to a 50% chance of splicing, and whether the operating point BT4
thresholds at is the right one. This module is the separate gate for that question
(``docs/DESIGN_splice_cnn_calibration.md`` Part B).

**Two tasks, two case types, never mixed.** BT4 exposes two different quantities to
a user decision, and they need different ground truth:

* :class:`SpliceSiteCase` -- *site prediction*. Is position ``i`` a splice site?
  Ground truth is annotated donor/acceptor positions. This is what
  ``DEFAULT_SITE_PROBABILITY`` thresholds, and what the localization count in BT4
  Studio is built from.
* :class:`SpliceVariantCase` -- *variant effect*. Does this sequence change disrupt
  splicing? Ground truth is a measured assay. This is what ``delta_splicing``
  answers.

Keeping them apart is not fussiness. The **exonic / intronic** stratification that
matters most for BT4 (median prAUC 0.419 exonic vs 0.773 intronic, Smith & Kitzman
2023) is a *variant-effect* result. Applying it to site prediction would be
near-degenerate, because annotated sites sit at exon/intron boundaries by
construction. So ``region`` lives only on the variant case.

**What this gate deliberately does not use.**

*Spearman*, the expression gate's primary metric, is excluded. On a binary label it
is an exact affine function of ROC-AUC -- it carries no information ROC-AUC lacks --
and at splice prevalence it is unusable as a threshold: measured, a **perfect**
classifier scores Spearman 0.055 at 0.1% prevalence, far under the expression gate's
``min_spearman=0.30``. Transplanting that bar would fail a flawless model.

*Bare accuracy*, per CLAUDE.md section 6, because calling "not a site" everywhere
scores >99% on this task.

**What it reports, and why each is there.**

* **PR-AUC** is the pass axis, because it is the metric the published anchors use
  (Pangolin AUPRC 0.85, SpliceAI 0.77) and because it is sensitive to the false
  positives that matter at this prevalence. But raw PR-AUC is **not comparable
  across panels**: its floor is the prevalence, which is a *construction choice*.
  So the report carries ``pr_auc_skill`` (0 at no-skill, 1 at perfect, for any
  prevalence) and records the panel's negative construction, without which a
  threshold is passable by thinning negatives.
* **ROC-AUC** is reported because it is stable under negative-sampling where PR-AUC
  is not -- not as a skill claim, and never as a pass axis: at splice prevalence a
  model with hopeless precision still scores above 0.9.
* **top-k accuracy** for comparability with the published anchors, in the pooled
  Zeng & Li construction (see :func:`~bt4.biomodels._stats.top_k_accuracy`).
* **Brier + ECE + reliability** for calibration proper, and **Brier skill** because
  ECE alone cannot tell an informative model from a vacuous one: a predictor
  emitting the base rate everywhere is perfectly calibrated and useless.
* **Per-stratum results**, never one blended number. A pooled figure lets a backend
  certify on intronic strength while being weak exactly where BT4 operates.

**This gate flips nothing.** It returns an honest report; a maintainer promotes a
backend only on the evidence, and the thresholds are inputs set at gate time so this
module never blesses a default a weak backend could clear.

Pure standard library; depends only on :mod:`bt4.biomodels._stats`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from bt4.biomodels._stats import (
    brier_score,
    brier_skill_score,
    expected_calibration_error,
    mcc,
    pr_auc,
    pr_auc_skill,
    reliability_bins,
    roc_auc,
    top_k_accuracy,
)
from bt4.biomodels.splice.base import DEFAULT_SITE_PROBABILITY

__all__ = [
    "SITE_PREDICTION",
    "VARIANT_EFFECT",
    "SpliceGateReport",
    "SpliceSiteCase",
    "SpliceStratumReport",
    "SpliceVariantCase",
    "verify_splice_gate",
]

SITE_PREDICTION = "site_prediction"
"""Task id: per-position "is this a splice site?" against annotated positions."""

VARIANT_EFFECT = "variant_effect"
"""Task id: "does this change disrupt splicing?" against a measured assay."""


@dataclass(frozen=True, slots=True)
class SpliceSiteCase:
    """One scored position for the **site-prediction** task.

    Attributes:
        predicted: The backend's per-position score at this position.
        label: ``1`` if the position is an annotated splice site, else ``0``.
        kind: ``"donor"`` or ``"acceptor"`` -- the stratum for this task. (Pangolin
            emits one combined track and cannot distinguish them; such a panel uses
            a single ``"splice"`` kind rather than pretending to a split.)
        group: The leakage-control group, a **chromosome**. Cases sharing a group
            are never split across folds.
    """

    predicted: float
    label: int
    kind: str
    group: str


@dataclass(frozen=True, slots=True)
class SpliceVariantCase:
    """One scored sequence change for the **variant-effect** task.

    Attributes:
        predicted: The backend's score for this variant -- typically ``abs`` of a
            delta. Note BT4's ``delta_splicing`` is a difference of whole-sequence
            pooled top-k log-odds, **not** the published SpliceAI/Pangolin delta, so
            published cutoffs do not transfer to it.
        label: ``1`` if the assay called the variant splice-disrupting, else ``0``.
        region: ``"exonic"`` or ``"intronic"`` -- the stratum that matters for BT4,
            which designs coding sequence and therefore lives in the weaker half.
        group: The leakage-control group -- a gene or assay. **Read the panel's
            label definition before pooling groups:** a benchmark may aggregate
            several assays' differing criteria under one boolean, in which case a
            single pooled figure is a composite, not a measurement.
    """

    predicted: float
    label: int
    region: str
    group: str


@dataclass(frozen=True, slots=True)
class SpliceStratumReport:
    """Metrics for one stratum (a ``kind`` or a ``region``), scored on its own.

    Attributes:
        name: The stratum's label.
        n_cases: Cases in this stratum.
        n_positive: Positives in this stratum.
        prevalence: ``n_positive / n_cases`` -- the floor of :attr:`pr_auc`, carried
            beside it so the number is never separated from the question it answers.
        pr_auc: Average precision within the stratum.
        pr_auc_skill: :attr:`pr_auc` rescaled to ``[no-skill = 0, perfect = 1]``, the
            only PR figure comparable across panels of differing prevalence.
        roc_auc: Prevalence-stable discrimination; not a skill claim here.
        top_k_accuracy: Pooled top-k in the Zeng & Li construction.
        mcc: Matthews correlation at :attr:`SpliceGateReport.threshold`.
        brier: Mean squared error of the probabilities.
        brier_skill: Brier against the base-rate predictor; ``0.0`` means no better.
        ece: Expected calibration error.
        reliability: ``(mean_probability, observed_rate, count)`` per non-empty bin.
    """

    name: str
    n_cases: int
    n_positive: int
    prevalence: float
    pr_auc: float
    pr_auc_skill: float
    roc_auc: float
    top_k_accuracy: float
    mcc: float
    brier: float
    brier_skill: float
    ece: float
    reliability: tuple[tuple[float, float, int], ...]


@dataclass(frozen=True, slots=True)
class SpliceGateReport:
    """The honest result of a splice acceptance gate. Nothing is flipped.

    Attributes:
        task: :data:`SITE_PREDICTION` or :data:`VARIANT_EFFECT`.
        passed: ``True`` only if **every** stratum cleared its own thresholds. A
            pooled pass is never sufficient: it would let a backend certify on
            intronic strength while failing where BT4 operates.
        reasons: Why the gate did not pass, one plain sentence each. Empty on a pass.
        n_cases / n_positive / n_groups: Panel size, positives, leakage groups.
        threshold: The operating point at which :attr:`SpliceStratumReport.mcc` was
            computed. A gate input, not a blessed value -- deriving it is one of the
            things this gate exists to inform.
        negative_construction: How the panel's negatives were built, recorded
            verbatim. Without it a PR-AUC threshold is passable by thinning
            negatives, so a report that omits it cannot be adjudicated.
        panel_note: Free text for provenance the numbers depend on -- annotation
            release, assay, label definition.
        overall: Metrics pooled across strata, **explicitly a comparability
            statistic and never the pass axis**.
        strata: Per-stratum metrics, sorted by name. These carry the verdict.
        min_pr_auc / min_pr_auc_skill / max_ece: The thresholds each stratum was
            required to clear, as supplied.
    """

    task: str
    passed: bool
    reasons: tuple[str, ...]
    n_cases: int
    n_positive: int
    n_groups: int
    threshold: float
    negative_construction: str
    panel_note: str
    overall: SpliceStratumReport
    strata: tuple[SpliceStratumReport, ...]
    min_pr_auc: float
    min_pr_auc_skill: float
    max_ece: float
    reliability_bins_used: int = field(default=10)


def _stratum(
    name: str,
    labels: Sequence[int],
    scores: Sequence[float],
    *,
    threshold: float,
    n_bins: int,
) -> SpliceStratumReport:
    """Score one stratum. Assumes a non-empty case list."""
    n = len(labels)
    n_pos = sum(1 for v in labels if v == 1)
    hard = [1 if s > threshold else 0 for s in scores]
    return SpliceStratumReport(
        name=name,
        n_cases=n,
        n_positive=n_pos,
        prevalence=n_pos / n,
        pr_auc=pr_auc(labels, scores),
        pr_auc_skill=pr_auc_skill(labels, scores),
        roc_auc=roc_auc(labels, scores),
        top_k_accuracy=top_k_accuracy(labels, scores),
        mcc=mcc(labels, hard),
        brier=brier_score(labels, scores),
        brier_skill=brier_skill_score(labels, scores),
        ece=expected_calibration_error(labels, scores, n_bins=n_bins),
        reliability=tuple(reliability_bins(labels, scores, n_bins=n_bins)),
    )


def verify_splice_gate(
    cases: Sequence[SpliceSiteCase | SpliceVariantCase],
    *,
    negative_construction: str,
    panel_note: str = "",
    threshold: float = DEFAULT_SITE_PROBABILITY,
    min_pr_auc: float = 0.0,
    min_pr_auc_skill: float = 0.0,
    max_ece: float = 1.0,
    n_bins: int = 10,
) -> SpliceGateReport:
    """Score a splice backend's predictions and report honestly.

    Strata are the ``kind`` of a :class:`SpliceSiteCase` or the ``region`` of a
    :class:`SpliceVariantCase`; the task is inferred from the case type, and mixing
    the two is refused.

    ``passed`` requires **every** stratum to clear every threshold. The defaults are
    deliberately permissive (``min_pr_auc=0.0``, ``max_ece=1.0``) so that calling
    this without arguments produces a *report*, not a pass -- the thresholds are the
    maintainer's to set at gate time, on the evidence, and this module does not bless
    a bar a weak backend could clear (CLAUDE.md sections 6 and 10.6).

    Args:
        cases: Held-out cases, **all of one type**. Typed as a sequence of the union
            rather than a union of sequences so element access stays checkable;
            mixing the two types is refused at runtime, not by the annotation.
            Must be non-empty.
        negative_construction: How the panel's negatives were built. Required,
            because a PR-AUC threshold is meaningless without a pinned denominator.
        panel_note: Any further provenance the numbers depend on.
        threshold: Operating point for the MCC and the hard call.
        min_pr_auc: Minimum per-stratum average precision to pass.
        min_pr_auc_skill: Minimum per-stratum prevalence-normalized PR skill.
        max_ece: Maximum per-stratum expected calibration error.
        n_bins: Reliability bin count.

    Returns:
        A :class:`SpliceGateReport`. This function never flips ``calibrated``.

    Raises:
        ValueError: If ``cases`` is empty, mixes case types, ``negative_construction``
            is blank, there are fewer than two leakage groups, or a threshold
            argument is out of range.
    """
    if not cases:
        raise ValueError("splice gate needs at least one case")
    if not negative_construction.strip():
        raise ValueError(
            "negative_construction is required: without a pinned negative denominator "
            "a PR-AUC threshold is passable by thinning negatives"
        )
    if not 0.0 < threshold < 1.0:
        raise ValueError(f"threshold must be in (0, 1), got {threshold}")

    kinds = {type(c) for c in cases}
    if len(kinds) != 1:
        raise ValueError(
            "cases must all be one type: site prediction and variant effect are "
            "different tasks against different ground truth and cannot be pooled"
        )
    case_type = kinds.pop()
    if case_type is SpliceSiteCase:
        task = SITE_PREDICTION
        stratum_field = "kind"
    elif case_type is SpliceVariantCase:
        task = VARIANT_EFFECT
        stratum_field = "region"
    else:
        raise ValueError(f"unknown case type {case_type.__name__}")

    def strata_of(case: SpliceSiteCase | SpliceVariantCase) -> str:
        """The stratum a case belongs to: its ``kind`` or its ``region``."""
        return str(getattr(case, stratum_field))

    groups = {c.group for c in cases}
    if len(groups) < 2:
        raise ValueError(
            "splice gate needs at least two leakage-control groups (chromosomes for "
            f"site prediction, genes/assays for variant effect), got {len(groups)}"
        )

    labels = [c.label for c in cases]
    scores = [c.predicted for c in cases]
    overall = _stratum("overall", labels, scores, threshold=threshold, n_bins=n_bins)

    # One pass to bucket indices by stratum. (`cases.index(...)` would be O(n^2)
    # and would mis-assign duplicate cases to the first match.)
    by_stratum: dict[str, list[int]] = {}
    for i, case in enumerate(cases):
        by_stratum.setdefault(strata_of(case), []).append(i)

    strata: list[SpliceStratumReport] = []
    reasons: list[str] = []
    for name in sorted(by_stratum):
        idx = by_stratum[name]
        s_labels = [labels[i] for i in idx]
        s_scores = [scores[i] for i in idx]
        n_pos = sum(1 for v in s_labels if v == 1)
        if n_pos == 0 or n_pos == len(s_labels):
            # A single-class stratum cannot be scored, in either direction. With no
            # positives `pr_auc` returns its honest 0.0, which would read as measured
            # incompetence; with no negatives it returns 1.0, which would read as a
            # perfect result. Both are artifacts of the panel, not findings about the
            # model, so the stratum is refused rather than reported.
            missing = "positives" if n_pos == 0 else "negatives"
            reasons.append(
                f"stratum {name!r} has no {missing} ({n_pos}/{len(s_labels)} positive), "
                "so it cannot be scored"
            )
            continue
        report = _stratum(name, s_labels, s_scores, threshold=threshold, n_bins=n_bins)
        strata.append(report)
        if report.pr_auc < min_pr_auc:
            reasons.append(
                f"stratum {name!r} PR-AUC {report.pr_auc:.4f} < required {min_pr_auc:.4f}"
            )
        if report.pr_auc_skill < min_pr_auc_skill:
            reasons.append(
                f"stratum {name!r} PR-AUC skill {report.pr_auc_skill:.4f} "
                f"< required {min_pr_auc_skill:.4f}"
            )
        if report.ece > max_ece:
            reasons.append(
                f"stratum {name!r} ECE {report.ece:.4f} > allowed {max_ece:.4f}"
            )

    if not strata:
        reasons.append("no stratum could be scored")

    return SpliceGateReport(
        task=task,
        passed=not reasons,
        reasons=tuple(reasons),
        n_cases=len(cases),
        n_positive=sum(1 for v in labels if v == 1),
        n_groups=len(groups),
        threshold=threshold,
        negative_construction=negative_construction.strip(),
        panel_note=panel_note.strip(),
        overall=overall,
        strata=tuple(strata),
        min_pr_auc=min_pr_auc,
        min_pr_auc_skill=min_pr_auc_skill,
        max_ece=max_ece,
        reliability_bins_used=n_bins,
    )
