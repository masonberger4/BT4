"""Run the expression acceptance gate over a measured panel, against fixed baselines.

:func:`bt4.biomodels.expression.verify_expression_gate` answers "does this head clear
these thresholds?". That is necessary and nowhere near sufficient, because a threshold
is only meaningful next to what a *dumb* predictor scores on the same panel. This module
is the orchestration that makes the answer usable: it scores a
:class:`~bt4.biomodels.expression.panel.ExpressionPanel` with a chosen backend, runs the
gate, runs **the same gate on every baseline**, and reports a verdict built from all
three of the conditions that have to hold at once.

**Why the baselines are not optional.** A within-protein Spearman of 0.3 is worthless if
plain CAI scores 0.35 -- BT4 already computes CAI, directly, inside the optimizer loop,
so a head that cannot beat it has added nothing. And because split conformal is valid
for *any* score function, a **constant predictor** achieves exactly valid coverage; it
is a permanent baseline precisely so that its "pass" on the coverage axis is visible in
the same table rather than being a trap the reader must remember.

**The verdict.** ``promotable`` requires (a) the gate's own thresholds, (b) the head's
cluster-bootstrap CI lower bound above **every** baseline's point estimate, and (c) an
interval narrower than the spread of the labels. All three are reported separately, so a
failure says *which* condition failed.

**Nothing here flips a flag.** ``promotable`` means "the pre-registered conditions held
on this panel". Promotion is a separate, deliberate, recorded step (CLAUDE.md
§6/§8/§10.6).
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.codon.tables import CodonUsageTable, load_table
from bt4.biomodels.expression import (
    BatchExpressionPredictor,
    ExpressionPanel,
    PanelRow,
    resolve_backend,
)
from bt4.biomodels.expression.gate import (
    ExpressionEvalCase,
    ExpressionGateReport,
    verify_expression_gate,
)

__all__ = [
    "BASELINES",
    "GateComparison",
    "GateSettings",
    "baseline_scores",
    "run_panel_gate",
    "score_panel",
]

BASELINES: tuple[str, ...] = ("permutation", "cai", "gc3", "length", "constant")
"""Baselines a head must beat. Kept permanently: a control that disappears once it is
inconvenient was never a control."""


@dataclass(frozen=True, slots=True)
class GateSettings:
    """The gate's thresholds and modes, gathered so they can be pre-registered as one.

    Attributes:
        within_group: Score inside each protein -- the strict bar, and BT4's regime.
        recalibrate: Fit the affine link on the calibration fold before residuals.
        target_coverage: Conformal level (1 - alpha).
        coverage_tolerance: Allowed absolute gap from the target, two-sided.
        min_spearman: Absolute rank threshold. A **pre-commitment**, not a community
            standard -- no such standard exists -- so it is recorded, never presented as
            authoritative.
        calibration_fraction: Fraction of *groups* in the calibration fold.
        bootstrap_resamples: Cluster-bootstrap resamples for the CI.
        seed: Seed for the bootstrap and the permutation baseline (invariant #7).
    """

    within_group: bool = False
    recalibrate: bool = False
    target_coverage: float = 0.90
    coverage_tolerance: float = 0.05
    min_spearman: float = 0.30
    calibration_fraction: float = 0.50
    bootstrap_resamples: int = 1000
    seed: int = 0


@dataclass(frozen=True, slots=True)
class GateComparison:
    """A head's gate report, every baseline's, and the verdict built from them.

    Attributes:
        panel_hash: The panel's content hash, so a result is bound to exact bytes.
        backend: The backend's ``name``.
        backend_calibrated: The backend's honesty flag at the time of the run.
        settings: The thresholds and modes used.
        head: The head's gate report.
        baselines: ``(name, report)`` for each baseline, in :data:`BASELINES` order.
        best_baseline: The baseline with the highest primary metric.
        best_baseline_spearman: That baseline's primary metric.
        beats_every_baseline: Head's CI lower bound above the best baseline's estimate.
        interval_is_informative: Median interval width below the label IQR.
        promotable: All three conditions held **on this panel**. Not a promotion.
        notes: Human-readable notes (e.g. how many backend invocations were needed).
    """

    panel_hash: str
    backend: str
    backend_calibrated: bool
    settings: GateSettings
    head: ExpressionGateReport
    baselines: tuple[tuple[str, ExpressionGateReport], ...]
    best_baseline: str
    best_baseline_spearman: float
    beats_every_baseline: bool
    interval_is_informative: bool
    promotable: bool
    notes: tuple[str, ...]


def _gc3(dna: str) -> float:
    """GC fraction at codon third positions -- the axis synonymous recoding moves most."""
    thirds = dna[2::3]
    return sum(1 for ch in thirds if ch in "GC") / len(thirds) if thirds else 0.0


def score_panel(
    panel: ExpressionPanel,
    backend: str,
    *,
    species: str = "human",
    cell_types: tuple[str, ...] = (),
    top_k: int = 5,
    batch_size: int = 64,
    num_workers: int = 0,
) -> tuple[list[float], list[str]]:
    """Score every panel row, building **one predictor per distinct UTR context**.

    A predictor carries its fixed ``utr5`` / ``utr3`` on the model rather than per call,
    so a panel spanning transcripts with different UTRs genuinely needs one predictor
    each -- and no more. Within a context the whole bucket goes through a single batched
    invocation, because a heavy backend's cost is dominated by fixed per-invocation
    overhead (for RiboNN: hashing 90 weight files and loading 50 models).

    Returns:
        ``(scores in panel order, notes)``. Order is restored after bucketing, so a
        measurement can never end up paired with a different sequence's score.
    """
    position = {row.variant_id: i for i, row in enumerate(panel.rows)}
    scores: list[float] = [0.0] * len(panel.rows)
    contexts = panel.contexts()
    for (utr5, utr3), rows in contexts.items():
        predictor = resolve_backend(
            backend,
            species=species,
            utr5=utr5,
            utr3=utr3,
            top_k=top_k,
            batch_size=batch_size,
            num_workers=num_workers,
            cell_types=cell_types,
        )
        dnas = [row.cds for row in rows]
        if isinstance(predictor, BatchExpressionPredictor):
            values = [result.score for result in predictor.score_many(dnas)]
        else:
            values = [predictor.score_sequence(dna).score for dna in dnas]
        for row, value in zip(rows, values, strict=True):
            scores[position[row.variant_id]] = value
    notes = [
        f"{len(contexts)} UTR context(s) => {len(contexts)} backend invocation(s) "
        f"for {len(panel.rows)} rows"
    ]
    return scores, notes


def baseline_scores(
    name: str,
    rows: Sequence[PanelRow],
    head_scores: Sequence[float],
    table: CodonUsageTable,
    seed: int,
) -> list[float]:
    """Return one baseline's predictions over the panel, in panel order.

    ``permutation`` is the null: the head's own predictions against a deterministic
    shuffle, so it measures what this panel yields from no real relationship at all. The
    rest are cheap features the head must beat to have added anything -- ``cai`` above
    all, since BT4 optimizes it directly and for free.

    Raises:
        ValueError: If ``name`` is not in :data:`BASELINES`.
    """
    if name == "permutation":
        shuffled = list(head_scores)
        random.Random(seed).shuffle(shuffled)
        return shuffled
    if name == "cai":
        return [table.cai(row.cds) for row in rows]
    if name == "gc3":
        return [_gc3(row.cds) for row in rows]
    if name == "length":
        return [float(len(row.cds)) for row in rows]
    if name == "constant":
        return [0.0] * len(rows)
    raise ValueError(f"unknown baseline {name!r}; choose from {list(BASELINES)}")


def _gate(
    rows: Sequence[PanelRow],
    predictions: Sequence[float],
    settings: GateSettings,
) -> ExpressionGateReport:
    """Run the gate over ``predictions`` with the shared settings."""
    cases = [
        ExpressionEvalCase(predicted=score, measured=row.measured, group=row.group)
        for row, score in zip(rows, predictions, strict=True)
    ]
    return verify_expression_gate(
        cases,
        target_coverage=settings.target_coverage,
        coverage_tolerance=settings.coverage_tolerance,
        min_spearman=settings.min_spearman,
        calibration_fraction=settings.calibration_fraction,
        within_group=settings.within_group,
        recalibrate=settings.recalibrate,
        bootstrap_resamples=settings.bootstrap_resamples,
        bootstrap_seed=settings.seed,
    )


def run_panel_gate(
    panel: ExpressionPanel,
    backend: str = "ribonn",
    *,
    settings: GateSettings | None = None,
    baselines: Sequence[str] = BASELINES,
    species: str = "human",
    cell_types: tuple[str, ...] = (),
    top_k: int = 5,
    batch_size: int = 64,
    num_workers: int = 0,
    organism: str = "homo_sapiens",
    reference_set: str | None = None,
    head_scores: Sequence[float] | None = None,
) -> GateComparison:
    """Score ``panel`` with ``backend``, gate it, gate every baseline, and compare.

    Args:
        panel: The measured CDS-variant panel.
        backend: Expression backend name (``"ribonn"`` / ``"null"``).
        settings: Thresholds and modes; defaults are the pooled, unlinked ones, so a
            caller must opt in to the strict bar deliberately.
        baselines: Which baselines to run. Defaults to all of :data:`BASELINES`.
        species / cell_types / top_k / batch_size / num_workers: Backend configuration.
            ``cell_types`` should name the panel's own cell line: averaging 78 tissues
            against a single-cell-line measurement is a scope error.
        organism / reference_set: Codon table for the CAI baseline, so that number
            always names the question it answers (CLAUDE.md §8).
        head_scores: Pre-computed scores, to gate a panel without re-running the model
            (or to evaluate a head this registry cannot construct). When ``None`` the
            panel is scored via :func:`score_panel`.

    Returns:
        A :class:`GateComparison`. **This function flips nothing.**

    Raises:
        ValueError: On an unknown backend or baseline, a length mismatch in
            ``head_scores``, or any refusal from the gate itself (too few groups, an
            empty fold, no rankable group in within-group mode).
    """
    settings = settings or GateSettings()
    unknown = [name for name in baselines if name not in BASELINES]
    if unknown:
        raise ValueError(f"unknown baseline(s) {unknown}; choose from {list(BASELINES)}")

    if head_scores is None:
        scores, notes = score_panel(
            panel,
            backend,
            species=species,
            cell_types=cell_types,
            top_k=top_k,
            batch_size=batch_size,
            num_workers=num_workers,
        )
    else:
        if len(head_scores) != len(panel.rows):
            raise ValueError(
                f"head_scores has {len(head_scores)} values for {len(panel.rows)} rows"
            )
        scores, notes = list(head_scores), ["scores supplied by the caller"]

    probe = resolve_backend(backend, species=species)
    table = load_table(organism, reference_set=reference_set)

    head = _gate(panel.rows, scores, settings)
    reports = tuple(
        (
            name,
            _gate(
                panel.rows,
                baseline_scores(name, panel.rows, scores, table, settings.seed),
                settings,
            ),
        )
        for name in baselines
    )

    best_name, best_rho = "none", float("-inf")
    for name, report in reports:
        if report.spearman > best_rho:
            best_name, best_rho = name, report.spearman

    beats = head.spearman_ci_low > best_rho
    informative = head.width_over_iqr < 1.0
    return GateComparison(
        panel_hash=panel.content_hash(),
        backend=probe.name,
        backend_calibrated=probe.calibrated,
        settings=settings,
        head=head,
        baselines=reports,
        best_baseline=best_name,
        best_baseline_spearman=best_rho,
        beats_every_baseline=beats,
        interval_is_informative=informative,
        promotable=bool(head.passed and beats and informative),
        notes=tuple(notes),
    )
