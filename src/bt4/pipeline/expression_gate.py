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
    utr_context_sha256,
)
from bt4.biomodels.expression.gate import (
    ExpressionEvalCase,
    ExpressionGateReport,
    verify_expression_gate,
)

__all__ = [
    "BASELINES",
    "GateComparison",
    "GateScope",
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
class GateScope:
    """How the run was actually configured -- the record that makes a claim checkable.

    :class:`GateSettings` says what the gate *demanded*; this says what it *scored*. The
    distinction is load-bearing, because everything here changes the number and none of
    it used to survive the run: an attestation's ``species`` and ``cell_types`` were
    caller-declared free text, and the JSON record omitted the cell types, ``top_k`` and
    the UTR context entirely, so a finished run could not be reconstructed from its own
    output. :func:`~bt4.biomodels.expression.attest_expression` now derives its scope from
    this object and refuses any declaration that disagrees with it.

    The ``panel_*`` fields are what the *panel* declares about itself, kept separate from
    what the run was configured with precisely so the two can be compared rather than
    conflated.

    Attributes:
        species: The weight set scored.
        cell_types: The cell-type selection averaged, sorted. Empty means every one of
            them, which is a different quantity from any single cell line.
        top_k: Cross-validation runs ensembled.
        batch_size: Inference batch size. Recorded for reconstructability only -- it
            cannot change a score (RiboNN pads to a fixed width, ``shuffle=False``), so
            no attestation binds it.
        num_workers: DataLoader workers. Recorded, not bound, for the same reason.
        readout: The panel's own declared readout when it declares exactly one, else
            ``""`` (zero or several, so the caller must name it).
        utr_context_sha256: One hash per distinct ``(utr5, utr3)`` context in the panel,
            sorted -- the transcript contexts the measurement was made in.
        panel_species: Distinct non-empty ``species`` values the panel declares.
        panel_cell_types: Distinct non-empty ``cell_type`` values the panel declares.
        panel_readouts: Distinct non-empty ``readout`` values the panel declares.
        scoring_source: ``"gate"`` when this module invoked the backend, or
            ``"caller_supplied"`` when ``head_scores`` were handed in. The point at which
            the link between the named backend and the numbers stops being mechanical, so
            it is recorded rather than assumed away.
    """

    species: str
    cell_types: tuple[str, ...]
    top_k: int
    batch_size: int
    num_workers: int
    readout: str
    utr_context_sha256: tuple[str, ...]
    panel_species: tuple[str, ...]
    panel_cell_types: tuple[str, ...]
    panel_readouts: tuple[str, ...]
    scoring_source: str


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
        scope: How the run was configured (:class:`GateScope`), so the result is
            reconstructable from its own output and a later scope declaration can be
            checked against it rather than believed.
        notes: Human-readable notes (e.g. how many backend invocations were needed).
    """

    panel_hash: str
    backend: str
    backend_calibrated: bool
    settings: GateSettings
    scope: GateScope
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
            # Never honour a standing $BT4_EXPRESSION_USE_ATTESTED here: the gate is
            # what *decides* promotion, so scoring with an already-promoted head would
            # let a prior attestation colour the run that judges the next one.
            use_attested=False,
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
            ``head_scores``, a configuration that contradicts what the panel declares it
            measured (see :func:`_refuse_scope_mismatch` -- checked **before** the
            scoring pass, because this is a run-once procedure), or any refusal from the
            gate itself (too few groups, an empty fold, no rankable group in within-group
            mode).
    """
    settings = settings or GateSettings()
    unknown = [name for name in baselines if name not in BASELINES]
    if unknown:
        raise ValueError(f"unknown baseline(s) {unknown}; choose from {list(BASELINES)}")
    _refuse_scope_mismatch(panel, species=species, cell_types=cell_types)

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

    probe = resolve_backend(backend, species=species, use_attested=False)
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
        scope=_scope(
            panel,
            species=species,
            cell_types=cell_types,
            top_k=top_k,
            batch_size=batch_size,
            num_workers=num_workers,
            scoring_source="caller_supplied" if head_scores is not None else "gate",
        ),
        head=head,
        baselines=reports,
        best_baseline=best_name,
        best_baseline_spearman=best_rho,
        beats_every_baseline=beats,
        interval_is_informative=informative,
        promotable=bool(head.passed and beats and informative),
        notes=tuple(notes),
    )


def _panel_column(panel: ExpressionPanel, attribute: str) -> tuple[str, ...]:
    """Return the distinct non-empty values of one optional panel column, sorted."""
    return tuple(
        sorted({value for value in (getattr(row, attribute) for row in panel.rows) if value})
    )


def _refuse_scope_mismatch(
    panel: ExpressionPanel, *, species: str, cell_types: tuple[str, ...]
) -> None:
    """Refuse a run whose configuration contradicts what the panel says it measured.

    The trap this closes is silent and expensive: leave ``--cell-type`` off and RiboNN
    averages all 78 human cell types, which against a single-cell-line panel is a
    different quantity entirely -- no error, no warning, and a clean run to a wrong
    verdict. The gate is meant to be run **once**, so the check belongs here, before the
    scoring pass, not at promotion time after the budget is spent.

    Only what the panel actually declares is checked; a panel without ``species`` /
    ``cell_type`` columns declares nothing and is run as configured (with the gap
    recorded in :class:`GateScope`, not hidden). A maintainer who deliberately wants the
    all-cell-type average against a single-line panel drops the column -- which changes
    the panel hash, so the record stays honest about being a different panel.

    Raises:
        ValueError: On a declared/configured mismatch, naming what to pass instead.
    """
    declared_species = _panel_column(panel, "species")
    if declared_species and declared_species != (species,):
        raise ValueError(
            f"panel declares species {list(declared_species)} but the head is "
            f"configured for {species!r}; RiboNN's human and mouse weights are "
            "different models. Pass the matching --species."
        )
    declared_cells = _panel_column(panel, "cell_type")
    selection = tuple(sorted(cell_types))
    if declared_cells and declared_cells != selection:
        shown = list(selection) if selection else "every cell type (no selection)"
        raise ValueError(
            f"panel was measured in {list(declared_cells)} but the head would score "
            f"{shown}. Averaging every cell type against a single-cell-line measurement "
            f"is a scope error, not a rounding one: pass "
            + " ".join(f"--cell-type {name}" for name in declared_cells)
            + " (or drop the panel's cell_type column if the all-tissue average really "
            "is what you mean)."
        )


def _scope(
    panel: ExpressionPanel,
    *,
    species: str,
    cell_types: tuple[str, ...],
    top_k: int,
    batch_size: int,
    num_workers: int,
    scoring_source: str,
) -> GateScope:
    """Build the :class:`GateScope` recording how this run was configured."""
    readouts = _panel_column(panel, "readout")
    return GateScope(
        species=species,
        cell_types=tuple(sorted(cell_types)),
        top_k=top_k,
        batch_size=batch_size,
        num_workers=num_workers,
        # One readout is unambiguous and can be carried; zero or several means the
        # caller has to name it, so the attestation never guesses the assay.
        readout=readouts[0] if len(readouts) == 1 else "",
        utr_context_sha256=tuple(
            sorted(utr_context_sha256(utr5, utr3) for utr5, utr3 in panel.contexts())
        ),
        panel_species=_panel_column(panel, "species"),
        panel_cell_types=_panel_column(panel, "cell_type"),
        panel_readouts=readouts,
        scoring_source=scoring_source,
    )
