"""Run the splice acceptance gate over an annotated panel, against fixed baselines.

:func:`bt4.biomodels.splice.verify_splice_gate` answers "does this backend clear these
thresholds?". That is necessary and nowhere near sufficient, because a threshold is
only meaningful next to what a *dumb* predictor scores on the same panel -- and on this
task the dumb predictors are unusually strong. This module is the orchestration that
makes the answer usable: it scores a
:class:`~bt4.biomodels.splice.panel.SplicePanel` with a chosen backend, runs the gate,
runs **the same gate on every baseline**, and reports a verdict built from all three
conditions that have to hold at once.

**Why the baselines are not optional here in particular.** Roughly 99% of human introns
open ``GT`` and close ``AG``, so "is the canonical dinucleotide here?" is already a
usable splice-site detector -- and BT4 ships
:class:`~bt4.biomodels.splice.baseline.ConsensusPwmSplicePredictor`, which reads that
motif plus its consensus context, for free and with no licence. A wrapped CNN that
cannot beat those two has not earned a PyTorch dependency, a hash-pinned weight set, or
a non-commercial licence term, however impressive its absolute PR-AUC looks. The
``constant`` baseline is permanent for the mirror-image reason: it is perfectly
calibrated and completely useless, so its excellent ECE is visible in the same table
rather than being a trap the reader has to remember.

**Two alignment traps, both reported rather than assumed.**

* *Anchor offset.* A backend anchors its per-position score somewhere -- BT4's PWM
  baseline on the intronic dinucleotide, another model perhaps on the exonic boundary.
  One base of disagreement turns a good model into a hopeless one, silently. So
  ``anchor_offset`` is an explicit input and the report carries an
  :class:`AlignmentDiagnostic` showing where the backend's score actually peaked around
  each true site. A modal offset that is not zero is a wiring bug, and the note says so.
* *Combined tracks.* Pangolin emits **one** ``P(splice)`` track and leaves ``acceptor``
  all-zero, because its head is a binary softmax that never separated the two. Scoring
  that panel with a donor/acceptor split would report it as perfectly hopeless at
  acceptors -- an artifact of the wrapper, not a finding about the model. So a combined
  track collapses to a single ``"splice"`` stratum, detected from the track itself and
  recorded in the report.

**Nothing here flips a flag.** ``promotable`` means "the pre-registered conditions held
on this panel" -- and it requires that conditions were *actually* pre-registered, since
the gate's own defaults are permissive on purpose and a bar nobody set is not a bar that
held. Promotion is a separate, deliberate, recorded step (CLAUDE.md
sections 6/8/10.6), and for a splice backend it additionally requires an integration
fidelity attestation -- a different gate answering a different question.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.base import DEFAULT_SITE_PROBABILITY, SpliceResult
from bt4.biomodels.splice.gate import (
    SpliceGateReport,
    SpliceSiteCase,
    verify_splice_gate,
)
from bt4.biomodels.splice.panel import SplicePanel, SpliceWindow, canonical_motif_at
from bt4.pipeline.splice_crosscheck import resolve_splice_backend

__all__ = [
    "COMBINED_TRACK_EPSILON",
    "PEAK_SEARCH",
    "SPLICE_BASELINES",
    "AlignmentDiagnostic",
    "SpliceGateComparison",
    "SpliceGateSettings",
    "baseline_predictions",
    "run_splice_panel_gate",
    "score_splice_panel",
]

SPLICE_BASELINES: tuple[str, ...] = ("permutation", "gt_ag", "pwm", "constant")
"""Baselines a splice backend must beat. Kept permanently: a control that disappears
once it is inconvenient was never a control."""

PEAK_SEARCH: tuple[int, ...] = (-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
"""Offsets probed by the alignment diagnostic. Reported, never fitted -- the anchor is
the caller's declared input, and this only says whether the declaration looks right."""


@dataclass(frozen=True, slots=True)
class SpliceGateSettings:
    """The gate's thresholds and modes, gathered so they can be pre-registered as one.

    Attributes:
        threshold: The operating point the MCC is computed at. Defaults to BT4's own
            shared operating point, so the gate measures the cutoff BT4 actually uses
            rather than an ad-hoc one.
        min_pr_auc: Per-stratum average-precision floor. A **pre-commitment** recorded
            at gate time, never a bar this module blesses.
        min_pr_auc_skill: Per-stratum floor on average precision rescaled so no-skill
            is 0 and perfect is 1 at any prevalence. Note this is **not** a
            cross-panel-comparable quantity -- see
            :func:`~bt4.biomodels._stats.pr_auc_skill` -- so a bar set here is only
            meaningful against a pinned ``negative_construction``.
        max_ece: Per-stratum expected-calibration-error ceiling.
        n_bins: Reliability bins.
        seed: Seed for the permutation baseline (invariant #7).
    """

    threshold: float = DEFAULT_SITE_PROBABILITY
    min_pr_auc: float = 0.0
    min_pr_auc_skill: float = 0.0
    max_ece: float = 1.0
    n_bins: int = 10
    seed: int = 0


@dataclass(frozen=True, slots=True)
class AlignmentDiagnostic:
    """Where a backend's score actually peaked, relative to each declared site.

    The runner's counterpart to the panel's motif check: the panel proves its positions
    mean what the format says, and this proves the *backend* agrees about where a site
    is anchored. Both failures look identical in the metrics -- a competent model
    scoring near zero -- and neither is visible without being measured.

    Attributes:
        counts: ``(offset, number of sites peaking there)`` over :data:`PEAK_SEARCH`.
            These are **residual** shifts, measured on the already-aligned track, so
            ``0`` means the declared ``anchor_offset`` is right.
        n_sites: Sites the probe ran over.
        n_flat: Sites whose probe window was flat, so no peak existed to align. Tracked
            separately because the tie-break would otherwise resolve them to ``0`` and
            report a silent backend as perfectly aligned.
        applied_offset: The ``anchor_offset`` that was in force during the probe.
            Carried so :meth:`note` can name the **absolute** value to declare rather
            than the residual -- a user told to "pass anchor_offset=+1" while already
            passing +2 would end up further from the truth, not closer.
    """

    counts: tuple[tuple[int, int], ...]
    n_sites: int
    applied_offset: int = 0
    n_flat: int = 0

    @property
    def modal_offset(self) -> int:
        """The most common peak offset, preferring ``0`` on a tie."""
        return max(self.counts, key=lambda item: (item[1], item[0] == 0), default=(0, 0))[0]

    @property
    def recommended_offset(self) -> int:
        """The absolute ``anchor_offset`` to declare: what was applied, plus the residual."""
        return self.applied_offset + self.modal_offset

    @property
    def fraction_at_zero(self) -> float:
        """Fraction of sites whose peak landed exactly on the declared position."""
        if not self.n_sites:
            return 0.0
        return dict(self.counts).get(0, 0) / self.n_sites

    def note(self) -> str:
        """Return a one-line verdict on whether the anchors agree."""
        if not self.n_sites:
            return "alignment: no sites to probe"
        # A completely flat probe window has no argmax, and the tie-break resolves it to
        # 0 -- so a backend emitting no signal at all was being credited with peaks
        # landing exactly on every declared site. That is a positive alignment claim
        # from an absence of evidence, in the one diagnostic whose job is to separate
        # "misaligned" from "scoring near zero".
        if self.n_flat == self.n_sites:
            return (
                "alignment: the backend's scores are FLAT around every declared site, so "
                "there is no peak to align. This says nothing about the anchors -- it "
                "says the backend produced no signal on this panel"
            )
        if self.n_flat:
            return (
                f"alignment: {self.n_flat}/{self.n_sites} sites have flat scores with no "
                f"peak to align; of the rest, "
            ) + self._peak_note()
        return "alignment: " + self._peak_note()

    def _peak_note(self) -> str:
        """The alignment verdict over the sites that actually carried a peak."""
        if self.modal_offset:
            return (
                f"the backend's score peaks {self.modal_offset:+d} from the "
                f"declared position for most sites (only {self.fraction_at_zero:.1%} "
                f"peak exactly on it) -- the anchors DISAGREE. Pass "
                f"anchor_offset={self.recommended_offset:+d} (currently "
                f"{self.applied_offset:+d}) after confirming the backend's convention, "
                "rather than reading these metrics as model quality"
            )
        return (
            f"peaks land on the declared position for "
            f"{self.fraction_at_zero:.1%} of sites with anchor_offset="
            f"{self.applied_offset:+d} -- anchors agree"
        )


@dataclass(frozen=True, slots=True)
class SpliceGateComparison:
    """A backend's gate report, every baseline's, and the verdict built from them.

    Attributes:
        panel_hash: The panel's content hash, so a result is bound to exact bytes.
        backend: The backend's ``name``.
        backend_calibrated: Its honesty flag at the time of the run. For a splice
            backend this reflects *integration fidelity*; this gate is about the
            separate statistical question, so the two are reported side by side and
            never conflated.
        settings: The thresholds and modes used.
        strata: The stratum names scored, in report order.
        head: The backend's gate report.
        baselines: ``(name, report)`` for each baseline, in :data:`SPLICE_BASELINES`
            order.
        best_baseline: Per stratum, ``(stratum, baseline name, its skill)`` for the
            strongest baseline there.
        beats_every_baseline: The backend's ``pr_auc_skill`` exceeds every baseline's
            **in every stratum** -- the per-stratum rule applied to the comparison, so
            beating the motif on donors cannot excuse losing to it on acceptors.
        held_out: The panel contains no chromosome either model trained on.
        thresholds_declared: Whether the caller actually set a bar. The gate ships
            permissive defaults on purpose -- a threshold this module blessed would be
            one a weak backend could be pointed at -- but that leaves ``passed``
            trivially true for a bare call, which would make "the pre-registered
            conditions held" vacuous. So :attr:`promotable` additionally requires that
            conditions were pre-registered at all.
        alignment: Where the backend's peaks landed relative to the declared sites.
        promotable: All three conditions held **on this panel**. Not a promotion.
        notes: Human-readable notes carried into the report.
    """

    panel_hash: str
    backend: str
    backend_calibrated: bool
    settings: SpliceGateSettings
    strata: tuple[str, ...]
    head: SpliceGateReport
    baselines: tuple[tuple[str, SpliceGateReport], ...]
    best_baseline: tuple[tuple[str, str, float], ...]
    beats_every_baseline: bool
    held_out: bool
    thresholds_declared: bool
    alignment: AlignmentDiagnostic
    promotable: bool
    notes: tuple[str, ...]


def _tracks(result: SpliceResult, combined: bool) -> dict[str, tuple[float, ...]]:
    """Split a :class:`SpliceResult` into the per-stratum tracks to be scored.

    A combined-track backend collapses to one ``"splice"`` stratum rather than being
    credited with an all-zero acceptor track it never claimed to produce.

    **The collapse unions the two tracks rather than taking the donor one.** For a
    genuinely combined backend the two are equivalent -- Pangolin's acceptor track is
    identically zero, so the maximum is its donor track -- but the *baselines* are
    scored through this same function, and BT4's PWM baseline is a real two-track
    predictor. Taking its donor track alone would make every acceptor site invisible
    to the control, understating it by ~0.31 skill on a mixed panel and making a
    combined-track head that much easier to beat. Per position, "is this a splice site
    of either kind" is the stronger of the two claims, which is also how ``gt_ag``
    unions its two motifs and how ``pooled_risk`` unions the tracks.
    """
    if combined:
        return {
            "splice": tuple(
                max(donor, acceptor)
                for donor, acceptor in zip(result.donor, result.acceptor, strict=True)
            )
        }
    return {"donor": result.donor, "acceptor": result.acceptor}


COMBINED_TRACK_EPSILON: float = 1e-9
"""Below this, an acceptor track counts as absent rather than merely small.

Exact ``== 0.0`` was too strict: a backend whose acceptor channel carries float32
softmax dust (or any caller-supplied ``results``, a path the API explicitly offers)
would be scored with a donor/acceptor split, producing exactly the artifact the
detection exists to prevent -- and silently, since the explanatory note only appears
when the collapse happens. No real per-position probability is meaningfully non-zero
at this scale."""


def _is_combined(results: Sequence[SpliceResult]) -> bool:
    """Return whether the backend emits one combined track rather than two.

    Detected from the output: an all-but-zero acceptor track over the whole panel is
    either Pangolin's binary head (which never separated the two) or a backend whose
    acceptor scores are unusable anyway. Either way it must not be scored as an acceptor
    prediction. Callers can override the detection when they know better.
    """
    return all(
        abs(value) < COMBINED_TRACK_EPSILON for result in results for value in result.acceptor
    )


def score_splice_panel(panel: SplicePanel, backend: str) -> list[SpliceResult]:
    """Score every window with ``backend``, in panel order.

    Args:
        panel: The annotated panel.
        backend: A splice backend name (``"pwm"``, ``"pangolin"``, ``"spliceai"``).
            Resolved through :func:`~bt4.pipeline.splice_crosscheck.resolve_splice_backend`,
            so a committed attestation is honored only under the standing opt-in.

    Returns:
        One :class:`~bt4.biomodels.splice.base.SpliceResult` per window.

    Raises:
        ValueError: On an unknown backend name.
    """
    predictor = resolve_splice_backend(backend)
    return [predictor.score_sequence(window.sequence) for window in panel.windows]


def _labels(window: SpliceWindow, kind: str) -> tuple[int, ...]:
    """Return per-position labels for ``kind``, merging both kinds for ``"splice"``."""
    if kind != "splice":
        return window.labels(kind)
    sites = set(window.donors) | set(window.acceptors)
    return tuple(1 if i in sites else 0 for i in range(len(window.sequence)))


def _aligned(track: Sequence[float], position: int, anchor_offset: int) -> float:
    """Return the backend's score for ``position`` under the declared anchor offset.

    A position whose aligned index falls outside the track scores ``0.0`` -- the same
    convention the backends themselves use for positions without full flanking context.
    """
    index = position + anchor_offset
    return track[index] if 0 <= index < len(track) else 0.0


def _build_cases(
    panel: SplicePanel,
    results: Sequence[SpliceResult],
    combined: bool,
    anchor_offset: int,
) -> tuple[list[SpliceSiteCase], list[tuple[int, int]]]:
    """Build one case per (position, stratum), plus a back-reference to the sequence.

    Returns:
        ``(cases, refs)`` where ``refs[j]`` is the ``(window index, position)`` case
        ``j`` came from -- what the sequence-derived baselines need in order to score
        the same positions without re-deriving the alignment.
    """
    cases: list[SpliceSiteCase] = []
    refs: list[tuple[int, int]] = []
    for index, (window, result) in enumerate(zip(panel.windows, results, strict=True)):
        for kind, track in _tracks(result, combined).items():
            labels = _labels(window, kind)
            for position, label in enumerate(labels):
                cases.append(
                    SpliceSiteCase(
                        predicted=_aligned(track, position, anchor_offset),
                        label=label,
                        kind=kind,
                        group=window.group,
                    )
                )
                refs.append((index, position))
    return cases, refs


def _alignment(
    panel: SplicePanel,
    results: Sequence[SpliceResult],
    combined: bool,
    anchor_offset: int,
) -> AlignmentDiagnostic:
    """Probe where the backend's score peaks around each declared site.

    Computed on the **aligned** track, so a modal offset of zero means the declared
    ``anchor_offset`` is right and any other value is the residual correction needed.
    """
    counts = dict.fromkeys(PEAK_SEARCH, 0)
    n_sites = 0
    n_flat = 0
    for window, result in zip(panel.windows, results, strict=True):
        tracks = _tracks(result, combined)
        for position, kind in window.sites():
            track = tracks.get("splice" if combined else kind)
            if track is None:
                continue
            n_sites += 1
            probed = [_aligned(track, position + d, anchor_offset) for d in PEAK_SEARCH]
            if max(probed) - min(probed) < COMBINED_TRACK_EPSILON:
                # No peak exists here. Counting the tie-break's 0 would manufacture an
                # "anchors agree" claim out of a backend that produced nothing.
                n_flat += 1
                continue
            # Ties resolve to the offset nearest zero, then the lower one, so the probe
            # is deterministic (invariant #7).
            best = max(
                zip(PEAK_SEARCH, probed, strict=True),
                key=lambda item: (item[1], -abs(item[0]), -item[0]),
            )[0]
            counts[best] += 1
    return AlignmentDiagnostic(
        counts=tuple(sorted(counts.items())),
        n_sites=n_sites,
        applied_offset=anchor_offset,
        n_flat=n_flat,
    )


def baseline_predictions(
    name: str,
    panel: SplicePanel,
    cases: Sequence[SpliceSiteCase],
    refs: Sequence[tuple[int, int]],
    head_predictions: Sequence[float],
    seed: int,
) -> list[float]:
    """Return one baseline's predictions over the same cases, in case order.

    ``permutation`` is the null: the backend's own predictions shuffled **within each
    stratum**, so it preserves that stratum's score distribution exactly and measures
    what this panel yields from no relationship at all. The rest are things BT4 already
    has for free -- ``gt_ag``, the canonical dinucleotide rule that ~99% of human
    introns follow, and ``pwm``, the consensus baseline
    :func:`~bt4.biomodels.splice.default` already returns -- plus ``constant``, the
    per-stratum base rate, which is perfectly calibrated and carries no information.

    Raises:
        ValueError: If ``name`` is not in :data:`SPLICE_BASELINES`.
    """
    if name == "permutation":
        rng = random.Random(seed)
        by_kind: dict[str, list[int]] = {}
        for index, case in enumerate(cases):
            by_kind.setdefault(case.kind, []).append(index)
        shuffled = list(head_predictions)
        for kind in sorted(by_kind):
            indices = by_kind[kind]
            values = [head_predictions[i] for i in indices]
            rng.shuffle(values)
            for index, value in zip(indices, values, strict=True):
                shuffled[index] = value
        return shuffled
    if name == "gt_ag":
        return [
            1.0
            if canonical_motif_at(panel.windows[window].sequence, position, case.kind)
            else 0.0
            for case, (window, position) in zip(cases, refs, strict=True)
        ]
    if name == "pwm":
        # The PWM baseline anchors on the same intronic dinucleotide the panel format
        # pins, so it is scored at offset 0 by construction rather than by assumption.
        # In combined mode `_tracks` unions its donor and acceptor tracks, so the
        # control keeps its full strength instead of going blind to every acceptor.
        results = score_splice_panel(panel, "pwm")
        combined = any(case.kind == "splice" for case in cases)
        tracks = [_tracks(result, combined) for result in results]
        return [
            _aligned(tracks[window][case.kind], position, 0)
            for case, (window, position) in zip(cases, refs, strict=True)
        ]
    if name == "constant":
        totals: dict[str, list[int]] = {}
        for case in cases:
            bucket = totals.setdefault(case.kind, [0, 0])
            bucket[0] += case.label
            bucket[1] += 1
        rates = {kind: hits / total for kind, (hits, total) in totals.items()}
        return [rates[case.kind] for case in cases]
    raise ValueError(f"unknown baseline {name!r}; choose from {list(SPLICE_BASELINES)}")


def _gate(
    cases: Sequence[SpliceSiteCase],
    predictions: Sequence[float],
    panel: SplicePanel,
    settings: SpliceGateSettings,
) -> SpliceGateReport:
    """Run the gate over ``predictions``, re-labelling the shared cases."""
    relabelled = [
        SpliceSiteCase(
            predicted=value, label=case.label, kind=case.kind, group=case.group
        )
        for case, value in zip(cases, predictions, strict=True)
    ]
    return verify_splice_gate(
        relabelled,
        negative_construction=panel.negative_construction,
        panel_note=panel.annotation,
        threshold=settings.threshold,
        min_pr_auc=settings.min_pr_auc,
        min_pr_auc_skill=settings.min_pr_auc_skill,
        max_ece=settings.max_ece,
        n_bins=settings.n_bins,
    )


def run_splice_panel_gate(
    panel: SplicePanel,
    backend: str = "pwm",
    *,
    settings: SpliceGateSettings | None = None,
    baselines: Sequence[str] = SPLICE_BASELINES,
    anchor_offset: int = 0,
    combined_track: bool | None = None,
    results: Sequence[SpliceResult] | None = None,
) -> SpliceGateComparison:
    """Score ``panel`` with ``backend``, gate it, gate every baseline, and compare.

    Args:
        panel: The annotated splice panel.
        backend: A splice backend name (``"pwm"``, ``"pangolin"``, ``"spliceai"``).
        settings: Thresholds and modes. The defaults set no bar, so a bare call
            produces a report and can never certify anything.
        baselines: Which baselines to run. Defaults to all of
            :data:`SPLICE_BASELINES`; dropping one is a deliberate, visible choice.
        anchor_offset: Where the backend anchors its score relative to the panel's
            convention (see the module docstring). Declared by the caller, never
            fitted -- :attr:`SpliceGateComparison.alignment` reports whether the
            declaration looks right.
        combined_track: Force single-``"splice"``-stratum scoring on (``True``) or off
            (``False``). ``None`` detects it from the backend's own output.
        results: Pre-computed per-window scores, to gate a panel without re-running a
            heavy model (or to evaluate a backend this registry cannot construct).

    Returns:
        A :class:`SpliceGateComparison`. **This function flips nothing.**

    Raises:
        ValueError: On an unknown backend or baseline, a ``results`` length mismatch,
            or any refusal from the gate itself.
    """
    settings = settings or SpliceGateSettings()
    unknown = [name for name in baselines if name not in SPLICE_BASELINES]
    if unknown:
        raise ValueError(
            f"unknown baseline(s) {unknown}; choose from {list(SPLICE_BASELINES)}"
        )

    notes: list[str] = []
    if results is None:
        scored = score_splice_panel(panel, backend)
        probe = resolve_splice_backend(backend)
        name, calibrated = probe.name, probe.calibrated
    else:
        if len(results) != len(panel.windows):
            raise ValueError(
                f"results has {len(results)} entries for {len(panel.windows)} windows"
            )
        scored = list(results)
        name = scored[0].model_name if scored else backend
        # `all`, not `any`: one calibrated window must never certify a whole run, and a
        # calibration flag is the last place to take the generous reading (section 10.6).
        # The `bool(scored)` guard is load-bearing -- `all(())` is True.
        calibrated = bool(scored) and all(result.calibrated for result in scored)
        notes.append("scores supplied by the caller")

    combined = _is_combined(scored) if combined_track is None else combined_track
    if combined:
        notes.append(
            "the backend emits one combined P(splice) track, so donor and acceptor are "
            "scored as a single 'splice' stratum rather than crediting it with an "
            "acceptor prediction it does not make"
        )

    cases, refs = _build_cases(panel, scored, combined, anchor_offset)
    head_predictions = [case.predicted for case in cases]
    alignment = _alignment(panel, scored, combined, anchor_offset)
    notes.append(alignment.note())

    head = _gate(cases, head_predictions, panel, settings)
    reports = tuple(
        (
            baseline,
            _gate(
                cases,
                baseline_predictions(
                    baseline, panel, cases, refs, head_predictions, settings.seed
                ),
                panel,
                settings,
            ),
        )
        for baseline in baselines
    )

    head_skill = {stratum.name: stratum.pr_auc_skill for stratum in head.strata}
    best: list[tuple[str, str, float]] = []
    beats = bool(head_skill)
    for stratum in sorted(head_skill):
        winner, top = "none", float("-inf")
        for baseline, report in reports:
            for scored_stratum in report.strata:
                if scored_stratum.name == stratum and scored_stratum.pr_auc_skill > top:
                    winner, top = baseline, scored_stratum.pr_auc_skill
        # A stratum no baseline scored is an UNCONTESTED stratum, not a won one. Left as
        # the initial `-inf` it read as a win against a control that never ran -- and it
        # reported `-inf` as though it were that control's skill.
        if winner == "none":
            beats = False
            best.append((stratum, "none (uncontested)", float("nan")))
            notes.append(
                f"no baseline scored the {stratum!r} stratum, so the backend was not "
                "contested there; an uncontested stratum is never counted as beaten"
            )
            continue
        best.append((stratum, winner, top))
        if head_skill[stratum] <= top:
            beats = False

    # `SPLICE_BASELINES` says the controls are "kept permanently: a control that
    # disappears once it is inconvenient was never a control". Nothing enforced that:
    # dropping `pwm` alone let BT4's own shipped default certify itself with one keyword.
    # Running a subset is still allowed -- it is useful for a quick look -- it just
    # cannot produce a recommendation.
    missing_baselines = tuple(b for b in SPLICE_BASELINES if b not in tuple(baselines))
    if missing_baselines:
        notes.append(
            f"baseline(s) {list(missing_baselines)} were not run, so this comparison is "
            "incomplete and cannot be promotable. Dropping the control a backend would "
            "lose to is exactly how a weak result gets certified"
        )

    # A bar the caller never set is not a bar that held. The defaults are permissive by
    # design (this module must not bless a threshold), so the honest consequence is that
    # a bare call reports numbers and cannot recommend anything.
    thresholds_declared = (
        settings.min_pr_auc > 0.0
        or settings.min_pr_auc_skill > 0.0
        or settings.max_ece < 1.0
    )
    if not thresholds_declared:
        notes.append(
            "no threshold was declared (min_pr_auc / min_pr_auc_skill / max_ece are all "
            "at their permissive defaults), so 'gate passed' is vacuous here and this "
            "run cannot be promotable. Set the bar deliberately, before the run"
        )

    overlap = panel.training_overlap
    unclassified = panel.unclassified_groups
    held_out = not overlap and not unclassified
    if overlap:
        notes.append(
            f"NOT HELD OUT: {list(overlap)} are chromosomes both models trained on, so "
            "these metrics are optimistic and cannot support promotion"
        )
    if unclassified:
        # "I do not recognise this name" must never read as "held out". A GENCODE/Ensembl
        # panel names chromosomes `2`, `4`, `X` -- bare, no `chr` -- so a panel drawn
        # ENTIRELY from training chromosomes used to report held_out=True.
        notes.append(
            f"group(s) {list(unclassified)} are not recognisable as human chromosomes, so "
            "whether this panel is held out could not be established. Name groups as "
            "chromosomes (`chr1` or `1`); an unrecognised group is never assumed clean"
        )

    return SpliceGateComparison(
        panel_hash=panel.content_hash(),
        backend=name,
        backend_calibrated=calibrated,
        settings=settings,
        strata=tuple(sorted(head_skill)),
        head=head,
        baselines=reports,
        best_baseline=tuple(best),
        beats_every_baseline=beats,
        held_out=held_out,
        thresholds_declared=thresholds_declared,
        alignment=alignment,
        promotable=bool(
            head.passed
            and beats
            and held_out
            and thresholds_declared
            and not missing_baselines
        ),
        notes=tuple(notes),
    )
