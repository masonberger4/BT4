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
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.agreement import SiteCallAgreement, site_call_agreement
from bt4.biomodels.splice.base import DEFAULT_SITE_PROBABILITY, SpliceResult
from bt4.biomodels.splice.designed_panel import DesignedCdsPanel
from bt4.biomodels.splice.gate import (
    SpliceGateReport,
    SpliceSiteCase,
    verify_splice_gate,
)
from bt4.biomodels.splice.panel import (
    DEFAULT_EDGE_MARGIN,
    SplicePanel,
    SpliceWindow,
    canonical_motif_at,
)
from bt4.pipeline.splice_crosscheck import resolve_splice_backend

__all__ = [
    "CNN_ANCHOR_OFFSETS",
    "COMBINED_TRACK_EPSILON",
    "PEAK_SEARCH",
    "SPLICE_BASELINES",
    "AlignmentDiagnostic",
    "DesignedCdsProbe",
    "KindAlignment",
    "ProgressCallback",
    "SpliceGateComparison",
    "SpliceGateSettings",
    "baseline_predictions",
    "probe_designed_cds",
    "run_panel_backend_agreement",
    "run_splice_panel_gate",
    "score_splice_panel",
]

ProgressCallback = Callable[[int, int, str, int], None]
"""``(index, total, window_id, length)`` -- reported before each window is scored.

``index`` is 1-based. ``length`` is the window's nucleotide count, which is what the
wait is actually proportional to: windows are whole gene spans and vary by more than
an order of magnitude, so "12/20" alone does not predict the remaining time and the
length is what tells a reader a 2 Mb gene is next."""


SPLICE_BASELINES: tuple[str, ...] = ("permutation", "gt_ag", "pwm", "constant")
"""Baselines a splice backend must beat. Kept permanently: a control that disappears
once it is inconvenient was never a control."""

_ECE_TIE = 1e-9
"""How close two ECEs must be to count as matching, for the note below.

Far above float-summation noise (a base-rate baseline computes 1.1e-16 where an exact
oracle computes 0.0) and far below any difference that would mean something."""


CNN_ANCHOR_OFFSETS: dict[str, int] = {"donor": -1, "acceptor": 1}
"""Where SpliceAI and Pangolin put a site's score, relative to BT4's panel convention.

**Both models anchor on the exonic boundary base**; BT4's panel anchors on the intronic
dinucleotide. The gap is one base **in opposite directions for the two kinds**: a donor's
score sits on the last *exonic* base (BT4 position - 1) and an acceptor's on the first
*exonic* base (BT4 position + 1). The two backends agree with each other; donor and
acceptor disagree within each.

Established three ways: SpliceAI's training-label construction (``Y0[c-tx_start] = 2`` at
exon *ends* for donors, ``= 1`` at exon *starts* for acceptors), Pangolin's CLI using
gffutils' first/last exonic base as the sites, and direct measurement against the
hash-verified weights (34 sites, unanimous, both strands).

There is deliberately **no ``"splice"`` entry**. A combined track cannot be shifted by one
value without breaking the other kind, so offsets are keyed by *site* kind everywhere and
applied per site **before** the union (see :func:`_positive_indices`). An accepted-but-
inert ``"splice"`` key was worse than no key: it scored identically whatever it was set
to while the diagnostic reported it as the offset in force."""

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
class KindAlignment:
    """Where one site kind's scores peaked, relative to its declared positions.

    Attributes:
        kind: ``"donor"`` or ``"acceptor"`` -- the **site** kind, which is reported
            separately even when the backend's track is combined, because the two were
            shifted by different amounts.
        counts: ``(residual offset, sites peaking there)`` over :data:`PEAK_SEARCH`,
            measured **after** this kind's ``anchor_offset`` is applied, so ``0`` means
            the declared value is right.
        n_sites: Sites probed for this kind.
        n_flat: Sites whose probe window was flat, so no peak existed to align.
        applied_offset: The offset in force for this kind during the probe.
    """

    kind: str
    counts: tuple[tuple[int, int], ...]
    n_sites: int
    n_flat: int
    applied_offset: int

    @property
    def modal_offset(self) -> int:
        """The most common residual offset, preferring ``0`` on a tie."""
        return max(self.counts, key=lambda item: (item[1], item[0] == 0), default=(0, 0))[0]

    @property
    def recommended_offset(self) -> int:
        """The absolute offset to declare for this kind: applied plus residual."""
        return self.applied_offset + self.modal_offset

    @property
    def fraction_at_zero(self) -> float:
        """Fraction of this kind's sites peaking exactly on the declared position."""
        if not self.n_sites:
            return 0.0
        return dict(self.counts).get(0, 0) / self.n_sites

    @property
    def aligned(self) -> bool:
        """Whether this kind's declared anchor looks right.

        Requires a **majority** of sites to peak on the declared position, not merely a
        plurality among the few that carried any signal. ``n_flat < n_sites`` alone let a
        single peak outvote any number of flat sites: a kind where 7 of 8 sites produced
        nothing was reported as aligned at 12% agreement -- the same over-claim the
        per-kind rewrite exists to prevent, one level down.
        """
        return (
            bool(self.n_sites)
            and self.modal_offset == 0
            and self.fraction_at_zero > 0.5
        )


@dataclass(frozen=True, slots=True)
class AlignmentDiagnostic:
    """Where a backend's score actually peaked, relative to each declared site.

    The runner's counterpart to the panel's motif check: the panel proves its positions
    mean what the format says, and this proves the *backend* agrees about where a site is
    anchored. Both failures look identical in the metrics -- a competent model scoring
    near zero -- and neither is visible without being measured.

    **It is per site kind, and that is load-bearing.** SpliceAI and Pangolin both anchor
    on the exonic boundary base, which is one base *before* BT4's donor position and one
    base *after* its acceptor position -- opposite directions. A pooled probe under
    either single value sees half the sites aligned and half two bases off, and the modal
    tie-break resolves that to ``0``: it reported "anchors agree" at 50% alignment while
    endorsing a setting that left an entire stratum scoring ~0. Reporting per kind makes
    that state unrepresentable.

    Attributes:
        kinds: One :class:`KindAlignment` per site kind present, sorted.
    """

    kinds: tuple[KindAlignment, ...] = ()

    @property
    def n_sites(self) -> int:
        """Sites probed across every kind."""
        return sum(kind.n_sites for kind in self.kinds)

    @property
    def n_flat(self) -> int:
        """Sites with no peak to align, across every kind."""
        return sum(kind.n_flat for kind in self.kinds)

    @property
    def aligned(self) -> bool:
        """Whether **every** kind's declared anchor looks right."""
        return bool(self.kinds) and all(kind.aligned for kind in self.kinds)

    @property
    def recommended_offsets(self) -> dict[str, int]:
        """The offset to declare per kind -- ready to pass back as ``anchor_offset``."""
        return {kind.kind: kind.recommended_offset for kind in self.kinds}

    @property
    def fraction_at_zero(self) -> float:
        """Fraction of all probed sites peaking exactly on their declared position."""
        if not self.n_sites:
            return 0.0
        return sum(dict(k.counts).get(0, 0) for k in self.kinds) / self.n_sites

    def note(self) -> str:
        """Return a verdict naming each kind separately."""
        if not self.kinds:
            return "alignment: no sites to probe"
        parts: list[str] = []
        for kind in self.kinds:
            if not kind.n_sites:
                continue
            if kind.n_flat == kind.n_sites:
                parts.append(
                    f"{kind.kind}: scores are FLAT around every site, so there is no "
                    "peak to align -- this says nothing about the anchor"
                )
            elif kind.modal_offset:
                parts.append(
                    f"{kind.kind}: peaks {kind.modal_offset:+d} from the declared "
                    f"position ({kind.fraction_at_zero:.0%} on it) -- DISAGREES; "
                    f"declare {kind.recommended_offset:+d} (currently "
                    f"{kind.applied_offset:+d})"
                )
            else:
                parts.append(
                    f"{kind.kind}: peaks on the declared position for "
                    f"{kind.fraction_at_zero:.0%} of sites at "
                    f"{kind.applied_offset:+d} -- agrees"
                )
        verdict = "anchors agree" if self.aligned else "the anchors DISAGREE"
        return f"alignment ({verdict}) -- " + "; ".join(parts)


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


def _resolve_offsets(anchor_offset: int | Mapping[str, int]) -> dict[str, int]:
    """Return a per-kind offset map from a scalar or a mapping.

    A single scalar was the original design and it is **not sufficient for a real
    backend**: SpliceAI and Pangolin need ``-1`` for donors and ``+1`` for acceptors, so
    no one value is right for a mixed panel. A scalar still works -- it applies to every
    kind, which is correct for a kind-separated panel -- but a mapping is what a mixed
    one needs.

    Raises:
        ValueError: On an unknown site kind, which would otherwise be silently ignored
            and leave that stratum misaligned.
    """
    if isinstance(anchor_offset, int):
        return dict.fromkeys(("donor", "acceptor"), anchor_offset)
    unknown = sorted(set(anchor_offset) - {"donor", "acceptor"})
    if unknown:
        raise ValueError(
            f"unknown site kind(s) {unknown} in anchor_offset; expected 'donor' and/or "
            "'acceptor'. There is no 'splice' key: a combined track is still made of "
            "donors and acceptors, and they are shifted separately before the union"
        )
    return {kind: anchor_offset.get(kind, 0) for kind in ("donor", "acceptor")}


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


def score_splice_panel(
    panel: SplicePanel,
    backend: str,
    *,
    progress: ProgressCallback | None = None,
) -> list[SpliceResult]:
    """Score every window with ``backend``, in panel order.

    Args:
        panel: The annotated panel.
        backend: A splice backend name (``"pwm"``, ``"pangolin"``, ``"spliceai"``).
            Resolved through :func:`~bt4.pipeline.splice_crosscheck.resolve_splice_backend`,
            so a committed attestation is honored only under the standing opt-in.
        progress: Called **before** each window is scored, as
            ``(index, total, window_id, length)`` with a 1-based ``index``. A real panel
            takes tens of minutes on a CPU -- a wrapped CNN reads ~10 kb of context per
            position -- and reporting nothing for that long is indistinguishable from
            hanging. Reporting *before* rather than after means the slow window is named
            while it is the one being waited on. ``None`` is silent, so the API default
            stays print-free (section 3: only ``cli`` prints).

    Returns:
        One :class:`~bt4.biomodels.splice.base.SpliceResult` per window.

    Raises:
        ValueError: On an unknown backend name.
    """
    predictor = resolve_splice_backend(backend)
    total = len(panel.windows)
    results: list[SpliceResult] = []
    for index, window in enumerate(panel.windows, start=1):
        if progress is not None:
            progress(index, total, window.window_id, len(window.sequence))
        results.append(predictor.score_sequence(window.sequence))
    return results


def _positive_indices(
    window: SpliceWindow, kind: str, offsets: Mapping[str, int]
) -> set[int]:
    """Return the track indices where a site of ``kind`` is expected.

    The reformulation the per-kind offsets require: rather than shifting the *track* by
    one scalar, each site's positive is placed where its own kind's offset says the
    backend scores it. That is the only way a **combined** ``"splice"`` stratum can be
    handled at all -- its union drops the kind, so the shift has to happen per site,
    before the union.

    A site whose shifted index falls outside the window is **clamped back into it**
    rather than dropped. Dropping it would quietly remove a positive, lowering the
    prevalence and *raising* every metric -- so a backend that structurally cannot find an
    annotated site would be rewarded for it. Under the previous track-shifting
    formulation such a site stayed ``label=1`` with ``predicted=0.0``, a forced miss that
    depressed the score, and that is the honest behaviour to preserve. Clamping keeps the
    positive at the nearest scoreable index, where the backend's score is whatever it is.
    """
    n = len(window.sequence)
    indices: set[int] = set()
    for position, site_kind in window.sites():
        if kind != "splice" and site_kind != kind:
            continue
        shifted = position + offsets.get(site_kind, 0)
        indices.add(min(max(shifted, 0), n - 1))
    return indices


def _build_cases(
    panel: SplicePanel,
    results: Sequence[SpliceResult],
    combined: bool,
    offsets: Mapping[str, int],
) -> tuple[list[SpliceSiteCase], list[tuple[int, int]]]:
    """Build one case per (position, stratum), plus a back-reference to the sequence.

    Returns:
        ``(cases, refs)`` where ``refs[j]`` is the ``(window index, position)`` case
        ``j`` came from -- what the sequence-derived baselines need in order to score the
        same positions. Baselines are scored at their own anchors, not the backend's, so
        ``refs`` carries the panel's positions rather than the shifted ones.
    """
    cases: list[SpliceSiteCase] = []
    refs: list[tuple[int, int]] = []
    for index, (window, result) in enumerate(zip(panel.windows, results, strict=True)):
        for kind, track in _tracks(result, combined).items():
            positives = _positive_indices(window, kind, offsets)
            for position in range(len(window.sequence)):
                cases.append(
                    SpliceSiteCase(
                        predicted=track[position] if position < len(track) else 0.0,
                        label=1 if position in positives else 0,
                        kind=kind,
                        group=window.group,
                    )
                )
                refs.append((index, position))
    return cases, refs


def _track_at(track: Sequence[float], index: int) -> float:
    """Return ``track[index]``, or ``0.0`` outside it.

    The same convention the backends themselves use for positions without full flanking
    context.
    """
    return track[index] if 0 <= index < len(track) else 0.0


def _alignment(
    panel: SplicePanel,
    results: Sequence[SpliceResult],
    combined: bool,
    offsets: Mapping[str, int],
) -> AlignmentDiagnostic:
    """Probe, **per site kind**, where the backend's score peaks around each site.

    Per kind rather than pooled, and that is the whole point. SpliceAI and Pangolin need
    opposite offsets for donors and acceptors, so a pooled probe sees half the sites
    aligned and half two bases off under *either* value -- and the modal tie-break
    resolves that to ``0``, reporting "anchors agree" at 50% alignment while endorsing a
    setting that leaves an entire stratum unscored.
    """
    per_kind: dict[str, dict[int, int]] = {}
    totals: dict[str, int] = {}
    flats: dict[str, int] = {}
    for window, result in zip(panel.windows, results, strict=True):
        tracks = _tracks(result, combined)
        for position, site_kind in window.sites():
            # Keyed by SITE kind even when the stratum is combined: the stratum is one
            # track, but the two kinds were shifted by different amounts, so reporting a
            # single "splice" anchor would name a value that was never applied.
            track = tracks.get("splice" if combined else site_kind)
            if track is None:
                continue
            counts = per_kind.setdefault(site_kind, dict.fromkeys(PEAK_SEARCH, 0))
            totals[site_kind] = totals.get(site_kind, 0) + 1
            applied = offsets.get(site_kind, 0)
            probed = [_track_at(track, position + applied + d) for d in PEAK_SEARCH]
            if max(probed) - min(probed) < COMBINED_TRACK_EPSILON:
                flats[site_kind] = flats.get(site_kind, 0) + 1
                continue
            best = max(
                zip(PEAK_SEARCH, probed, strict=True),
                key=lambda item: (item[1], -abs(item[0]), -item[0]),
            )[0]
            counts[best] += 1
    return AlignmentDiagnostic(
        kinds=tuple(
            KindAlignment(
                kind=kind,
                counts=tuple(sorted(counts.items())),
                n_sites=totals.get(kind, 0),
                n_flat=flats.get(kind, 0),
                applied_offset=offsets.get(kind, 0),
            )
            for kind, counts in sorted(per_kind.items())
        )
    )


def _panel_position(index: int, kind: str, offsets: Mapping[str, int]) -> dict[str, int]:
    """Map a track index back to the panel position(s) it could describe, per site kind.

    Case ``j`` sits at *track* index ``i`` and is positive when a site of kind ``k`` sits
    at panel position ``i - offsets[k]``. The sequence-derived baselines anchor on the
    **panel's** convention, so they must be read there and not at ``i`` -- otherwise a
    declared ``anchor_offset`` silently moves every control one base off its own anchor
    while leaving the head correctly aligned.

    A combined ``"splice"`` stratum has no single inverse (its two kinds were shifted in
    opposite directions), so both are returned and the caller unions them -- exactly as
    the label was built.
    """
    kinds = ("donor", "acceptor") if kind == "splice" else (kind,)
    return {k: index - offsets.get(k, 0) for k in kinds}


def baseline_predictions(
    name: str,
    panel: SplicePanel,
    cases: Sequence[SpliceSiteCase],
    refs: Sequence[tuple[int, int]],
    head_predictions: Sequence[float],
    seed: int,
    offsets: Mapping[str, int] | None = None,
) -> list[float]:
    """Return one baseline's predictions over the same cases, in case order.

    ``permutation`` is the null: the backend's own predictions shuffled **within each
    stratum**, so it preserves that stratum's score distribution exactly and measures
    what this panel yields from no relationship at all. The rest are things BT4 already
    has for free -- ``gt_ag``, the canonical dinucleotide rule that ~99% of human
    introns follow, and ``pwm``, the consensus baseline
    :func:`~bt4.biomodels.splice.default` already returns -- plus ``constant``, the
    per-stratum base rate, which is perfectly calibrated and carries no information.

    Args:
        offsets: The per-kind anchor offsets in force for the head. The baselines are
            read at the **panel's** anchor, so this is used to invert the shift; ``None``
            means no shift was applied.

    Raises:
        ValueError: If ``name`` is not in :data:`SPLICE_BASELINES`.
    """
    shift = dict(offsets or {})
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
        motif: list[float] = []
        for case, (window, position) in zip(cases, refs, strict=True):
            sequence = panel.windows[window].sequence
            sites = _panel_position(position, case.kind, shift)
            hit = any(canonical_motif_at(sequence, at, k) for k, at in sites.items())
            motif.append(1.0 if hit else 0.0)
        return motif
    if name == "pwm":
        # The PWM baseline anchors on the same intronic dinucleotide the panel format
        # pins, so it is scored at offset 0 by construction rather than by assumption.
        # In combined mode `_tracks` unions its donor and acceptor tracks, so the
        # control keeps its full strength instead of going blind to every acceptor.
        results = score_splice_panel(panel, "pwm")
        # The PWM is a genuine two-track predictor scored at the panel's own anchor, so
        # its tracks are read per kind and unioned only where the head's stratum is.
        pwm_tracks = [
            {"donor": result.donor, "acceptor": result.acceptor} for result in results
        ]
        scored: list[float] = []
        for case, (window, position) in zip(cases, refs, strict=True):
            sites = _panel_position(position, case.kind, shift)
            scored.append(
                max(
                    (_track_at(pwm_tracks[window][k], at) for k, at in sites.items()),
                    default=0.0,
                )
            )
        return scored
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
    anchor_offset: int | Mapping[str, int] = 0,
    combined_track: bool | None = None,
    results: Sequence[SpliceResult] | None = None,
    progress: ProgressCallback | None = None,
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
            convention -- a scalar, or a **per-kind mapping** such as
            :data:`CNN_ANCHOR_OFFSETS` (``{"donor": -1, "acceptor": 1}``). A scalar is
            not sufficient for a mixed panel scored by a real backend: SpliceAI and
            Pangolin both anchor on the exonic boundary base, which is one base before
            BT4's donor position and one base *after* its acceptor position. Declared by
            the caller, never fitted -- :attr:`SpliceGateComparison.alignment` reports
            per kind whether the declaration looks right.
        combined_track: Force single-``"splice"``-stratum scoring on (``True``) or off
            (``False``). ``None`` detects it from the backend's own output.
        results: Pre-computed per-window scores, to gate a panel without re-running a
            heavy model (or to evaluate a backend this registry cannot construct).
        progress: Forwarded to :func:`score_splice_panel`; ignored when ``results`` is
            supplied, since nothing is scored then. ``None`` is silent.

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
        scored = score_splice_panel(panel, backend, progress=progress)
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
    else:
        # The mirror of the note above, and the one a reader actually needs. A
        # kind-separated run solves a HARDER problem than a combined one: in the donor
        # stratum an acceptor site is a negative, so the backend must locate the site
        # AND get its kind right. Setting a separated backend's number beside a
        # combined backend's single figure therefore understates the separated one --
        # and nothing in the table says so, which is exactly the comparison two
        # consecutive gate runs invite.
        notes.append(
            "this backend emits separate donor and acceptor tracks, so each kind is its "
            "own stratum -- and in the donor stratum an ACCEPTOR site is a negative, "
            "which is a harder task than a single combined 'is this a site at all'. "
            "These figures are therefore NOT comparable with a combined-track backend's "
            "single number; pass combined_track=True (CLI: --combined-track on) to score "
            "this backend on the same task and get a like-for-like figure"
        )

    edge = panel.edge_sites()
    if edge:
        # A forced miss is not a model failure, and a metric that silently absorbs one
        # is misleading in the direction that matters least obviously: downward.
        notes.append(
            f"{len(edge)} annotated site(s) sit within {DEFAULT_EDGE_MARGIN} nt of a "
            "window edge, where a backend has no flanking sequence and scores 0.0. They "
            "are FORCED MISSES that depress every metric here; extend those windows "
            f"rather than reading the result as model quality (first: {edge[0]})"
        )

    offsets = _resolve_offsets(anchor_offset)
    cases, refs = _build_cases(panel, scored, combined, offsets)
    head_predictions = [case.predicted for case in cases]
    alignment = _alignment(panel, scored, combined, offsets)
    notes.append(alignment.note())

    head = _gate(cases, head_predictions, panel, settings)
    reports = tuple(
        (
            baseline,
            _gate(
                cases,
                baseline_predictions(
                    baseline, panel, cases, refs, head_predictions, settings.seed, offsets
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
    # `max_ece` deliberately does NOT count. At splice prevalence ECE is not a bar a
    # backend can fail: the `constant` baseline predicts the base rate, so its ECE is
    # 0 BY CONSTRUCTION -- measured at this panel's prevalence, exactly the score a
    # PERFECT classifier gets -- and an all-zero predictor lands at 0.0003. A run whose
    # only declared threshold is on ECE has pre-registered a condition no predictor can
    # fail, which is the same hole as dropping the baseline a backend would lose to.
    # Promotion needs a DISCRIMINATION bar.
    discrimination_declared = settings.min_pr_auc > 0.0 or settings.min_pr_auc_skill > 0.0
    thresholds_declared = discrimination_declared
    if not thresholds_declared:
        if settings.max_ece < 1.0:
            notes.append(
                f"the only declared threshold is max_ece={settings.max_ece:g}, which is "
                "NOT a bar: at splice prevalence a base-rate predictor scores ECE 0.0 -- "
                "the same as a perfect classifier -- so no backend can fail it. Declare "
                "min_pr_auc_skill (or min_pr_auc) as well; this run cannot be promotable"
            )
        else:
            notes.append(
                "no threshold was declared (min_pr_auc / min_pr_auc_skill are at their "
                "permissive defaults), so 'gate passed' is vacuous here and this run "
                "cannot be promotable. Set the bar deliberately, before the run"
            )

    # And say it wherever the ECE column could be read as evidence, not only when it was
    # declared as a threshold: a baseline matching or beating the head on ECE is the
    # measured demonstration that the column is not discriminating on this panel.
    head_ece = {stratum.name: stratum.ece for stratum in head.strata}
    outscored = sorted(
        {
            f"{baseline_name}/{stratum.name}"
            for baseline_name, report in reports
            for stratum in report.strata
            # A tolerance, not `<=`: the note claims "match or beat", and two
            # equally-calibrated predictors differ by float noise (a base-rate
            # baseline lands at 1.1e-16 where an exact oracle lands at 0.0).
            if stratum.ece <= head_ece.get(stratum.name, float("inf")) + _ECE_TIE
        }
    )
    if outscored:
        notes.append(
            f"baseline(s) {outscored} match or beat the backend's ECE. ECE rewards "
            "predicting the base rate, which every one of these baselines does better "
            "than a model that commits, so read it as a description of the score "
            "distribution and NEVER as evidence of quality -- the skill column carries "
            "the verdict"
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


def run_panel_backend_agreement(
    panel: SplicePanel,
    backends: tuple[str, str] = ("pangolin", "spliceai"),
    *,
    anchor_offset: int | Mapping[str, int] = CNN_ANCHOR_OFFSETS,
    progress: ProgressCallback | None = None,
    results: Mapping[str, Sequence[SpliceResult]] | None = None,
) -> SiteCallAgreement:
    """Report whether two backends call the **same positions** on a site panel.

    Section 6 names cross-backend agreement a first-class uncertainty signal, and
    :func:`~bt4.biomodels.splice.agreement.backend_agreement` provides it for the
    *design* flow -- ranking candidate CDSs by Delta-splicing. It cannot answer the
    positional question a site panel poses, because a site panel has no candidates
    to rank. This does.

    **Comparing the two gate reports is not a substitute.** Two backends can each
    reach skill 0.98 on the same panel while being confident about different bases;
    the metrics would agree and the models would not.

    Both wrapped CNNs anchor identically (donor -1, acceptor +1), so no shift is
    applied *between* them; ``anchor_offset`` is used only to place the annotated
    sites in the backends' shared frame, which is what makes ``both`` /
    ``only_first`` / ``only_second`` meaningful.

    A backend emitting one combined track is collapsed exactly as the gate collapses
    it -- unioning donor and acceptor -- so Pangolin's single ``P(splice)`` series is
    compared against the stronger of SpliceAI's two per position, rather than
    against a track it never claimed to produce.

    Args:
        panel: The annotated site panel.
        backends: Exactly two backend names.
        anchor_offset: Where the backends score a site relative to the panel's
            convention. Defaults to :data:`CNN_ANCHOR_OFFSETS`, correct for both
            wrapped CNNs; a scalar applies to every kind.
        progress: Forwarded to :func:`score_splice_panel` for each backend in turn.
        results: Pre-computed per-window scores keyed by backend name, to avoid
            re-scoring a panel that takes tens of minutes per backend.

    Returns:
        A :class:`~bt4.biomodels.splice.agreement.SiteCallAgreement`.

    Raises:
        ValueError: If ``backends`` is not two distinct names, or supplied
            ``results`` do not cover both of them at the panel's length.
    """
    if len(backends) != 2 or backends[0] == backends[1]:
        raise ValueError(f"need two distinct backend names, got {backends!r}")

    scored: dict[str, Sequence[SpliceResult]] = {}
    for name in backends:
        if results is not None and name in results:
            supplied = results[name]
            if len(supplied) != len(panel.windows):
                raise ValueError(
                    f"results[{name!r}] has {len(supplied)} entries for "
                    f"{len(panel.windows)} windows"
                )
            scored[name] = supplied
        else:
            scored[name] = score_splice_panel(panel, name, progress=progress)

    # One backend may be combined and the other not -- exactly the Pangolin/SpliceAI
    # case -- but both collapse the same way: per position, "a site of either kind",
    # the union. For a combined backend the acceptor track is identically zero, so the
    # union IS its single track; for a two-track backend it is the stronger claim. So
    # there is no branch here, and adding one would imply a distinction that does not
    # exist.
    offsets = _resolve_offsets(anchor_offset)
    tracks = {
        name: [_tracks(result, combined=True)["splice"] for result in per_window]
        for name, per_window in scored.items()
    }

    site_indices = [
        sorted(_positive_indices(window, "splice", offsets)) for window in panel.windows
    ]
    return site_call_agreement(
        tracks[backends[0]], tracks[backends[1]], site_indices, backends
    )


@dataclass(frozen=True, slots=True)
class DesignedCdsProbe:
    """What splice backends do on designed synonymous CDS. **Not a gate.**

    There is no field named ``passed`` or ``promotable``, and there will not be: the
    panel carries no splice labels because designed coding sequence has no splice ground
    truth. Everything here is label-free, which is precisely why it is measurable in the
    one regime BT4 actually operates in.

    Attributes:
        group: The protein. Every measurement is *within* it -- see
            ``designed_panel``'s note on why pooling across proteins answers the wrong
            question.
        n_designs: Designs scored against the group's native reference.
        delta_spread: Per backend, ``max - min`` of Δsplicing across the designs.
            Synonymous positions are the only thing BT4 changes, so a backend whose
            spread is ~0 **and whose scores actually reached the pooling background**
            cannot distinguish the candidates BT4 produces. Read it with
            :attr:`sub_background` or it will be misread: a spread of zero is the
            *expected* output here whatever the model does, because the pooling hinge
            floors every sub-background score to zero.
        delta_range: Per backend, the ``(min, max)`` Δ, so the spread has a location --
            a spread of 0.1 around zero and one around -3 mean different things.
        rank_correlations: Pairwise Spearman agreement of the backends' Δ *rankings*.
            If BT4 ranked candidates by splice Δ, this is whether the winner would
            depend on which backend ran -- **but only when the Δs are not degenerate**.
            Over a backend whose Δs are all zero this is a correlation of constants.
        sign_agreement: Fraction of designs where all backends agree on the Δ's sign --
            better or worse than native. Same caveat as ``rank_correlations``.
        response_spread: Per backend, ``max - min`` of the **background-free** Δ
            (:func:`~bt4.biomodels.splice.base.pool_top_k_logit`). This is the field
            that answers question 1 below when the risk Δ has been floored: it is
            monotone in the model's scores at every score, so it separates candidates
            that ``delta_spread`` cannot. It is a *ranking* statistic and **not a
            risk** -- it has no calibrated zero and says nothing about how
            spliceogenic any sequence is.
        response_range: Per backend, the ``(min, max)`` background-free Δ.
        response_rank_correlations: ``rank_correlations`` recomputed over the
            background-free Δs -- the agreement number that survives a floored risk.
        response_sign_agreement: ``sign_agreement`` over the background-free Δs.
        sub_background: Per backend, how many of the ``n_designs + 1`` scored sequences
            had **no** position above the pooling background. Equal to
            ``n_designs + 1`` means the backend's every risk Δ is zero *by
            construction* and the risk-based fields above carry no information.
            Measured with the hash-verified Pangolin weights, this is what happens:
            all of them.
        max_score: Per backend, the highest per-position score seen anywhere in the
            group. The magnitude the hinge discarded -- ``0.44`` and ``0.001`` both
            floor to a risk of zero and mean entirely different things.
        backends: The backends' **own** names, in input order -- ``pwm`` resolves to
            ``consensus-pwm-baseline`` and Pangolin to a name carrying its tissue set.
            Those identify what actually ran, which a registry alias does not.
    """

    group: str
    n_designs: int
    backends: tuple[str, ...]
    delta_spread: dict[str, float]
    delta_range: dict[str, tuple[float, float]]
    rank_correlations: dict[tuple[str, str], float]
    sign_agreement: float
    response_spread: dict[str, float]
    response_range: dict[str, tuple[float, float]]
    response_rank_correlations: dict[tuple[str, str], float]
    response_sign_agreement: float
    sub_background: dict[str, int]
    max_score: dict[str, float]

    def degenerate(self, backend: str) -> bool:
        """``True`` when ``backend``'s risk Δs are all zero because pooling floored them.

        The distinction the probe's whole conclusion rests on: "this backend does not
        respond to synonymous change" versus "this backend responds and BT4 threw the
        response away". Without it a reader takes a zero spread for the former, which
        is what happened when this probe was first run.
        """
        return self.sub_background.get(backend, 0) == self.n_designs + 1


def probe_designed_cds(
    panel: DesignedCdsPanel,
    backends: Sequence[str] = ("pangolin", "spliceai"),
    *,
    progress: ProgressCallback | None = None,
) -> tuple[DesignedCdsProbe, ...]:
    """Measure splice backends on designed synonymous CDS, without any labels.

    Two questions, both answerable with no ground truth, and both unanswerable from the
    natural-sequence panels BT4 has measured so far:

    1. **Does a synonymous change move the score at all?** Every design in a group
       encodes the same protein, so any Δ between them is attributable to codon choice
       alone. A backend with ~0 spread is insensitive to exactly the axis BT4 varies.
    2. **Would two backends pick the same candidate?** Spearman over their Δ rankings.
       If BT4 used splice Δ to select among candidates, this says whether the selection
       is a property of the sequence or of the backend.

    **The agreement number here is NOT the site panel's Jaccard**, and the two must not
    be set side by side: that one is positional overlap of called sites on genomic
    sequence, this one is rank agreement over candidate Δs. Different statistics on
    different inputs.

    Args:
        panel: The designed-CDS panel.
        backends: Backend names to compare. One is allowed -- the spread is still
            meaningful -- but agreement needs two.
        progress: Called per sequence scored, as :data:`ProgressCallback`.

    Returns:
        One :class:`DesignedCdsProbe` per group, in the panel's group order.

    Raises:
        ValueError: On an unknown or duplicated backend name.
    """
    names = tuple(backends)
    if len(set(names)) != len(names):
        raise ValueError(f"backend names must be distinct, got {names!r}")

    from bt4.biomodels.splice.agreement import agreement_from_deltas, backend_agreement

    predictors = [resolve_splice_backend(name) for name in names]
    # Report under the predictors' OWN names, not the registry aliases asked for.
    # `backend_agreement` keys by `predictor.name`, and those names carry configuration
    # a bare alias loses -- "pangolin[heart+liver+brain+testis]" versus "pangolin" -- so
    # they are the honest provenance as well as the working key.
    resolved = tuple(predictor.name for predictor in predictors)

    probes: list[DesignedCdsProbe] = []
    for group in panel.groups:
        native = panel.native(group)
        designs = panel.designed(group)
        if progress is not None:
            progress(
                panel.groups.index(group) + 1, len(panel.groups), group, len(native.cds)
            )
        report = backend_agreement(
            predictors, [m.cds for m in designs], native.cds
        )
        # The background-free deltas get the identical rank/sign treatment, through the
        # same function -- so the two agreement numbers are the same statistic on two
        # poolings, and neither can be the flattering one by construction.
        response_report = agreement_from_deltas(
            report.response_by_backend, sign_epsilon=report.sign_epsilon
        )

        def _spread(
            by_backend: dict[str, tuple[float, ...]],
        ) -> tuple[dict[str, float], dict[str, tuple[float, float]]]:
            spread: dict[str, float] = {}
            ranges: dict[str, tuple[float, float]] = {}
            for name in resolved:
                values = by_backend[name]
                lo, hi = (min(values), max(values)) if values else (0.0, 0.0)
                spread[name] = hi - lo
                ranges[name] = (lo, hi)
            return spread, ranges

        delta_spread, delta_range = _spread(report.delta_by_backend)
        response_spread, response_range = _spread(report.response_by_backend)

        probes.append(
            DesignedCdsProbe(
                group=group,
                n_designs=len(designs),
                backends=resolved,
                delta_spread=delta_spread,
                delta_range=delta_range,
                rank_correlations=dict(report.rank_correlations),
                sign_agreement=report.sign_agreement,
                response_spread=response_spread,
                response_range=response_range,
                response_rank_correlations=dict(response_report.rank_correlations),
                response_sign_agreement=response_report.sign_agreement,
                sub_background=dict(report.n_sub_background),
                max_score=dict(report.max_score_by_backend),
            )
        )
    return tuple(probes)
