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
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.base import DEFAULT_SITE_PROBABILITY, SpliceResult
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
    "KindAlignment",
    "SpliceGateComparison",
    "SpliceGateSettings",
    "baseline_predictions",
    "run_splice_panel_gate",
    "score_splice_panel",
]

SPLICE_BASELINES: tuple[str, ...] = ("permutation", "gt_ag", "pwm", "constant")
"""Baselines a splice backend must beat. Kept permanently: a control that disappears
once it is inconvenient was never a control."""

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
