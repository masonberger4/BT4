"""Opt-in, out-of-loop splice cross-check on a delivered sequence (graceful).

This is the *final audit / validation pass* end of BT4's splice story (CLAUDE.md
section 6): given one already-delivered coding sequence and a chosen splice
backend, it reports that backend's predicted cryptic-splice sites and pooled risk.
It is what ``bt4 validate --splice-backend ...`` and ``bt4 optimize --check-splice
...`` call.

Two properties make it honest where BT3's in-loop ASSP scrape was not:

* **Never blocking.** Any backend failure -- an ASSP outage / garbled response
  (:class:`~bt4.biomodels.splice.assp.AsspError`), or a wrapped CNN's missing deps
  / weights -- is caught and turned into an ``available is False`` report with the
  reason, so a cross-check can **never** fail an optimization (CLAUDE.md section
  10.15). A genuinely invalid sequence (non-ACGT) still raises up front -- that is
  a caller error, not a service outage.
* **Network-derived numbers are labeled and non-reproducible.** The report carries
  the backend's :attr:`network_derived` flag; when it is ``True`` (ASSP) the caller
  keeps those numbers out of the reproducible-from-manifest guarantee (they are
  reported as a separate advisory section, never folded into a
  :class:`~bt4.domain.result.Result` audit or provenance manifest).

This module composes ``biomodels.splice`` (the backends) with ``domain`` -- a
legal ``pipeline`` dependency (CLAUDE.md section 3). It never edits the sequence;
it annotates.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.assp import AsspError, AsspSplicePredictor
from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SplicePredictor,
    SpliceResult,
    pooled_risk,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.pangolin import PangolinSplicePredictor
from bt4.biomodels.splice.spliceai import SpliceAiSplicePredictor
from bt4.domain.sequence import validate_dna

__all__ = [
    "DEFAULT_CROSSCHECK_THRESHOLD",
    "CrossCheckSite",
    "SpliceCrossCheck",
    "resolve_splice_backend",
    "run_splice_crosscheck",
]

# Backend failures the cross-check turns into a graceful "unavailable" report
# instead of raising: an ASSP outage / garbled response (AsspError), or a wrapped
# CNN's missing deps / weights. A ValueError (invalid sequence) is a CALLER error
# and is deliberately absent -- it is raised up front by validate_dna.
_DEGRADE_ERRORS: tuple[type[BaseException], ...] = (
    AsspError,
    ModuleNotFoundError,
    RuntimeError,
    FileNotFoundError,
    KeyError,
    OSError,
)

DEFAULT_CROSSCHECK_THRESHOLD: float = 0.5
"""Per-position score above which a dense-backend track is localized into a site.

A heuristic display knob (every cross-check backend is ``calibrated is False``),
used only for backends that report a dense per-position track (the PWM baseline and
the CNNs). ASSP already reports discrete classified sites, so its own site list is
used directly and this threshold does not apply to it.
"""

# The cross-check backends the CLI/API name-selects. ASSP is the opt-in network
# validator; the others let the same surface run an offline / installed backend.
_BACKENDS: dict[str, str] = {
    "assp": "assp",
    "pwm": "pwm",
    "baseline": "pwm",
    "consensus": "pwm",
    "pangolin": "pangolin",
    "spliceai": "spliceai",
}


@dataclass(frozen=True, slots=True)
class CrossCheckSite:
    """One localized splice site the cross-check backend flagged on the sequence.

    Attributes:
        position: 0-based anchor index of the site in the sequence.
        kind: ``"donor"``, ``"acceptor"``, or ``"splice"`` (a combined-track
            backend, e.g. Pangolin, cannot distinguish donor from acceptor).
        score: The backend's per-position score at ``position`` -- an uncalibrated
            pseudo-score (ASSP's confidence, the PWM baseline's pseudo-probability,
            or a CNN probability uncalibrated in BT4's ledger); never a calibrated
            splice probability.
        site_class: ASSP's classification (``constitutive`` / ``alternative`` /
            ``cryptic`` / ``unknown``) when available, else ``""``.
    """

    position: int
    kind: str
    score: float
    site_class: str = ""


@dataclass(frozen=True, slots=True)
class SpliceCrossCheck:
    """The result of a single-sequence splice cross-check (advisory, never edits).

    Attributes:
        dna: The audited coding sequence (returned unchanged).
        backend: The backend's :attr:`SplicePredictor.name`.
        available: Whether the backend actually produced a result. ``False`` marks a
            graceful degradation (see :attr:`reason`); the run was not failed.
        reason: Why the backend was unavailable (``None`` when ``available``).
        calibrated: The backend's ``calibrated`` flag -- ``False`` for every
            cross-check backend today, so the whole report is advisory.
        network_derived: ``True`` for ASSP -- its numbers are network-derived and
            excluded from the reproducible-from-manifest guarantee. ``False`` for the
            offline / hash-pinned backends.
        threshold: The localization threshold used for dense-track backends.
        top_k: The pooling depth used for :attr:`pooled_risk`.
        pooled_risk: The whole-sequence pooled top-k splice risk (``0.0`` when
            unavailable). Larger means more predicted risk.
        sites: The localized / classified splice sites (empty when unavailable),
            ordered by ``(position, kind)``.
    """

    dna: str
    backend: str
    available: bool
    reason: str | None
    calibrated: bool
    network_derived: bool
    threshold: float
    top_k: int
    pooled_risk: float
    sites: tuple[CrossCheckSite, ...]


def resolve_splice_backend(name: str) -> SplicePredictor:
    """Construct a splice cross-check backend by name.

    Args:
        name: One of ``"assp"`` (the opt-in online cross-check), ``"pwm"`` (aliases
            ``"baseline"`` / ``"consensus"``; the offline PWM baseline),
            ``"pangolin"``, or ``"spliceai"`` (the wrapped CNNs, when installed).

    Returns:
        A :class:`~bt4.biomodels.splice.base.SplicePredictor`.

    Raises:
        ValueError: If ``name`` is not a known backend.
    """
    key = _BACKENDS.get(name.strip().lower())
    if key is None:
        raise ValueError(
            f"unknown splice backend {name!r}; choose from {sorted(set(_BACKENDS))}"
        )
    if key == "assp":
        return AsspSplicePredictor()
    if key == "pwm":
        return ConsensusPwmSplicePredictor()
    if key == "pangolin":
        return PangolinSplicePredictor()
    return SpliceAiSplicePredictor()


def _localize_track(scores: Sequence[float], kind: str, threshold: float) -> list[CrossCheckSite]:
    """Localize a dense per-position track into one site per above-threshold run.

    Non-maximal suppression: a broad site spanning several positions yields exactly
    one flag at the run's peak (first index on ties, for determinism), rather than a
    cluster. Mirrors the audit's localization but stays a pipeline-level display
    helper (no cross-layer private import).
    """
    sites: list[CrossCheckSite] = []
    n = len(scores)
    i = 0
    while i < n:
        if scores[i] > threshold:
            best = i
            j = i + 1
            while j < n and scores[j] > threshold:
                if scores[j] > scores[best]:
                    best = j
                j += 1
            sites.append(CrossCheckSite(position=best, kind=kind, score=scores[best]))
            i = j
        else:
            i += 1
    return sites


def _sites_from_predictor(
    predictor: SplicePredictor, result: SpliceResult, dna: str, threshold: float
) -> tuple[CrossCheckSite, ...]:
    """Extract cross-check sites, preferring a backend's own discrete site list.

    ASSP reports discrete, classified sites via a ``sites(dna)`` method, which is
    richer than thresholding its sparse track (it keeps constitutive sites below the
    display threshold and their classification). Any other backend reports a dense
    per-position track, localized here at ``threshold``. A combined-track backend
    (Pangolin, all-zero acceptor) has its donor-track sites labelled ``"splice"``.
    """
    sites_method = getattr(predictor, "sites", None)
    if callable(sites_method):
        return tuple(
            CrossCheckSite(
                position=s.position, kind=s.kind, score=s.confidence, site_class=s.site_class
            )
            for s in sites_method(dna)
        )
    combined = not any(result.acceptor)
    out: list[CrossCheckSite] = []
    out += _localize_track(result.donor, "splice" if combined else "donor", threshold)
    if not combined:
        out += _localize_track(result.acceptor, "acceptor", threshold)
    out.sort(key=lambda s: (s.position, s.kind))
    return tuple(out)


def run_splice_crosscheck(
    dna: str,
    *,
    backend: str = "assp",
    predictor: SplicePredictor | None = None,
    threshold: float = DEFAULT_CROSSCHECK_THRESHOLD,
    top_k: int = DEFAULT_TOP_K,
) -> SpliceCrossCheck:
    """Cross-check one delivered sequence with a splice backend (never fails a run).

    Runs ``predictor`` (or the backend named by ``backend``) over ``dna`` and
    reports its predicted sites and pooled risk. Any backend failure -- an ASSP
    outage / garbled response, or a wrapped CNN's missing deps / weights -- is
    caught and reported as ``available is False`` with the reason; only a genuinely
    invalid sequence (non-ACGT) raises.

    Args:
        dna: The delivered coding sequence to audit.
        backend: Which backend to construct when ``predictor`` is not given (see
            :func:`resolve_splice_backend`).
        predictor: An explicit backend, overriding ``backend`` (e.g. a
            fixture-backed ASSP predictor in tests).
        threshold: Localization threshold for dense-track backends.
        top_k: Pooling depth for :attr:`SpliceCrossCheck.pooled_risk`.

    Returns:
        A :class:`SpliceCrossCheck`.

    Raises:
        ValueError: If ``dna`` is empty / non-ACGT, or ``backend`` is unknown.
    """
    seq = validate_dna(dna)
    pred = predictor if predictor is not None else resolve_splice_backend(backend)
    network_derived = bool(getattr(pred, "network_derived", False))
    try:
        result = pred.score_sequence(seq)
        pooled = pooled_risk(result, top_k)
        sites = _sites_from_predictor(pred, result, seq, threshold)
    except _DEGRADE_ERRORS as exc:
        return SpliceCrossCheck(
            dna=seq,
            backend=pred.name,
            available=False,
            reason=str(exc),
            calibrated=pred.calibrated,
            network_derived=network_derived,
            threshold=threshold,
            top_k=top_k,
            pooled_risk=0.0,
            sites=(),
        )
    return SpliceCrossCheck(
        dna=seq,
        backend=pred.name,
        available=True,
        reason=None,
        calibrated=pred.calibrated,
        network_derived=network_derived,
        threshold=threshold,
        top_k=top_k,
        pooled_risk=pooled,
        sites=sites,
    )
