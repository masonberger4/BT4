"""Localize-and-flag cryptic-splice audit over a candidate set (design step 4).

This is the out-of-loop splice audit of the expression/splice design flow
(``docs/DESIGN_expression_splice_flow.md`` Stage C / step 4). Given a set of
candidate coding sequences and one reference, it runs each available
:class:`~bt4.biomodels.splice.base.SplicePredictor` over every candidate and:

* **localizes** residual cryptic splice sites the in-loop motif constraint missed
  -- one flag per contiguous above-threshold run, at its peak (non-maximal
  suppression), each with its position, kind, score, and the *intra-backend*
  added risk vs the reference at that position;
* attaches the whole-panel **backend agreement** (pooled rank + sign, via
  :func:`~bt4.biomodels.splice.agreement.backend_agreement`) as the authoritative
  cross-backend confidence signal.

**It annotates; it never edits.** The candidate sequences are returned unchanged
-- the audit is advisory only. A targeted synonymous *auto-edit* at flagged loci
is a deliberately deferred, calibrated-gated future step (step 6): it will unlock
only when a backend passes its own integration-fidelity gate
(``calibrated is True``), and even then it may claim only *reduced predicted
cryptic-splice risk*, never expression gain. Today **every shipped backend is
``calibrated is False``** (the PWM baseline structurally; the wrapped SpliceAI /
Pangolin CNNs until their BT4 fidelity gate is recorded), so:

* :attr:`SpliceAuditReport.all_calibrated` is ``False`` and every
  :class:`SpliceFlag` carries its emitting backend's ``calibrated`` flag -- a flag
  is an advisory annotation, **never** a calibrated risk assertion (CLAUDE.md §6,
  §10.6);
* the site ``threshold`` is a **display / localization knob, heuristic for every
  backend** whose ``calibrated`` is ``False`` -- not a validated cutoff. The
  baseline's per-position ``score`` is an uncalibrated PWM pseudo-score in
  *arbitrary units* (never a probability); the CNNs' are real upstream
  probabilities but uncalibrated *in BT4's ledger*.

**Cross-backend honesty.** The backends disagree on convention -- the baseline
separates donor / acceptor, Pangolin reports one combined ``P(splice)`` (so it can
never *disagree* on kind), SpliceAI splits cleanly -- so per-site cross-backend
matching is inherently approximate. Each flag's :attr:`SpliceFlag.also_flagged_by`
is therefore a **raw positional co-occurrence** (other backends that localized
*any* site within ``match_window`` nt), explicitly **not** a kind-level agreement
and **not** a probability. The authoritative cross-backend signal is the pooled
:class:`~bt4.biomodels.splice.agreement.AgreementReport` (rank + sign of
Delta-splicing), which is backend-agnostic by construction.

Determinism (#7): pure, seedless, no wall-clock -- identical inputs yield
identical reports. This module depends only on :mod:`bt4.domain` and the splice
backends; it never imports ``pipeline`` (the ``CandidateSet`` adapter lives in
:mod:`bt4.pipeline.splice_audit`).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.agreement import AgreementReport, agreement_from_deltas
from bt4.biomodels.splice.base import (
    DEFAULT_SITE_PROBABILITY,
    DEFAULT_TOP_K,
    SplicePredictor,
    SpliceResult,
    pooled_risk_detail,
)

__all__ = [
    "DEFAULT_MATCH_WINDOW",
    "DEFAULT_SITE_THRESHOLD",
    "BackendCandidateAudit",
    "CandidateSpliceAudit",
    "SpliceAuditReport",
    "SpliceFlag",
    "audit_splice",
]

DEFAULT_SITE_THRESHOLD: float = DEFAULT_SITE_PROBABILITY
"""Default per-position score above which a site is localized.

A **display / localization knob**, not a calibrated cutoff: every shipped backend
is ``calibrated is False``, so this threshold is heuristic for all of them (and
doubly so for the PWM baseline, whose scores are arbitrary-units pseudo-scores).

Aliases :data:`~bt4.biomodels.splice.base.DEFAULT_SITE_PROBABILITY` rather than
repeating the literal, because the same operating point also sets the *pooling
background* inside :func:`~bt4.biomodels.splice.base.pool_log_odds`. Two copies of
``0.5`` could be moved independently, and then a site flagged here would contribute
zero pooled risk there.
"""

DEFAULT_MATCH_WINDOW: int = 3
"""Default +/- nt window for the approximate cross-backend co-occurrence.

Sized to absorb the backends' differing anchor conventions (the baseline anchors
the donor on the ``GT`` ``G`` and the acceptor on the ``AG`` ``G``; SpliceAI /
Pangolin anchor at the junction). Matching within this window is deliberately
approximate -- the authoritative cross-backend signal is the pooled agreement.
"""


@dataclass(frozen=True, slots=True)
class SpliceFlag:
    """One localized cryptic splice site flagged by one backend on one candidate.

    Attributes:
        position: The site's anchor index in the coding sequence (the peak of a
            contiguous above-threshold run).
        kind: ``"donor"`` or ``"acceptor"`` for a backend that separates the two,
            or ``"splice"`` for a combined-track backend (e.g. Pangolin, which
            reports one ``P(splice)`` and cannot distinguish donor from acceptor).
        score: The backend's per-position score at ``position``. An uncalibrated
            PWM pseudo-score in *arbitrary units* for the baseline; a real (but
            BT4-uncalibrated) probability for the CNNs. Never read as a calibrated
            splice probability when :attr:`calibrated` is ``False``.
        added_risk_vs_reference: ``score - reference_score`` at the same position
            **for the same backend** -- *positive means the redesign added risk*
            here (the opposite sign to the panel-level ``delta_splicing``, which is
            larger-is-better). Strictly intra-backend: cross-backend position
            deltas are meaningless given the anchor mismatch.
        backend: The emitting backend's name.
        calibrated: The emitting backend's ``calibrated`` flag (``False`` for every
            shipped backend today) -- an advisory annotation, not a validated risk.
        also_flagged_by: Names of *other* backends that localized any site within
            ``match_window`` nt of ``position`` on this same candidate. A raw
            **positional co-occurrence**, explicitly not a kind-level agreement or a
            probability (see the module docstring).
    """

    position: int
    kind: str
    score: float
    added_risk_vs_reference: float
    backend: str
    calibrated: bool
    also_flagged_by: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BackendCandidateAudit:
    """One backend's audit of one candidate sequence.

    Attributes:
        backend: The backend's name.
        calibrated: The backend's ``calibrated`` flag (``False`` today).
        flags: The localized splice-site flags, ordered by position.
        pooled_risk: The whole-sequence pooled top-k splice risk (larger = more
            predicted risk), from :func:`~bt4.biomodels.splice.base.pooled_risk`.
        delta_splicing: ``pooled_risk(reference) - pooled_risk(candidate)`` for this
            backend -- **larger is better** (positive = the redesign lowered pooled
            risk vs the reference). Note the opposite sign to each flag's
            :attr:`SpliceFlag.added_risk_vs_reference`. Read it with
            :attr:`risk_floored`: on designed coding sequence a zero here is usually
            not a measurement.
        risk_floored: ``True`` when **no** position of this candidate exceeded the
            pooling background, so :attr:`pooled_risk` is zero *by construction* and
            :attr:`delta_splicing` carries no information about this candidate.
            Measured against the hash-verified Pangolin weights on designed CDS, this
            was true of every sequence -- peak scores of 0.128 to 0.445 all pooled to
            zero. A consumer must not present that zero as "no predicted risk".
        max_site_score: The highest per-position score on this candidate. What
            :attr:`risk_floored` discarded, so a floored zero can be reported with the
            magnitude behind it -- ``0.44`` and ``0.001`` pool identically and mean
            entirely different things.
    """

    backend: str
    calibrated: bool
    flags: tuple[SpliceFlag, ...]
    pooled_risk: float
    delta_splicing: float
    risk_floored: bool = False
    max_site_score: float = 0.0


@dataclass(frozen=True, slots=True)
class CandidateSpliceAudit:
    """The per-backend splice audit of one candidate sequence.

    Attributes:
        index: The candidate's index in the audited set.
        dna: The candidate coding sequence (returned unchanged -- no editing).
        by_backend: One :class:`BackendCandidateAudit` per backend, in input order.
    """

    index: int
    dna: str
    by_backend: tuple[BackendCandidateAudit, ...]


@dataclass(frozen=True, slots=True)
class SpliceAuditReport:
    """The whole localize-and-flag audit over a candidate set.

    Attributes:
        backends: The backend names, in input order.
        all_calibrated: ``True`` only if *every* backend is calibrated -- ``False``
            today, marking the entire report advisory (annotate, never assert).
        threshold: The site-localization threshold used (a heuristic display knob).
        match_window: The +/- nt window used for the approximate cross-backend
            co-occurrence in :attr:`SpliceFlag.also_flagged_by`.
        top_k: The pooling depth used for ``pooled_risk`` / ``delta_splicing``.
        candidates: The per-candidate audits, in input order.
        agreement: The whole-panel pooled backend agreement (rank + sign of
            Delta-splicing) -- the authoritative, backend-agnostic cross-backend
            confidence signal.
    """

    backends: tuple[str, ...]
    all_calibrated: bool
    threshold: float
    match_window: int
    top_k: int
    candidates: tuple[CandidateSpliceAudit, ...]
    agreement: AgreementReport


def _localize(scores: Sequence[float], kind: str, threshold: float) -> list[tuple[int, str, float]]:
    """Localize sites as the peak of each contiguous above-threshold run (NMS).

    A broad site spanning several positions yields exactly one flag -- the argmax
    of its run (the first index on ties, for determinism) -- rather than a cluster
    of adjacent flags that would inflate counts and cross-backend co-occurrence.
    """
    flags: list[tuple[int, str, float]] = []
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
            flags.append((best, kind, scores[best]))
            i = j
        else:
            i += 1
    return flags


def _reference_score(ref: SpliceResult, position: int, kind: str) -> float:
    """Return the same-backend reference score at ``position`` for ``kind`` (0 if OOB)."""
    track = ref.acceptor if kind == "acceptor" else ref.donor
    return track[position] if 0 <= position < len(track) else 0.0


def _raw_flags(
    result: SpliceResult, combined: bool, threshold: float
) -> list[tuple[int, str, float]]:
    """Localize a SpliceResult's donor/acceptor tracks into raw ``(pos, kind, score)``.

    ``combined`` is a fixed property of the *backend* (see
    :func:`_backend_is_combined`), not inferred from this one sequence: a combined
    backend (e.g. Pangolin, one ``P(splice)`` in ``donor`` with ``acceptor``
    all-zero) has its donor-track flags labelled ``"splice"`` -- never claimed as
    donor-specific -- and its (empty) acceptor track is not localized.
    """
    donor_kind = "splice" if combined else "donor"
    raw = _localize(result.donor, donor_kind, threshold)
    if not combined:
        raw += _localize(result.acceptor, "acceptor", threshold)
    raw.sort(key=lambda f: (f[0], f[1]))
    return raw


def _backend_is_combined(results: Sequence[SpliceResult]) -> bool:
    """Return whether a backend reports a *combined* track (acceptor always all-zero).

    Combined-vs-separated is a property of the backend, so it is decided across
    **all** of a backend's scored sequences (the reference and every candidate),
    not per-sequence: a separated backend is treated as combined only if its
    acceptor track is entirely zero on *every* sequence -- so it is never misread as
    combined just because one short or degenerate sequence produced no acceptor
    sites (whereas Pangolin's acceptor track is all-zero by construction, always).
    """
    return not any(any(res.acceptor) for res in results)


def audit_splice(
    predictors: Sequence[SplicePredictor],
    candidates: Sequence[str],
    reference: str,
    *,
    threshold: float = DEFAULT_SITE_THRESHOLD,
    match_window: int = DEFAULT_MATCH_WINDOW,
    top_k: int = DEFAULT_TOP_K,
) -> SpliceAuditReport:
    """Localize-and-flag cryptic splice sites across ``candidates`` (no editing).

    For each backend and candidate, localizes above-``threshold`` sites (peak per
    run), records each with its intra-backend added risk vs ``reference``, and
    annotates cross-backend positional co-occurrence within ``match_window``. Also
    attaches the whole-panel pooled :func:`backend_agreement`. The sequences are
    never modified -- the audit is advisory (see the module docstring).

    Args:
        predictors: The splice backends to run. Names must be distinct. Every
            shipped backend is ``calibrated is False`` today, so the report is
            advisory.
        candidates: Candidate coding sequences (e.g. a step-3 candidate set).
        reference: The reference sequence each candidate's added risk is measured
            against (e.g. the delivered/wild-type sequence). All candidates should
            encode the same protein as the reference so positions align.
        threshold: Site-localization threshold (a heuristic display knob).
        match_window: +/- nt window for the approximate cross-backend co-occurrence.
        top_k: Pooling depth for ``pooled_risk`` / ``delta_splicing``.

    Returns:
        A :class:`SpliceAuditReport`.

    Raises:
        ValueError: If ``predictors`` or ``candidates`` is empty, two predictors
            share a name, ``match_window`` is negative, or ``top_k`` is not
            positive.
    """
    if not predictors:
        raise ValueError("need at least one splice backend to audit")
    if not candidates:
        raise ValueError("need at least one candidate sequence to audit")
    names = tuple(p.name for p in predictors)
    if len(set(names)) != len(names):
        raise ValueError(f"backend names must be distinct, got {names}")
    if match_window < 0:
        raise ValueError(f"match_window must be >= 0, got {match_window}")
    if top_k <= 0:
        raise ValueError(f"top_k must be a positive integer, got {top_k}")

    # Score the reference once per backend (reused for every candidate's added-risk
    # and pooled delta -- never re-run per candidate, CLAUDE.md §7).
    ref_results = {p.name: p.score_sequence(reference) for p in predictors}
    ref_pooled = {
        name: pooled_risk_detail(res, top_k, background=threshold)
        for name, res in ref_results.items()
    }
    # Score every candidate once per backend, retaining the results so combined-vs-
    # separated can be decided per backend across the whole panel (not per sequence).
    cand_results: list[dict[str, SpliceResult]] = [
        {p.name: p.score_sequence(cand) for p in predictors} for cand in candidates
    ]
    combined = {
        name: _backend_is_combined([ref_results[name], *(cr[name] for cr in cand_results)])
        for name in names
    }

    audits: list[CandidateSpliceAudit] = []
    # Accumulate each backend's Delta-splicing vector as we go, so the whole-panel
    # agreement is built from these *already-computed* values rather than re-running
    # every backend over the candidate set a second time (CLAUDE.md §7 -- doubling a
    # ~10 kb-context CNN's forward passes would defeat the point of the out-of-loop
    # audit).
    delta_by_backend: dict[str, list[float]] = {name: [] for name in names}
    for index in range(len(candidates)):
        # First pass: raw flags per backend (positions only), so cross-backend
        # co-occurrence can be computed before the frozen flags are built.
        results = cand_results[index]
        raw_by_backend = {
            name: _raw_flags(results[name], combined[name], threshold) for name in names
        }

        by_backend: list[BackendCandidateAudit] = []
        for predictor in predictors:
            name = predictor.name
            ref_res = ref_results[name]
            flags: list[SpliceFlag] = []
            for pos, kind, score in raw_by_backend[name]:
                also = tuple(
                    other
                    for other in names
                    if other != name
                    and any(abs(pos - opos) <= match_window for opos, _, _ in raw_by_backend[other])
                )
                flags.append(
                    SpliceFlag(
                        position=pos,
                        kind=kind,
                        score=score,
                        added_risk_vs_reference=score - _reference_score(ref_res, pos, kind),
                        backend=name,
                        calibrated=predictor.calibrated,
                        also_flagged_by=also,
                    )
                )
            # `threshold` is passed as the pooling background deliberately: it is the
            # same operating point the flags above are localized at, and the audit
            # would otherwise flag a site at the caller's threshold while pooling it
            # against a different one. `pooled_risk_detail` returns the same number
            # `pooled_risk` does, plus what says whether that number is a measurement.
            cand_pooled = pooled_risk_detail(results[name], top_k, background=threshold)
            delta = ref_pooled[name].risk - cand_pooled.risk
            delta_by_backend[name].append(delta)
            by_backend.append(
                BackendCandidateAudit(
                    backend=name,
                    calibrated=predictor.calibrated,
                    flags=tuple(flags),
                    pooled_risk=cand_pooled.risk,
                    delta_splicing=delta,
                    risk_floored=cand_pooled.below_background,
                    max_site_score=cand_pooled.max_score,
                )
            )
        audits.append(
            CandidateSpliceAudit(index=index, dna=candidates[index], by_backend=tuple(by_backend))
        )

    # Reuse the deltas already computed above (no second scoring pass); this equals
    # backend_agreement(predictors, candidates, reference, top_k=top_k).
    agreement = agreement_from_deltas({name: tuple(delta_by_backend[name]) for name in names})
    return SpliceAuditReport(
        backends=names,
        all_calibrated=all(p.calibrated for p in predictors),
        threshold=threshold,
        match_window=match_window,
        top_k=top_k,
        candidates=tuple(audits),
        agreement=agreement,
    )
