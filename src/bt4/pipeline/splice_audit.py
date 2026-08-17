"""Splice audit over a step-3 candidate set (design step 4, pipeline adapter).

:func:`~bt4.biomodels.splice.audit.audit_splice` operates on raw sequences so it
can live in ``biomodels`` (which may import only ``domain``). This thin adapter
bridges it to the step-3 :class:`~bt4.pipeline.candidates.CandidateSet` -- a
``pipeline`` type -- and picks the backends. It composes ``pipeline`` +
``biomodels``, so it belongs here, not in ``biomodels`` (CLAUDE.md §3, strict
acyclic layering).

The audit is **advisory and never edits**: it annotates the candidates with
localized cryptic-splice flags and the whole-panel backend agreement, and returns
them unchanged. Every shipped backend is ``calibrated is False`` today, so
``SpliceAuditReport.all_calibrated`` is ``False`` -- the flags advise, they do not
assert (CLAUDE.md §6, §10.6).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.splice import SpliceResult, pooled_risk, score_in_context
from bt4.biomodels.splice import default as splice_default
from bt4.biomodels.splice.audit import (
    DEFAULT_MATCH_WINDOW,
    DEFAULT_SITE_THRESHOLD,
    SpliceAuditReport,
    audit_splice,
)
from bt4.biomodels.splice.base import DEFAULT_TOP_K, SplicePredictor
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.pangolin import PangolinSplicePredictor
from bt4.biomodels.splice.spliceai import SpliceAiSplicePredictor
from bt4.domain import ConstructContext
from bt4.pipeline.candidates import CandidateSet

__all__ = ["audit_candidate_set", "available_splice_backends"]


def available_splice_backends() -> list[SplicePredictor]:
    """Return the splice backends that can actually run here, most-informative first.

    Always includes the honest PWM baseline; adds the wrapped **Pangolin** and
    **SpliceAI** CNNs only when the user's own install and weights are present
    (their ``available()`` gate). Running *all* available backends is what makes
    the cross-backend agreement signal meaningful -- but every one is
    ``calibrated is False`` today, so the audit stays advisory.
    """
    backends: list[SplicePredictor] = [ConsensusPwmSplicePredictor()]
    for cnn in (PangolinSplicePredictor(), SpliceAiSplicePredictor()):
        if cnn.available():
            backends.append(cnn)
    return backends


@dataclass(frozen=True, slots=True)
class _FlankedPredictor:
    """A splice backend that scores each candidate inside its real flanks.

    Wrapping rather than changing the backends is what keeps this safe: every
    adapter (baseline, SpliceAI, Pangolin, ASSP) is used unmodified, and because
    :func:`~bt4.biomodels.splice.base.score_in_context` slices its result back to
    the coding sequence, every downstream coordinate -- localization, pooling,
    Delta-splicing -- is unchanged.
    """

    inner: SplicePredictor
    upstream: str
    downstream: str

    @property
    def name(self) -> str:
        """The wrapped backend's name (unchanged: this is the same model)."""
        return self.inner.name

    @property
    def calibrated(self) -> bool:
        """The wrapped backend's calibration flag -- better input is not calibration."""
        return self.inner.calibrated

    def score_sequence(self, dna: str) -> SpliceResult:
        """Score ``dna`` in context, returning CDS-aligned per-position scores."""
        return score_in_context(self.inner, dna, self.upstream, self.downstream)

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Added-risk objective with **both** sequences scored in the same flanks.

        Scoring the two in different contexts would conflate a CDS change with a
        context change, so the flanks are held fixed across the pair -- which is
        also what keeps ``delta_splicing(seq, seq) == 0.0`` exactly.
        """
        designed = pooled_risk(self.score_sequence(designed_dna))
        reference = pooled_risk(self.score_sequence(reference_dna))
        # Negated to obey the fixed larger-is-better orientation: positive means
        # the redesign REMOVED risk.
        return reference - designed


def audit_candidate_set(
    candidate_set: CandidateSet,
    *,
    reference: str | None = None,
    predictors: Sequence[SplicePredictor] | None = None,
    threshold: float = DEFAULT_SITE_THRESHOLD,
    match_window: int = DEFAULT_MATCH_WINDOW,
    top_k: int = DEFAULT_TOP_K,
    context: ConstructContext | None = None,
) -> SpliceAuditReport:
    """Localize-and-flag cryptic splice sites across a candidate set (no editing).

    Extracts the candidate sequences from ``candidate_set`` and runs
    :func:`~bt4.biomodels.splice.audit.audit_splice` over them. The sequences are
    returned unchanged -- this is an advisory annotation pass (Stage C /
    ``docs/DESIGN_expression_splice_flow.md`` step 4), never an edit.

    Args:
        candidate_set: A step-3 :class:`~bt4.pipeline.candidates.CandidateSet`.
        reference: The reference each candidate's added risk is measured against.
            Defaults to the delivered (``chosen``) candidate's sequence.
        predictors: Splice backends to run; defaults to the honest baseline only
            (``bt4.biomodels.splice.default()``). Pass
            :func:`available_splice_backends` to include the wrapped SpliceAI /
            Pangolin CNNs when they are installed.
        threshold: Site-localization threshold (a heuristic display knob).
        match_window: +/- nt window for the approximate cross-backend co-occurrence.
        top_k: Pooling depth for ``pooled_risk`` / ``delta_splicing``.
        context: Known sequence around the CDS. When given, each candidate is
            scored inside its real flanks instead of alone -- which matters because
            the wrapped CNNs otherwise pad their wide context with literal ``N``,
            and both pad identically, so their agreement cannot detect the shared
            artifact. Scores are sliced back to the CDS, so every reported position
            still indexes the candidate. Better input is **not** calibration: a
            backend reporting ``calibrated=False`` still does.

    Returns:
        A :class:`~bt4.biomodels.splice.audit.SpliceAuditReport`.

    Raises:
        ValueError: If the candidate set is empty, or per
            :func:`~bt4.biomodels.splice.audit.audit_splice`.
    """
    delivered = candidate_set.delivered()
    if delivered is None:
        raise ValueError("candidate set is empty; nothing to audit")
    dnas = [c.result.dna for c in candidate_set.candidates]
    ref = reference if reference is not None else delivered.result.dna
    preds = list(predictors) if predictors is not None else [splice_default()]
    if context is not None and not context.is_empty:
        preds = [
            _FlankedPredictor(p, context.upstream, context.downstream) for p in preds
        ]
    return audit_splice(
        preds, dnas, ref, threshold=threshold, match_window=match_window, top_k=top_k
    )
