"""The ``SplicePredictor`` contract for cryptic-splice-site risk.

BT4 treats cryptic splicing as a genuinely non-local objective that lives in the
refinement / validation layer (CLAUDE.md sections 6 and 4.3). A splice backend
reduces a coding sequence to *per-position* donor / acceptor site scores, and
reframes the design objective as **Delta-splicing** (CLAUDE.md section 6):

    P(site | designed) - P(site | reference)

i.e. the *added* splice risk a synonymous redesign introduces relative to a
reference. Two honesty rules shape this module, mirroring
:mod:`bt4.biomodels.folding`:

* **``calibrated`` is a first-class flag.** A backend may only claim
  ``calibrated is True`` when its scores are real, calibrated probabilities from
  a validated, hash-pinned model that has passed its held-out acceptance gate --
  a SpliceAI / Pangolin-class per-nucleotide CNN (the slot documented below).
  Anything else -- the consensus / PWM baseline shipped today -- must report
  ``calibrated is False`` and never be presented as a real probability
  (CLAUDE.md sections 6 and 10.6, "no placeholder model presented as a
  feature").
* **Pooling is top-k / log-odds, never saturating noisy-OR.** The per-position
  site scores are pooled into a whole-sequence risk with :func:`pool_log_odds`
  -- a sum of the top-k logits. This is additive: two strong sites contribute
  roughly twice one strong site, so real risk keeps accumulating. BT3's
  saturating noisy-OR (``1 - prod(1 - p_i)``) pegged at 1.0 and hid every site
  after the first (CLAUDE.md section 10.14); BT4 does not repeat it.

**Orientation is fixed and documented.** :meth:`SplicePredictor.delta_splicing`
follows BT4's convention that **larger is better**. The raw Delta-splicing
``P(site|designed) - P(site|reference)`` is *added risk*, where a **negative**
value (the redesign removed splice risk) is the good outcome. So
:meth:`delta_splicing` returns the **negative** of the raw added risk --
``pooled_risk(reference) - pooled_risk(designed)`` -- which is larger exactly
when the designed sequence carries *less* splice risk than the reference.
Maximizing it therefore drives cryptic-splice risk down.

**The calibrated backend slot.** No calibrated model ships yet, so
:func:`bt4.biomodels.splice.default` returns the labeled baseline. A future
SpliceAI / Pangolin-class CNN backend would live beside this module (e.g.
``cnn.py``) behind a lazy import -- exactly as :class:`ViennaFoldingModel`
lazily imports ViennaRNA -- report ``calibrated is True`` only after passing its
held-out-chromosome gate, and be selected by ``default()`` ahead of the
baseline. Until then no such stub exists: an unbuilt model with fake weights
would be the dishonest placeholder CLAUDE.md forbids.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

__all__ = [
    "DEFAULT_SITE_PROBABILITY",
    "DEFAULT_TOP_K",
    "PooledRisk",
    "SplicePredictor",
    "SpliceResult",
    "logit",
    "pool_log_odds",
    "pool_top_k_logit",
    "pooled_risk",
    "pooled_risk_detail",
    "score_in_context",
]

DEFAULT_SITE_PROBABILITY: float = 0.5
"""The probability above which a position counts as a *site* rather than background.

One operating point, referenced by every consumer, because it appears in two
structurally different roles that must not drift apart:

* as the **localization** cutoff -- which positions get flagged and counted
  (``audit.DEFAULT_SITE_THRESHOLD``, ``splice_crosscheck``, ``api.splice_audit`` /
  ``api.splice_crosscheck``);
* as the **background reference** inside :func:`pool_log_odds` -- the probability
  whose log-odds is zero, below which a position contributes no pooled risk.

The second role was previously implicit: pooling took ``max(0.0, logit(p))``, whose
zero crossing is *exactly* ``p = 0.5``, unparameterized and undocumented as a
threshold. That made the two roles silently coupled, so moving the visible cutoff
without moving the hidden one would have made the two disagree -- a flagged site at
``p = 0.35`` contributing exactly zero risk, and a variant introducing five such
sites reporting ``delta_splicing == 0.0``, indistinguishable from introducing none.

**This is a display / localization knob, not a calibrated cutoff.** Every shipped
backend's per-position score is an uncalibrated pseudo-probability, so ``0.5`` is a
convention, not evidence. Published work points elsewhere for the *delta* quantity
(Walker et al., AJHG 2023 calibrate SpliceAI at 0.2, concluding 0.5 "may be
calibrated too high"), but those cutoffs describe a different functional than BT4's
pooled top-k log-odds and cannot be imported directly. Deriving BT4's own operating
point is the job of the statistical-calibration gate
(``docs/DESIGN_splice_cnn_calibration.md`` Part B), and *this* constant is what it
would move.

**Measured consequence, no longer hypothetical.** Run against the hash-verified
Pangolin weights on designed coding sequence -- BT4's own regime -- **no position on
any sequence exceeded 0.5**, across a native CDS and thirty synonymous redesigns of
each of three proteins. Peak scores ran 0.128 to 0.445 and differed more than twofold
between the native and its designs, and every one of them pooled to a risk of exactly
``0.0``. So ``delta_splicing`` was identically zero for every candidate, the rank
agreements computed from those deltas were Spearman correlations of constants, and
none of it was visible in the output.

Lowering this constant is **not** the fix, and must not be done to make the signal
reappear: that is the same uncalibrated knob pointed somewhere more flattering, and
this docstring already says the number is a convention. The fix is to stop reporting a
floored zero as a measured one -- see :class:`PooledRisk`, which carries the counts
that tell the two apart, and :func:`pool_top_k_logit`, the background-free statistic
that still separates sequences when this hinge has flattened them.
"""

DEFAULT_TOP_K: int = 3
"""Default number of strongest sites summed by :func:`pool_log_odds`.

Top-k log-odds pooling counts the ``k`` highest-scoring sites additively rather
than collapsing every site into a single saturating probability. ``k = 3`` lets
several genuine cryptic sites accumulate risk while ignoring the long tail of
background positions.
"""

_LOGIT_EPS: float = 1e-6
"""Clamp applied to probabilities before :func:`logit` to keep it finite."""


def logit(p: float) -> float:
    """Return the log-odds ``ln(p / (1 - p))`` of a probability.

    The input is clamped to ``[_LOGIT_EPS, 1 - _LOGIT_EPS]`` so the result is
    always finite (a calibrated ``p`` of exactly 0 or 1 would otherwise map to
    +/- infinity).

    Args:
        p: A probability-like score. Values outside ``[0, 1]`` are clamped.

    Returns:
        The finite log-odds of ``p``.
    """
    q = min(max(p, _LOGIT_EPS), 1.0 - _LOGIT_EPS)
    return math.log(q / (1.0 - q))


def pool_log_odds(
    probs: Iterable[float],
    top_k: int = DEFAULT_TOP_K,
    *,
    background: float = DEFAULT_SITE_PROBABILITY,
) -> float:
    """Pool per-position site probabilities into one risk via top-k log-odds.

    The positive part of the ``top_k`` largest *background-relative* logits is
    summed. Only positions scoring above ``background`` count as sites and
    contribute risk; the background tail contributes nothing. This pooling is **additive
    and does not saturate**: two strong sites contribute roughly twice the risk
    of one, unlike the noisy-OR aggregation ``1 - prod(1 - p_i)`` that pegs at
    1.0 and masks every site after the first (CLAUDE.md section 10.14).

    Args:
        probs: Per-position site probabilities (donor and/or acceptor). May be
            empty.
        top_k: How many of the highest-scoring sites to sum. Must be positive.
            Fewer than ``top_k`` positive sites simply sum what is present.
        background: The probability treated as background -- positions at or below
            it contribute nothing. Defaults to :data:`DEFAULT_SITE_PROBABILITY`,
            which reproduces the previous ``max(0.0, logit(p))`` behavior exactly
            (``logit(0.5) == 0``). Made explicit so the pooling background and the
            localization cutoff can be moved together rather than drifting apart.

    Returns:
        A non-negative pooled risk: the summed positive part of the ``top_k``
        largest logits, or ``0.0`` when no position is above background (or the
        input is empty). Larger means more pooled splice risk.

    Raises:
        ValueError: If ``top_k`` is not positive.
    """
    if top_k <= 0:
        raise ValueError(f"top_k must be a positive integer, got {top_k}")
    if not 0.0 < background < 1.0:
        raise ValueError(f"background must be in (0, 1), got {background}")
    offset = logit(background)
    logits = sorted((logit(p) - offset for p in probs), reverse=True)
    return math.fsum(max(0.0, value) for value in logits[:top_k])


def pool_top_k_logit(probs: Iterable[float], top_k: int = DEFAULT_TOP_K) -> float:
    """Pool per-position scores into a **background-free** top-k logit sum.

    The same top-k / log-odds shape as :func:`pool_log_odds` with the hinge and the
    background removed: the sum of the ``top_k`` largest ``logit(p)``, whatever
    those are. It takes **no background parameter**, so unlike a pooled *risk* it
    carries no operating point at all.

    **This is not a risk, and it is not a probability.** It is negative whenever the
    strongest positions score below 0.5, which on designed coding sequence is the
    normal case rather than the exception. It has no zero point that means "no risk",
    and it must never be presented as one.

    **What it is for.** It is *monotone* in the per-position scores at every score,
    which pooled risk is not: :func:`pool_log_odds` floors each of its ``top_k``
    contributions at zero, so once every position falls below background the pooled
    risk is a **constant zero** and the differences between sequences are gone.
    Measured on Pangolin over designed CDS, that is not a corner case -- no position
    on any sequence reached 0.5, so every pooled risk and every ``delta_splicing``
    was identically zero while the underlying scores varied by more than twofold.
    This statistic answers the question that survives there: **does the model respond
    to the change at all, and do two backends rank the candidates the same way?**

    Being background-free is what makes it honest to use in that regime. Lowering the
    background to whatever happens to make the signal reappear would be the same
    uncalibrated knob pointed somewhere more flattering; removing it means the
    ranking rests on the model's own scores and on no threshold of BT4's choosing.

    Not comparable across different ``top_k``, and not comparable to
    :func:`pool_log_odds` -- a different functional, reported alongside it, never in
    place of it.

    Args:
        probs: Per-position site probabilities (donor and/or acceptor). May be empty.
        top_k: How many of the highest-scoring positions to sum. Must be positive.

    Returns:
        The sum of the ``top_k`` largest logits, or ``0.0`` for empty input. May be
        negative. Larger means the model scored the sequence's strongest positions
        higher -- **not** that the sequence is riskier in any calibrated sense.

    Raises:
        ValueError: If ``top_k`` is not positive.
    """
    if top_k <= 0:
        raise ValueError(f"top_k must be a positive integer, got {top_k}")
    logits = sorted((logit(p) for p in probs), reverse=True)
    return math.fsum(logits[:top_k])


@dataclass(frozen=True, slots=True)
class SpliceResult:
    """Per-position donor / acceptor site scores for one sequence (immutable).

    The two arrays are aligned to the coding sequence: ``donor[i]`` and
    ``acceptor[i]`` describe the site anchored at nucleotide ``i`` (a backend
    defines its own anchor convention; the baseline anchors on the invariant
    ``G`` of the ``GT`` donor and the ``G`` of the ``AG`` acceptor). Positions
    with insufficient flanking context score ``0.0``.

    Attributes:
        donor: Per-position 5' splice-site (donor) scores, one per nucleotide.
            Calibrated probabilities in ``[0, 1]`` when ``calibrated`` is
            ``True``; uncalibrated pseudo-probabilities otherwise.
        acceptor: Per-position 3' splice-site (acceptor) scores, one per
            nucleotide, with the same calibration caveat as ``donor``.
        model_name: The :attr:`SplicePredictor.name` of the producing backend.
        calibrated: Mirror of :attr:`SplicePredictor.calibrated` for the
            producing backend -- ``False`` marks the scores as uncalibrated.
    """

    donor: tuple[float, ...]
    acceptor: tuple[float, ...]
    model_name: str
    calibrated: bool


def pooled_risk(
    result: SpliceResult,
    top_k: int = DEFAULT_TOP_K,
    *,
    background: float = DEFAULT_SITE_PROBABILITY,
) -> float:
    """Pool a :class:`SpliceResult` into one whole-sequence splice risk.

    Combines the donor and acceptor per-position scores and pools them with
    :func:`pool_log_odds`. Larger means more cryptic-splice risk.

    Args:
        result: Per-position scores from :meth:`SplicePredictor.score_sequence`.
        top_k: How many of the strongest sites to sum (see :func:`pool_log_odds`).
        background: The probability below which a position contributes no risk.
            Threaded through so a caller that moves the localization cutoff can
            move the pooling background with it; see
            :data:`DEFAULT_SITE_PROBABILITY` for why the two must agree.

    Returns:
        The pooled top-k log-odds splice risk of the sequence.
    """
    return pool_log_odds((*result.donor, *result.acceptor), top_k, background=background)


@dataclass(frozen=True, slots=True)
class PooledRisk:
    """A pooled splice risk **together with what would make its value a lie.**

    :func:`pooled_risk` returns a bare float, and that float has two completely
    different meanings when it is ``0.0``:

    * the sequence carries no predicted splice risk, or
    * every position scored at or below :attr:`background`, so the pooling hinge
      floored the whole sequence to zero and the risk is zero **by construction**.

    Nothing distinguished them, and in BT4's own regime the second is the universal
    case. Measured with the hash-verified Pangolin weights on designed coding
    sequence, **no position on any sequence exceeded the 0.5 background** -- so every
    pooled risk was ``0.0``, every ``delta_splicing`` was ``0.0``, and the model's
    scores (which varied more than twofold between the native and the designs) were
    discarded in full without a word. The rank agreements computed over those deltas
    were Spearman correlations of constants.

    This type is the fix: the same number, plus the evidence needed to attribute it.
    A consumer that reports a risk should report :attr:`below_background` with it, and
    must not present a floored zero as a measured one.

    Attributes:
        risk: Exactly :func:`pooled_risk` -- the hinged, background-relative top-k
            log-odds. Unchanged, non-negative, and the only field that is a *risk*.
        response: :func:`pool_top_k_logit` -- the same pooling with no hinge and no
            background. Monotone in the scores everywhere, so it still separates
            sequences that :attr:`risk` has flattened. **Not a risk**: it is negative
            whenever the top positions score below 0.5, and it has no meaningful zero.
        background: The background :attr:`risk` was pooled against.
        top_k: The pooling depth both statistics used.
        n_scored: Positions pooled (donor plus acceptor).
        n_above_background: How many of them exceeded :attr:`background`. Zero is what
            makes :attr:`risk` uninformative rather than merely small.
        max_score: The highest per-position score seen, so a floored zero can be
            reported with the magnitude that was thrown away -- ``0.44`` and ``0.001``
            both pool to ``0.0`` and mean very different things.
    """

    risk: float
    response: float
    background: float
    top_k: int
    n_scored: int
    n_above_background: int
    max_score: float

    @property
    def below_background(self) -> bool:
        """``True`` when :attr:`risk` is zero *by construction*, not by measurement.

        Every scored position fell at or below :attr:`background`, so the hinge in
        :func:`pool_log_odds` floored all of them and no difference between sequences
        can survive into :attr:`risk` or into any Delta-splicing computed from it.
        A consumer seeing this must say so rather than report the zero alone.
        """
        return self.n_scored > 0 and self.n_above_background == 0


def pooled_risk_detail(
    result: SpliceResult,
    top_k: int = DEFAULT_TOP_K,
    *,
    background: float = DEFAULT_SITE_PROBABILITY,
) -> PooledRisk:
    """Pool a :class:`SpliceResult` and report **why** the pooled value is what it is.

    Same pooling as :func:`pooled_risk`, and :attr:`PooledRisk.risk` is equal to it
    exactly -- this adds the attribution, it does not change the number. See
    :class:`PooledRisk` for why a bare pooled risk is not enough.

    Args:
        result: Per-position scores from :meth:`SplicePredictor.score_sequence`.
        top_k: How many of the strongest positions to pool (see :func:`pool_log_odds`).
        background: The probability below which a position contributes no risk.

    Returns:
        A :class:`PooledRisk` carrying the risk, the background-free response, and
        the counts that say whether the risk was floored.

    Raises:
        ValueError: If ``top_k`` is not positive or ``background`` is outside ``(0, 1)``.
    """
    scores = (*result.donor, *result.acceptor)
    return PooledRisk(
        risk=pool_log_odds(scores, top_k, background=background),
        response=pool_top_k_logit(scores, top_k),
        background=background,
        top_k=top_k,
        n_scored=len(scores),
        n_above_background=sum(1 for value in scores if value > background),
        max_score=max(scores, default=0.0),
    )


@runtime_checkable
class SplicePredictor(Protocol):
    """A backend that scores cryptic donor / acceptor splice sites.

    Implementations are swappable behind this contract (CLAUDE.md section 4.3):
    consumers depend only on the protocol and never on a concrete backend. Every
    implementation must satisfy the orientation and honesty rules documented at
    the module level.
    """

    @property
    def name(self) -> str:
        """Stable identifier for the backend (read-only).

        Declared as a read-only property so concrete backends may be frozen
        dataclasses exposing ``name`` as a property.
        """
        ...

    @property
    def calibrated(self) -> bool:
        """Honesty flag: ``True`` only for a validated, hash-pinned model.

        A backend returns ``True`` here **only** when its scores are real,
        calibrated probabilities from a model that has passed its held-out
        acceptance gate (a SpliceAI / Pangolin-class CNN). The consensus / PWM
        baseline returns ``False`` so its pseudo-probabilities are never mistaken
        for calibrated splice probabilities (CLAUDE.md sections 6 and 10.6).
        """
        ...

    def score_sequence(self, dna: str) -> SpliceResult:
        """Return per-position donor / acceptor site scores for ``dna``.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            A :class:`SpliceResult` whose ``donor`` and ``acceptor`` arrays are
            aligned to the sequence. Scores are calibrated probabilities when
            :attr:`calibrated` is ``True`` and uncalibrated pseudo-probabilities
            otherwise.
        """
        ...

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Return the added-splice-risk objective of ``designed`` vs ``reference``.

        Reframes the CLAUDE.md section 6 objective
        ``P(site|designed) - P(site|reference)`` using top-k / log-odds pooling
        (:func:`pool_log_odds`), then **negates it** so the result obeys BT4's
        fixed *larger-is-better* orientation: the value is
        ``pooled_risk(reference) - pooled_risk(designed)``, which is positive
        when the redesign *removed* splice risk and negative when it *added*
        risk. ``delta_splicing(seq, seq)`` is therefore exactly ``0.0``.

        Args:
            designed_dna: The candidate coding sequence.
            reference_dna: The reference (e.g. the input / wild-type) sequence.

        Returns:
            The negated added splice risk; larger means less added risk.
        """
        ...


def score_in_context(
    predictor: SplicePredictor,
    cds: str,
    upstream: str = "",
    downstream: str = "",
) -> SpliceResult:
    """Score ``cds`` **in its real flanking sequence**, aligned back to the CDS.

    Splice models read a wide window around each position, so a coding sequence
    scored on its own is scored in a nucleotide vacuum: the wrapped SpliceAI and
    Pangolin adapters pad their ~10 kb context with literal ``N``, which is a
    documented boundary-artifact regime, and -- because *both* pad identically --
    their cross-backend agreement cannot detect the shared artifact. Real flanking
    sequence, when the user knows it, replaces that padding for the region that
    matters.

    The returned scores are sliced back to the CDS, so ``donor[i]`` still describes
    coding position ``i``. That is what keeps this a drop-in: every downstream
    consumer (localization, pooling, Delta-splicing, the per-site track) keeps its
    coordinates and needs no remapping, and no backend has to change.

    Honest scope: better *input* is not calibration. A backend that reports
    ``calibrated=False`` still does after this -- and a fidelity attestation
    captured on the N-padded path does **not** transfer to a flanked one, because
    it is a different input regime.

    Args:
        predictor: Any backend implementing :class:`SplicePredictor`.
        cds: The coding sequence whose per-position scores are wanted.
        upstream: Known sequence immediately 5' of ``cds`` (may be empty).
        downstream: Known sequence immediately 3' of ``cds`` (may be empty).

    Returns:
        A :class:`SpliceResult` covering exactly ``cds``.
    """
    if not upstream and not downstream:
        return predictor.score_sequence(cds)
    assembled = f"{upstream}{cds}{downstream}".upper()
    scored = predictor.score_sequence(assembled)
    lo = len(upstream)
    hi = lo + len(cds)
    return SpliceResult(
        donor=scored.donor[lo:hi],
        acceptor=scored.acceptor[lo:hi],
        model_name=scored.model_name,
        # An attestation is REGIME-SCOPED, and this is a different regime. The
        # integration-fidelity gate is captured on the bare-CDS path, where the
        # adapter pads its ~10 kb window with literal `N`; a real-flank score is
        # outside what that attestation covers. Forwarding the inner flag here
        # would let better *input* silently masquerade as calibration -- exactly
        # what this function's own docstring forbids -- so a flanked result is
        # always reported uncalibrated, whatever the backend claims.
        calibrated=False,
    )
