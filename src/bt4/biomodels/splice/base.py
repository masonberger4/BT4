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
    "DEFAULT_TOP_K",
    "SplicePredictor",
    "SpliceResult",
    "logit",
    "pool_log_odds",
    "pooled_risk",
]

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


def pool_log_odds(probs: Iterable[float], top_k: int = DEFAULT_TOP_K) -> float:
    """Pool per-position site probabilities into one risk via top-k log-odds.

    The positive part of the ``top_k`` largest logits is summed. Only
    above-background positions (log-odds > 0) count as sites and contribute
    risk; the background tail contributes nothing. This pooling is **additive
    and does not saturate**: two strong sites contribute roughly twice the risk
    of one, unlike the noisy-OR aggregation ``1 - prod(1 - p_i)`` that pegs at
    1.0 and masks every site after the first (CLAUDE.md section 10.14).

    Args:
        probs: Per-position site probabilities (donor and/or acceptor). May be
            empty.
        top_k: How many of the highest-scoring sites to sum. Must be positive.
            Fewer than ``top_k`` positive sites simply sum what is present.

    Returns:
        A non-negative pooled risk: the summed positive part of the ``top_k``
        largest logits, or ``0.0`` when no position is above background (or the
        input is empty). Larger means more pooled splice risk.

    Raises:
        ValueError: If ``top_k`` is not positive.
    """
    if top_k <= 0:
        raise ValueError(f"top_k must be a positive integer, got {top_k}")
    logits = sorted((logit(p) for p in probs), reverse=True)
    return math.fsum(max(0.0, value) for value in logits[:top_k])


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


def pooled_risk(result: SpliceResult, top_k: int = DEFAULT_TOP_K) -> float:
    """Pool a :class:`SpliceResult` into one whole-sequence splice risk.

    Combines the donor and acceptor per-position scores and pools them with
    :func:`pool_log_odds`. Larger means more cryptic-splice risk.

    Args:
        result: Per-position scores from :meth:`SplicePredictor.score_sequence`.
        top_k: How many of the strongest sites to sum (see :func:`pool_log_odds`).

    Returns:
        The pooled top-k log-odds splice risk of the sequence.
    """
    return pool_log_odds((*result.donor, *result.acceptor), top_k)


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
