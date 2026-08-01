"""An honestly-labeled, uncalibrated consensus / PWM splice baseline -- NOT SpliceAI.

:class:`ConsensusPwmSplicePredictor` exists so that
:func:`bt4.biomodels.splice.default` can always return a working
:class:`~bt4.biomodels.splice.base.SplicePredictor`, since no calibrated model
ships yet (CLAUDE.md section 4.3: "else a safe baseline; NEVER crashes"). It
scores each position with small **position weight matrices** encoding the
canonical human U2 (GT-AG) splice-site consensus:

* donor / 5' splice site: ``MAG | GT RAGT`` (the invariant ``GT`` opens the
  intron), scored over a 9-nt window;
* acceptor / 3' splice site: a polypyrimidine tract followed by ``... YAG | G``
  (the invariant ``AG`` closes the intron), scored over a 15-nt window.

Honesty (CLAUDE.md sections 6 and 10.6, "no placeholder model presented as a
feature"; section 6, "loudly labels ... consensus-baseline only"):

* :attr:`ConsensusPwmSplicePredictor.calibrated` is ``False`` and its
  :attr:`~ConsensusPwmSplicePredictor.name` is ``"consensus-pwm-baseline"`` --
  the name screams baseline so no consumer can mistake it for a validated model.
* This is **NOT** a validated SpliceAI / Pangolin-class model. Its per-position
  scores are the logistic of a log-odds PWM match -- an uncalibrated
  pseudo-probability, **not** a calibrated probability of splicing. The PWM
  weights encode the *canonical textbook consensus* (graded preferences around
  the invariant ``GT`` / ``AG``); they are illustrative consensus weights, not
  measured frequencies extracted from a held-out annotated dataset. Until a real
  model passes its acceptance gate, these numbers must never be presented as
  calibrated splice probabilities.

The baseline is dependency-free (standard-library ``math`` only) and fully
deterministic: no RNG, no wall-clock, so repeated calls agree exactly
(invariant #7). It exists to give the refinement / validation layer an honest,
GT-AG-aware structural signal and to keep :func:`default` non-crashing.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SpliceResult,
    pooled_risk,
)
from bt4.domain.sequence import validate_dna

__all__ = ["ConsensusPwmSplicePredictor"]

_BACKGROUND: float = 0.25
"""Uniform per-base background probability used for the PWM log-odds ratio."""

_PSEUDOCOUNT: float = 1e-3
"""Additive smoothing so every PWM cell is strictly positive (finite log-odds)."""

# Canonical U2 donor (5' splice site) consensus MAG|GTRAGT over positions
# -3,-2,-1,+1,+2,+3,+4,+5,+6. The invariant GT sits at +1,+2. Values are
# illustrative consensus weights (graded textbook preferences), NOT measured
# frequencies -- see the module docstring's honesty note.
_DONOR_RAW: tuple[dict[str, float], ...] = (
    {"A": 0.34, "C": 0.36, "G": 0.18, "T": 0.12},  # -3  M (A/C)
    {"A": 0.60, "C": 0.13, "G": 0.14, "T": 0.13},  # -2  A
    {"A": 0.09, "C": 0.03, "G": 0.80, "T": 0.08},  # -1  G
    {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},  # +1  G invariant
    {"A": 0.00, "C": 0.00, "G": 0.00, "T": 1.00},  # +2  T invariant
    {"A": 0.52, "C": 0.03, "G": 0.42, "T": 0.03},  # +3  R (A/G)
    {"A": 0.71, "C": 0.08, "G": 0.12, "T": 0.09},  # +4  A
    {"A": 0.07, "C": 0.05, "G": 0.84, "T": 0.04},  # +5  G
    {"A": 0.16, "C": 0.18, "G": 0.20, "T": 0.46},  # +6  T (weak)
)

# Canonical U2 acceptor (3' splice site) consensus: a polypyrimidine tract, then
# ... Y A G | G. The invariant AG closes the intron at window indices 12,13; the
# first exonic base is index 14. Illustrative consensus weights, NOT measured
# frequencies -- see the module docstring's honesty note.
_ACCEPTOR_RAW: tuple[dict[str, float], ...] = (
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},  # polypyrimidine tract
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},
    {"A": 0.14, "C": 0.35, "G": 0.14, "T": 0.37},  # -4
    {"A": 0.16, "C": 0.42, "G": 0.12, "T": 0.30},  # -3  Y (C/T), often C
    {"A": 1.00, "C": 0.00, "G": 0.00, "T": 0.00},  # -2  A invariant
    {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},  # -1  G invariant
    {"A": 0.22, "C": 0.16, "G": 0.50, "T": 0.12},  # +1  exon, G preferred
)

_DONOR_ANCHOR: int = 3
"""Index within the donor window of the anchor base (the invariant ``G`` of GT)."""

_ACCEPTOR_ANCHOR: int = 13
"""Index within the acceptor window of the anchor base (the invariant ``G`` of AG)."""


def _build_pwm(raw: tuple[dict[str, float], ...]) -> tuple[dict[str, float], ...]:
    """Return ``raw`` smoothed by a pseudocount and normalized per column.

    Args:
        raw: Per-column base weights (need not be normalized; may contain zeros).

    Returns:
        A PWM whose every column sums to 1.0 with strictly positive cells, so
        the log-odds ratio against :data:`_BACKGROUND` is always finite.
    """
    pwm: list[dict[str, float]] = []
    for col in raw:
        smoothed = {base: col.get(base, 0.0) + _PSEUDOCOUNT for base in "ACGT"}
        total = math.fsum(smoothed.values())
        pwm.append({base: value / total for base, value in smoothed.items()})
    return tuple(pwm)


_DONOR_PWM: tuple[dict[str, float], ...] = _build_pwm(_DONOR_RAW)
_ACCEPTOR_PWM: tuple[dict[str, float], ...] = _build_pwm(_ACCEPTOR_RAW)


def _sigmoid(x: float) -> float:
    """Return the logistic ``1 / (1 + e^-x)`` without overflow for large ``|x|``."""
    if x >= 0.0:
        return 1.0 / (1.0 + math.exp(-x))
    ex = math.exp(x)
    return ex / (1.0 + ex)


def _window_score(window: str, pwm: tuple[dict[str, float], ...]) -> float:
    """Score one fixed-length window against a PWM, as a pseudo-probability.

    Computes the summed log-odds ``sum_k ln(pwm[k][base] / background)`` and maps
    it through the logistic. A perfect consensus match scores near 1.0; a window
    lacking the invariant ``GT`` / ``AG`` scores near 0.0.

    Args:
        window: A ``{A,C,G,T}`` string whose length equals ``len(pwm)``.
        pwm: A normalized PWM from :func:`_build_pwm`.

    Returns:
        An uncalibrated pseudo-probability in ``(0, 1)``.
    """
    log_odds = math.fsum(math.log(pwm[k][base] / _BACKGROUND) for k, base in enumerate(window))
    return _sigmoid(log_odds)


@dataclass(frozen=True, slots=True)
class ConsensusPwmSplicePredictor:
    """Uncalibrated consensus / PWM splice baseline (NOT a validated model).

    A safe, dependency-free fallback
    :class:`~bt4.biomodels.splice.base.SplicePredictor`. It scores each position
    with the canonical GT-AG donor and acceptor PWMs; the scores are
    uncalibrated pseudo-probabilities, **never** calibrated splice probabilities
    (see the module docstring).

    Attributes:
        top_k: Number of strongest sites summed by :meth:`delta_splicing`'s
            top-k / log-odds pooling. Defaults to
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`.
    """

    top_k: int = DEFAULT_TOP_K

    def __post_init__(self) -> None:
        """Validate the pooling depth.

        Raises:
            ValueError: If ``top_k`` is not a positive integer.
        """
        if self.top_k <= 0:
            raise ValueError(f"top_k must be a positive integer, got {self.top_k}")

    @property
    def name(self) -> str:
        """Backend identifier -- deliberately screams that this is a baseline."""
        return "consensus-pwm-baseline"

    @property
    def calibrated(self) -> bool:
        """Always ``False``: PWM pseudo-probabilities are not calibrated (honesty)."""
        return False

    def score_sequence(self, dna: str) -> SpliceResult:
        """Score every position's donor and acceptor consensus match.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            A :class:`~bt4.biomodels.splice.base.SpliceResult` whose ``donor[i]``
            is the donor PWM match anchored at the ``GT`` ``G`` at position ``i``
            and whose ``acceptor[i]`` is the acceptor PWM match anchored at the
            ``AG`` ``G`` at position ``i``. Positions without full flanking
            context score ``0.0``. Scores are uncalibrated pseudo-probabilities.

        Raises:
            ValueError: If ``dna`` is empty or contains non-ACGT characters.
        """
        seq = validate_dna(dna)
        n = len(seq)
        donor = [0.0] * n
        acceptor = [0.0] * n

        donor_len = len(_DONOR_PWM)
        for start in range(n - donor_len + 1):
            anchor = start + _DONOR_ANCHOR
            donor[anchor] = _window_score(seq[start : start + donor_len], _DONOR_PWM)

        acceptor_len = len(_ACCEPTOR_PWM)
        for start in range(n - acceptor_len + 1):
            anchor = start + _ACCEPTOR_ANCHOR
            acceptor[anchor] = _window_score(seq[start : start + acceptor_len], _ACCEPTOR_PWM)

        return SpliceResult(
            donor=tuple(donor),
            acceptor=tuple(acceptor),
            model_name=self.name,
            calibrated=False,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Return the negated added splice risk of ``designed`` vs ``reference``.

        See :meth:`bt4.biomodels.splice.base.SplicePredictor.delta_splicing` for
        the fixed orientation and pooling. Concretely this returns
        ``pooled_risk(reference) - pooled_risk(designed)`` using top-k / log-odds
        pooling, so it is ``0.0`` for identical sequences, positive when the
        redesign lowers splice risk, and negative when it raises it.

        Args:
            designed_dna: The candidate coding sequence.
            reference_dna: The reference (e.g. input / wild-type) sequence.

        Returns:
            The negated added splice risk (uncalibrated); larger is better.

        Raises:
            ValueError: If either sequence is empty or contains non-ACGT
                characters.
        """
        designed_risk = pooled_risk(self.score_sequence(designed_dna), self.top_k)
        reference_risk = pooled_risk(self.score_sequence(reference_dna), self.top_k)
        return reference_risk - designed_risk
