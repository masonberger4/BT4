"""Cryptic-splice-risk models behind the ``SplicePredictor`` contract.

This package provides the splice half of BT4's non-local biology (CLAUDE.md
sections 6 and 4.3): a swappable :class:`SplicePredictor` contract, its top-k /
log-odds pooling primitives, and an honestly-labeled, dependency-free consensus
/ PWM baseline (:class:`ConsensusPwmSplicePredictor`).

No calibrated SpliceAI / Pangolin-class model ships yet, so :func:`default`
returns the labeled baseline. When a validated CNN backend lands (behind a lazy
import, reporting ``calibrated is True`` only after passing its held-out gate),
``default`` will prefer it -- exactly as :func:`bt4.biomodels.folding.default`
prefers the calibrated ViennaRNA backend over its baseline. The slot is
documented in :mod:`bt4.biomodels.splice.base`; no fake-weight stub ships in the
meantime (CLAUDE.md sections 6 and 10.6).

This package depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SplicePredictor,
    SpliceResult,
    logit,
    pool_log_odds,
    pooled_risk,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor

__all__ = [
    "DEFAULT_TOP_K",
    "ConsensusPwmSplicePredictor",
    "SplicePredictor",
    "SpliceResult",
    "default",
    "logit",
    "pool_log_odds",
    "pooled_risk",
]


def default() -> SplicePredictor:
    """Return the best available splice predictor, never crashing.

    Returns:
        An uncalibrated :class:`ConsensusPwmSplicePredictor`. No calibrated
        SpliceAI / Pangolin-class model ships yet, so the baseline is the only
        backend today; it is honestly labeled (``calibrated is False``) so its
        pseudo-probabilities are never mistaken for calibrated splice
        probabilities (CLAUDE.md sections 4.3, 6, and 10.6). When a validated CNN
        backend lands, this function will prefer it.
    """
    return ConsensusPwmSplicePredictor()
