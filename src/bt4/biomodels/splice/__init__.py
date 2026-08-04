"""Cryptic-splice-risk models behind the ``SplicePredictor`` contract.

This package provides the splice half of BT4's non-local biology (CLAUDE.md
sections 6 and 4.3): a swappable :class:`SplicePredictor` contract, its top-k /
log-odds pooling primitives, an honestly-labeled, dependency-free consensus / PWM
baseline (:class:`ConsensusPwmSplicePredictor`), a wrapped published **Pangolin**
CNN backend (:class:`PangolinSplicePredictor`, lazily importing the user's own
GPL ``pangolin`` install), and a cross-backend agreement report
(:func:`backend_agreement`).

No *calibrated* SpliceAI / Pangolin-class model ships yet, so :func:`default`
returns the labeled baseline. The Pangolin backend is wrapped but reports
``calibrated is False`` until it passes its integration-fidelity gate
(:func:`verify_pangolin_fidelity`), so ``default`` keeps returning the honest
baseline -- exactly as :func:`bt4.biomodels.folding.default` returns its baseline
until the calibrated ViennaRNA backend is available. No fake-weight stub ships in
the meantime (CLAUDE.md sections 6 and 10.6).

Importing this package stays lightweight: the Pangolin adapter imports ``torch``
and the GPL ``pangolin`` package **only inside its methods**, never at load
(CLAUDE.md section 3). This package therefore depends only on :mod:`bt4.domain`
and the standard library at import time.
"""

from __future__ import annotations

from bt4.biomodels.splice.agreement import AgreementReport, backend_agreement, spearman
from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SplicePredictor,
    SpliceResult,
    logit,
    pool_log_odds,
    pooled_risk,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.pangolin import (
    DEFAULT_TISSUES,
    FidelityCase,
    FidelityReport,
    PangolinSplicePredictor,
    verify_pangolin_fidelity,
)

__all__ = [
    "DEFAULT_TISSUES",
    "DEFAULT_TOP_K",
    "AgreementReport",
    "ConsensusPwmSplicePredictor",
    "FidelityCase",
    "FidelityReport",
    "PangolinSplicePredictor",
    "SplicePredictor",
    "SpliceResult",
    "backend_agreement",
    "default",
    "logit",
    "pool_log_odds",
    "pooled_risk",
    "spearman",
    "verify_pangolin_fidelity",
]


def default() -> SplicePredictor:
    """Return the best available *calibrated* splice predictor, never crashing.

    Returns:
        An uncalibrated :class:`ConsensusPwmSplicePredictor`. The Pangolin backend
        is wrapped (:class:`PangolinSplicePredictor`) but reports
        ``calibrated is False`` until its integration-fidelity gate passes
        (:func:`verify_pangolin_fidelity`), and no calibrated model ships, so the
        honestly-labeled baseline is the default today. It is labeled
        ``calibrated is False`` so its pseudo-probabilities are never mistaken for
        calibrated splice probabilities (CLAUDE.md sections 4.3, 6, and 10.6).
        When a calibrated backend lands, this function will prefer it.
    """
    return ConsensusPwmSplicePredictor()
