"""Cryptic-splice-risk models behind the ``SplicePredictor`` contract.

This package provides the splice half of BT4's non-local biology (CLAUDE.md
sections 6 and 4.3): a swappable :class:`SplicePredictor` contract, its top-k /
log-odds pooling primitives, an honestly-labeled, dependency-free consensus / PWM
baseline (:class:`ConsensusPwmSplicePredictor`), two wrapped published CNN
backends -- **Pangolin** (:class:`PangolinSplicePredictor`, lazily importing the
user's own GPL ``pangolin`` install) and **SpliceAI**
(:class:`SpliceAiSplicePredictor`, lazily importing the user's own CC BY-NC
``spliceai`` install) -- and a cross-backend agreement report
(:func:`backend_agreement`) that turns running both into an uncertainty signal.

No *calibrated* SpliceAI / Pangolin-class model ships yet, so :func:`default`
returns the labeled baseline. Both CNN backends are wrapped but report
``calibrated is False`` until they pass their integration-fidelity gates
(:func:`verify_pangolin_fidelity` / :func:`verify_spliceai_fidelity`), so
``default`` keeps returning the honest baseline -- exactly as
:func:`bt4.biomodels.folding.default` returns its baseline until the calibrated
ViennaRNA backend is available. No fake-weight stub ships in the meantime
(CLAUDE.md sections 6 and 10.6).

Importing this package stays lightweight: the CNN adapters import ``torch`` /
``tensorflow`` and the GPL ``pangolin`` / CC BY-NC ``spliceai`` packages **only
inside their methods**, never at load (CLAUDE.md section 3). This package
therefore depends only on :mod:`bt4.domain` and the standard library at import
time.
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
from bt4.biomodels.splice.spliceai import (
    SpliceAiFidelityCase,
    SpliceAiFidelityReport,
    SpliceAiSplicePredictor,
    verify_spliceai_fidelity,
)

__all__ = [
    "DEFAULT_TISSUES",
    "DEFAULT_TOP_K",
    "AgreementReport",
    "ConsensusPwmSplicePredictor",
    "FidelityCase",
    "FidelityReport",
    "PangolinSplicePredictor",
    "SpliceAiFidelityCase",
    "SpliceAiFidelityReport",
    "SpliceAiSplicePredictor",
    "SplicePredictor",
    "SpliceResult",
    "backend_agreement",
    "default",
    "logit",
    "pool_log_odds",
    "pooled_risk",
    "spearman",
    "verify_pangolin_fidelity",
    "verify_spliceai_fidelity",
]


def default() -> SplicePredictor:
    """Return the best available *calibrated* splice predictor, never crashing.

    Returns:
        An uncalibrated :class:`ConsensusPwmSplicePredictor`. The Pangolin and
        SpliceAI backends are wrapped (:class:`PangolinSplicePredictor` /
        :class:`SpliceAiSplicePredictor`) but report ``calibrated is False`` until
        their integration-fidelity gates pass (:func:`verify_pangolin_fidelity` /
        :func:`verify_spliceai_fidelity`), and no calibrated model ships, so the
        honestly-labeled baseline is the default today. It is labeled
        ``calibrated is False`` so its pseudo-probabilities are never mistaken for
        calibrated splice probabilities (CLAUDE.md sections 4.3, 6, and 10.6).
        When a calibrated backend lands, this function will prefer it.
    """
    return ConsensusPwmSplicePredictor()
