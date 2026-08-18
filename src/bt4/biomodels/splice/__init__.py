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
from bt4.biomodels.splice.assp import (
    ASSP_ENDPOINT,
    FIXTURE_DIR_ENV_VAR,
    AsspError,
    AsspReportError,
    AsspSite,
    AsspSplicePredictor,
    AsspTransport,
    AsspUnavailableError,
    CachingAsspTransport,
    FixtureAsspTransport,
    HttpAsspTransport,
    cache_key,
    default_assp_transport,
    parse_assp_report,
)
from bt4.biomodels.splice.attestation import (
    MAX_ATTESTATION_TOLERANCE,
    AttestationError,
    FidelityAttestation,
    attest_backend,
    load_attestation,
    verified_predictor,
)
from bt4.biomodels.splice.attestations import (
    USE_ATTESTED_ENV_VAR,
    attested_promotion_enabled,
    bundled_attestation,
    promote_if_attested,
)
from bt4.biomodels.splice.audit import (
    DEFAULT_MATCH_WINDOW,
    DEFAULT_SITE_THRESHOLD,
    BackendCandidateAudit,
    CandidateSpliceAudit,
    SpliceAuditReport,
    SpliceFlag,
    audit_splice,
)
from bt4.biomodels.splice.base import (
    DEFAULT_SITE_PROBABILITY,
    DEFAULT_TOP_K,
    SplicePredictor,
    SpliceResult,
    logit,
    pool_log_odds,
    pooled_risk,
    score_in_context,
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

# Names carrying the `splice` qualifier, for re-export through `bt4.api` where
# "attestation" alone would be ambiguous against the expression head's.
USE_ATTESTED_SPLICE_ENV_VAR = USE_ATTESTED_ENV_VAR
attested_splice_promotion_enabled = attested_promotion_enabled
bundled_splice_attestation = bundled_attestation

__all__ = [
    "ASSP_ENDPOINT",
    "DEFAULT_MATCH_WINDOW",
    "DEFAULT_SITE_PROBABILITY",
    "DEFAULT_SITE_THRESHOLD",
    "DEFAULT_TISSUES",
    "DEFAULT_TOP_K",
    "FIXTURE_DIR_ENV_VAR",
    "MAX_ATTESTATION_TOLERANCE",
    "USE_ATTESTED_ENV_VAR",
    "USE_ATTESTED_SPLICE_ENV_VAR",
    "AgreementReport",
    "AsspError",
    "AsspReportError",
    "AsspSite",
    "AsspSplicePredictor",
    "AsspTransport",
    "AsspUnavailableError",
    "AttestationError",
    "BackendCandidateAudit",
    "CachingAsspTransport",
    "CandidateSpliceAudit",
    "ConsensusPwmSplicePredictor",
    "FidelityAttestation",
    "FidelityCase",
    "FidelityReport",
    "FixtureAsspTransport",
    "HttpAsspTransport",
    "PangolinSplicePredictor",
    "SpliceAiFidelityCase",
    "SpliceAiFidelityReport",
    "SpliceAiSplicePredictor",
    "SpliceAuditReport",
    "SpliceFlag",
    "SplicePredictor",
    "SpliceResult",
    "attest_backend",
    "attested_promotion_enabled",
    "attested_splice_promotion_enabled",
    "audit_splice",
    "backend_agreement",
    "bundled_attestation",
    "bundled_splice_attestation",
    "cache_key",
    "default",
    "default_assp_transport",
    "load_attestation",
    "logit",
    "parse_assp_report",
    "pool_log_odds",
    "pooled_risk",
    "promote_if_attested",
    "score_in_context",
    "spearman",
    "verified_predictor",
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

        A calibrated CNN backend is obtained explicitly, not here: construct the
        backend with your own installed weights and promote it with
        :func:`verified_predictor` against a passing
        :class:`FidelityAttestation` (see :mod:`bt4.biomodels.splice.attestation`).
        ``default`` needs no per-user weight/tissue configuration and so keeps
        returning the honest baseline.
    """
    return ConsensusPwmSplicePredictor()
