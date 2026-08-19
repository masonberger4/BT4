"""Tests for the committed Pangolin attestation and the opt-in promotion seam.

A fidelity attestation now ships in the tree: BT4's Pangolin adapter was shown to
reproduce upstream's per-position scores **exactly** (max abs deviation 0.0) over
an 18-case panel, on a maintainer machine holding the GPL weights. This is the
first wrapped model in BT4 to pass its gate.

These tests pin what that does and does not change.

**Does:** with an explicit opt-in, a Pangolin predictor built by BT4's own
pipeline reports ``calibrated=True``.

**Does not:** anything by default. The attestation records that BT4's *wrapper* is
faithful, not that the model's scores are calibrated probabilities for designed
coding sequence -- a separate, unmet gate, and one where these models are measured
weakest (median prAUC 0.419 exonic vs 0.773 intronic, Smith & Kitzman 2023). So
promotion stays opt-in, ``default()`` still returns the PWM baseline, and a
real-flank score still reports uncalibrated even from a promoted backend.

None of this needs torch or the licensed weights: the attestation is eight scalars
plus public weight hashes, and promotion only sets a flag.
"""

from __future__ import annotations

import json

import pytest

import bt4
from bt4 import api
from bt4.biomodels.splice import (
    USE_ATTESTED_ENV_VAR,
    AttestationError,
    ConsensusPwmSplicePredictor,
    PangolinSplicePredictor,
    SpliceAiSplicePredictor,
    attested_promotion_enabled,
    bundled_attestation,
    default,
    promote_if_attested,
)
from bt4.biomodels.splice.attestation import MAX_ATTESTATION_TOLERANCE
from bt4.biomodels.splice.attestations import bundled_attestation_path
from bt4.biomodels.splice.pangolin import DEFAULT_TISSUES, PINNED_WEIGHT_SHA256
from bt4.pipeline.splice_audit import _FlankedPredictor

# The content hash of the committed attestation, reproduced independently from the
# gate's scalars. Pinning it here means an edit to the file -- even one that still
# parses -- fails CI rather than silently changing what BT4 claims.
_EXPECTED_CONTENT_HASH = "5176032cabca1eecc68a9b313a9594378316b562190f0c1f2267fa7fd3164e76"


# --------------------------------------------------------------------------
# The committed file


def test_pangolin_attestation_ships_and_passed() -> None:
    """A passing Pangolin attestation is bundled."""
    att = bundled_attestation("pangolin")
    assert att is not None
    assert att.backend == "pangolin"
    assert att.passed is True
    assert att.n_cases == 18


def test_attestation_records_an_exact_reproduction() -> None:
    """The gate did not merely pass -- it matched bit-for-bit.

    A deviation of exactly zero is a stronger claim than "within tolerance", and
    worth pinning: a future adapter change that introduced float drift would still
    pass the gate but should not silently inherit this attestation.
    """
    att = bundled_attestation("pangolin")
    assert att is not None
    assert att.max_abs_deviation == 0.0
    assert att.tolerance <= MAX_ATTESTATION_TOLERANCE


def test_attestation_content_hash_is_pinned() -> None:
    """The committed attestation's content hash must not drift unnoticed."""
    att = bundled_attestation("pangolin")
    assert att is not None
    assert att.content_hash() == _EXPECTED_CONTENT_HASH


def test_attestation_covers_the_adapter_pinned_weights_exactly() -> None:
    """It claims the adapter's full pinned map -- a subset could never promote."""
    att = bundled_attestation("pangolin")
    assert att is not None
    assert dict(att.weight_sha256) == dict(PINNED_WEIGHT_SHA256)
    assert len(att.weight_sha256) == 12


def test_committed_file_carries_no_raw_scores() -> None:
    """The licensed per-position outputs must never reach the repository."""
    raw = json.loads(bundled_attestation_path("pangolin").read_text(encoding="utf-8"))
    assert set(raw) == {
        "backend",
        "passed",
        "max_abs_deviation",
        "n_cases",
        "tolerance",
        "weight_sha256",
        "bt4_version",
        "schema_version",
    }
    blob = json.dumps(raw)
    for banned in ("sequence", "expected", "acceptor", "donor", "site_scores"):
        assert banned not in blob


def test_no_spliceai_attestation_ships() -> None:
    """SpliceAI's gate has not been run, so nothing is claimed for it."""
    assert bundled_attestation("spliceai") is None


def test_attestation_records_the_producing_version() -> None:
    """Provenance: the attestation names the BT4 that produced it."""
    att = bundled_attestation("pangolin")
    assert att is not None
    assert att.bt4_version == bt4.__version__


# --------------------------------------------------------------------------
# Promotion is opt-in


def test_promotion_is_off_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    """Without the opt-in, a bundled attestation changes nothing."""
    monkeypatch.delenv(USE_ATTESTED_ENV_VAR, raising=False)
    assert attested_promotion_enabled() is False
    assert promote_if_attested(PangolinSplicePredictor()).calibrated is False


def test_explicit_opt_in_promotes(monkeypatch: pytest.MonkeyPatch) -> None:
    """With the opt-in, Pangolin reports calibrated."""
    monkeypatch.delenv(USE_ATTESTED_ENV_VAR, raising=False)
    promoted = promote_if_attested(PangolinSplicePredictor(), enabled=True)
    assert promoted.calibrated is True
    assert isinstance(promoted, PangolinSplicePredictor)
    assert promoted.tissues == DEFAULT_TISSUES


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
def test_env_var_truthy_values(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    """The documented truthy spellings all opt in."""
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, value)
    assert attested_promotion_enabled() is True


@pytest.mark.parametrize("value", ["", "0", "false", "no", "off", "maybe"])
def test_env_var_falsy_values(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    """Anything else leaves promotion off -- the safe direction."""
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, value)
    assert attested_promotion_enabled() is False


def test_spliceai_is_not_promoted_even_when_opted_in() -> None:
    """No attestation ships for SpliceAI, so the opt-in cannot promote it."""
    assert promote_if_attested(SpliceAiSplicePredictor(), enabled=True).calibrated is False


def test_baseline_passes_through_unharmed() -> None:
    """The seam is safe to call on any backend, including non-attestable ones."""
    baseline = ConsensusPwmSplicePredictor()
    assert promote_if_attested(baseline, enabled=True) is baseline


def test_promotion_refuses_a_configuration_the_gate_did_not_cover() -> None:
    """A partial tissue set loads only part of the attested weights -- refuse it."""
    with pytest.raises(AttestationError, match="tissue"):
        promote_if_attested(PangolinSplicePredictor(tissues=("brain",)), enabled=True)


# --------------------------------------------------------------------------
# What must NOT change


def test_default_still_returns_the_uncalibrated_baseline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Even opted in, ``default()`` is the honest PWM baseline.

    ``default()`` needs no per-user weight configuration, so it cannot assume the
    licensed model is installed or wanted. Promotion is requested explicitly.
    """
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, "1")
    model = default()
    assert isinstance(model, ConsensusPwmSplicePredictor)
    assert model.calibrated is False


def test_a_promoted_backend_still_reports_uncalibrated_with_flanks() -> None:
    """Regime scoping survives promotion.

    The gate was captured on the bare-CDS, ``N``-padded path. A real-flank score is
    a regime it never exercised, so the flag must clear -- promotion does not widen
    what the attestation covers.
    """
    promoted = promote_if_attested(PangolinSplicePredictor(), enabled=True)
    assert promoted.calibrated is True
    assert _FlankedPredictor(promoted, "ACGT", "").calibrated is False
    assert _FlankedPredictor(promoted, "", "").calibrated is True


def test_api_exposes_the_opt_in_surface() -> None:
    """cli/app reach promotion through ``api``, never by importing biomodels."""
    assert api.USE_ATTESTED_SPLICE_ENV_VAR == USE_ATTESTED_ENV_VAR
    assert api.bundled_splice_attestation("pangolin") is not None
    assert callable(api.attested_splice_promotion_enabled)


# --------------------------------------------------------------------------
# The promotion opt-in is reachable per call, not only per process


def test_available_backends_take_the_opt_in_as_an_argument() -> None:
    """A GUI must be able to offer the choice without mutating os.environ.

    Reading only ``$BT4_SPLICE_USE_ATTESTED`` meant the switch was set before launch or
    not at all, which is not a control a user can see.
    """
    import inspect

    from bt4.api import available_splice_backends

    parameter = inspect.signature(available_splice_backends).parameters["use_attested"]
    assert parameter.default is None  # None = read the standing env-var opt-in
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY


def test_the_default_stays_uncalibrated_however_it_is_called() -> None:
    """Promotion is opt-in in every spelling; the baseline is never promoted."""
    from bt4.api import available_splice_backends

    for kwargs in ({}, {"use_attested": False}, {"use_attested": True}):
        backends = available_splice_backends(**kwargs)  # type: ignore[arg-type]
        baseline = next(b for b in backends if b.name == "consensus-pwm-baseline")
        assert baseline.calibrated is False, kwargs


def test_the_attested_probe_reports_only_installed_and_attested_backends() -> None:
    """What a caller needs to know whether offering the opt-in is meaningful.

    A control that silently does nothing is worse than one that says why it cannot.
    In CI no licensed weights are installed, so this is empty and a GUI disables the
    switch with an explanation.
    """
    from bt4.api import attested_backends_available

    reported = attested_backends_available()
    assert isinstance(reported, tuple)
    assert set(reported) <= {"pangolin", "spliceai"}
    # Every reported backend must actually carry a committed attestation.
    from bt4.biomodels.splice import bundled_attestation

    assert all(bundled_attestation(name) is not None for name in reported)
