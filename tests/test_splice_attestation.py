"""Tests for the license-clean splice fidelity-attestation layer.

These pin the honesty guarantees of :mod:`bt4.biomodels.splice.attestation`:

* an attestation carries only license-clean scalars + public weight SHAs -- never
  a raw per-position model score (the licensed output);
* a failing or too-loose gate is never recorded as an attestation;
* a backend is promoted to ``calibrated=True`` only against an attestation that
  structurally matches the adapter's pinned weights, and refused otherwise;
* ``default()`` still returns the honest baseline (no attestation ships).

None of this needs torch / tensorflow / the licensed weights: the predictors are
frozen dataclasses that construct with defaults, and the reports are built by hand.
"""

from __future__ import annotations

import json

import pytest

from bt4.biomodels.splice import (
    AttestationError,
    ConsensusPwmSplicePredictor,
    FidelityAttestation,
    FidelityReport,
    PangolinSplicePredictor,
    SpliceAiSplicePredictor,
    attest_backend,
    default,
    load_attestation,
    verified_predictor,
)
from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256 as PANGOLIN_PINS
from bt4.biomodels.splice.spliceai import PINNED_WEIGHT_SHA256 as SPLICEAI_PINS

_VERSION = "0.4.0"


def _passing_pangolin_attestation() -> FidelityAttestation:
    report = FidelityReport(passed=True, max_abs_deviation=1e-6, n_cases=5, tolerance=1e-3)
    return attest_backend("pangolin", report, PANGOLIN_PINS, bt4_version=_VERSION)


def test_attestation_fields_are_license_clean() -> None:
    # The dataclass shape must be exactly the allowed scalar set -- no field that
    # could ever hold a raw per-position model score (the licensed output).
    from dataclasses import fields

    from bt4.biomodels.splice.attestation import _ALLOWED_FIELDS

    assert {f.name for f in fields(FidelityAttestation)} == _ALLOWED_FIELDS
    forbidden = ("score", "expected", "acceptor", "donor", "per_position", "sequence")
    for name in _ALLOWED_FIELDS:
        assert not any(tok in name for tok in forbidden), name


def test_serialized_form_has_only_license_clean_keys() -> None:
    att = _passing_pangolin_attestation()
    d = att.to_dict()
    assert set(d) <= {
        "backend", "passed", "max_abs_deviation", "n_cases", "tolerance",
        "weight_sha256", "bt4_version", "schema_version",
    }
    # weight_sha256 is a name->hexdigest map of PUBLIC weight-file hashes, not scores.
    assert all(isinstance(k, str) and isinstance(v, str) for k, v in d["weight_sha256"].items())


def test_attest_and_roundtrip() -> None:
    att = _passing_pangolin_attestation()
    assert att.backend == "pangolin"
    assert att.passed is True
    again = FidelityAttestation.from_dict(att.to_dict())
    assert again == att
    assert again.content_hash() == att.content_hash()


def test_content_hash_is_deterministic_and_weight_sensitive() -> None:
    a = _passing_pangolin_attestation()
    b = _passing_pangolin_attestation()
    assert a.content_hash() == b.content_hash()  # timestamp-free, deterministic
    # A different weight-hash set must change the stamp (provenance honesty).
    tampered = FidelityAttestation.from_dict(
        {**a.to_dict(), "weight_sha256": {"pangolin.v2": "deadbeef"}}
    )
    assert tampered.content_hash() != a.content_hash()


def test_attest_refuses_failing_gate() -> None:
    report = FidelityReport(passed=False, max_abs_deviation=0.5, n_cases=3, tolerance=1e-3)
    with pytest.raises(AttestationError, match="failing"):
        attest_backend("pangolin", report, PANGOLIN_PINS, bt4_version=_VERSION)


def test_attest_refuses_too_loose_tolerance() -> None:
    report = FidelityReport(passed=True, max_abs_deviation=1e-6, n_cases=3, tolerance=1e-1)
    with pytest.raises(AttestationError, match=r"floor|looser"):
        attest_backend("pangolin", report, PANGOLIN_PINS, bt4_version=_VERSION)


def test_attest_refuses_unknown_backend() -> None:
    report = FidelityReport(passed=True, max_abs_deviation=0.0, n_cases=1, tolerance=1e-3)
    with pytest.raises(AttestationError, match="unknown backend"):
        attest_backend("nope", report, PANGOLIN_PINS, bt4_version=_VERSION)


def test_from_dict_rejects_smuggled_raw_score_field() -> None:
    att = _passing_pangolin_attestation()
    poisoned = {**att.to_dict(), "expected_site_scores": [0.1, 0.9, 0.2]}
    with pytest.raises(AttestationError, match="unexpected"):
        FidelityAttestation.from_dict(poisoned)


def test_verified_predictor_flips_calibrated_on_match() -> None:
    pred = PangolinSplicePredictor()
    assert pred.calibrated is False
    promoted = verified_predictor(pred, _passing_pangolin_attestation())
    assert promoted.calibrated is True
    assert promoted.fidelity_verified is True
    # Original is untouched (frozen dataclass, replace returns a new instance).
    assert pred.calibrated is False


def test_verified_predictor_refuses_weight_mismatch() -> None:
    pred = PangolinSplicePredictor()
    bad = FidelityAttestation.from_dict(
        {**_passing_pangolin_attestation().to_dict(), "weight_sha256": {"pangolin.v2": "00"}}
    )
    with pytest.raises(AttestationError, match=r"weight SHA-256|different weights"):
        verified_predictor(pred, bad)


def test_verified_predictor_refuses_backend_mismatch() -> None:
    # A pangolin attestation must not calibrate a SpliceAI predictor.
    with pytest.raises(AttestationError, match=r"is for 'pangolin'|spliceai"):
        verified_predictor(SpliceAiSplicePredictor(), _passing_pangolin_attestation())


def test_verified_predictor_refuses_non_backend_predictor() -> None:
    with pytest.raises(AttestationError, match="not an attestable"):
        verified_predictor(ConsensusPwmSplicePredictor(), _passing_pangolin_attestation())


def test_spliceai_attestation_matches_its_pins() -> None:
    report = FidelityReport(passed=True, max_abs_deviation=1e-6, n_cases=4, tolerance=1e-3)
    att = attest_backend("spliceai", report, SPLICEAI_PINS, bt4_version=_VERSION)
    promoted = verified_predictor(SpliceAiSplicePredictor(), att)
    assert promoted.calibrated is True


def test_load_attestation_roundtrip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    att = _passing_pangolin_attestation()
    path = tmp_path / "pangolin.attestation.json"
    path.write_text(json.dumps(att.to_dict()), encoding="utf-8")
    assert load_attestation(path) == att


def test_load_attestation_missing_file_raises(tmp_path) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(AttestationError, match="cannot load"):
        load_attestation(tmp_path / "nope.json")


def test_default_returns_uncalibrated_baseline() -> None:
    # No attestation ships, so the honest PWM baseline remains the default.
    pred = default()
    assert pred.calibrated is False
    assert type(pred).__name__ == "ConsensusPwmSplicePredictor"
