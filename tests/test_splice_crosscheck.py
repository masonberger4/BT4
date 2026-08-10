"""Tests for the opt-in, out-of-loop splice cross-check (pipeline + API + CLI).

Covers the graceful-degradation wrapper around a :class:`SplicePredictor`:

* an available backend reports sites + pooled risk, an unavailable one degrades to
  ``available is False`` (never raising) -- for ASSP outages *and* for a wrapped
  CNN's missing deps;
* ASSP results are stamped ``network_derived`` and kept out of the reproducible
  stdout artifact (the cross-check prints to stderr; ``optimize --json`` stdout
  stays a clean, ASSP-free manifest);
* backend name resolution and the ``bt4 validate --splice-backend`` /
  ``bt4 optimize --check-splice`` CLI wiring, driven by the committed offline
  fixtures.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bt4 import api
from bt4.biomodels.splice.assp import (
    AsspSplicePredictor,
    CachingAsspTransport,
    FixtureAsspTransport,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.pangolin import PangolinSplicePredictor
from bt4.biomodels.splice.spliceai import SpliceAiSplicePredictor
from bt4.cli.__main__ import main
from bt4.pipeline.splice_crosscheck import (
    SpliceCrossCheck,
    resolve_splice_backend,
    run_splice_crosscheck,
)

FIXTURES = Path(__file__).parent / "fixtures" / "assp"
SEQ_WITH_SITES = "ATGGCCGGCGATCGATCGATCGTAA"
SEQ_NO_FIXTURE = "ATGCACCACCACCACCACCACTGA"


def _fixture_assp() -> AsspSplicePredictor:
    return AsspSplicePredictor(transport=CachingAsspTransport(FixtureAsspTransport(str(FIXTURES))))


# --------------------------------------------------------------------------- #
# Pipeline / API
# --------------------------------------------------------------------------- #


def test_crosscheck_available_with_fixture() -> None:
    cc = run_splice_crosscheck(SEQ_WITH_SITES, predictor=_fixture_assp())
    assert isinstance(cc, SpliceCrossCheck)
    assert cc.available is True
    assert cc.reason is None
    assert cc.backend == "assp"
    assert cc.network_derived is True
    assert cc.calibrated is False
    assert cc.pooled_risk > 0.0
    assert len(cc.sites) == 3
    assert {s.kind for s in cc.sites} == {"donor", "acceptor"}
    assert any(s.site_class == "cryptic" for s in cc.sites)


def test_crosscheck_degrades_on_missing_fixture() -> None:
    # No committed fixture for this sequence -> the wrapper degrades, never raises.
    cc = run_splice_crosscheck(SEQ_NO_FIXTURE, predictor=_fixture_assp())
    assert cc.available is False
    assert cc.reason is not None
    assert cc.pooled_risk == 0.0
    assert cc.sites == ()
    # Still honestly labeled network-derived even when unavailable.
    assert cc.network_derived is True


def test_crosscheck_pwm_is_offline_and_local() -> None:
    cc = run_splice_crosscheck(SEQ_WITH_SITES, backend="pwm")
    assert cc.available is True
    assert cc.backend == "consensus-pwm-baseline"
    assert cc.network_derived is False


def test_crosscheck_cnn_degrades_without_weights() -> None:
    # A wrapped CNN with no installed deps/weights must DEGRADE in the cross-check
    # (available False), not raise -- same graceful path as an ASSP outage.
    if PangolinSplicePredictor().available():
        pytest.skip("pangolin is installed with weights; the missing path is not exercised")
    cc = run_splice_crosscheck(SEQ_WITH_SITES, backend="pangolin")
    assert cc.available is False
    assert cc.reason is not None
    assert cc.network_derived is False


class _WeightMismatchPredictor:
    """A fake CNN backend that raises ValueError like a failed SHA-256 weight pin."""

    name = "fake-cnn"
    calibrated = False

    def score_sequence(self, dna: str) -> object:
        raise ValueError("weight file failed its SHA-256 pin; refusing to load")

    def delta_splicing(self, designed: str, reference: str) -> float:
        return 0.0


def test_crosscheck_degrades_on_weight_integrity_valueerror() -> None:
    # A wrapped CNN's hash-pin refusal raises ValueError; the cross-check must
    # DEGRADE (never fail the run), not propagate -- the guarantee that a corrupted
    # weight file can never turn a successful optimize into a failure (§10.15).
    cc = run_splice_crosscheck(SEQ_WITH_SITES, predictor=_WeightMismatchPredictor())  # type: ignore[arg-type]
    assert cc.available is False
    assert cc.reason is not None
    assert "SHA-256" in cc.reason


def test_resolve_splice_backend_names() -> None:
    assert isinstance(resolve_splice_backend("assp"), AsspSplicePredictor)
    assert isinstance(resolve_splice_backend("pwm"), ConsensusPwmSplicePredictor)
    assert isinstance(resolve_splice_backend("baseline"), ConsensusPwmSplicePredictor)
    assert isinstance(resolve_splice_backend("consensus"), ConsensusPwmSplicePredictor)
    assert isinstance(resolve_splice_backend("pangolin"), PangolinSplicePredictor)
    assert isinstance(resolve_splice_backend("spliceai"), SpliceAiSplicePredictor)
    assert isinstance(resolve_splice_backend("ASSP"), AsspSplicePredictor)  # case-insensitive


def test_resolve_unknown_backend_raises() -> None:
    with pytest.raises(ValueError):
        resolve_splice_backend("nope")


def test_crosscheck_invalid_dna_raises() -> None:
    # A bad sequence is a CALLER error, not a service outage -> it raises.
    with pytest.raises(ValueError):
        run_splice_crosscheck("ATGZ", backend="pwm")


def test_crosscheck_bad_top_k_raises_not_degrades() -> None:
    # top_k is a CALLER error and must surface up front -- NOT be swallowed by the
    # degrade path (which also catches ValueError, for a CNN weight-hash refusal).
    with pytest.raises(ValueError):
        run_splice_crosscheck(SEQ_WITH_SITES, backend="pwm", top_k=0)


def test_api_splice_crosscheck_delegates() -> None:
    cc = api.splice_crosscheck(SEQ_WITH_SITES, predictor=_fixture_assp())
    assert cc.available is True
    assert cc.backend == "assp"


# --------------------------------------------------------------------------- #
# CLI wiring
# --------------------------------------------------------------------------- #


def test_cli_validate_assp_available(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(FIXTURES))
    assert main(["validate", SEQ_WITH_SITES, "--splice-backend", "assp"]) == 0
    captured = capsys.readouterr()
    # The audit stays on stdout; the cross-check is an out-of-band advisory on stderr.
    assert "feasible" in captured.out
    assert "splice cross-check [assp]" in captured.err
    assert "network-derived" in captured.err
    assert "cryptic" in captured.err
    # An AVAILABLE network-derived backend's numbers must NOT leak into stdout
    # (they are excluded from the reproducible artifact / manifest).
    assert "network-derived" not in captured.out
    assert "cryptic" not in captured.out
    assert "pooled risk" not in captured.out


def test_cli_validate_assp_degrades_rc0(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(FIXTURES))
    # A sequence with no fixture -> unavailable, but the run still succeeds (rc 0).
    assert main(["validate", SEQ_NO_FIXTURE, "--splice-backend", "assp"]) == 0
    assert "unavailable" in capsys.readouterr().err


def test_cli_optimize_check_splice_pwm(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["optimize", "MAAL", "--check-splice", "pwm"]) == 0
    captured = capsys.readouterr()
    assert "proven_optimal" in captured.out  # normal summary on stdout
    assert "splice cross-check [consensus-pwm-baseline]" in captured.err


def test_cli_optimize_json_crosscheck_is_out_of_band(capsys: pytest.CaptureFixture[str]) -> None:
    # The cross-check must NOT contaminate the reproducible stdout JSON / manifest.
    assert main(["optimize", "MAAL", "--json", "--check-splice", "pwm"]) == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)  # stdout is still clean, valid JSON
    assert "splice cross-check" not in captured.out
    assert "crosscheck" not in json.dumps(payload).lower()
    assert "splice cross-check" in captured.err  # the advisory went to stderr


def test_cli_optimize_assp_degrades_rc0(capsys: pytest.CaptureFixture[str]) -> None:
    # No httpx and no fixture env -> ASSP unavailable, but optimize still succeeds.
    if AsspSplicePredictor().available():
        pytest.skip("ASSP transport is available here; the degrade path is not exercised")
    assert main(["optimize", "MAAL", "--check-splice", "assp"]) == 0
    captured = capsys.readouterr()
    assert "ATG" in captured.out  # optimize output unaffected on stdout
    assert "unavailable" in captured.err
