"""Tests for the wrapped-SpliceAI splice backend.

The SpliceAI adapter drives an *installed* CC BY-NC ``spliceai`` package and its
Keras weights, neither of which ships with BT4 (or with CI). These tests cover
everything that does not need the real model:

* the contract surface and honesty flags (``calibrated is False`` until the
  fidelity gate; name is ``"spliceai"``);
* configuration validation (``top_k``);
* ``available()`` never raising and reporting ``False`` without the deps;
* graceful *refusal* -- scoring raises rather than fabricating scores when
  tensorflow / spliceai / weights are absent;
* the pure encoding (``_one_hot_rows``, position-major) and the acceptor/donor
  channel split (``_split_acceptor_donor``);
* the SHA-256 weight-pin logic and the pinned-registry shape (5 ``.h5`` files);
* the integration-fidelity gate machinery (both tracks), with a fake predictor;
* the guarantee that importing the module does not pull ``tensorflow``.

The real inference path is gated behind ``importorskip`` so CI (no TF, no
spliceai) still passes, mirroring ``tests/test_splice_pangolin.py``.
"""

from __future__ import annotations

import dataclasses
import subprocess
import sys

import pytest

from bt4.biomodels.splice import (
    SpliceAiFidelityCase,
    SpliceAiSplicePredictor,
    SplicePredictor,
    default,
    verify_spliceai_fidelity,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.spliceai import (
    ACCEPTOR_CHANNEL,
    DONOR_CHANNEL,
    PINNED_WEIGHT_SHA256,
    SPLICEAI_CONTEXT,
    SPLICEAI_ENSEMBLE_SIZE,
    SPLICEAI_FLANK,
    _one_hot_rows,
    _sha256_file,
    _split_acceptor_donor,
    _verify_weight_file,
)


def test_adapter_is_a_splice_predictor() -> None:
    model = SpliceAiSplicePredictor()
    assert isinstance(model, SplicePredictor)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)


def test_adapter_is_uncalibrated_by_default() -> None:
    model = SpliceAiSplicePredictor()
    assert model.calibrated is False
    assert model.fidelity_verified is False
    # default() must NOT prefer the uncalibrated adapter -- it stays the baseline.
    assert isinstance(default(), ConsensusPwmSplicePredictor)
    assert default().calibrated is False


def test_fidelity_verified_flag_flips_calibrated() -> None:
    model = SpliceAiSplicePredictor()
    promoted = dataclasses.replace(model, fidelity_verified=True)
    assert promoted.calibrated is True
    assert promoted.top_k == model.top_k


def test_name_is_spliceai() -> None:
    assert SpliceAiSplicePredictor().name == "spliceai"


def test_config_validation() -> None:
    with pytest.raises(ValueError):
        SpliceAiSplicePredictor(top_k=0)
    with pytest.raises(ValueError):
        SpliceAiSplicePredictor(top_k=-3)


def test_context_constants() -> None:
    assert SPLICEAI_CONTEXT == 10_000
    assert SPLICEAI_FLANK == 5_000
    assert SPLICEAI_ENSEMBLE_SIZE == 5


def test_channel_constants() -> None:
    # SpliceAI's 3-way softmax: 0 = null, 1 = acceptor, 2 = donor.
    assert ACCEPTOR_CHANNEL == 1
    assert DONOR_CHANNEL == 2


def test_available_never_raises_and_is_false_without_deps() -> None:
    result = SpliceAiSplicePredictor().available()
    assert isinstance(result, bool)
    # CI has neither tensorflow nor the spliceai package.
    try:
        import tensorflow  # type: ignore[import-not-found]  # noqa: F401
    except ImportError:
        try:
            import keras  # type: ignore[import-not-found]  # noqa: F401
        except ImportError:
            assert result is False


def test_scoring_refuses_rather_than_fabricates_without_deps() -> None:
    model = SpliceAiSplicePredictor()
    if model.available():
        pytest.skip("SpliceAI is installed; the refusal path is not exercised here")
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.site_scores("ATGGCCGGCTAA")
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.score_sequence("ATGGCCGGCTAA")
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.delta_splicing("ATGGCCGGCTAA", "ATGGCTGGGTAA")


def test_one_hot_rows_encoding() -> None:
    # Position-major [L][4], channel order A, C, G, T; N is an all-zero row.
    rows = _one_hot_rows("ACGTN")
    assert rows == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    # Each real base contributes exactly one 1.0; N contributes none.
    for row in rows[:4]:
        assert sum(row) == 1.0
    assert sum(rows[4]) == 0.0


def test_split_acceptor_donor_maps_channels() -> None:
    # Rows are [null, acceptor, donor]; the split pulls channels 1 and 2.
    per_position = [
        [0.98, 0.01, 0.01],
        [0.10, 0.85, 0.05],
        [0.20, 0.05, 0.75],
    ]
    acceptor, donor = _split_acceptor_donor(per_position)
    assert acceptor == (0.01, 0.85, 0.05)
    assert donor == (0.01, 0.05, 0.75)


def test_pinned_registry_shape() -> None:
    expected = {f"spliceai{index}.h5" for index in range(1, 6)}
    assert set(PINNED_WEIGHT_SHA256) == expected
    for name, digest in PINNED_WEIGHT_SHA256.items():
        assert len(digest) == 64, name
        assert all(ch in "0123456789abcdef" for ch in digest), name


def test_weight_verification(tmp_path: object) -> None:
    import hashlib
    from pathlib import Path

    assert isinstance(tmp_path, Path)
    pinned = tmp_path / "spliceai1.h5"
    pinned.write_bytes(b"not the real weights")
    assert _sha256_file(pinned) == hashlib.sha256(b"not the real weights").hexdigest()

    # A pinned file name whose content does not match its pin -> refused.
    with pytest.raises(ValueError):
        _verify_weight_file(pinned)

    # A file whose name is not pinned -> KeyError.
    unknown = tmp_path / "spliceai9.h5"
    unknown.write_bytes(b"x")
    with pytest.raises(KeyError):
        _verify_weight_file(unknown)

    # A missing file -> FileNotFoundError.
    with pytest.raises(FileNotFoundError):
        _verify_weight_file(tmp_path / "spliceai2.h5")


class _FakeSites:
    """A stand-in exposing only ``site_scores`` for the fidelity-gate tests."""

    def __init__(self, mapping: dict[str, tuple[tuple[float, ...], tuple[float, ...]]]) -> None:
        self._mapping = mapping

    def site_scores(self, dna: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
        return self._mapping[dna.upper()]


def test_fidelity_gate_pass_and_fail() -> None:
    seq = "ATGGCCGGCTAA"
    n = len(seq)
    acc = (0.1,) * n
    don = (0.2,) * n
    exact = _FakeSites({seq: (acc, don)})
    panel = [SpliceAiFidelityCase(sequence=seq, expected_acceptor=acc, expected_donor=don)]
    report = verify_spliceai_fidelity(exact, panel)  # type: ignore[arg-type]
    assert report.passed is True
    assert report.n_cases == 1
    assert report.max_abs_deviation == pytest.approx(0.0)

    # Deviation in the donor track beyond tolerance -> fail, deviation reported.
    off = _FakeSites({seq: (acc, (0.9,) * n)})
    report2 = verify_spliceai_fidelity(off, panel, tolerance=1e-3)  # type: ignore[arg-type]
    assert report2.passed is False
    assert report2.max_abs_deviation == pytest.approx(0.7)


def test_fidelity_gate_validation() -> None:
    with pytest.raises(ValueError):
        verify_spliceai_fidelity(SpliceAiSplicePredictor(), [])
    seq = "ATGGCC"
    bad = [
        SpliceAiFidelityCase(
            sequence=seq,
            expected_acceptor=(0.1, 0.2),  # wrong length
            expected_donor=(0.1,) * 6,
        )
    ]
    with pytest.raises(ValueError):
        verify_spliceai_fidelity(
            _FakeSites({seq: ((0.1,) * 6, (0.1,) * 6)}),  # type: ignore[arg-type]
            bad,
        )


def test_importing_module_does_not_load_tensorflow() -> None:
    code = (
        "import bt4.biomodels.splice, bt4.biomodels.splice.spliceai, sys;"
        "bad=[m for m in ('tensorflow','keras','tf_keras','spliceai') if m in sys.modules];"
        "print(bad); assert not bad, bad"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
