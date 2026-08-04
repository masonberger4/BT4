"""Tests for the wrapped-Pangolin splice backend.

The Pangolin adapter drives an *installed* GPL ``pangolin`` package and its
weights, neither of which ships with BT4 (or with CI). These tests therefore
cover everything that does not need the real model:

* the contract surface and honesty flags (``calibrated is False`` until the
  fidelity gate; name carries the tissue set);
* configuration validation (tissues, ``top_k``);
* ``available()`` never raising and reporting ``False`` without the deps;
* graceful *refusal* -- ``score_sequence`` / ``site_scores`` raise a clear error
  rather than fabricating scores when torch / pangolin / weights are absent;
* the pure encoding (``_one_hot_channels``) and the SHA-256 weight-pin logic;
* the pinned-hash registry shape (20 P(splice) files, all valid digests);
* the integration-fidelity gate machinery, exercised with a fake predictor; and
* the guarantee that importing the module does not pull ``torch``.

The real inference path is gated behind ``importorskip`` so CI (no torch, no
pangolin) still passes, mirroring ``tests/test_folding.py``'s ViennaRNA gate.
"""

from __future__ import annotations

import dataclasses
import subprocess
import sys

import pytest

from bt4.biomodels.splice import (
    DEFAULT_TISSUES,
    FidelityCase,
    PangolinSplicePredictor,
    SplicePredictor,
    default,
    verify_pangolin_fidelity,
)
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.pangolin import (
    PANGOLIN_CONTEXT,
    PANGOLIN_FLANK,
    PINNED_WEIGHT_SHA256,
    TISSUE_OUTPUTS,
    _one_hot_channels,
    _sha256_file,
    _verify_weight_file,
)


def test_adapter_is_a_splice_predictor() -> None:
    model = PangolinSplicePredictor()
    assert isinstance(model, SplicePredictor)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)


def test_adapter_is_uncalibrated_by_default() -> None:
    model = PangolinSplicePredictor()
    # Wrapping a validated model is not the same as verifying the wrapping.
    assert model.calibrated is False
    assert model.fidelity_verified is False
    # default() must NOT prefer the uncalibrated adapter -- it stays the baseline.
    assert isinstance(default(), ConsensusPwmSplicePredictor)
    assert default().calibrated is False


def test_fidelity_verified_flag_flips_calibrated() -> None:
    model = PangolinSplicePredictor()
    promoted = dataclasses.replace(model, fidelity_verified=True)
    assert promoted.calibrated is True
    # The flag is the ONLY thing that changed.
    assert promoted.tissues == model.tissues


def test_name_carries_tissue_set() -> None:
    assert PangolinSplicePredictor().name == "pangolin[heart+liver+brain+testis]"
    assert PangolinSplicePredictor(tissues=("brain",)).name == "pangolin[brain]"


def test_config_validation() -> None:
    with pytest.raises(ValueError):
        PangolinSplicePredictor(tissues=())
    with pytest.raises(ValueError):
        PangolinSplicePredictor(tissues=("spleen",))
    with pytest.raises(ValueError):
        PangolinSplicePredictor(top_k=0)


def test_default_tissues_constant() -> None:
    assert DEFAULT_TISSUES == ("heart", "liver", "brain", "testis")
    assert PangolinSplicePredictor().tissues == DEFAULT_TISSUES


def test_context_constants() -> None:
    assert PANGOLIN_CONTEXT == 10_000
    assert PANGOLIN_FLANK == 5_000


def test_available_never_raises_and_is_false_without_deps() -> None:
    # Must return a bool regardless of whether torch/pangolin are installed.
    result = PangolinSplicePredictor().available()
    assert isinstance(result, bool)
    # CI has neither torch nor the pangolin package.
    if "torch" not in sys.modules:
        try:
            import torch  # type: ignore[import-not-found]  # noqa: F401
        except ImportError:
            assert result is False


def test_scoring_refuses_rather_than_fabricates_without_deps() -> None:
    model = PangolinSplicePredictor()
    if model.available():
        pytest.skip("Pangolin is installed; the refusal path is not exercised here")
    # Never returns fake scores when it cannot run -- it raises.
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.site_scores("ATGGCCGGCTAA")
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.score_sequence("ATGGCCGGCTAA")
    with pytest.raises((ModuleNotFoundError, RuntimeError, FileNotFoundError)):
        model.delta_splicing("ATGGCCGGCTAA", "ATGGCTGGGTAA")


def test_one_hot_channels_encoding() -> None:
    # Channel order is A, C, G, T; N is an all-zero column (matches Pangolin IN_MAP).
    channels = _one_hot_channels("ACGTN")
    assert len(channels) == 4
    assert all(len(c) == 5 for c in channels)
    assert channels[0] == [1.0, 0.0, 0.0, 0.0, 0.0]  # A row
    assert channels[1] == [0.0, 1.0, 0.0, 0.0, 0.0]  # C row
    assert channels[2] == [0.0, 0.0, 1.0, 0.0, 0.0]  # G row
    assert channels[3] == [0.0, 0.0, 0.0, 1.0, 0.0]  # T row
    # Each real base contributes exactly one 1.0 across the channels.
    for pos in range(4):
        assert sum(channels[c][pos] for c in range(4)) == 1.0
    # The N column is all zero.
    assert sum(channels[c][4] for c in range(4)) == 0.0


def test_tissue_output_channels() -> None:
    # Pangolin's INDEX_MAP for the four P(splice) heads.
    assert TISSUE_OUTPUTS["heart"].weight_index == 0
    assert TISSUE_OUTPUTS["heart"].output_channel == 1
    assert TISSUE_OUTPUTS["liver"].output_channel == 4
    assert TISSUE_OUTPUTS["brain"].output_channel == 7
    assert TISSUE_OUTPUTS["testis"].output_channel == 10


def test_pinned_registry_shape() -> None:
    # Exactly the 12 production P(splice) files: weight_index in {0,2,4,6} x the
    # .v2 folds 1..3 that the pangolin CLI loads (NOT the older .3 / 5-fold demo).
    expected = {
        f"final.{fold}.{index}.3.v2"
        for index in (0, 2, 4, 6)
        for fold in (1, 2, 3)
    }
    assert set(PINNED_WEIGHT_SHA256) == expected
    # Every pinned name is a production .v2 file over folds 1..3.
    assert all(name.endswith(".3.v2") for name in PINNED_WEIGHT_SHA256)
    # Every pin is a lowercase 64-hex SHA-256 digest.
    for name, digest in PINNED_WEIGHT_SHA256.items():
        assert len(digest) == 64, name
        assert all(ch in "0123456789abcdef" for ch in digest), name


def test_weight_verification(tmp_path: object) -> None:
    import hashlib
    from pathlib import Path

    assert isinstance(tmp_path, Path)
    # _sha256_file streams the exact bytes on disk.
    pinned = tmp_path / "final.1.0.3.v2"
    pinned.write_bytes(b"not the real weights")
    assert _sha256_file(pinned) == hashlib.sha256(b"not the real weights").hexdigest()

    # A pinned file name whose content does not match its pin -> refused.
    with pytest.raises(ValueError):
        _verify_weight_file(pinned)

    # A file whose name is not pinned -> KeyError (e.g. a usage head, or a
    # pre-.v2 demo weight the adapter deliberately does not use).
    unknown = tmp_path / "final.1.1.3.v2"  # a usage head, not a pinned P(splice) file
    unknown.write_bytes(b"x")
    with pytest.raises(KeyError):
        _verify_weight_file(unknown)
    legacy = tmp_path / "final.1.0.3"  # the older .3 demo weight is not pinned
    legacy.write_bytes(b"x")
    with pytest.raises(KeyError):
        _verify_weight_file(legacy)

    # A missing file -> FileNotFoundError.
    with pytest.raises(FileNotFoundError):
        _verify_weight_file(tmp_path / "final.2.0.3.v2")


class _FakeSites:
    """A stand-in exposing only ``site_scores`` for the fidelity-gate tests."""

    def __init__(self, mapping: dict[str, tuple[float, ...]]) -> None:
        self._mapping = mapping

    def site_scores(self, dna: str) -> tuple[float, ...]:
        return self._mapping[dna.upper()]


def test_fidelity_gate_pass_and_fail() -> None:
    seq = "ATGGCCGGCTAA"
    # Adapter reproduces the reference exactly -> gate passes.
    exact = _FakeSites({seq: (0.1,) * len(seq)})
    panel = [FidelityCase(sequence=seq, expected_site_scores=(0.1,) * len(seq))]
    report = verify_pangolin_fidelity(exact, panel)  # type: ignore[arg-type]
    assert report.passed is True
    assert report.n_cases == 1
    assert report.max_abs_deviation == pytest.approx(0.0)

    # Adapter deviates beyond tolerance -> gate fails, deviation reported.
    off = _FakeSites({seq: (0.5,) * len(seq)})
    report2 = verify_pangolin_fidelity(off, panel, tolerance=1e-3)  # type: ignore[arg-type]
    assert report2.passed is False
    assert report2.max_abs_deviation == pytest.approx(0.4)


def test_fidelity_gate_validation() -> None:
    with pytest.raises(ValueError):
        verify_pangolin_fidelity(PangolinSplicePredictor(), [])
    seq = "ATGGCC"
    bad = [FidelityCase(sequence=seq, expected_site_scores=(0.1, 0.2))]  # wrong length
    with pytest.raises(ValueError):
        verify_pangolin_fidelity(_FakeSites({seq: (0.1,) * len(seq)}), bad)  # type: ignore[arg-type]


def test_importing_module_does_not_load_torch() -> None:
    # Guard against a regression that adds a top-level torch/pangolin import,
    # which would make `import bt4` heavy (CLAUDE.md section 3).
    code = (
        "import bt4.biomodels.splice, bt4.biomodels.splice.pangolin, sys;"
        "bad=[m for m in ('torch','pangolin') if m in sys.modules];"
        "print(bad); assert not bad, bad"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
