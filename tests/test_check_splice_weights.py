"""Tests for the standalone splice-weight hash-check script.

The script lives under ``scripts/`` (outside the ``bt4`` import graph), so it is
loaded by file path here, exactly as ``test_sensitivity.py`` does.

These tests pin the honesty behaviors that make step A3 of
``docs/DESIGN_splice_cnn_calibration.md`` worth running at all:

* an **absent** backend is reported as "not installed", never as a pass -- a
  missing install must not read like a verified one;
* a **mismatched** or **missing** weight file is reported per file and fails;
* only a **complete, exactly-matching** set passes, mirroring
  :func:`~bt4.biomodels.splice.verified_predictor`'s all-or-nothing comparison
  against ``PINNED_WEIGHT_SHA256``;
* the exit status distinguishes "nothing to check" from "checked and wrong".

No real weights are needed: the tests synthesize files whose SHA-256 is made to
match by pinning the *content*, so CI exercises every path without any licensed
bytes (and without torch or TensorFlow).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

_SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "check_splice_weights.py"


def _load_script() -> ModuleType:
    """Load ``scripts/check_splice_weights.py`` as a module by file path."""
    spec = importlib.util.spec_from_file_location("bt4_check_splice_weights", _SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def script() -> ModuleType:
    """The loaded hash-check script module."""
    return _load_script()


def _write_matching_weights(directory: Path, pins: dict[str, str]) -> dict[str, str]:
    """Write files whose real SHA-256 matches ``pins``, by *searching* for content.

    The pinned digests are of real licensed weights we do not have, so instead of
    faking the hash we fake the *pin*: return a pin map computed from content we
    wrote. Callers monkeypatch the adapter's ``PINNED_WEIGHT_SHA256`` with the
    returned map, which exercises the comparison logic honestly.
    """
    computed: dict[str, str] = {}
    for name in pins:
        payload = f"synthetic-weight-bytes-for-{name}".encode()
        (directory / name).write_bytes(payload)
        computed[name] = hashlib.sha256(payload).hexdigest()
    return computed


def test_absent_backend_is_not_a_pass(script: ModuleType, tmp_path: Path) -> None:
    """A directory that does not resolve reports 'not checked', never 'passed'."""
    check = script.check_backend("pangolin", str(tmp_path / "does-not-exist"))
    assert check.checked is False
    assert check.passed is False  # the load-bearing distinction
    assert check.reason is not None
    assert check.files == []


def test_mismatch_and_missing_are_reported_per_file(
    script: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Wrong bytes report 'mismatch'; absent files report 'missing'; neither passes."""
    from bt4.biomodels.splice import spliceai as mod

    pins = dict(mod.PINNED_WEIGHT_SHA256)
    names = sorted(pins)
    # Two files present with the wrong bytes, the rest absent entirely.
    (tmp_path / names[0]).write_bytes(b"not the real weights")
    (tmp_path / names[1]).write_bytes(b"x")

    check = script.check_backend("spliceai", str(tmp_path))

    assert check.checked is True
    assert check.passed is False
    by_name = {f.name: f.status for f in check.files}
    assert by_name[names[0]] == "mismatch"
    assert by_name[names[1]] == "mismatch"
    assert all(by_name[n] == "missing" for n in names[2:])
    # A mismatch records both digests so a maintainer can tell which release it is.
    mismatched = next(f for f in check.files if f.name == names[0])
    assert mismatched.actual is not None
    assert mismatched.actual != mismatched.expected
    # A missing file has no computed digest to report.
    absent = next(f for f in check.files if f.name == names[2])
    assert absent.actual is None


def test_complete_matching_set_passes(
    script: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only a complete, exactly-matching weight set is a pass."""
    from bt4.biomodels.splice import spliceai as mod

    real_pins = dict(mod.PINNED_WEIGHT_SHA256)
    synthetic = _write_matching_weights(tmp_path, real_pins)
    monkeypatch.setattr(mod, "PINNED_WEIGHT_SHA256", synthetic)

    check = script.check_backend("spliceai", str(tmp_path))

    assert check.checked is True
    assert check.passed is True
    assert all(f.status == "ok" for f in check.files)
    assert len(check.files) == len(synthetic)


def test_a_subset_never_passes(
    script: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Removing one file from an otherwise-matching set fails the whole backend.

    This mirrors ``verified_predictor``: it compares the *full* sorted tuple, so a
    subset can never satisfy an attestation. The check must agree.
    """
    from bt4.biomodels.splice import spliceai as mod

    synthetic = _write_matching_weights(tmp_path, dict(mod.PINNED_WEIGHT_SHA256))
    monkeypatch.setattr(mod, "PINNED_WEIGHT_SHA256", synthetic)
    (tmp_path / sorted(synthetic)[0]).unlink()

    check = script.check_backend("spliceai", str(tmp_path))

    assert check.passed is False
    assert sum(1 for f in check.files if f.status == "missing") == 1


def test_exit_status_separates_nothing_to_check_from_checked_and_wrong(
    script: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """No install exits 0 (not a failure); a resolved-but-wrong backend exits 1."""
    missing = tmp_path / "nope"
    assert script.main(["--backend", "spliceai", "--spliceai-dir", str(missing)]) == 0
    assert "NOT INSTALLED" in capsys.readouterr().out

    (tmp_path / "spliceai1.h5").write_bytes(b"wrong")
    assert script.main(["--backend", "spliceai", "--spliceai-dir", str(tmp_path)]) == 1
    out = capsys.readouterr().out
    assert "DOES NOT MATCH" in out
    # The refusal must tell the maintainer what NOT to do (CLAUDE.md 10.15).
    assert "Do NOT edit PINNED_WEIGHT_SHA256" in out


def test_json_output_is_valid_and_machine_readable(
    script: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """``--json`` emits parseable JSON carrying the per-file verdicts."""
    (tmp_path / "spliceai1.h5").write_bytes(b"wrong")
    script.main(["--backend", "spliceai", "--spliceai-dir", str(tmp_path), "--json"])

    payload = json.loads(capsys.readouterr().out)
    backend = payload["backends"][0]
    assert backend["backend"] == "spliceai"
    assert backend["checked"] is True
    assert backend["passed"] is False
    assert backend["n_ok"] == 0
    assert {f["status"] for f in backend["files"]} == {"mismatch", "missing"}


def test_weights_dir_is_public_on_both_adapters() -> None:
    """The script must reach the resolver through a public API (CLAUDE.md 10.9)."""
    from bt4.biomodels.splice import PangolinSplicePredictor, SpliceAiSplicePredictor

    for predictor in (PangolinSplicePredictor(), SpliceAiSplicePredictor()):
        resolver = predictor.weights_dir
        assert callable(resolver)
        assert not resolver.__name__.startswith("_")
        # Never raises, even with nothing installed.
        assert predictor.weights_dir() is None or isinstance(predictor.weights_dir(), Path)


def test_weights_dir_honours_an_explicit_directory(tmp_path: Path) -> None:
    """An explicit ``model_dir`` wins, and a non-directory resolves to ``None``."""
    from bt4.biomodels.splice import PangolinSplicePredictor

    assert PangolinSplicePredictor(model_dir=str(tmp_path)).weights_dir() == tmp_path
    assert PangolinSplicePredictor(model_dir=str(tmp_path / "nope")).weights_dir() is None
