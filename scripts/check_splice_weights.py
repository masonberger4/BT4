"""Verify installed SpliceAI / Pangolin weights against BT4's pinned SHA-256s.

Step **A3** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md):
before capturing a reference panel or running an integration-fidelity gate, prove
the weights on this machine are byte-identical to the ones BT4 pins. This matters
because :func:`bt4.biomodels.splice.verified_predictor` -- the single seam that can
flip a wrapped CNN to ``calibrated=True`` -- compares the attestation's weight-hash
map against the adapter's :data:`PINNED_WEIGHT_SHA256` with ``!=`` on the **full
sorted tuple**. A missing or altered file, or a subset, can never satisfy it, so a
mismatch found here saves a wasted panel capture later.

The adapters already hash-verify each file immediately before loading it (so
unverified bytes are never unpickled). This script front-loads the same check for
*every* pinned file at once, and reports rather than raises, so a maintainer sees
the whole picture in one run.

**A mismatch is not something to fix by editing the pins.** Either the download came
from a mirror or a fork (re-download from the official source), or upstream
re-released the weights -- in which case updating
:data:`PINNED_WEIGHT_SHA256` is a deliberate, separately-described commit naming
the upstream release, never a silent edit to make an error go away (CLAUDE.md
sections 6, 10.15).

Weights are **not** bundled with BT4 (Pangolin is GPL-3.0; SpliceAI code is
PolyForm Strict 1.0.0 and its weights are CC BY-NC 4.0, noncommercial). Point this
script at your own install with ``--pangolin-dir`` / ``--spliceai-dir``, or set
``$BT4_PANGOLIN_MODEL_DIR`` / ``$BT4_SPLICEAI_MODEL_DIR``; with neither, each
backend's own resolver is used (which falls back to the installed package's
``models`` directory).

Run it directly::

    python scripts/check_splice_weights.py                    # both, auto-resolved
    python scripts/check_splice_weights.py --backend pangolin # just one
    python scripts/check_splice_weights.py --json             # machine-readable

Exit status is ``0`` when every *checked* backend matched, ``1`` when any backend
resolved but failed, and ``0`` when a backend simply is not installed (not being
installed is not a failure).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Sequence
from pathlib import Path

__all__ = ["BackendCheck", "FileCheck", "check_backend", "main"]

_READ_BLOCK = 1 << 20
"""Streaming read size (1 MiB) -- the weight files are millions of bytes each."""

_BACKENDS = ("pangolin", "spliceai")


class FileCheck:
    """One weight file's verdict.

    Attributes:
        name: The weight file's base name (the key in ``PINNED_WEIGHT_SHA256``).
        status: ``"ok"``, ``"mismatch"``, or ``"missing"``.
        expected: The pinned SHA-256.
        actual: The computed SHA-256, or ``None`` when the file is missing.
    """

    __slots__ = ("actual", "expected", "name", "status")

    def __init__(self, name: str, status: str, expected: str, actual: str | None) -> None:
        self.name = name
        self.status = status
        self.expected = expected
        self.actual = actual

    def to_dict(self) -> dict[str, object]:
        """Return the JSON-serializable form."""
        return {
            "name": self.name,
            "status": self.status,
            "expected": self.expected,
            "actual": self.actual,
        }


class BackendCheck:
    """One backend's verdict over its whole pinned weight set.

    Attributes:
        backend: ``"pangolin"`` or ``"spliceai"``.
        resolved_dir: The weights directory that was checked, or ``None`` when the
            backend is not installed / no directory resolved.
        reason: Why nothing was checked (``None`` when ``resolved_dir`` is set).
        files: Per-file verdicts, sorted by name.
    """

    __slots__ = ("backend", "files", "reason", "resolved_dir")

    def __init__(
        self,
        backend: str,
        resolved_dir: Path | None,
        reason: str | None,
        files: list[FileCheck],
    ) -> None:
        self.backend = backend
        self.resolved_dir = resolved_dir
        self.reason = reason
        self.files = files

    @property
    def checked(self) -> bool:
        """Whether this backend's weights were actually found and hashed."""
        return self.resolved_dir is not None

    @property
    def passed(self) -> bool:
        """Whether every pinned file was present and matched.

        An unchecked backend is **not** a pass -- callers must consult
        :attr:`checked` first. This keeps "not installed" and "verified" distinct
        rather than letting an absent backend read as a clean result.
        """
        return self.checked and bool(self.files) and all(f.status == "ok" for f in self.files)

    def to_dict(self) -> dict[str, object]:
        """Return the JSON-serializable form."""
        return {
            "backend": self.backend,
            "checked": self.checked,
            "passed": self.passed,
            "resolved_dir": str(self.resolved_dir) if self.resolved_dir else None,
            "reason": self.reason,
            "n_pinned": len(self.files),
            "n_ok": sum(1 for f in self.files if f.status == "ok"),
            "files": [f.to_dict() for f in self.files],
        }


def _sha256_file(path: Path) -> str:
    """Return the lowercase hex SHA-256 of ``path`` (streamed, constant memory)."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_READ_BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()


def _pins_and_dir(backend: str, explicit_dir: str | None) -> tuple[dict[str, str], Path | None]:
    """Return ``(pinned_sha256, resolved_weights_dir)`` for ``backend``.

    Uses each adapter's **public** :meth:`weights_dir` accessor rather than its
    private resolver, so this script stays free of cross-layer private imports
    (CLAUDE.md section 10.9). Imports are lazy so ``--help`` works regardless of
    which optional extras are installed.
    """
    if backend == "pangolin":
        from bt4.biomodels.splice import PangolinSplicePredictor
        from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256

        predictor = PangolinSplicePredictor(model_dir=explicit_dir)
        return dict(PINNED_WEIGHT_SHA256), predictor.weights_dir()
    if backend == "spliceai":
        from bt4.biomodels.splice import SpliceAiSplicePredictor
        from bt4.biomodels.splice.spliceai import PINNED_WEIGHT_SHA256

        predictor = SpliceAiSplicePredictor(model_dir=explicit_dir)
        return dict(PINNED_WEIGHT_SHA256), predictor.weights_dir()
    raise ValueError(f"unknown backend {backend!r}")  # pragma: no cover - argparse guards


def check_backend(backend: str, explicit_dir: str | None = None) -> BackendCheck:
    """Hash every pinned weight file for ``backend`` and report the verdict.

    Args:
        backend: ``"pangolin"`` or ``"spliceai"``.
        explicit_dir: A weights directory to check. When ``None`` the backend's own
            resolver is used (explicit arg -> env var -> installed package's
            ``models`` directory).

    Returns:
        A :class:`BackendCheck`. Never raises for a missing install or a bad hash --
        both are reported, so one absent backend cannot mask the other's result.
    """
    pins, model_dir = _pins_and_dir(backend, explicit_dir)
    if model_dir is None:
        where = explicit_dir or f"$BT4_{backend.upper()}_MODEL_DIR or the installed package"
        return BackendCheck(backend, None, f"no weights directory resolved ({where})", [])

    files: list[FileCheck] = []
    for name, expected in sorted(pins.items()):
        path = Path(model_dir) / name
        if not path.is_file():
            files.append(FileCheck(name, "missing", expected, None))
            continue
        actual = _sha256_file(path)
        status = "ok" if actual == expected else "mismatch"
        files.append(FileCheck(name, status, expected, actual))
    return BackendCheck(backend, Path(model_dir), None, files)


def _print_backend(check: BackendCheck) -> None:
    """Render one backend's verdict as a human-readable block."""
    print(f"\n=== {check.backend} ===")
    if not check.checked:
        print(f"  not checked: {check.reason}")
        print(f"  -> {check.backend}: NOT INSTALLED (this is not a failure)")
        return

    print(f"  weights dir: {check.resolved_dir}")
    for f in check.files:
        if f.status == "ok":
            print(f"  ok        {f.name}")
        elif f.status == "missing":
            print(f"  MISSING   {f.name}")
        else:
            print(f"  MISMATCH  {f.name}")
            print(f"      pinned {f.expected}")
            print(f"      actual {f.actual}")

    n_ok = sum(1 for f in check.files if f.status == "ok")
    total = len(check.files)
    if check.passed:
        print(f"  -> {check.backend}: ALL {total} PINS MATCH")
    else:
        print(f"  -> {check.backend}: {n_ok}/{total} matched - DOES NOT MATCH, STOP")
        print("     Do NOT edit PINNED_WEIGHT_SHA256 to silence this. See the")
        print("     module docstring and docs/DESIGN_splice_cnn_calibration.md step A3.")


def main(argv: Sequence[str] | None = None) -> int:
    """Verify pinned splice-CNN weights and report per backend.

    Returns:
        ``0`` if every backend that resolved matched (or none were installed),
        ``1`` if any resolved backend had a missing or mismatched file.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--backend",
        choices=(*_BACKENDS, "both"),
        default="both",
        help="Which backend to check (default: both).",
    )
    parser.add_argument("--pangolin-dir", default=None, help="Explicit Pangolin weights directory.")
    parser.add_argument("--spliceai-dir", default=None, help="Explicit SpliceAI weights directory.")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a table.")
    args = parser.parse_args(argv)

    wanted = _BACKENDS if args.backend == "both" else (args.backend,)
    explicit = {"pangolin": args.pangolin_dir, "spliceai": args.spliceai_dir}
    checks = [check_backend(name, explicit[name]) for name in wanted]

    if args.json:
        print(json.dumps({"backends": [c.to_dict() for c in checks]}, indent=2, sort_keys=True))
    else:
        for check in checks:
            _print_backend(check)
        installed = [c for c in checks if c.checked]
        print()
        if not installed:
            print("No splice CNN weights found. Install Pangolin and/or SpliceAI first")
            print("(docs/DESIGN_splice_cnn_calibration.md step A1), then re-run.")
        elif all(c.passed for c in installed):
            names = ", ".join(c.backend for c in installed)
            print(f"RESULT: all pins match for {names}.")
            print("Safe to proceed to panel capture (step A4).")
        else:
            failed = ", ".join(c.backend for c in installed if not c.passed)
            print(f"RESULT: {failed} did NOT match its pins. Stop and diagnose before step A4.")

    return 0 if all(c.passed for c in checks if c.checked) else 1


if __name__ == "__main__":
    raise SystemExit(main())
