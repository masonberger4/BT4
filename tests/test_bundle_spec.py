"""The PyInstaller spec must collect every data file BT4 ships.

This is the packaging half of invariant-style coverage: `packaging/bt4-studio.spec`
names *glob patterns* for the non-Python files that go into the frozen BT4 Studio
bundle, and those patterns can silently fall behind the tree. They did: the spec
collected ``**/*.provenance.json``, and ``biomodels/expression/data/ribonn_sha256.json``
(read at import time by :mod:`bt4.biomodels.expression.ribonn`, which ``bt4.api``
imports) matched none of the patterns -- so the packaged app raised
``FileNotFoundError`` before its first window appeared, while every from-source test
passed. The committed splice fidelity attestations were missing from the bundle for
the same reason, and failed quieter: an absent attestation file reads as "none ships".

The check mirrors PyInstaller's own mechanism exactly -- ``collect_data_files``
resolves its ``includes`` with :meth:`pathlib.Path.glob` against the package
directory -- so a pattern that passes here collects the same files there.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import bt4

_SPEC_PATH = Path(__file__).resolve().parents[1] / "packaging" / "bt4-studio.spec"
_PACKAGE_DIR = Path(bt4.__file__).resolve().parent


def _spec_include_patterns() -> list[str]:
    """Return the ``includes`` patterns the spec passes to ``collect_data_files("bt4")``.

    Parsed from the spec's source with :mod:`ast` rather than executed: the spec
    imports PyInstaller and is evaluated in PyInstaller's own namespace (``Analysis``,
    ``EXE``, ``PYZ`` are injected), so it is not importable from a test run.
    """
    tree = ast.parse(_SPEC_PATH.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        if node.func.id != "collect_data_files":
            continue
        first = node.args[0] if node.args else None
        if not (isinstance(first, ast.Constant) and first.value == "bt4"):
            continue
        for kw in node.keywords:
            if kw.arg == "includes":
                return [ast.literal_eval(elt) for elt in kw.value.elts]  # type: ignore[attr-defined]
    raise AssertionError(f"no collect_data_files('bt4', includes=...) call in {_SPEC_PATH}")


def _shipped_data_files() -> set[Path]:
    """Return every non-Python file inside the installed ``bt4`` package."""
    return {
        path
        for path in _PACKAGE_DIR.rglob("*")
        if path.is_file()
        and path.suffix not in {".py", ".pyc", ".pyo"}
        and "__pycache__" not in path.parts
    }


def test_spec_declares_include_patterns() -> None:
    """The spec still has the call this module's coverage check reads."""
    patterns = _spec_include_patterns()
    assert patterns, "the spec collects no bt4 data files at all"


def test_bundle_collects_every_shipped_data_file() -> None:
    """No packaged data file is left out of the frozen bundle.

    A file that ships in the wheel but not in the bundle is a defect the from-source
    test suite structurally cannot see -- it only appears once the app is frozen and
    launched, i.e. for the user.
    """
    collected = {
        path
        for pattern in _spec_include_patterns()
        for path in _PACKAGE_DIR.glob(pattern)
        if path.is_file()
    }
    missing = sorted(str(p.relative_to(_PACKAGE_DIR)) for p in _shipped_data_files() - collected)
    assert not missing, (
        "packaging/bt4-studio.spec does not collect these data files, so they are "
        f"absent from the frozen BT4 Studio bundle: {missing}"
    )


def test_check_would_catch_a_dropped_pattern() -> None:
    """The check fails when a pattern stops covering the tree (it is not vacuous).

    Pins the guard itself: with only the TSV pattern, the JSON sidecars must be
    reported missing. A coverage test that passes no matter what is worse than none.
    """
    collected = {path for path in _PACKAGE_DIR.glob("**/*.tsv") if path.is_file()}
    missing = _shipped_data_files() - collected
    assert any(p.suffix == ".json" for p in missing)


@pytest.mark.parametrize(
    "relative",
    [
        "biomodels/expression/data/ribonn_sha256.json",
        "biomodels/splice/data/pangolin.attestation.json",
        "biomodels/splice/data/spliceai.attestation.json",
    ],
)
def test_the_files_that_broke_the_bundle_are_collected(relative: str) -> None:
    """The three files the old pattern missed are named, so a regression is legible."""
    target = _PACKAGE_DIR / relative
    assert target.is_file(), f"{relative} is no longer shipped; update this test"
    collected = {
        path
        for pattern in _spec_include_patterns()
        for path in _PACKAGE_DIR.glob(pattern)
        if path.is_file()
    }
    assert target in collected
