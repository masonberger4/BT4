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

The check mirrors PyInstaller's own mechanism: ``collect_data_files`` resolves its
``includes`` with :meth:`pathlib.Path.glob` against the package directory, then
resolves ``excludes`` the same way and *discards* what they match (a matching
directory takes everything under it). Both passes are reproduced here, because a
guard that read only ``includes`` would keep passing while an ``excludes`` entry
quietly dropped a file from the bundle -- the same defect one level down.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import bt4

_SPEC_PATH = Path(__file__).resolve().parents[1] / "packaging" / "bt4-studio.spec"
_PACKAGE_DIR = Path(bt4.__file__).resolve().parent


def _spec_clude_patterns() -> tuple[list[str], list[str]]:
    """Return the ``(includes, excludes)`` the spec passes to ``collect_data_files("bt4")``.

    Parsed from the spec's source with :mod:`ast` rather than executed: the spec imports
    PyInstaller and is evaluated in PyInstaller's own namespace (``Analysis``, ``EXE``,
    ``PYZ`` are injected), so it is not importable from a test run.

    ``excludes`` is read even though the spec does not currently pass one, because the
    day it does, a guard that ignored it would certify a coverage the bundle no longer
    has.
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
        cludes: dict[str, list[str]] = {"includes": [], "excludes": []}
        for kw in node.keywords:
            if kw.arg in cludes and isinstance(kw.value, ast.List):
                cludes[kw.arg] = [ast.literal_eval(elt) for elt in kw.value.elts]
        if not cludes["includes"]:
            continue
        return cludes["includes"], cludes["excludes"]
    raise AssertionError(f"no collect_data_files('bt4', includes=...) call in {_SPEC_PATH}")


def _spec_include_patterns() -> list[str]:
    """Return just the ``includes`` patterns (see :func:`_spec_clude_patterns`)."""
    return _spec_clude_patterns()[0]


def _collected_by_the_spec() -> set[Path]:
    """Return the files the spec's patterns actually collect, includes minus excludes.

    Reproduces ``collect_data_files``'s two passes: glob the includes and add, then glob
    the excludes and discard -- where a pattern matching a *directory* takes every file
    under it, which is how PyInstaller treats it.
    """
    includes, excludes = _spec_clude_patterns()

    def _walk(patterns: list[str]) -> set[Path]:
        found: set[Path] = set()
        for pattern in patterns:
            for path in _PACKAGE_DIR.glob(pattern):
                if path.is_dir():
                    found |= {p for p in path.rglob("*") if p.is_file()}
                elif path.is_file():
                    found.add(path)
        return found

    return _walk(includes) - _walk(excludes)


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
    collected = _collected_by_the_spec()
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
    assert target in _collected_by_the_spec()


def test_an_excludes_entry_cannot_hide_a_dropped_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A spec that excludes a shipped file fails the guard.

    The first version of this module read only ``includes``, so adding
    ``excludes=["**/*.attestation.json"]`` to the spec would have dropped both committed
    splice attestations from the bundle while every test here still passed -- the same
    defect this module exists to catch, one level down. Rebuilt against a real spec file
    rather than trusting a reading of PyInstaller's source.
    """
    spec = tmp_path / "excluding.spec"
    spec.write_text(
        'datas = collect_data_files(\n'
        '    "bt4", includes=["**/*.tsv", "**/*.json", "py.typed"],\n'
        '    excludes=["**/*.attestation.json"],\n'
        ')\n',
        encoding="utf-8",
    )
    monkeypatch.setitem(globals(), "_SPEC_PATH", spec)

    includes, excludes = _spec_clude_patterns()
    assert includes and excludes == ["**/*.attestation.json"]
    dropped = _shipped_data_files() - _collected_by_the_spec()
    assert any(p.name.endswith(".attestation.json") for p in dropped)

    with pytest.raises(AssertionError, match=r"attestation\.json"):
        test_bundle_collects_every_shipped_data_file()
