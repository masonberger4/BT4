"""Tests for the ``bt4`` command-line interface (the print-only shell)."""

from __future__ import annotations

import json

import pytest

from bt4.cli.__main__ import main


def test_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert "bt4" in capsys.readouterr().out


def test_organisms(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["organisms"]) == 0
    assert "homo_sapiens" in capsys.readouterr().out


def test_optimize_summary(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["optimize", "MAAL", "--max-homopolymer", "5"]) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "CAI" in out


def test_optimize_fasta(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["optimize", "MAAL", "--fasta", "--header", "demo"]) == 0
    assert capsys.readouterr().out.startswith(">demo")


def test_optimize_json(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["optimize", "MAAL", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["schema_version"] == "1"
    assert payload["protein"] == "MAAL"
    assert "manifest" in payload["audit"]


def test_validate_flags_homopolymer(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["validate", "AAAAAA", "--max-homopolymer", "3"]) == 0
    out = capsys.readouterr().out
    assert "homopolymer" in out
    assert "False" in out


def test_invalid_protein_exits_2(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["optimize", "MAZX"]) == 2
    assert "error" in capsys.readouterr().err
