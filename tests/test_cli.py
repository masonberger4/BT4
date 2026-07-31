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


def test_enzymes(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["enzymes"]) == 0
    assert "EcoRI" in capsys.readouterr().out


def test_optimize_with_enzyme(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECF", "--enzyme", "EcoRI", "--max-homopolymer", "5"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "0 hard" in out


def test_build_table(tmp_path: pytest.TempPathFactory, capsys: pytest.CaptureFixture[str]) -> None:
    from pathlib import Path

    from bt4.biomodels.codon.tables import load_table_from_file

    out = Path(str(tmp_path))
    cds = out / "cds.fasta"
    cds.write_text(
        ">g1\nATGGCCGCCCTGAAGCACGAGACCCAGTGGTAA\n>g2\nATGGCTGCACTGAAACATGAAACGCAATGGTAG\n",
        encoding="utf-8",
    )
    code = main(["build-table", str(cds), "--organism", "demo", "--out", str(out)])
    assert code == 0
    assert "wrote" in capsys.readouterr().out
    # The smoothed table must load straight back and cover every amino acid.
    table = load_table_from_file(out / "demo.tsv")
    assert table.organism == "demo"
    assert 0.0 < table.cai("ATGGCC") <= 1.0
