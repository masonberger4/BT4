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


def test_optimize_with_cpg(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--cpg-weight", "3", "--cpg-mode", "deplete"]
    assert main(argv) == 0
    assert "proven_optimal" in capsys.readouterr().out


def test_optimize_gc_budget(capsys: pytest.CaptureFixture[str]) -> None:
    pytest.importorskip("ortools")
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--gc-min", "36", "--gc-max", "40",
            "--max-homopolymer", "0"]
    assert main(argv) == 0
    assert "cpsat" in capsys.readouterr().out


def test_optimize_cpg_budget(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MARPGARSTKLE", "--cpg-max", "3"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    # The dinucleotide budget uses the exact bucketed DP -> proven optimal, and the
    # recomputed CpG count is printed and within the cap.
    assert "proven_optimal" in out
    assert "lagrangian" in out
    assert "CG count" in out


def test_optimize_upa_budget(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MYLIVFYLIV", "--upa-min", "1", "--upa-max", "4"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "TA count" in out


def test_optimize_both_dinuc_families_errors() -> None:
    # Only one dinucleotide budget at a time: mixing --cpg-* and --upa-* is an error.
    argv = ["optimize", "MAAL", "--cpg-max", "3", "--upa-max", "3"]
    assert main(argv) == 2


def test_optimize_with_minmax(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--minmax-weight", "2", "--minmax-direction", "min"]
    assert main(argv) == 0
    assert "proven_optimal" in capsys.readouterr().out


def test_optimize_with_repeats(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--tandem-unit", "3", "--inverted-stem", "4"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "0 hard" in out


def test_optimize_avoid_internal_start(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--avoid-internal-start"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "0 hard" in out


def test_optimize_refine(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--refine", "--refine-iterations", "200"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "heuristic" in out
    assert "folding" in out


def test_optimize_with_tai(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["optimize", "MAALKHETQWSNDECFGR", "--tai-weight", "2"]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "proven_optimal" in out
    assert "tAI" in out


def test_library_fasta(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["library", "MAALKHETQWSNDECFGR", "--n", "4", "--max-homopolymer", "5"]
    assert main(argv) == 0
    out = capsys.readouterr()
    # Four FASTA records on stdout; the honest SAMPLED note goes to stderr.
    assert out.out.count(">") == 4
    assert ">bt4_lib_1" in out.out
    assert "SAMPLED" in out.err
    assert "not optimized" in out.err.lower() or "not optimized" in out.err


def test_library_json(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["library", "MAALKHETQW", "--n", "3", "--seed", "7", "--json"]
    assert main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["n"] == 3
    assert payload["certificate"] == "sampled"
    assert len(payload["sequences"]) == 3
    assert "manifest" in payload
    for seq in payload["sequences"]:
        assert seq["certificate"]["status"] == "sampled"


def test_library_deterministic_cli(capsys: pytest.CaptureFixture[str]) -> None:
    argv = ["library", "MAALKHETQW", "--n", "3", "--seed", "11", "--json"]
    assert main(argv) == 0
    first = json.loads(capsys.readouterr().out)
    assert main(argv) == 0
    second = json.loads(capsys.readouterr().out)
    assert [s["dna"] for s in first["sequences"]] == [s["dna"] for s in second["sequences"]]


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
