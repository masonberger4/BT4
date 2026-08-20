"""Tests for the ``bt4`` command-line interface (the print-only shell)."""

from __future__ import annotations

import json
import os

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


# The subcommands `bt4 --help` lists. Kept explicit rather than scraped from the
# parser so that adding a command without a rendering help string fails here.
_SUBCOMMANDS = (
    "optimize",
    "library",
    "validate",
    "organisms",
    "enzymes",
    "presets",
    "tracks",
    "build-table",
)


@pytest.mark.parametrize("command", ["", *_SUBCOMMANDS])
def test_help_renders_for_every_command(
    command: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """Every ``--help`` must actually render.

    argparse %-formats help strings *while formatting help*, so a literal
    percent sign written as ``%M`` instead of ``%%M`` raises ValueError from
    inside the formatter. That is exactly what happened to ``bt4 tracks`` -- and,
    because a subparser's one-line help is reprinted in the top-level listing, to
    ``bt4 --help`` itself. Nothing else in the suite invokes ``--help``, so the
    first thing a new user types was broken and every other test passed.
    """
    argv = [command, "--help"] if command else ["--help"]
    with pytest.raises(SystemExit) as excinfo:
        main(argv)
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip()


def test_reference_set_flag_changes_the_delivered_sequence(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """``--reference-set`` must reach the engine, and the run must report it."""
    protein = "MKTAYIAKQRQISFVKSHFSRQ"
    assert main(["optimize", protein, "--organism", "escherichia_coli", "--json"]) == 0
    default = json.loads(capsys.readouterr().out)
    assert main(
        [
            "optimize", protein, "--organism", "escherichia_coli",
            "--reference-set", "genome_wide", "--json",
        ]
    ) == 0
    genome_wide = json.loads(capsys.readouterr().out)

    assert default["audit"]["codon_reference_set"] == "highly_expressed"
    assert genome_wide["audit"]["codon_reference_set"] == "genome_wide"
    assert default["dna"] != genome_wide["dna"]


def test_reference_set_refuses_rather_than_substituting(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Asking for a table an organism lacks must fail, not silently swap.

    A. thaliana ships only the genome-wide table (PaxDb identifies its proteins
    by UniProt accession the pinned annotation does not carry). Quietly handing
    that table back would make the run's CAI answer a different question than
    the one asked, while still exiting 0.
    """
    code = main(
        [
            "optimize", "MAAL", "--organism", "arabidopsis_thaliana",
            "--reference-set", "highly_expressed",
        ]
    )
    assert code != 0
    assert "no 'highly_expressed' codon table" in capsys.readouterr().err


def test_organisms_lists_reference_sets(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["organisms"]) == 0
    out = capsys.readouterr().out
    assert "highly_expressed" in out
    # A. thaliana has only the genome-wide table, and the listing must say so
    # rather than implying every organism offers both.
    arabidopsis = next(
        line for line in out.splitlines() if line.startswith("arabidopsis_thaliana")
    )
    assert "highly_expressed" not in arabidopsis


# --------------------------------------------------------------------------- #
# A user-facing flag must actually do something.
# --------------------------------------------------------------------------- #


def test_use_attested_splice_flag_is_wired(monkeypatch: pytest.MonkeyPatch) -> None:
    """``--use-attested-splice`` must have an effect, not merely parse.

    It shipped dead: `_enable_attested_splice` was defined and `--use-attested-splice`
    was advertised in `--help` on both `optimize` and `validate`, but nothing ever
    called the helper, so the flag silently did nothing and the cross-check kept
    printing UNCALIBRATED. A control that parses and no-ops is worse than an absent
    one -- it tells the user they opted in when they did not.

    Asserted through the observable effect (the env var the library actually
    consults), not by checking that a function was called, so the test survives the
    wiring being done a different way.
    """
    from bt4 import api
    from bt4.cli.__main__ import _parser, main

    monkeypatch.delenv(api.USE_ATTESTED_SPLICE_ENV_VAR, raising=False)

    for command in ("optimize", "validate"):
        monkeypatch.delenv(api.USE_ATTESTED_SPLICE_ENV_VAR, raising=False)
        argv = (
            ["optimize", "MAAK", "--use-attested-splice"]
            if command == "optimize"
            else ["validate", "ATGGCCGCCAAATAA", "--use-attested-splice"]
        )
        # Parse-and-dispatch through `main`, which is where the wiring has to live
        # for BOTH subcommands; a per-command call site would pass one and fail the
        # other, which is the shape of the original defect.
        assert main(argv) == 0, command
        assert os.environ.get(api.USE_ATTESTED_SPLICE_ENV_VAR) == "1", command

    # Absent the flag, the switch stays untouched -- opting in must remain explicit.
    monkeypatch.delenv(api.USE_ATTESTED_SPLICE_ENV_VAR, raising=False)
    assert main(["validate", "ATGGCCGCCAAATAA"]) == 0
    assert api.USE_ATTESTED_SPLICE_ENV_VAR not in os.environ

    # And the flag is genuinely reachable on both parsers (guards a silent rename).
    parser = _parser()
    for command in ("optimize", "validate"):
        args = parser.parse_args(
            [command, "MAAK" if command == "optimize" else "ATGGCCGCCAAATAA",
             "--use-attested-splice"]
        )
        assert args.use_attested_splice is True, command


def test_crosscheck_never_prints_the_bare_word_calibrated() -> None:
    """A fidelity attestation must not be reported as statistical calibration.

    Wiring `--use-attested-splice` made this reachable for the first time: with the
    flag dead, a promoted backend's tag could never be printed from the CLI. The tag
    read simply "calibrated", which is the stronger of the two claims BT4 keeps
    apart -- the flag is set by a *fidelity* attestation (the adapter reproduces the
    published model bit-for-bit) and asserts nothing about whether a score of 0.5
    means a 50% chance of splicing. Naming the weaker claim explicitly is the whole
    point; a reader who sees one word will take the stronger reading.
    """
    from pathlib import Path as _Path

    source = (
        _Path(__file__).resolve().parent.parent / "src" / "bt4" / "cli" / "__main__.py"
    ).read_text(encoding="utf-8")
    assert '"calibrated" if cc.calibrated' not in source
    assert "fidelity-attested (reproduces upstream; NOT statistically calibrated)" in source
