"""Tests for the gate runner (``scripts/run_expression_gate.py``).

What matters here is not the arithmetic -- the gate itself is tested in
``test_expression_gate.py`` -- but the *judgement* the runner layers on top: a head is
only promotable if it passes the thresholds, **beats every baseline**, and produces an
interval that is actually informative. The baselines are the point: a within-protein
Spearman of 0.3 means nothing if plain CAI scores 0.35, because BT4 already optimizes
CAI directly and for free.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import pathlib
import sys
from types import ModuleType

import pytest

from bt4.biomodels.expression import PanelRow, panel_from_rows
from bt4.pipeline.expression_gate import (
    BASELINES,
    GateSettings,
    baseline_scores,
    run_panel_gate,
    score_panel,
)

_SCRIPTS = pathlib.Path(__file__).resolve().parent.parent / "scripts"


def _load_script(name: str) -> ModuleType:
    path = _SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"bt4_script_{name}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rg = _load_script("run_expression_gate")

_SETTINGS = GateSettings(
    within_group=True,
    recalibrate=True,
    coverage_tolerance=0.20,  # loose: these panels are small, coverage is granular
    bootstrap_resamples=200,
)


def _panel(n_groups: int = 8, n_variants: int = 6) -> object:
    """Panel whose within-protein signal is an index, uncorrelated with codon content."""
    rows = []
    for g in range(n_groups):
        for v in range(n_variants):
            # Same codon content in every variant of a group, permuted -- so CAI and GC3
            # carry no within-protein information at all, and only a head that knows the
            # actual driver can score.
            body = ("AAA" * (v + 1)) + ("AAG" * (n_variants - v))
            rows.append(
                PanelRow(
                    group=f"P{g:02d}",
                    variant_id=f"g{g}v{v}",
                    cds="ATG" + body + "TAA",
                    measured=100.0 * g + v,  # protein baseline + variant rank
                    utr5="GCCACC",
                    utr3="GCTAAT",
                )
            )
    return panel_from_rows(rows)


# --- baselines ----------------------------------------------------------------


def test_every_baseline_returns_one_score_per_row() -> None:
    from bt4.biomodels.codon.tables import load_table

    panel = _panel()
    table = load_table("homo_sapiens")
    head = [float(i) for i in range(len(panel.rows))]
    for name in BASELINES:
        scores = baseline_scores(name, panel.rows, head, table, seed=0)
        assert len(scores) == len(panel.rows)


def test_permutation_baseline_is_a_deterministic_shuffle_of_the_head() -> None:
    from bt4.biomodels.codon.tables import load_table

    panel = _panel()
    table = load_table("homo_sapiens")
    head = [float(i) for i in range(len(panel.rows))]

    first = baseline_scores("permutation", panel.rows, head, table, seed=3)
    again = baseline_scores("permutation", panel.rows, head, table, seed=3)

    assert first == again  # invariant #7
    assert sorted(first) == sorted(head)  # a permutation, not new numbers
    assert first != head  # and actually shuffled


def test_constant_baseline_is_constant_and_unknown_names_raise() -> None:
    from bt4.biomodels.codon.tables import load_table

    panel = _panel()
    table = load_table("homo_sapiens")
    assert set(baseline_scores("constant", panel.rows, [], table, 0)) == {0.0}
    with pytest.raises(ValueError, match="unknown baseline"):
        baseline_scores("nope", panel.rows, [], table, 0)


# --- the verdict --------------------------------------------------------------


def _report(panel: object, head: list[float]) -> dict[str, object]:
    """Run the pipeline gate with caller-supplied scores, then shape it like the CLI."""
    comparison = run_panel_gate(
        panel,  # type: ignore[arg-type]
        "null",
        settings=_SETTINGS,
        head_scores=head,
    )
    return rg.build_report(comparison)


def test_an_oracle_head_beats_every_baseline() -> None:
    panel = _panel()
    oracle = [row.measured for row in panel.rows]  # knows the answer exactly
    verdict = _report(panel, oracle)["verdict"]
    assert isinstance(verdict, dict)
    assert verdict["beats_every_baseline"] is True
    assert verdict["promotable"] is True


def test_a_head_that_merely_reproduces_cai_does_not_beat_the_cai_baseline() -> None:
    # The judgement that stops a head being credited for what BT4 already computes for
    # free, inside the optimizer loop.
    from bt4.biomodels.codon.tables import load_table

    panel = _panel()
    table = load_table("homo_sapiens")
    as_cai = [table.cai(row.cds) for row in panel.rows]

    verdict = _report(panel, as_cai)["verdict"]
    assert isinstance(verdict, dict)
    assert verdict["beats_every_baseline"] is False
    assert verdict["promotable"] is False


def test_a_blind_head_is_not_promotable() -> None:
    panel = _panel()
    verdict = _report(panel, [1.0] * len(panel.rows))["verdict"]
    assert isinstance(verdict, dict)
    assert verdict["promotable"] is False


def test_the_constant_baseline_shows_its_coverage_pass_in_the_table() -> None:
    # Split conformal is valid for ANY score function, so a constant predictor gets
    # correct coverage. It is in the report precisely so that fact is visible next to
    # the head's, rather than being a trap the reader has to remember.
    panel = _panel()
    report = _report(panel, [row.measured for row in panel.rows])
    baselines = report["baselines"]
    assert isinstance(baselines, dict)
    constant = baselines["constant"]
    assert constant["spearman"] == pytest.approx(0.0)  # no information
    assert constant["width_over_iqr"] > 1.0  # and a uselessly wide interval


def test_promotable_requires_an_informative_interval() -> None:
    panel = _panel()
    report = _report(panel, [row.measured for row in panel.rows])
    verdict = report["verdict"]
    head = report["head"]
    assert isinstance(verdict, dict) and isinstance(head, dict)
    # The three conditions are reported separately, so a reader can see WHICH one failed.
    assert set(verdict) >= {
        "gate_passed", "beats_every_baseline", "interval_is_informative", "promotable"
    }
    assert verdict["promotable"] == (
        verdict["gate_passed"]
        and verdict["beats_every_baseline"]
        and verdict["interval_is_informative"]
    )


# --- panel scoring ------------------------------------------------------------


def test_score_panel_uses_one_invocation_per_utr_context() -> None:
    # A predictor carries its UTR context on the model, so a multi-transcript panel
    # genuinely needs one predictor each -- but no more than that.
    rows = [
        PanelRow(group="P1", variant_id="a", cds="ATGAAATAA", measured=1.0,
                 utr5="GCCACC", utr3="GCTAAT"),
        PanelRow(group="P1", variant_id="b", cds="ATGAAGTAA", measured=2.0,
                 utr5="GCCACC", utr3="GCTAAT"),
        PanelRow(group="P2", variant_id="c", cds="ATGGATTAA", measured=3.0,
                 utr5="AAACCC", utr3="TTTGGG"),
    ]
    panel = panel_from_rows(rows)
    scores, notes = score_panel(
        panel, "null", species="human", cell_types=(), top_k=5,
        batch_size=64, num_workers=0,
    )
    assert len(scores) == 3
    assert "2 UTR context(s) => 2 backend invocation(s)" in notes[0]


def test_score_panel_returns_scores_in_panel_order() -> None:
    # Bucketing by context reorders internally; the caller must get panel order back or
    # every measurement would be paired with the wrong sequence.
    rows = [
        PanelRow(group="P1", variant_id=f"v{i}", cds="ATGAAATAA", measured=float(i),
                 utr5="GCCACC" if i % 2 else "AAACCC", utr3="GCTAAT")
        for i in range(6)
    ]
    panel = panel_from_rows(rows)
    scores, _notes = score_panel(
        panel, "null", species="human", cell_types=(), top_k=5,
        batch_size=64, num_workers=0,
    )
    assert scores == [0.0] * 6  # the null placeholder, but one per row and in order


# --- CLI ----------------------------------------------------------------------


def test_cli_runs_end_to_end_on_the_null_backend(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(4):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    code = rg.main(
        ["--panel", str(path), "--backend", "null", "--within-group",
         "--bootstrap-resamples", "50", "--json"]
    )
    assert code == 0
    report = json.loads(capsys.readouterr().out)
    assert report["panel_summary"]["n_rows"] == 16
    assert len(report["panel_hash"]) == 64
    assert report["verdict"]["promotable"] is False  # the placeholder cannot promote
    assert report["backend"]["calibrated"] is False
    assert set(report["baselines"]) == set(BASELINES)


def test_cli_warns_when_run_in_pooled_mode(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Pooled mode credits between-protein skill, which is not BT4's regime. Running it
    # is allowed (it is a useful contrast) but must never be silent.
    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(3):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    code = rg.main(["--panel", str(path), "--backend", "null", "--bootstrap-resamples", "0"])
    assert code == 0
    assert "NOT the regime BT4 deploys in" in capsys.readouterr().err


def test_cli_records_the_scope_it_ran_with(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # The record used to omit the cell-type selection, top_k and the UTR context
    # entirely, so a finished run could not be reconstructed from its own output -- and
    # an attestation could later declare a scope the record had no way to contradict.
    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(4):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    code = rg.main(
        ["--panel", str(path), "--backend", "null", "--within-group",
         "--cell-type", "HEK293T", "--top-k", "3",
         "--bootstrap-resamples", "50", "--json"]
    )
    assert code == 0
    scope = json.loads(capsys.readouterr().out)["scope"]
    assert scope["cell_types"] == ["HEK293T"]
    assert scope["top_k"] == 3
    assert len(scope["utr_context_sha256"]) == 1  # one (utr5, utr3) pair in this panel
    assert scope["scoring_source"] == "gate"


def test_attest_refuses_a_run_that_is_not_promotable(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # --attest is not a way to record a verdict you did not get.
    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(4):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out = tmp_path / "attestation.json"

    code = rg.main(
        ["--panel", str(path), "--backend", "null", "--within-group",
         "--bootstrap-resamples", "50", "--attest", str(out)]
    )
    assert code == 3
    assert "not promotable" in capsys.readouterr().err
    assert not out.exists()


def test_attest_writes_the_scope_the_comparison_actually_used(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The success path of --attest, driven without licensed weights.

    The CLI cannot reach it (the placeholder backend can never be promotable), so the
    helper is called directly against an oracle comparison. What this pins is the point
    of the whole change: the record's scope comes from the comparison, so it cannot be
    typed into disagreeing with the run.
    """
    panel = _panel(n_groups=10)
    comparison = run_panel_gate(
        panel,  # type: ignore[arg-type]
        "null",
        # The module-wide _SETTINGS run at coverage_tolerance=0.20, which is looser than
        # the attestation layer's own floor -- so a record built from them is refused.
        # That refusal is correct and load-bearing, hence a tighter tolerance here.
        settings=GateSettings(
            within_group=True, recalibrate=True, coverage_tolerance=0.10,
            bootstrap_resamples=200,
        ),
        cell_types=("HEK293T",),
        top_k=3,
        head_scores=[row.measured for row in panel.rows],  # type: ignore[attr-defined]
    )
    assert comparison.promotable is True

    # A gate run against the neutral placeholder cannot be filed as a RiboNN result --
    # the record's backend comes from the head the gate constructed, not from a label.
    # The message is checked, not just the exit code: this run also has supplied scores,
    # which is a *different* refusal, and an exit code alone cannot tell them apart.
    assert rg._attest(comparison, str(tmp_path / "no.json"), readout="mrl") == 3
    assert "not an attestable expression head" in capsys.readouterr().err
    assert not (tmp_path / "no.json").exists()

    # A test double standing in for a gate-scored RiboNN run (no licensed weights in CI).
    comparison = dataclasses.replace(
        comparison,
        backend="ribonn[human]",
        scope=dataclasses.replace(comparison.scope, scoring_source="gate"),
    )

    out = tmp_path / "attestation.json"
    assert rg._attest(comparison, str(out), readout="mean_ribosome_load") == 0
    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["cell_types"] == ["HEK293T"]
    assert record["top_k"] == 3
    assert record["panel_sha256"] == comparison.panel_hash
    assert record["scoring_source"] == "gate"
    assert record["backend"] == "ribonn"
    assert "$BT4_EXPRESSION_ATTESTATION" in capsys.readouterr().err
    # This panel declares no readout column, so the readout was taken on the
    # maintainer's word -- and the record SAYS so rather than implying it was checked.
    assert record["verified_against_panel"] == []

    # Give the panel its own readout column and a disagreeing declaration is refused
    # rather than written.
    declared = panel_from_rows(
        [
            PanelRow(
                group=row.group, variant_id=row.variant_id, cds=row.cds,
                measured=row.measured, utr5=row.utr5, utr3=row.utr3,
                readout="mean_ribosome_load",
            )
            for row in panel.rows  # type: ignore[attr-defined]
        ]
    )
    checked = dataclasses.replace(
        run_panel_gate(
            declared,
            "null",
            settings=GateSettings(
                within_group=True, recalibrate=True, coverage_tolerance=0.10,
                bootstrap_resamples=200,
            ),
            cell_types=("HEK293T",),
            head_scores=[row.measured for row in declared.rows],
        ),
        backend="ribonn[human]",
    )
    checked = dataclasses.replace(
        checked, scope=dataclasses.replace(checked.scope, scoring_source="gate")
    )
    other = tmp_path / "lie.json"
    assert rg._attest(checked, str(other), readout="something_else") == 3
    assert not other.exists()
    assert rg._attest(checked, str(other), readout=None) == 0
    assert json.loads(other.read_text(encoding="utf-8"))["verified_against_panel"] == [
        "readout"
    ]


def test_cli_refuses_an_unknown_baseline(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "panel.tsv"
    path.write_text(
        "group\tvariant_id\tcds\tmeasured\tutr5\tutr3\n"
        "P1\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTAAT\n",
        encoding="utf-8",
    )
    # Validation lives in the pipeline now, so this exits 2 with a message rather than
    # raising out of argparse -- the same shape as every other bad input.
    code = rg.main(["--panel", str(path), "--backend", "null", "--baselines", "cai,bogus"])
    assert code == 2


def test_cli_reports_a_bad_panel_without_a_traceback(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "panel.tsv"
    path.write_text(
        "group\tvariant_id\tcds\tmeasured\tutr5\tutr3\n"
        "P1\tv1\tATGAAATA\t1.0\tGCCACC\tGCTAAT\n",  # not length-3N
        encoding="utf-8",
    )
    code = rg.main(["--panel", str(path), "--backend", "null"])
    assert code == 2
    assert "not a multiple of 3" in capsys.readouterr().err


# --- the CLI subcommand -------------------------------------------------------


def test_bt4_expression_gate_subcommand_runs(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # The gate must be reachable without writing Python, and the CLI and the script
    # must agree about what a result means -- which they do because both render the
    # same GateComparison from bt4.pipeline.
    from bt4.cli.__main__ import main as cli_main

    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(4):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    code = cli_main(
        ["expression-gate", str(path), "--backend", "null", "--within-group",
         "--bootstrap-resamples", "50"]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "PROMOTABLE on this panel : False" in out
    assert "UNCALIBRATED" in out
    assert "flips nothing" in out


def test_bt4_expression_gate_warns_in_pooled_mode(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from bt4.cli.__main__ import main as cli_main

    path = tmp_path / "panel.tsv"
    lines = ["group\tvariant_id\tcds\tmeasured\tutr5\tutr3"]
    for g in range(4):
        for v in range(3):
            lines.append(
                f"P{g}\tg{g}v{v}\tATG{'AAA' * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    assert cli_main(["expression-gate", str(path), "--backend", "null",
                     "--bootstrap-resamples", "0"]) == 0
    assert "NOT the regime BT4 deploys in" in capsys.readouterr().err


def test_bt4_expression_gate_reports_a_bad_panel_cleanly(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from bt4.cli.__main__ import main as cli_main

    path = tmp_path / "panel.tsv"
    path.write_text(
        "group\tvariant_id\tcds\tmeasured\tutr5\tutr3\n"
        "P1\tv1\tATGAAATA\t1.0\tGCCACC\tGCTAAT\n",
        encoding="utf-8",
    )
    assert cli_main(["expression-gate", str(path), "--backend", "null"]) == 2
    assert "not a multiple of 3" in capsys.readouterr().err


def test_run_panel_gate_refuses_mismatched_head_scores() -> None:
    panel = _panel(n_groups=2, n_variants=2)
    with pytest.raises(ValueError, match="head_scores has 3 values for 4 rows"):
        run_panel_gate(panel, "null", head_scores=[1.0, 2.0, 3.0])


def test_run_panel_gate_refuses_an_unknown_baseline() -> None:
    panel = _panel(n_groups=2, n_variants=2)
    with pytest.raises(ValueError, match="unknown baseline"):
        run_panel_gate(panel, "null", baselines=["cai", "bogus"], head_scores=[1.0] * 4)
