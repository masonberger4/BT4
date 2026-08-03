"""Tests for the standalone ``scripts/compare_reproducibility.py`` report.

Like ``tests/test_compare_tools.py`` this loads the off-graph script by file
path and pins the honesty-critical contract: the full 93-record Tab 4 panel
loads over three human proteins, each anonymized algorithm aggregates its ten
repeat runs (n_runs == 10) while Native and BT4 are single deterministic rows
(n_runs == 1, zero run-to-run span), at least one algorithm shows genuine
run-to-run variability, metrics are recomputed from the sequences, and ``--json``
emits valid JSON.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
from types import ModuleType

import pytest


def _load() -> ModuleType:
    path = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "compare_reproducibility.py"
    spec = importlib.util.spec_from_file_location("bt4_compare_reproducibility", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


repro = _load()


def test_panel_loads_three_proteins() -> None:
    records = repro.load_panel()
    assert len(records) == 93  # 3 proteins x (Native + 3 algorithms x 10 runs)
    proteins = {repro.parse_header(h)[0] for h, _ in records}
    assert proteins == {"KRas4B", "Beclin1", "PDE3A"}


def test_header_parsing() -> None:
    assert repro.parse_header("KRas4B|Native|acc=x") == ("KRas4B", "Native", None)
    assert repro.parse_header("KRas4B|Algorithm2|run7|acc=x") == ("KRas4B", "Algorithm2", 7)


def test_algorithms_aggregate_ten_runs_native_and_bt4_are_single() -> None:
    rows = repro.reproducibility(repro.load_panel())
    algo_rows = [r for r in rows if str(r["source"]).startswith("Algorithm")]
    assert algo_rows and all(r["n_runs"] == 10 for r in algo_rows)
    single = [r for r in rows if r["source"] in ("Native", "BT4")]
    assert single and all(r["n_runs"] == 1 for r in single)
    # Every protein contributes a BT4 reference row.
    bt4_proteins = {r["protein"] for r in rows if r["source"] == "BT4"}
    assert bt4_proteins == {"KRas4B", "Beclin1", "PDE3A"}


def test_bt4_is_deterministic_zero_span() -> None:
    rows = repro.reproducibility(repro.load_panel())
    for row in rows:
        if row["source"] == "BT4":
            for metric in ("cai", "tai", "gc_pct", "cpg", "max_homo"):
                assert row[metric]["span"] == pytest.approx(0.0), (row["protein"], metric)


def test_algorithms_show_real_run_to_run_variability() -> None:
    # The point of the panel: at least one stochastic algorithm varies run-to-run.
    rows = repro.reproducibility(repro.load_panel())
    spans = [row["cai"]["span"] for row in rows if str(row["source"]).startswith("Algorithm")]
    assert max(spans) > 0.0


def test_metrics_recomputed_from_sequence() -> None:
    from bt4.biomodels.codon.tables import load_table

    records = repro.load_panel()
    header, dna = records[0]  # KRas4B|Native
    assert repro.parse_header(header)[1] == "Native"
    rows = repro.reproducibility(records)
    native = next(
        r for r in rows if r["protein"] == "KRas4B" and r["source"] == "Native"
    )
    # n=1 so the aggregate mean is the single recomputed value.
    assert native["cai"]["mean"] == pytest.approx(load_table("homo_sapiens").cai(dna))


def test_json_output_is_valid(capsys: pytest.CaptureFixture[str]) -> None:
    assert repro.main(["--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert isinstance(payload, list) and len(payload) == 15  # 3 proteins x 5 sources
    assert all("protein" in row and "source" in row and "cai" in row for row in payload)
