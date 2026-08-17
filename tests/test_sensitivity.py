"""Tests for the standalone sensitivity/uncertainty analysis script.

The script lives under ``scripts/`` (outside the ``bt4`` import graph), so it is
loaded by file path here. These tests pin the load-bearing honesty behaviors:
the analysis runs for a small protein; changing the organism actually moves at
least one delivered metric (a spread, not a constant); a beam=1 budget honestly
degrades the certificate to ``beam_truncated`` while the exact DP stays
``proven_optimal``; and ``--json`` emits valid JSON.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

_SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "sensitivity.py"

# A short protein with genuinely degenerate residues, so beams prune and tables
# disagree.
_PROTEIN = "MKTAYIAKQR"


def _load_script() -> ModuleType:
    """Load ``scripts/sensitivity.py`` as a module by file path."""
    spec = importlib.util.spec_from_file_location("bt4_sensitivity_script", _SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def script() -> ModuleType:
    """The loaded sensitivity script module."""
    return _load_script()


def test_analyze_runs_for_small_protein(script: ModuleType) -> None:
    """The full analysis produces both axes with at least one usable organism."""
    report = script.analyze(_PROTEIN)
    org = report["organism_sensitivity"]
    bud = report["budget_sensitivity"]
    assert report["protein"] == _PROTEIN
    assert org["rows"], "expected at least one usable codon-usage table"
    assert bud["rows"], "expected at least the exact solve plus beams"
    # Every organism row recomputes CAI and GC% from the delivered sequence.
    for row in org["rows"]:
        assert 0.0 <= float(row["cai"]) <= 1.0
        assert 0.0 <= float(row["gc_percent"]) <= 100.0


def test_organism_choice_changes_a_metric(script: ModuleType) -> None:
    """Sweeping organisms moves at least one delivered metric (a real spread)."""
    org = script.organism_sensitivity(_PROTEIN)
    rows = org["rows"]
    assert len(rows) >= 2, "need >= 2 usable tables to observe a spread"
    ranges = [
        float(org["spread"][metric]["range"])
        for metric in ("cai", "gc_percent", "tai_logw")
        if org["spread"][metric]["range"] is not None
    ]
    assert any(r > 0.0 for r in ranges), "organism choice changed no metric at all"


def test_skipped_entries_recorded_not_dropped(script: ModuleType) -> None:
    """Non-usable tables (e.g. tRNA-only *.trna) are recorded with a reason."""
    org = script.organism_sensitivity(_PROTEIN)
    # The bundled organism discovery surfaces *.trna tables that are not codon
    # usage tables; they must be reported as skipped, never silently dropped.
    for entry in org["skipped"]:
        assert entry["organism"]
        assert entry["reason"]


def test_beam_one_truncates_exact_proven(script: ModuleType) -> None:
    """beam=1 yields a beam_truncated certificate; exact yields proven_optimal."""
    bud = script.budget_sensitivity(_PROTEIN, beams=(1,))
    by_budget = {row["budget"]: row for row in bud["rows"]}
    assert by_budget["exact"]["certificate"] == "proven_optimal"
    assert by_budget["beam=1"]["certificate"] == "beam_truncated"
    assert bud["certificate_degrades"] is True


def test_tai_logw_available_for_human(script: ModuleType) -> None:
    """tAI is recomputed through the API for an organism shipping a tRNA table."""
    result = script.api.optimize(_PROTEIN, script.api.OptimizeConfig(organism="homo_sapiens"))
    value = script.tai_logw(result.dna, "homo_sapiens")
    assert value is not None
    # E. coli now ships a tRNA table too (it was the last gap), so tAI resolves
    # there as well; only an organism with no bundled table returns None.
    assert script.tai_logw(result.dna, "escherichia_coli") is not None
    assert script.tai_logw(result.dna, "nonexistent_organism") is None


def test_json_output_is_valid(script: ModuleType, capsys: pytest.CaptureFixture[str]) -> None:
    """``--json`` prints valid, round-trippable JSON with both axes present."""
    rc = script.main(["--protein", _PROTEIN, "--beams", "1,2", "--json"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["protein"] == _PROTEIN
    assert payload["organism_sensitivity"]["axis"] == "organism"
    assert payload["budget_sensitivity"]["axis"] == "budget"
    assert "honesty_note" in payload


def test_text_output_runs(script: ModuleType, capsys: pytest.CaptureFixture[str]) -> None:
    """The default (non-JSON) render prints both sections and the honesty note."""
    rc = script.main(["--protein", _PROTEIN])
    assert rc == 0
    out = capsys.readouterr().out
    assert "organism sensitivity" in out
    assert "solver-budget sensitivity" in out
    assert "SPREAD" in out
