"""Tests for the standalone ``scripts/benchmark.py`` report harness.

``scripts/`` is not an importable package (no ``__init__.py`` and not on the
``bt4`` import graph), so the module is loaded by file path via
:mod:`importlib.util`. These tests pin the enriched harness's contract: the row
schema (including the CpG / tandem / %MinMax columns), the round-trip honesty of
the naive baseline, independent recomputation of the CpG count, the ``main``
table / JSON entry points, and determinism.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
from types import ModuleType

import pytest

from bt4 import api
from bt4.domain.genetic_code import translate


def _load_benchmark() -> ModuleType:
    """Load ``scripts/benchmark.py`` as a module by file path.

    Returns:
        The imported ``benchmark`` module.
    """
    path = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "benchmark.py"
    spec = importlib.util.spec_from_file_location("bt4_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


benchmark = _load_benchmark()

# The full set of keys every row must carry after the Phase 2 enrichment.
_EXPECTED_KEYS: frozenset[str] = frozenset(
    {
        "name",
        "len_nt",
        "naive_gc",
        "bt4_gc",
        "naive_cpg",
        "bt4_cpg",
        "naive_maxhomo",
        "bt4_maxhomo",
        "naive_maxtandem",
        "bt4_maxtandem",
        "naive_minmax_mean",
        "bt4_minmax_mean",
        "bt4_cai",
        "bt4_certificate",
    }
)


def test_benchmark_rows_have_expected_schema() -> None:
    """One row per panel protein, each carrying exactly the expected keys."""
    rows = benchmark.benchmark(benchmark.DEFAULT_PANEL)
    assert len(rows) == len(benchmark.DEFAULT_PANEL)
    names = [row["name"] for row in rows]
    assert names == list(benchmark.DEFAULT_PANEL)
    for row in rows:
        assert frozenset(row) == _EXPECTED_KEYS


@pytest.mark.parametrize("protein", ["MAAL", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIE"])
def test_naive_backtranslate_round_trips(protein: str) -> None:
    """The naive baseline translates back to the input protein plus a stop."""
    dna = benchmark.naive_backtranslate(protein)
    assert translate(dna) == protein + "*"


def test_cpg_column_matches_independent_recompute() -> None:
    """``bt4_cpg`` equals ``result.dna.count("CG")`` recomputed independently."""
    protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIE"
    row = benchmark.benchmark({"one": protein})[0]
    result = api.optimize(protein, api.OptimizeConfig())
    assert row["bt4_cpg"] == result.dna.count("CG")


def test_minmax_mean_none_when_shorter_than_window() -> None:
    """A sequence shorter than one %MinMax window reports ``None`` for the mean."""
    # "MAAL" is 5 codons with the stop -- far shorter than the 18-codon window.
    row = benchmark.benchmark({"short": "MAAL"})[0]
    assert row["naive_minmax_mean"] is None
    assert row["bt4_minmax_mean"] is None


def test_longest_tandem_span_definition() -> None:
    """The tandem helper matches its documented period-scan definition."""
    assert benchmark._longest_tandem_span("AAAA") == 4  # period-1 homopolymer run
    assert benchmark._longest_tandem_span("ATGATGATG") == 9  # 3x "ATG", period 3
    assert benchmark._longest_tandem_span("ACGT") == 0  # no >=2-copy repeat


def test_main_json_prints_valid_json(capsys: pytest.CaptureFixture[str]) -> None:
    """``main(["--json"])`` returns 0 and emits parseable JSON with the new keys."""
    code = benchmark.main(["--json"])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert isinstance(payload, list) and payload
    assert frozenset(payload[0]) == _EXPECTED_KEYS


def test_main_table_has_new_headers(capsys: pytest.CaptureFixture[str]) -> None:
    """``main([])`` returns 0 and prints a table including the new column headers."""
    code = benchmark.main([])
    assert code == 0
    out = capsys.readouterr().out
    for header in ("naive_cpg", "bt4_cpg", "naive_tand", "bt4_tand", "naive_mm", "bt4_mm"):
        assert header in out


def test_benchmark_is_deterministic() -> None:
    """Two runs over the same panel produce byte-identical rows."""
    first = benchmark.benchmark(benchmark.DEFAULT_PANEL)
    second = benchmark.benchmark(benchmark.DEFAULT_PANEL)
    assert first == second
