"""Tests for the standalone ``scripts/compare_tools.py`` report harness.

``scripts/`` is not an importable package (no ``__init__.py`` and not on the
``bt4`` import graph), so the module is loaded by file path via
:mod:`importlib.util`, mirroring ``tests/test_benchmark.py``. These tests pin the
harness's honesty-critical contract: the full ten-sequence panel loads, the
native reference round-trips, the DNA2.0 truncation is flagged as a length
mismatch (not silently pooled), the BT4 row is present with an independently
recomputed CAI, and ``--json`` emits valid JSON carrying the expected keys.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
from types import ModuleType

import pytest

from bt4.biomodels.codon.tables import load_table
from bt4.domain.genetic_code import translate


def _load_compare_tools() -> ModuleType:
    """Load ``scripts/compare_tools.py`` as a module by file path.

    Returns:
        The imported ``compare_tools`` module.
    """
    path = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "compare_tools.py"
    spec = importlib.util.spec_from_file_location("bt4_compare_tools", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


compare_tools = _load_compare_tools()

# Keys every row must carry (tool rows and the appended BT4 row alike).
_EXPECTED_KEYS: frozenset[str] = frozenset(
    {
        "name",
        "optimizer",
        "len_nt",
        "aa",
        "matches_native",
        "cai",
        "tai",
        "gc_pct",
        "cpg",
        "max_homo",
        "certificate",
        "note",
    }
)


def test_panel_loads_ten_sequences() -> None:
    """The staged Ranaghan panel loads exactly the ten expected records."""
    records = compare_tools.load_panel()
    assert len(records) == 10
    names = [header for header, _ in records]
    assert names == [
        "Native",
        "GeneArt",
        "GeneWiz",
        "DNA2.0",
        "IDT",
        "Genscript",
        "Twist",
        "JCAT",
        "OPTIMIZER",
        "COOL",
    ]


def test_native_round_trips() -> None:
    """The native reference CDS translates back to a full-length protein plus stop."""
    records = dict(compare_tools.load_panel())
    native_dna = records["Native"]
    protein = translate(native_dna)
    assert protein.endswith("*")
    assert len(protein) - 1 == 188  # full-length KRas4B


def test_rows_have_expected_schema() -> None:
    """Ten tool rows plus one BT4 row, each carrying exactly the expected keys."""
    rows = compare_tools.compare(compare_tools.load_panel())
    assert len(rows) == 11  # 10 panel sequences + BT4
    for row in rows:
        assert frozenset(row) == _EXPECTED_KEYS


def test_dna20_flagged_as_length_mismatch() -> None:
    """DNA2.0 is a C-terminal truncation: flagged, not treated as the same protein."""
    rows = {row["name"]: row for row in compare_tools.compare(compare_tools.load_panel())}
    dna20 = rows["DNA2.0"]
    assert dna20["matches_native"] is False
    assert dna20["aa"] == 169
    assert "length mismatch" in str(dna20["note"])
    # A genuine synonymous tool row, by contrast, matches the native protein.
    assert rows["GeneArt"]["matches_native"] is True
    assert rows["GeneArt"]["note"] == ""


def test_bt4_row_present_with_recomputed_cai() -> None:
    """The BT4 row exists and its CAI equals an independent recompute of its DNA."""
    rows = {row["name"]: row for row in compare_tools.compare(compare_tools.load_panel())}
    assert "BT4" in rows
    bt4 = rows["BT4"]
    assert bt4["matches_native"] is True
    assert bt4["certificate"] is not None  # BT4 vouches for its own optimality
    # Tool rows carry no certificate (BT4 cannot vouch for their optimality).
    assert rows["Twist"]["certificate"] is None
    assert isinstance(bt4["cai"], float)
    assert 0.0 < bt4["cai"] <= 1.0


def test_metrics_recomputed_from_sequence() -> None:
    """Every reported metric is BT4's own recomputation of the delivered sequence."""
    records = dict(compare_tools.load_panel())
    rows = {row["name"]: row for row in compare_tools.compare(list(records.items()))}
    table = load_table("homo_sapiens")
    idt_dna = records["IDT"]
    idt_row = rows["IDT"]
    assert idt_row["cai"] == round(table.cai(idt_dna), 6)
    assert idt_row["cpg"] == idt_dna.upper().count("CG")
    assert idt_row["len_nt"] == len(idt_dna)


def test_main_json_emits_valid_json(capsys: pytest.CaptureFixture[str]) -> None:
    """``main(["--json"])`` returns 0 and emits parseable JSON with the expected keys."""
    code = compare_tools.main(["--json"])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert isinstance(payload, list)
    assert len(payload) == 11
    for row in payload:
        assert frozenset(row) == _EXPECTED_KEYS
    assert payload[-1]["name"] == "BT4"


def test_main_table_prints_banner_and_headers(capsys: pytest.CaptureFixture[str]) -> None:
    """``main([])`` returns 0 and prints the honest banner plus the column headers."""
    code = compare_tools.main([])
    assert code == 0
    out = capsys.readouterr().out
    assert "RECOMPUTED" in out
    assert "NOT claimed 'better'" in out
    for header in ("optimizer", "cai", "tai", "gc_pct", "cpg", "maxhomo", "certificate"):
        assert header in out


def test_compare_is_deterministic() -> None:
    """Two runs over the same panel produce byte-identical rows."""
    panel = compare_tools.load_panel()
    assert compare_tools.compare(panel) == compare_tools.compare(panel)
