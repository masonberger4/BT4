"""Tests for the per-site risk tracks surface (api.tracks / bt4 tracks).

Pins that the tracks are honest sliding-window reporting profiles: correct
lengths and value ranges, the %MinMax track present only when codon-aligned, and
that the whole thing is deterministic and flows through the CLI.
"""

from __future__ import annotations

import json

import pytest

from bt4 import api
from bt4.cli.__main__ import main


def _delivered() -> str:
    return api.optimize("MAALKHETQWSNDECFGRPVIY").dna


def test_tracks_present_and_correct_lengths() -> None:
    dna = _delivered()
    result = api.tracks(dna, nt_window=30, codon_window=9)
    names = {t.name for t in result.tracks}
    assert {"gc_fraction", "cpg_density", "minmax"} <= names
    gc = result.get("gc_fraction")
    assert gc is not None
    assert len(gc.values) == len(dna) - 30 + 1  # sliding-window count
    minmax = result.get("minmax")
    assert minmax is not None
    assert len(minmax.values) == len(dna) // 3 - 9 + 1


def test_track_value_ranges() -> None:
    result = api.tracks(_delivered(), nt_window=30, codon_window=9)
    gc = result.get("gc_fraction")
    cpg = result.get("cpg_density")
    minmax = result.get("minmax")
    assert gc is not None and cpg is not None and minmax is not None
    assert all(0.0 <= v <= 1.0 for v in gc.values)
    assert all(0.0 <= v <= 1.0 for v in cpg.values)
    # %MinMax is in [-100, 100]; allow float-rounding epsilon at the extremes.
    assert all(-100.0 - 1e-9 <= v <= 100.0 + 1e-9 for v in minmax.values)


def test_minmax_track_omitted_when_not_codon_aligned() -> None:
    # A non-multiple-of-3 sequence still gets GC/CpG but no codon-based %MinMax.
    result = api.tracks("ACGTACGTAC", nt_window=4)
    names = {t.name for t in result.tracks}
    assert "gc_fraction" in names and "cpg_density" in names
    assert "minmax" not in names


def test_tracks_deterministic() -> None:
    dna = _delivered()
    first = api.tracks(dna)
    second = api.tracks(dna)
    assert [t.values for t in first.tracks] == [t.values for t in second.tracks]


def test_tracks_rejects_non_acgt() -> None:
    with pytest.raises(ValueError):
        api.tracks("ACGTX")


def test_cli_tracks_summary(capsys: pytest.CaptureFixture[str]) -> None:
    dna = _delivered()
    assert main(["tracks", dna, "--nt-window", "30"]) == 0
    out = capsys.readouterr().out
    assert "gc_fraction" in out
    assert "cpg_density" in out


def test_cli_tracks_json(capsys: pytest.CaptureFixture[str]) -> None:
    dna = _delivered()
    assert main(["tracks", dna, "--json", "--nt-window", "30"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["dna"] == dna
    names = {t["name"] for t in payload["tracks"]}
    assert "gc_fraction" in names
    gc = next(t for t in payload["tracks"] if t["name"] == "gc_fraction")
    assert gc["window"] == 30 and gc["window_unit"] == "nt"
    assert len(gc["values"]) == len(dna) - 30 + 1
