"""Integration tests for the BT3 ``almost-there`` parity features.

Covers the three knobs BT3's mature branch had that BT4 was missing, now wired
end-to-end through the pipeline: the "max GC length" (GC-run) local constraint,
the "max repeat length" non-local constraint (refinement-enforced, honestly
reported), and named forbidden-sequence presets. The load-bearing honesty checks:

* a GC-run limit is enforced exactly (proven-optimal, zero residual);
* a max-repeat limit is driven down by refinement without ever raising the
  hard-violation count (invariant #5) and its residual is reported truthfully;
* the certificate degrades to heuristic only when refinement actually runs.
"""

from __future__ import annotations

import pytest

from bt4 import api
from bt4.constraints.forbidden import available_forbidden_presets, resolve_forbidden_motifs
from bt4.constraints.max_repeat import MaxRepeatConstraint
from bt4.domain import translate
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.result import Severity
from bt4.pipeline.optimize import OptimizeConfig

_PROTEIN = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"


def _hard(dna: str, constraint: object) -> int:
    return sum(1 for v in constraint.validate(dna) if v.severity is Severity.HARD)  # type: ignore[attr-defined]


# --------------------------------------------------------------------------- #
# Max GC length (GcRunConstraint) -- local, exact DP.
# --------------------------------------------------------------------------- #


def test_max_gc_run_is_enforced_exactly() -> None:
    result = api.optimize(_PROTEIN, OptimizeConfig(max_gc_run=4))
    assert translate(result.dna) == _PROTEIN + "*"
    # No run of > 4 consecutive G/C bases anywhere.
    assert not any(v.constraint == "gc_run" for v in result.violations)
    # A pure local constraint keeps the proven-optimal certificate.
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_max_gc_run_actually_limits_runs() -> None:
    result = api.optimize(_PROTEIN, OptimizeConfig(max_gc_run=3))
    dna = result.dna
    longest = 0
    run = 0
    for ch in dna:
        run = run + 1 if ch in "GC" else 0
        longest = max(longest, run)
    assert longest <= 3


# --------------------------------------------------------------------------- #
# Max repeat length (MaxRepeatConstraint) -- non-local, refinement-enforced.
# --------------------------------------------------------------------------- #


def test_max_repeat_refinement_drives_repeats_down() -> None:
    limit = 8
    seed = api.optimize(_PROTEIN, OptimizeConfig())
    mr = MaxRepeatConstraint(limit)
    seed_hard = _hard(seed.dna, mr)
    seed_was_dirty = seed_hard > 0

    result = api.optimize(_PROTEIN, OptimizeConfig(max_repeat_length=limit))
    assert translate(result.dna) == _PROTEIN + "*"
    residual = _hard(result.dna, mr)
    # Refinement never makes repeats worse; on this protein it removes them all.
    assert residual <= seed_hard
    assert residual == 0
    assert result.audit["max_repeat_enforced"] == "clean"
    assert result.audit["max_repeat_residual"] == 0
    if seed_was_dirty:
        # Refinement ran, so the certificate is honestly heuristic.
        assert result.certificate.status is OptimalityStatus.HEURISTIC


def test_max_repeat_clean_seed_stays_proven_optimal() -> None:
    # A loose limit the exact-DP seed already satisfies needs no refinement.
    result = api.optimize(_PROTEIN, OptimizeConfig(max_repeat_length=20))
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert result.audit["max_repeat_enforced"] == "clean"


def test_max_repeat_partial_is_reported_honestly() -> None:
    # A highly repetitive protein forces repeats no synonymous choice can remove.
    repetitive = "MKAILVDEQTR" * 4
    result = api.optimize(repetitive, OptimizeConfig(max_repeat_length=5))
    assert translate(result.dna) == repetitive + "*"
    # Honesty: if repeats remain they are reported, not hidden.
    residual = result.audit["max_repeat_residual"]
    assert result.audit["max_repeat_enforced"] == ("clean" if residual == 0 else "partial")
    reported = sum(1 for v in result.violations if v.constraint == "max_repeat")
    assert reported == residual


def test_max_repeat_with_gc_budget_raises() -> None:
    with pytest.raises(ValueError, match="not supported together with a GC budget"):
        api.optimize("MKAILV", OptimizeConfig(max_repeat_length=6, gc_min=5, gc_max=9))


# --------------------------------------------------------------------------- #
# Forbidden-sequence presets.
# --------------------------------------------------------------------------- #


def test_presets_catalog_resolves_motifs() -> None:
    keys = [p.key for p in available_forbidden_presets()]
    assert "poly_a_signal" in keys and "tata_box" in keys
    motifs = resolve_forbidden_motifs(("poly_a_signal", "tata_box"))
    assert "AATAAA" in motifs and "ATTAAA" in motifs and "TATAAA" in motifs


def test_unknown_preset_raises() -> None:
    with pytest.raises(KeyError, match="unknown forbidden preset"):
        resolve_forbidden_motifs(("not_a_preset",))


def test_preset_motifs_are_banned_in_output() -> None:
    result = api.optimize(
        _PROTEIN, OptimizeConfig(forbidden_presets=("poly_a_signal", "tata_box"))
    )
    for motif in ("AATAAA", "ATTAAA", "TATAAA"):
        assert motif not in result.dna
    assert not any(v.constraint == "forbidden_motif" for v in result.violations)


def test_presets_merge_with_custom_motifs() -> None:
    result = api.optimize(
        _PROTEIN,
        OptimizeConfig(forbidden_motifs=("GAATTC",), forbidden_presets=("poly_a_signal",)),
    )
    assert "GAATTC" not in result.dna
    assert "AATAAA" not in result.dna


# --------------------------------------------------------------------------- #
# Determinism: the refinement path is seeded and reproducible.
# --------------------------------------------------------------------------- #


def test_max_repeat_run_is_deterministic() -> None:
    cfg = OptimizeConfig(max_repeat_length=8)
    a = api.optimize(_PROTEIN, cfg)
    b = api.optimize(_PROTEIN, cfg)
    assert a.dna == b.dna
