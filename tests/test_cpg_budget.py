"""Tests for the whole-sequence CpG/UpA dinucleotide count budget.

This is the last Phase 2 item (CLAUDE.md §9): an exact, whole-sequence
dinucleotide *count* budget enforced by the same amount-bucketed exact DP that
carries the GC budget. The load-bearing subtlety is that a dinucleotide count is
**not** per-codon decomposable -- a 2-mer can straddle a codon boundary -- so the
budget's per-codon ``amount`` is context-aware, reading the last base of the
prefix exactly as :meth:`DinucleotideTerm.delta` attributes each occurrence to
the codon holding its END base.

The tests pin the honesty invariants: the budget is honored on the recomputed
count (invariant #2), the solver's optimum matches an independent brute-force
enumeration and honestly reports ``PROVEN_OPTIMAL`` (invariant #6), the per-codon
amounts sum to the whole-sequence overlapping count, a two-sided UpA range binds,
the guards refuse an incompatible combination, and the solve is deterministic
(invariant #7).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from itertools import product

import pytest

from bt4 import api
from bt4.constraints.rules import HomopolymerConstraint
from bt4.domain import Severity
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import STOP, synonymous_codons, translate
from bt4.objectives.dinucleotide import DinucleotideTerm, dinucleotide_amount
from bt4.optimize.lagrangian import solve_lagrangian

ScalarDelta = Callable[[str, str, int], float]


def _overlapping_count(dna: str, dinuc: str) -> int:
    """Overlapping occurrences of a 2-mer across the whole sequence."""
    return sum(1 for i in range(len(dna) - 1) if dna[i : i + 2] == dinuc)


def _accumulate(scalar_delta: ScalarDelta, dna: str) -> float:
    """Accumulate ``scalar_delta`` over the codons of ``dna`` (as the DP would)."""
    acc = 0.0
    prefix = ""
    for i in range(0, len(dna), 3):
        codon = dna[i : i + 3]
        acc += scalar_delta(prefix, codon, i // 3)
        prefix += codon
    return acc


def _achievable_counts(residues: Sequence[str], dinuc: str) -> list[int]:
    """Sorted, de-duplicated set of achievable overlapping counts (brute force)."""
    seen = {
        _overlapping_count("".join(combo), dinuc)
        for combo in product(*[synonymous_codons(r) for r in residues])
    }
    return sorted(seen)


def _brute_force_optimum(
    residues: Sequence[str],
    scalar_delta: ScalarDelta,
    dinuc: str,
    constraints: Sequence[Constraint],
    *,
    budget_min: int | None,
    budget_max: int | None,
) -> tuple[float, str]:
    """Best (max-objective, lex-smallest DNA) assignment under budget + constraints.

    Mirrors the solver's objective (accumulated ``scalar_delta``), its budget (the
    overlapping ``dinuc`` count), and its constraint set (filtered by ``validate``,
    which agrees with ``ok_suffix`` by invariant #3), with the solver's own
    deterministic tie-break (lexicographically smaller DNA).
    """
    best_obj: float | None = None
    best_dna: str | None = None
    for combo in product(*[synonymous_codons(r) for r in residues]):
        dna = "".join(combo)
        cnt = _overlapping_count(dna, dinuc)
        if budget_min is not None and cnt < budget_min:
            continue
        if budget_max is not None and cnt > budget_max:
            continue
        if any(v.severity is Severity.HARD for c in constraints for v in c.validate(dna)):
            continue
        obj = _accumulate(scalar_delta, dna)
        if best_obj is None or obj > best_obj or (obj == best_obj and dna < (best_dna or "")):
            best_obj, best_dna = obj, dna
    assert best_dna is not None, "brute-force instance is infeasible"
    return best_obj if best_obj is not None else 0.0, best_dna


# --------------------------------------------------------------------------- #
# Budget honored end-to-end (invariant #2: reported == recomputed).
# --------------------------------------------------------------------------- #


def test_cpg_cap_honored_via_api() -> None:
    protein = "MARPGARSTKLE"
    cap = 3
    result = api.optimize(protein, api.OptimizeConfig(dinuc_budget="CG", dinuc_max=cap))

    count = _overlapping_count(result.dna, "CG")
    assert count <= cap
    assert translate(result.dna) == protein + STOP
    # The reported count is recomputed from the DNA and matches (invariant #2).
    assert result.audit["cg_count"] == count
    # A dinucleotide budget uses the amount-bucketed exact DP -> proven optimal.
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert result.certificate.solver == "lagrangian"
    assert result.metrics.hard_violations == 0


def test_cpg_floor_honored_via_api() -> None:
    # Arginine (CGN) and proline (CCN) give ready CpG, so a floor is reachable.
    protein = "MRPRPRPRA"
    floor = 3
    result = api.optimize(protein, api.OptimizeConfig(dinuc_budget="CG", dinuc_min=floor))

    count = _overlapping_count(result.dna, "CG")
    assert count >= floor
    assert translate(result.dna) == protein + STOP
    assert result.audit["cg_count"] == count
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


# --------------------------------------------------------------------------- #
# Proven-optimal + exactness (invariant #6): the solver matches brute force.
# --------------------------------------------------------------------------- #


def test_proven_optimal_matches_brute_force() -> None:
    # Small enough to enumerate every synonymous assignment. The objective rewards
    # CpG (elevate) while a CpG *cap* pushes back, so the budget genuinely binds; a
    # homopolymer rule is a local constraint the budget DP must keep honoring.
    residues = [*"MARPG", STOP]
    dinuc = "CG"
    scalar_delta = DinucleotideTerm(dinuc, "elevate").delta
    constraints: list[Constraint] = [HomopolymerConstraint(4)]
    budget_max = 2

    result = solve_lagrangian(
        residues,
        scalar_delta=scalar_delta,
        constraints=constraints,
        amount=dinucleotide_amount(dinuc),
        budget_max=budget_max,
        objective_context=1,
        budget_context=1,
        budget_name="cg_budget",
    )

    best_obj, best_dna = _brute_force_optimum(
        residues, scalar_delta, dinuc, constraints, budget_min=None, budget_max=budget_max
    )
    # Same optimal objective value and a budget-feasible, constraint-feasible seq.
    assert result.objective_scalar == pytest.approx(best_obj)
    assert _overlapping_count(result.dna, dinuc) <= budget_max
    assert translate(result.dna) == "MARPG" + STOP
    # The lex-smallest optimum is deterministic and matches brute force exactly.
    assert result.dna == best_dna
    # No beam -> honestly proven optimal (invariant #6).
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_proven_optimal_interior_band_matches_brute_force() -> None:
    # A two-sided band that excludes both extremes, with a flat objective so the
    # solver has no gradient -- landing in-band is purely the budget DP's doing.
    residues = [*"MRPRA", STOP]
    dinuc = "CG"
    counts = _achievable_counts(residues, dinuc)
    assert len(counts) >= 3, counts
    band_lo, band_hi = counts[1], counts[-2]

    def flat(prefix: str, codon: str, pos: int) -> float:
        return 0.0

    result = solve_lagrangian(
        residues,
        scalar_delta=flat,
        constraints=[],
        amount=dinucleotide_amount(dinuc),
        budget_min=band_lo,
        budget_max=band_hi,
        budget_context=1,
        budget_name="cg_budget",
    )

    assert band_lo <= _overlapping_count(result.dna, dinuc) <= band_hi
    assert translate(result.dna) == "MRPRA" + STOP
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


# --------------------------------------------------------------------------- #
# Sum-of-amounts == whole-sequence overlapping count (the counting invariant).
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("dinuc", ["CG", "TA"])
def test_sum_of_amounts_equals_whole_sequence_count(dinuc: str) -> None:
    # A hand-built sequence containing boundary-straddling occurrences.
    dna = "TGCGATCGCCGGTATACG"
    amount = dinucleotide_amount(dinuc)
    acc = 0
    prefix = ""
    for i in range(0, len(dna), 3):
        codon = dna[i : i + 3]
        acc += amount(prefix, codon)
        prefix += codon
    direct = _overlapping_count(dna, dinuc)
    assert acc == direct
    # And it matches the (unsigned) magnitude of the elevate term's score.
    assert float(acc) == DinucleotideTerm(dinuc, "elevate").score(dna)


# --------------------------------------------------------------------------- #
# UpA range (both bounds) via the public API.
# --------------------------------------------------------------------------- #


def test_upa_range_both_bounds() -> None:
    residues = [*"MYLIV", STOP]
    dinuc = "TA"
    counts = _achievable_counts(residues, dinuc)
    assert len(counts) >= 3, counts
    band_lo, band_hi = counts[1], counts[-2]

    def flat(prefix: str, codon: str, pos: int) -> float:
        return 0.0

    result = solve_lagrangian(
        residues,
        scalar_delta=flat,
        constraints=[],
        amount=dinucleotide_amount(dinuc),
        budget_min=band_lo,
        budget_max=band_hi,
        budget_context=1,
        budget_name="ta_budget",
    )

    assert band_lo <= _overlapping_count(result.dna, dinuc) <= band_hi
    assert translate(result.dna) == "MYLIV" + STOP
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_upa_range_via_api_reports_count() -> None:
    protein = "MYLIVFYLIV"
    result = api.optimize(protein, api.OptimizeConfig(dinuc_budget="TA", dinuc_min=1, dinuc_max=4))
    count = _overlapping_count(result.dna, "TA")
    assert 1 <= count <= 4
    assert result.audit["ta_count"] == count
    assert translate(result.dna) == protein + STOP
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


# --------------------------------------------------------------------------- #
# Guards: one budget at a time, and no budget with a refinement-enforced rule.
# --------------------------------------------------------------------------- #


def test_gc_and_dinuc_budget_are_mutually_exclusive() -> None:
    with pytest.raises(ValueError, match="cannot be combined"):
        api.optimize("MAAL", api.OptimizeConfig(gc_min=5, dinuc_budget="CG", dinuc_max=2))


def test_dinuc_budget_incompatible_with_refine() -> None:
    with pytest.raises(ValueError, match="refine"):
        api.optimize("MAAL", api.OptimizeConfig(dinuc_budget="CG", dinuc_max=2, refine=True))


def test_dinuc_budget_incompatible_with_max_repeat() -> None:
    config = api.OptimizeConfig(dinuc_budget="CG", dinuc_max=2, max_repeat_length=8)
    with pytest.raises(ValueError, match="not supported together with"):
        api.optimize("MAAL", config)


def test_dinuc_budget_incompatible_with_uorf() -> None:
    with pytest.raises(ValueError, match="not supported together with"):
        api.optimize("MAAL", api.OptimizeConfig(dinuc_budget="CG", dinuc_max=2, avoid_uorf=True))


@pytest.mark.parametrize("bad", ["C", "CGT", "", "CX", "XY", "CGA"])
def test_bad_dinuc_budget_raises(bad: str) -> None:
    with pytest.raises(ValueError):
        api.optimize("MAAL", api.OptimizeConfig(dinuc_budget=bad, dinuc_max=2))


def test_dinuc_bounds_without_budget_raise() -> None:
    with pytest.raises(ValueError, match="require dinuc_budget"):
        api.optimize("MAAL", api.OptimizeConfig(dinuc_min=1, dinuc_max=3))


# --------------------------------------------------------------------------- #
# Determinism (invariant #7).
# --------------------------------------------------------------------------- #


def test_determinism_same_config_twice() -> None:
    config = api.OptimizeConfig(dinuc_budget="CG", dinuc_max=4)
    first = api.optimize("MARPGARSTKLE", config)
    second = api.optimize("MARPGARSTKLE", config)
    assert first.dna == second.dna
    assert first.audit["cg_count"] == second.audit["cg_count"]
