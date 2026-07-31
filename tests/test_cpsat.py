"""Tests for the CP-SAT (OR-Tools) additive-plus-GC-budget solver backend.

The whole module is skipped when OR-Tools is not installed, so CI jobs without
the ``[ilp]`` extra stay green. Each behavioural claim is checked against a
brute-force enumeration over every codon assignment, which is tractable for the
tiny proteins used here and gives an independent optimum to compare against.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Sequence

import pytest

pytest.importorskip("ortools")

from bt4._accel import gc_count
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.genetic_code import synonymous_codons, translate
from bt4.optimize import InfeasibleError
from bt4.optimize.cpsat import solve_cpsat

CodonScore = Callable[[str, int], float]

# Distinct per-codon weights so optima are unambiguous for the test proteins.
WEIGHTS: dict[str, float] = {
    "ATG": 0.0,  # M - sole codon
    "GCT": 0.15, "GCC": 0.91, "GCA": 0.33, "GCG": 0.22,  # A
    "AAA": 0.42, "AAG": 0.88,  # K
    "TGT": 0.51, "TGC": 0.73,  # C
    "TAA": 0.64, "TAG": 0.27, "TGA": 0.19,  # stop
}


def make_codon_score(weights: dict[str, float]) -> CodonScore:
    """Build a context-free per-codon reward closure (larger is better)."""

    def score(codon: str, pos: int) -> float:
        return weights.get(codon, 0.0)

    return score


def _brute_force(
    residues: Sequence[str],
    score: CodonScore,
    *,
    gc_min: int | None = None,
    gc_max: int | None = None,
) -> tuple[float, set[str]]:
    """Enumerate every assignment; return (max score, argmax DNAs) under budget."""
    choices = [synonymous_codons(r) for r in residues]
    best = float("-inf")
    optima: set[str] = set()
    for combo in itertools.product(*choices):
        dna = "".join(combo)
        gc = gc_count(dna)
        if gc_min is not None and gc < gc_min:
            continue
        if gc_max is not None and gc > gc_max:
            continue
        total = sum(score(codon, pos) for pos, codon in enumerate(combo))
        if total > best:
            best = total
            optima = {dna}
        elif total == best:
            optima.add(dna)
    return best, optima


@pytest.mark.parametrize("protein", ["MA", "KK", "AK", "CC", "MAK"])
def test_brute_force_equivalence_no_budget(protein: str) -> None:
    residues = [*protein, "*"]
    score = make_codon_score(WEIGHTS)
    best, optima = _brute_force(residues, score)

    result = solve_cpsat(residues, codon_score=score)

    assert result.objective_scalar == pytest.approx(best)
    assert result.dna in optima
    assert translate(result.dna) == f"{protein}*"


def test_gc_budget_is_respected_and_optimal() -> None:
    residues = [*"AAK", "*"]
    score = make_codon_score(WEIGHTS)

    # The tightest still-feasible upper bound: the lowest attainable GC count.
    min_gc = sum(
        min(gc_count(codon) for codon in synonymous_codons(r)) for r in residues
    )
    best, optima = _brute_force(residues, score, gc_max=min_gc)

    result = solve_cpsat(residues, codon_score=score, gc_max=min_gc)

    assert gc_count(result.dna) <= min_gc
    assert result.objective_scalar == pytest.approx(best)
    assert result.dna in optima
    assert translate(result.dna) == "AAK*"


def test_gc_lower_bound_is_respected_and_optimal() -> None:
    residues = [*"AAK", "*"]
    score = make_codon_score(WEIGHTS)

    max_gc = sum(
        max(gc_count(codon) for codon in synonymous_codons(r)) for r in residues
    )
    best, optima = _brute_force(residues, score, gc_min=max_gc)

    result = solve_cpsat(residues, codon_score=score, gc_min=max_gc)

    assert gc_count(result.dna) >= max_gc
    assert result.objective_scalar == pytest.approx(best)
    assert result.dna in optima


def test_infeasible_gc_budget_raises() -> None:
    residues = [*"AAK", "*"]
    score = make_codon_score(WEIGHTS)
    max_gc = sum(
        max(gc_count(codon) for codon in synonymous_codons(r)) for r in residues
    )

    with pytest.raises(InfeasibleError) as excinfo:
        solve_cpsat(residues, codon_score=score, gc_min=max_gc + 1)

    assert "gc_min" in excinfo.value.constraints


def test_certificate_proven_optimal_for_small_instance() -> None:
    result = solve_cpsat(
        [*"MAK", "*"], codon_score=make_codon_score(WEIGHTS), max_time_s=30.0
    )
    assert result.certificate.is_proven_optimal
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert result.certificate.solver == "cpsat"
    assert result.certificate.gap == 0.0


def test_determinism_identical_solves_match() -> None:
    residues = [*"MAKCA", "*"]
    score = make_codon_score(WEIGHTS)
    first = solve_cpsat(residues, codon_score=score)
    second = solve_cpsat(residues, codon_score=score)
    assert first.dna == second.dna
    assert first.objective_scalar == second.objective_scalar
