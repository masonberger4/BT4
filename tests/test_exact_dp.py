"""Tests for the exact-DP codon-trellis solver.

These tests are deliberately self-contained: they define inline fake objective
and constraint objects rather than importing the concrete ``objectives``/
``constraints`` modules (built concurrently by another agent). The solver only
consumes a ``scalar_delta`` callable and a constraint's ``context_len`` /
``ok_suffix``, so the fakes stub the rest of each protocol trivially.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import pytest

from bt4.domain.certificate import OptimalityStatus
from bt4.domain.genetic_code import synonymous_codons, translate
from bt4.optimize import InfeasibleError, SolveResult, solve_exact

ScalarDelta = Callable[[str, str, int], float]


def _max_run(seq: str) -> int:
    """Return the length of the longest single-character run in ``seq``."""
    best = 0
    run = 0
    prev = ""
    for ch in seq:
        run = run + 1 if ch == prev else 1
        prev = ch
        best = max(best, run)
    return best


def make_scalar_delta(weights: dict[str, float]) -> ScalarDelta:
    """Build a per-codon reward closure (larger is better)."""

    def delta(prefix: str, codon: str, pos: int) -> float:
        return weights.get(codon, 0.0)

    return delta


@dataclass(frozen=True)
class NoConstraint:
    """A constraint that never vetoes anything (context-free)."""

    name: str = "none"

    def scope(self) -> str:
        return "local"

    def context_len(self) -> int:
        return 0

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        return 0.0

    def validate(self, dna: str) -> list[object]:
        return []


@dataclass(frozen=True)
class FakeHomopolymer:
    """Veto any extension that creates a single-base run longer than ``max_run``."""

    max_run: int
    name: str = field(default="homopolymer")

    def scope(self) -> str:
        return "local"

    def context_len(self) -> int:
        return self.max_run

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        window = prefix[-self.max_run :] + next_codon
        return _max_run(window) <= self.max_run

    def penalty(self, prefix: str, next_codon: str) -> float:
        return 0.0

    def validate(self, dna: str) -> list[object]:
        return []


@dataclass(frozen=True)
class AlwaysReject:
    """A constraint whose ``ok_suffix`` always vetoes - forces infeasibility."""

    name: str = "always_reject"

    def scope(self) -> str:
        return "local"

    def context_len(self) -> int:
        return 0

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        return False

    def penalty(self, prefix: str, next_codon: str) -> float:
        return 0.0

    def validate(self, dna: str) -> list[object]:
        return []


# A per-codon weight table covering every codon used by the tiny test proteins,
# with distinct values so optima are meaningful (and some ties are exercised).
WEIGHTS: dict[str, float] = {
    "ATG": 0.0,  # M - sole codon
    "GCT": 0.1, "GCC": 0.9, "GCA": 0.3, "GCG": 0.2,  # A
    "AAA": 0.4, "AAG": 0.8,  # K
    "TGT": 0.5, "TGC": 0.7,  # C
    "TAA": 0.6, "TAG": 0.2, "TGA": 0.1,  # stop
}


def _brute_force(
    residues: Sequence[str],
    scalar_delta: ScalarDelta,
    constraints: Sequence[object],
) -> tuple[float, set[str]]:
    """Enumerate every codon assignment; return (max scalar, feasible optima).

    Feasibility is judged incrementally with ``ok_suffix`` - exactly the rule the
    solver uses - so the two must agree.
    """
    choices = [synonymous_codons(r) for r in residues]
    best_score = float("-inf")
    optima: set[str] = set()
    for combo in itertools.product(*choices):
        prefix = ""
        total = 0.0
        feasible = True
        for pos, codon in enumerate(combo):
            if not all(c.ok_suffix(prefix, codon) for c in constraints):  # type: ignore[attr-defined]
                feasible = False
                break
            total += scalar_delta(prefix, codon, pos)
            prefix += codon
        if not feasible:
            continue
        if total > best_score:
            best_score = total
            optima = {prefix}
        elif total == best_score:
            optima.add(prefix)
    return best_score, optima


@pytest.mark.parametrize(
    "protein",
    ["MA", "AK", "CC", "KK", "MAK"],
)
@pytest.mark.parametrize(
    "constraints",
    [
        (NoConstraint(),),
        (FakeHomopolymer(3),),
        (FakeHomopolymer(2),),
    ],
)
def test_brute_force_equivalence(
    protein: str, constraints: tuple[object, ...]
) -> None:
    residues = [*protein, "*"]
    scalar_delta = make_scalar_delta(WEIGHTS)
    best_score, optima = _brute_force(residues, scalar_delta, constraints)
    result = solve_exact(
        residues,
        scalar_delta=scalar_delta,
        constraints=constraints,  # type: ignore[arg-type]
    )
    assert result.objective_scalar == pytest.approx(best_score)
    assert result.dna in optima


def test_roundtrip_includes_appended_stop() -> None:
    residues = [*"MA", "*"]
    result = solve_exact(
        residues,
        scalar_delta=make_scalar_delta(WEIGHTS),
        constraints=(NoConstraint(),),
    )
    assert translate(result.dna) == "MA*"
    assert result.dna[-3:] in synonymous_codons("*")


def test_certificate_proven_for_exact_solve() -> None:
    result = solve_exact(
        [*"MA", "*"],
        scalar_delta=make_scalar_delta(WEIGHTS),
        constraints=(NoConstraint(),),
    )
    assert isinstance(result, SolveResult)
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL
    assert result.certificate.solver == "exact_dp"


def test_certificate_beam_truncated_when_pruned() -> None:
    # FakeHomopolymer gives K=3, so the A layer holds 4 distinct states; beam=1
    # must drop three of them and record the loss honestly.
    result = solve_exact(
        [*"MA", "*"],
        scalar_delta=make_scalar_delta(WEIGHTS),
        constraints=(FakeHomopolymer(3),),
        beam=1,
    )
    assert result.certificate.status is OptimalityStatus.BEAM_TRUNCATED
    assert result.certificate.solver == "beam_dp"


def test_beam_large_enough_stays_proven() -> None:
    # A beam wider than any layer never prunes, so optimality is preserved.
    result = solve_exact(
        [*"MA", "*"],
        scalar_delta=make_scalar_delta(WEIGHTS),
        constraints=(FakeHomopolymer(3),),
        beam=100,
    )
    assert result.certificate.status is OptimalityStatus.PROVEN_OPTIMAL


def test_infeasible_raises_with_constraint_names() -> None:
    with pytest.raises(InfeasibleError) as excinfo:
        solve_exact(
            [*"MA", "*"],
            scalar_delta=make_scalar_delta(WEIGHTS),
            constraints=(AlwaysReject(),),
        )
    assert "always_reject" in excinfo.value.constraints
    assert "always_reject" in str(excinfo.value)


def test_determinism_identical_inputs_identical_output() -> None:
    scalar_delta = make_scalar_delta(WEIGHTS)
    constraints = (FakeHomopolymer(3),)
    first = solve_exact([*"MAK", "*"], scalar_delta=scalar_delta, constraints=constraints)
    second = solve_exact([*"MAK", "*"], scalar_delta=scalar_delta, constraints=constraints)
    assert first.dna == second.dna
    assert first.objective_scalar == second.objective_scalar
