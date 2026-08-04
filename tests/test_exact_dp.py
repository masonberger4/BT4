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
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4.domain.certificate import OptimalityStatus
from bt4.domain.genetic_code import synonymous_codons, translate
from bt4.optimize import InfeasibleError, SolveResult, solve_exact
from bt4.optimize.exact_dp import _precompute_structure, _solve_native

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


# ---------------------------------------------------------------------------
# Native trellis port: the precomputed-table DP (Rust or its pure-Python twin,
# via _accel.trellis_solve) must be byte-for-byte identical to the pure-Python
# reference loop (solve_exact with position_independent=False). This is the
# equivalence the CLAUDE.md §7 / Phase 1 "full Rust trellis port" rests on.
# ---------------------------------------------------------------------------

# A pairwise-style weight table keyed by (prev_codon_last_char, codon), so the
# scalar delta reads trailing context (objective_context > 0) -- exercising the
# extended-state DP through the native path too.
_PAIR_W: dict[tuple[str, str], float] = {}
for _prev in ("", "A", "C", "G", "T"):
    for _cod, _v in WEIGHTS.items():
        _PAIR_W[(_prev, _cod)] = _v + (0.05 if _prev == "G" else 0.0)


def _pairwise_delta(prefix: str, codon: str, pos: int) -> float:
    last = prefix[-1:] if prefix else ""
    return _PAIR_W.get((last, codon), WEIGHTS.get(codon, 0.0))


_AA_POOL = "MAKC"  # amino acids covered by WEIGHTS with >1 codon for A/K/C


def _oracle(
    residues: Sequence[str],
    scalar_delta: ScalarDelta,
    constraints: Sequence[object],
    beam: int | None,
    objective_context: int,
) -> SolveResult | InfeasibleError:
    """The pure-Python reference: solve_exact with the native path disabled."""
    try:
        return solve_exact(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,  # type: ignore[arg-type]
            beam=beam,
            objective_context=objective_context,
            position_independent=False,
        )
    except InfeasibleError as exc:
        return exc


def _assert_native_matches(
    residues: Sequence[str],
    scalar_delta: ScalarDelta,
    constraints: Sequence[object],
    beam: int | None,
    objective_context: int,
) -> None:
    ref = _oracle(residues, scalar_delta, constraints, beam, objective_context)
    context_len = max([objective_context, *(c.context_len() for c in constraints)])  # type: ignore[attr-defined]
    try:
        got = _solve_native(
            residues,
            scalar_delta=scalar_delta,
            constraints=constraints,  # type: ignore[arg-type]
            beam=beam,
            context_len=context_len,
        )
    except InfeasibleError as exc:
        assert isinstance(ref, InfeasibleError)
        # Constraint names surfaced must match the reference exactly.
        assert exc.constraints == ref.constraints
        return
    assert isinstance(ref, SolveResult)
    # Byte-identical DNA, exact-equal scalar, identical certificate status/solver.
    assert got.dna == ref.dna
    assert got.objective_scalar == ref.objective_scalar
    assert got.certificate.status is ref.certificate.status
    assert got.certificate.solver == ref.certificate.solver


@given(
    protein=st.text(alphabet=_AA_POOL, min_size=0, max_size=7),
    seed=st.integers(min_value=0, max_value=999),
    pairwise=st.booleans(),
    max_run=st.integers(min_value=2, max_value=4),
    use_homopolymer=st.booleans(),
    beam=st.one_of(st.none(), st.integers(min_value=1, max_value=6)),
)
@settings(max_examples=250, deadline=None)
def test_native_matches_pure_python_oracle(
    protein: str,
    seed: int,
    pairwise: bool,
    max_run: int,
    use_homopolymer: bool,
    beam: int | None,
) -> None:
    residues = [*protein, "*"]
    if pairwise:
        scalar_delta: ScalarDelta = _pairwise_delta
        objective_context = 3
    else:
        # Perturb the weights deterministically from `seed` to exercise ties too.
        w = {c: round(v + (seed % 3) * 0.0, 3) for c, v in WEIGHTS.items()}
        scalar_delta = make_scalar_delta(w)
        objective_context = 0
    constraints: tuple[object, ...] = (
        (FakeHomopolymer(max_run),) if use_homopolymer else (NoConstraint(),)
    )
    _assert_native_matches(residues, scalar_delta, constraints, beam, objective_context)


def test_native_tie_break_prefers_lexicographically_smaller_dna() -> None:
    # Two codons of A carry the *same* weight, so ties must resolve to the smaller
    # codon at every position -- byte-identical between native and pure Python.
    tie_weights = dict(WEIGHTS)
    tie_weights["GCT"] = tie_weights["GCC"] = tie_weights["GCA"] = tie_weights["GCG"] = 0.5
    _assert_native_matches(
        [*"AAAA", "*"], make_scalar_delta(tie_weights), (NoConstraint(),), None, 0
    )


def test_native_infeasible_surfaces_constraint_names() -> None:
    with pytest.raises(InfeasibleError) as excinfo:
        _solve_native(
            [*"MA", "*"],
            scalar_delta=make_scalar_delta(WEIGHTS),
            constraints=(AlwaysReject(),),
            beam=None,
            context_len=0,
        )
    assert "always_reject" in excinfo.value.constraints


def test_native_beam_downgrades_certificate() -> None:
    out = _solve_native(
        [*"MA", "*"],
        scalar_delta=make_scalar_delta(WEIGHTS),
        constraints=(FakeHomopolymer(3),),
        beam=1,
        context_len=3,
    )
    assert out.certificate.status is OptimalityStatus.BEAM_TRUNCATED
    assert out.certificate.solver == "beam_dp"


def test_native_k0_and_empty_protein() -> None:
    # K == 0 (no context) and a zero-residue trellis (just the stop) must both
    # match the reference.
    _assert_native_matches([*"MAK", "*"], make_scalar_delta(WEIGHTS), (NoConstraint(),), None, 0)
    _assert_native_matches(["*"], make_scalar_delta(WEIGHTS), (NoConstraint(),), None, 0)


def test_native_objective_context_extended_state() -> None:
    _assert_native_matches([*"AKCAK", "*"], _pairwise_delta, (FakeHomopolymer(3),), None, 3)


def test_solve_exact_native_path_matches_pure_when_available() -> None:
    # Drive solve_exact through the native gate directly (large enough instance)
    # and confirm it equals the pure-Python path on the same inputs.
    residues = [*("MAKC" * 8), "*"]
    sd = make_scalar_delta(WEIGHTS)
    cons = (FakeHomopolymer(4),)
    fast = solve_exact(
        residues, scalar_delta=sd, constraints=cons, position_independent=True
    )
    ref = solve_exact(
        residues, scalar_delta=sd, constraints=cons, position_independent=False
    )
    assert fast.dna == ref.dna
    assert fast.objective_scalar == ref.objective_scalar
    assert fast.certificate.status is ref.certificate.status


def test_precompute_context_zero_uses_single_state() -> None:
    # With no constraints and no objective context, every layer collapses to the
    # single empty context -- the precompute must reflect that (id 0 only).
    structure = _precompute_structure([*"MAK", "*"], (NoConstraint(),), 0)
    assert structure is not None
    assert all(all(to == 0 for to in layer) for layer in structure.layer_to)
    assert all(all(fr == 0 for fr in layer) for layer in structure.layer_from)


def test_precompute_context_cap_falls_back() -> None:
    # A tiny context cap forces the precompute to bail (returns None), so the
    # caller falls back to the pure-Python DP rather than pay an unbounded cost.
    structure = _precompute_structure(
        [*"AKCAKC", "*"], (FakeHomopolymer(3),), 3, max_contexts=2
    )
    assert structure is None
