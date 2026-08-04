"""Tests for the incremental simulated-annealing refinement engine.

These pin the load-bearing invariants of :func:`bt4.optimize.anneal_refine`:

* invariant #1 - synonymous-only moves keep the round-trip exact;
* invariant #5 - refinement from a feasible seed never adds a hard violation;
* invariant #7 - a fixed seed yields byte-identical output;
* quality - on a small additive-CAI instance, SA started from a poor seed reaches
  the exact-DP optimum (:func:`bt4.optimize.solve_exact`).

The objective is built from the real ``CaiTerm`` / codon-usage table so the
delta-scoring path is exercised exactly as callers use it.
"""

from __future__ import annotations

from collections.abc import Sequence

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4.biomodels.codon.tables import load_table
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint
from bt4.domain.certificate import OptimalityStatus
from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.objectives.terms import CaiTerm
from bt4.optimize import InfeasibleError, anneal_refine, solve_exact

# Degenerate amino acids only (exclude M, W) so proposals always have a choice.
_DEGENERATE = "ACDEFGHIKLNPQRSTVY"

_TABLE = load_table("human")
_CAI = CaiTerm(_TABLE.relative_adaptiveness())


def _cai_delta(prefix: str, codon: str, pos: int) -> float:
    """Additive CAI log-w delta (the exact-DP / SA scalar objective)."""
    return _CAI.delta(prefix, codon, pos)


def _cai_swap_delta(dna: str, pos: int, old_codon: str, new_codon: str) -> float:
    """Incremental score of one synonymous swap for the additive CAI term."""
    return _CAI.delta("", new_codon, pos) - _CAI.delta("", old_codon, pos)


def _worst_seed(residues: Sequence[str]) -> str:
    """Build the lowest-CAI synonymous back-translation (a deliberately poor seed)."""
    return "".join(
        min(synonymous_codons(res), key=lambda c: _CAI.delta("", c, 0)) for res in residues
    )


def _hard_violations(dna: str, constraints: Sequence[Constraint]) -> int:
    """Count HARD whole-sequence violations across ``constraints``."""
    return sum(
        1
        for constraint in constraints
        for v in constraint.validate(dna)
        if v.severity is Severity.HARD
    )


def test_roundtrip_preserved_after_refinement() -> None:
    protein = "MAKKLRSDEF"
    residues = [*protein, "*"]
    constraints = (HomopolymerConstraint(4),)
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)

    result = anneal_refine(
        seed.dna,
        residues,
        _CAI.score,
        constraints,
        iterations=500,
        seed=3,
        delta_score=_cai_swap_delta,
    )

    assert translate(result.dna) == protein + "*"
    assert result.dna[-3:] in synonymous_codons("*")


@settings(max_examples=60, deadline=None)
@given(
    protein=st.text(alphabet=_DEGENERATE, min_size=2, max_size=9),
    seed=st.integers(min_value=0, max_value=12),
    use_delta=st.booleans(),
)
def test_never_increases_hard_violations(protein: str, seed: int, use_delta: bool) -> None:
    residues = [*protein, "*"]
    constraints = (HomopolymerConstraint(4), ForbiddenMotifConstraint(("GAATTC",)))
    try:
        start = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)
    except InfeasibleError:
        return  # no feasible seed for this protein under these constraints

    assert _hard_violations(start.dna, constraints) == 0  # seed is feasible

    result = anneal_refine(
        start.dna,
        residues,
        _CAI.score,
        constraints,
        iterations=200,
        seed=seed,
        delta_score=_cai_swap_delta if use_delta else None,
    )

    # Invariant #5: refinement from a feasible seed never adds a hard violation.
    assert _hard_violations(result.dna, constraints) == 0
    # Invariant #1: still a valid back-translation.
    assert translate(result.dna) == protein + "*"


def test_determinism_same_seed_identical_output() -> None:
    protein = "MARKLESDEQVK"
    residues = [*protein, "*"]
    constraints = (HomopolymerConstraint(3),)
    start = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)

    def run() -> str:
        return anneal_refine(
            start.dna,
            residues,
            _CAI.score,
            constraints,
            iterations=400,
            seed=7,
            delta_score=_cai_swap_delta,
            temp0=0.8,
            cooling=0.99,
        ).dna

    assert run() == run()


def test_certificate_is_heuristic_never_proven() -> None:
    residues = [*"MAKL", "*"]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    result = anneal_refine(
        seed.dna, residues, _CAI.score, (), iterations=50, delta_score=_cai_swap_delta
    )
    assert result.certificate.status is OptimalityStatus.HEURISTIC
    assert result.certificate.solver == "anneal_refine"
    assert not result.certificate.is_proven_optimal


@pytest.mark.parametrize("use_delta", [True, False])
def test_quality_reaches_exact_dp_optimum(use_delta: bool) -> None:
    # Additive CAI objective, no hard constraints: the exact-DP optimum is the
    # per-position argmax, which SA hill-climbing must reach from a poor seed.
    protein = "AKLRSVDEFH"
    residues = [*protein, "*"]
    exact = solve_exact(residues, scalar_delta=_cai_delta, constraints=())

    seed_dna = _worst_seed(residues)
    seed_score = _CAI.score(seed_dna)
    assert seed_score < exact.objective_scalar  # the seed really is poor

    result = anneal_refine(
        seed_dna,
        residues,
        _CAI.score,
        (),
        iterations=4000,
        seed=1,
        delta_score=_cai_swap_delta if use_delta else None,
        temp0=0.5,
        cooling=0.997,
    )

    # Never worse than the seed, and it reaches the proven optimum.
    assert result.objective_scalar >= seed_score - 1e-9
    assert result.objective_scalar >= exact.objective_scalar - 1e-9
    assert result.objective_scalar == pytest.approx(exact.objective_scalar)


def test_zero_iterations_returns_feasible_seed() -> None:
    residues = [*"MAKL", "*"]
    constraints = (HomopolymerConstraint(3),)
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)
    result = anneal_refine(
        seed.dna, residues, _CAI.score, constraints, iterations=0, delta_score=_cai_swap_delta
    )
    assert result.dna == seed.dna
    assert result.certificate.status is OptimalityStatus.HEURISTIC


def test_seed_not_translating_to_residues_raises() -> None:
    residues = [*"MAK", "*"]
    # Correct length (4 codons) but translates to MAA* != MAK*, so the round-trip
    # check - not the length check - must fire.
    with pytest.raises(ValueError, match="translate"):
        anneal_refine("ATGGCTGCTTAA", residues, _CAI.score, ())


def test_seed_wrong_length_raises() -> None:
    residues = [*"MAK", "*"]
    with pytest.raises(ValueError, match="length"):
        anneal_refine("ATGGCT", residues, _CAI.score, ())


def test_infeasible_seed_raises() -> None:
    # Seed translates correctly (KK*) but "AAAAAA" is a run of 6 > max_run 3, so
    # the new seed-feasibility guard must reject it up front (invariant #5 made
    # unconditional), rather than silently refining from an infeasible start.
    residues = [*"KK", "*"]
    with pytest.raises(ValueError, match="feasible"):
        anneal_refine("AAAAAATAA", residues, _CAI.score, (HomopolymerConstraint(3),))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"iterations": -1}, "iterations"),
        ({"temp0": -1.0}, "temp0"),
        ({"cooling": 0.0}, "cooling"),
        ({"cooling": 1.5}, "cooling"),
    ],
)
def test_out_of_range_arguments_raise(kwargs: dict[str, float], match: str) -> None:
    residues = [*"MAKL", "*"]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    with pytest.raises(ValueError, match=match):
        anneal_refine(seed.dna, residues, _CAI.score, (), **kwargs)  # type: ignore[arg-type]


def test_greedy_temp0_zero_never_decreases_objective() -> None:
    # temp0 == 0 is greedy hill-climbing: only improving-or-flat moves accepted,
    # so the delivered objective must be >= the (poor) seed's.
    protein = "AKLRSVDEFH"
    residues = [*protein, "*"]
    seed_dna = _worst_seed(residues)
    result = anneal_refine(
        seed_dna,
        residues,
        _CAI.score,
        (),
        iterations=3000,
        seed=2,
        delta_score=_cai_swap_delta,
        temp0=0.0,
    )
    assert result.objective_scalar >= _CAI.score(seed_dna) - 1e-9


def test_default_new_kwargs_are_a_no_op_for_the_trajectory() -> None:
    """Passing the block/tempering knobs at their defaults changes nothing (#7).

    ``block_size=1``, ``block_prob=0.0``, ``replicas=1``, ``swap_every=0`` must
    reproduce the single-chain trajectory byte-for-byte, so no caller that never
    touches the new arguments can see a different result.
    """
    protein = "MARKLESDEQVK"
    residues = [*protein, "*"]
    constraints = (HomopolymerConstraint(3),)
    start = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)

    common = dict(
        iterations=400, seed=7, delta_score=_cai_swap_delta, temp0=0.8, cooling=0.99
    )
    baseline = anneal_refine(start.dna, residues, _CAI.score, constraints, **common)  # type: ignore[arg-type]
    explicit = anneal_refine(
        start.dna, residues, _CAI.score, constraints,
        block_size=1, block_prob=0.0, replicas=1, temps=None, swap_every=0,
        **common,  # type: ignore[arg-type]
    )
    assert explicit.dna == baseline.dna
    assert explicit.objective_scalar == baseline.objective_scalar


@settings(max_examples=40, deadline=None)
@given(
    protein=st.text(alphabet=_DEGENERATE, min_size=3, max_size=9),
    seed=st.integers(min_value=0, max_value=10),
)
def test_block_tempering_determinism_and_roundtrip(protein: str, seed: int) -> None:
    """Block moves + parallel tempering stay deterministic (#7) and round-trip (#1)."""
    residues = [*protein, "*"]
    constraints = (HomopolymerConstraint(4),)
    try:
        start = solve_exact(residues, scalar_delta=_cai_delta, constraints=constraints)
    except InfeasibleError:
        return

    def run() -> str:
        return anneal_refine(
            start.dna,
            residues,
            _CAI.score,
            constraints,
            iterations=150,
            seed=seed,
            block_size=3,
            block_prob=0.4,
            replicas=3,
            temps=(0.2, 0.6, 1.5),
            swap_every=10,
        ).dna

    first = run()
    assert first == run()  # byte-identical across runs (#7)
    assert translate(first) == protein + "*"  # still a valid back-translation (#1)


def test_block_tempering_never_raises_global_hard_count() -> None:
    """Invariant #5 under block moves + tempering: the global count never rises.

    Block moves swap several codons at once and hot replicas accept uphill moves,
    yet every replica passes the global gate against its own count and the delivered
    result is ranked lower-count-first -- so the delivered hard count is <= the
    seed's, and on this instance reaches zero (a coordinated barrier block moves are
    built to cross).
    """
    from bt4.constraints.max_repeat import MaxRepeatConstraint
    from bt4.domain import STOP

    protein = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLM"
    residues = [*protein, STOP]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    mr = MaxRepeatConstraint(8)

    def _hard(dna: str) -> int:
        return sum(1 for v in mr.validate(dna) if v.severity is Severity.HARD)

    seed_hard = _hard(seed.dna)

    def score(dna: str) -> float:
        return _CAI.score(dna) - 1e9 * _hard(dna)

    result = anneal_refine(
        seed.dna,
        residues,
        score,
        (),
        global_constraints=(mr,),
        iterations=1500,
        seed=1,
        block_size=4,
        block_prob=0.5,
        replicas=3,
        temps=(0.1, 0.5, 2.0),
        swap_every=8,
    )
    assert translate(result.dna) == "".join(residues)
    assert _hard(result.dna) <= seed_hard
    assert result.certificate.status is OptimalityStatus.HEURISTIC


def test_feasibility_floor_immovable_repeat_reported_not_hidden() -> None:
    """A repeat pinned to synonymously-immovable codons is an honest residual.

    Poly-methionine has only ``ATG`` per residue, so *no* synonymous scheme -- not a
    block move, not a hot replica -- can break the ``ATGATG...`` repeat. Refinement
    must leave the violation in place (count unchanged from the seed), never claim it
    clean.
    """
    from bt4.constraints.max_repeat import MaxRepeatConstraint
    from bt4.domain import STOP

    residues = [*("M" * 8), STOP]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    mr = MaxRepeatConstraint(4)

    def _hard(dna: str) -> int:
        return sum(1 for v in mr.validate(dna) if v.severity is Severity.HARD)

    seed_hard = _hard(seed.dna)
    assert seed_hard > 0  # the poly-Met repeat really is present

    def score(dna: str) -> float:
        return _CAI.score(dna) - 1e9 * _hard(dna)

    result = anneal_refine(
        seed.dna,
        residues,
        score,
        (),
        global_constraints=(mr,),
        iterations=800,
        seed=3,
        block_size=4,
        block_prob=0.6,
        replicas=2,
        temps=(0.3, 1.2),
        swap_every=5,
    )
    # Unremovable by construction: the residual is still reported, never hidden.
    assert _hard(result.dna) == seed_hard
    assert translate(result.dna) == "".join(residues)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"block_size": 0}, "block_size"),
        ({"block_prob": -0.1}, "block_prob"),
        ({"block_prob": 1.5}, "block_prob"),
        ({"replicas": 0}, "replicas"),
        ({"replicas": 2, "temps": (0.1,)}, "temps"),
        ({"replicas": 2, "temps": (0.1, -1.0)}, "temps"),
    ],
)
def test_block_tempering_argument_validation(kwargs: dict[str, object], match: str) -> None:
    residues = [*"MAKL", "*"]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    with pytest.raises(ValueError, match=match):
        anneal_refine(seed.dna, residues, _CAI.score, (), **kwargs)  # type: ignore[arg-type]


def test_global_constraints_never_raise_hard_count() -> None:
    """Invariant #5, global edition: a non-local constraint's hard count never rises.

    The dispersed :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint` cannot be
    checked by a bounded ``ok_suffix`` window, so ``anneal_refine`` re-counts its
    whole-sequence violations per move and rejects any that would increase them.
    Starting from a repeat-dirty seed and *rewarding* repeat removal, the count
    must fall or stay flat every step -- never climb.
    """
    from bt4.constraints.max_repeat import MaxRepeatConstraint
    from bt4.domain import STOP

    protein = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLM"
    residues = [*protein, STOP]
    seed = solve_exact(residues, scalar_delta=_cai_delta, constraints=())
    mr = MaxRepeatConstraint(8)

    def _hard(dna: str) -> int:
        return sum(1 for v in mr.validate(dna) if v.severity is Severity.HARD)

    seed_hard = _hard(seed.dna)

    def score(dna: str) -> float:
        return _CAI.score(dna) - 1e9 * _hard(dna)

    result = anneal_refine(
        seed.dna,
        residues,
        score,
        (),
        global_constraints=(mr,),
        iterations=4000,
        seed=1,
    )
    assert translate(result.dna) == "".join(residues)
    # The count can only fall or stay flat; on this instance it reaches zero.
    assert _hard(result.dna) <= seed_hard
    assert _hard(result.dna) == 0
    assert result.certificate.status is OptimalityStatus.HEURISTIC
