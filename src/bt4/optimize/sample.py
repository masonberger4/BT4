"""Deterministic constrained sampler over the codon trellis (library mode).

BT4's exact DP returns a single most-favored-codon optimum. Library /
degenerate-design mode (CLAUDE.md §9, Phase 5) instead draws a *library* of
sequences by **sampling** each residue's synonymous-codon distribution. This is
a stochastic sampler, **not** an optimizer: nothing here claims optimality, and
the caller stamps every sampled sequence with the ``SAMPLED`` certificate.

The honesty that survives from the exact core is feasibility: at each residue we
keep only the synonymous codons that pass every (LOCAL) constraint's
``ok_suffix`` veto, then sample one with probability proportional to
``weight ** (1 / temperature)``. So a sampled sequence respects the same local
hard constraints the DP does (invariant #3), while its non-local composition is
left to the caller to validate and report honestly.

Determinism (invariant #7) is total: identical ``(residues, codon_weights,
constraints, n, seed, temperature)`` yields a byte-identical library, driven by
stdlib :class:`random.Random` alone.

Layering: this module imports only :mod:`bt4.domain` (plus
:class:`~bt4.optimize.exact_dp.InfeasibleError` from its sibling), never
biomodels or pipeline.
"""

from __future__ import annotations

import random
from collections.abc import Mapping, Sequence

from bt4.domain.contracts import Constraint
from bt4.domain.genetic_code import synonymous_codons
from bt4.optimize.exact_dp import InfeasibleError

__all__ = ["sample_sequences"]

# How many times a single library member may restart from scratch before we
# declare the instance infeasible. A restart happens only when some residue has
# no feasible codon under the current prefix; because feasibility can depend on
# earlier (randomly chosen) codons, a fresh draw may succeed where the last one
# stalled. If none of this many independent draws completes, no feasible
# assignment is being found and we raise rather than loop forever.
_MAX_RESTARTS = 128


def _choose_codon(
    rng: random.Random,
    codons: Sequence[str],
    weights: Sequence[float],
    inv_temperature: float,
) -> str:
    """Sample one codon ∝ ``weight ** (1 / temperature)`` (numerically stable).

    ``codons`` are the feasible codons for a residue and ``weights`` their
    non-negative sampling weights (codon-usage frequencies). Weights are first
    divided by their maximum so the exponentiation stays in ``(0, 1]`` and cannot
    overflow: at ``temperature -> 0`` this approaches the argmax, ``temperature
    == 1`` reproduces the natural (weight-proportional) distribution, and large
    ``temperature`` approaches uniform. If every feasible weight is zero (or
    missing), fall back to a uniform choice so a feasible codon is never dropped.
    """
    mx = max(weights)
    if mx <= 0.0:
        return codons[rng.randrange(len(codons))]
    scaled = [(w / mx) ** inv_temperature for w in weights]
    total = sum(scaled)
    if total <= 0.0:
        return codons[rng.randrange(len(codons))]
    target = rng.random() * total
    upto = 0.0
    for codon, s in zip(codons, scaled, strict=True):
        upto += s
        if target < upto:
            return codon
    return codons[-1]  # pragma: no cover - float rounding safety net


def _draw_one(
    residues: Sequence[str],
    codon_weights: Mapping[str, float],
    constraints: Sequence[Constraint],
    rng: random.Random,
    inv_temperature: float,
) -> str | None:
    """Draw one full sequence, or ``None`` if a residue admitted no feasible codon.

    The ``rng`` stream advances whether or not the draw succeeds, so retrying
    with the same generator is exactly the "restart with a perturbed seed"
    behaviour -- deterministic, yet different from the stalled attempt.
    """
    dna_parts: list[str] = []
    prefix = ""
    for residue in residues:
        feasible: list[str] = []
        weights: list[float] = []
        for codon in synonymous_codons(residue):
            if all(c.ok_suffix(prefix, codon) for c in constraints):
                feasible.append(codon)
                weights.append(float(codon_weights.get(codon, 0.0)))
        if not feasible:
            return None
        codon = _choose_codon(rng, feasible, weights, inv_temperature)
        dna_parts.append(codon)
        prefix += codon
    return "".join(dna_parts)


def sample_sequences(
    residues: Sequence[str],
    *,
    codon_weights: Mapping[str, float],
    constraints: Sequence[Constraint],
    n: int,
    seed: int,
    temperature: float,
) -> list[str]:
    """Sample ``n`` feasible coding sequences from the codon distribution.

    Each residue is back-translated by sampling one of its synonymous codons that
    passes every constraint's ``ok_suffix`` (LOCAL hard feasibility), with
    probability proportional to ``codon_weights[codon] ** (1 / temperature)``.
    The result is a *library*, not an optimum: no optimality is claimed.

    Args:
        residues: Amino-acid letters to back-translate, **including** a trailing
            ``"*"`` as the final residue (the stop), mirroring
            :func:`~bt4.optimize.exact_dp.solve_exact`.
        codon_weights: Map codon -> non-negative sampling weight (typically the
            organism's per-codon usage frequencies). Missing codons weigh zero.
        constraints: LOCAL hard-feasibility rules; only codons whose ``ok_suffix``
            holds are ever sampled. Non-local (GLOBAL) rules are **not** enforced
            here -- the caller validates and reports those on each output.
        n: Number of sequences to draw (``>= 1``). Members use per-member derived
            seeds so the library is reproducible yet its members differ.
        seed: Master seed; identical inputs and seed reproduce the library
            byte-for-byte (invariant #7).
        temperature: Sampling temperature (``> 0``). ``-> 0`` approaches the
            per-residue argmax, ``1.0`` is the natural distribution, large values
            approach uniform.

    Returns:
        A list of ``n`` coding sequences, each translating back to ``residues``.

    Raises:
        ValueError: If ``n < 1`` or ``temperature <= 0``.
        InfeasibleError: If, after :data:`_MAX_RESTARTS` independent draws, some
            residue still admits no codon under the constraints (naming them).
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if temperature <= 0.0:
        raise ValueError(f"temperature must be > 0, got {temperature}")
    inv_temperature = 1.0 / temperature
    # A master generator spawns one child seed per member, so members are
    # independent and reproducible regardless of how many restarts each needs.
    master = random.Random(seed)
    child_seeds = [master.randrange(2**63) for _ in range(n)]

    library: list[str] = []
    for child_seed in child_seeds:
        rng = random.Random(child_seed)
        dna: str | None = None
        for _attempt in range(_MAX_RESTARTS):
            dna = _draw_one(residues, codon_weights, constraints, rng, inv_temperature)
            if dna is not None:
                break
        if dna is None:
            raise InfeasibleError([c.name for c in constraints])
        library.append(dna)
    return library
