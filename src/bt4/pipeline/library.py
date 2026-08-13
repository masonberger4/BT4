"""Library / degenerate-design mode: a sampled library, honestly stamped.

This composes the codon table, the LOCAL constraints, and the deterministic
:func:`~bt4.optimize.sample.sample_sequences` sampler into a *library* of coding
sequences for one protein (CLAUDE.md §9, Phase 5). It is the sampling counterpart
of :func:`~bt4.pipeline.optimize.run_optimize`:

* it reuses the exact-core machinery -- ``_build_constraints``,
  ``_active_terms``, and ``_make_result`` -- so every member round-trips
  (invariant #1) and carries metrics **recomputed from its own DNA** (invariant
  #2) and a content-addressed manifest (invariant #9);
* but the certificate is :attr:`~bt4.domain.OptimalityStatus.SAMPLED`: this is a
  stochastic draw, not an optimum, and it claims **neither** optimality **nor**
  an expression prediction (CLAUDE.md §1, §10.6);
* GLOBAL constraints (``max_repeat_length``, ``avoid_uorf``) are **not** enforced
  during sampling. Instead every output is validated against the full local +
  global set and any residual violation is reported honestly in that member's
  ``violations`` -- the library never silently claims a clean non-local audit.

It imports only :mod:`bt4.domain`, the biomodel/optimize layers, and helpers from
its sibling :mod:`bt4.pipeline.optimize`; nothing above ``pipeline``.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from bt4.biomodels.codon.tables import load_table
from bt4.domain import (
    STOP,
    OptimalityCertificate,
    OptimalityStatus,
    Result,
    validate_protein,
)
from bt4.optimize.sample import sample_sequences
from bt4.pipeline.optimize import (
    OptimizeConfig,
    _active_terms,
    _build_constraints,
    _build_global_constraints,
    _make_result,
    _manifest,
)
from bt4.provenance import Manifest

__all__ = ["LibraryResult", "run_library"]

# Default library size when a caller does not name one.
_DEFAULT_N = 8


@dataclass(frozen=True, slots=True)
class LibraryResult:
    """A sampled library of sequences plus its provenance and diversity stats.

    Attributes:
        results: The sampled sequences, each a full :class:`~bt4.domain.Result`
            with recomputed metrics, a ``SAMPLED`` certificate, and any residual
            (including non-local) violations reported honestly.
        manifest: The run's content-addressed provenance stamp. It is identical
            across members (they share config, seed, ``n``, and temperature), so
            the whole library is reproducible from this one stamp.
        distinct: Number of distinct sequences among ``results`` -- an honest
            diversity readout (a sampler may repeat a draw).
        mean_pairwise_hamming: Mean over all member pairs of the fraction of
            positions at which two sequences differ (``0.0`` for a single member).
    """

    results: tuple[Result, ...]
    manifest: Manifest
    distinct: int
    mean_pairwise_hamming: float


def _diversity(dnas: list[str]) -> tuple[int, float]:
    """Return ``(distinct count, mean pairwise Hamming fraction)`` for ``dnas``.

    All members share the same length (they encode the same protein), so the
    Hamming distance is well defined; the fraction normalizes it to ``[0, 1]``.
    With fewer than two members there are no pairs, so the mean is ``0.0``.
    """
    distinct = len(set(dnas))
    if len(dnas) < 2:
        return distinct, 0.0
    length = len(dnas[0])
    total = 0.0
    pairs = 0
    for i in range(len(dnas)):
        for j in range(i + 1, len(dnas)):
            diff = sum(1 for a, b in zip(dnas[i], dnas[j], strict=True) if a != b)
            total += (diff / length) if length else 0.0
            pairs += 1
    return distinct, total / pairs


def run_library(
    protein: str,
    config: OptimizeConfig | None = None,
    n: int = _DEFAULT_N,
    *,
    seed: int | None = None,
    temperature: float = 1.0,
) -> LibraryResult:
    """Sample a library of ``n`` sequences for ``protein`` (stochastic, not optimal).

    Each residue is back-translated by sampling its synonymous codons in
    proportion to the organism's usage frequencies (raised to ``1 / temperature``),
    keeping only codons that satisfy every LOCAL constraint's ``ok_suffix``. The
    delivered sequences are **sampled, not optimized**: their certificate is
    :attr:`~bt4.domain.OptimalityStatus.SAMPLED` and they carry no optimality or
    expression claim. GLOBAL constraints are not enforced during sampling but are
    validated and reported per member.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`. Objective
            weights do not steer the draw (this is a sampler); only the codon
            table and the LOCAL constraints shape it. ``gc_min``/``gc_max`` and
            ``refine`` do not apply to sampling and are ignored.
        n: Number of sequences to sample (``>= 1``).
        seed: Master sampling seed; when ``None`` the run uses ``config.seed``.
            The effective seed enters the manifest, so the library is reproducible
            from its stamp (invariants #7/#9).
        temperature: Sampling temperature (``> 0``). ``-> 0`` approaches the
            per-residue argmax, ``1.0`` is the natural distribution, large values
            approach uniform.

    Returns:
        A :class:`LibraryResult`: the sampled members, a shared manifest, and
        honest diversity statistics.

    Raises:
        ValueError: On an invalid protein, unknown organism, ``n < 1``, or
            ``temperature <= 0``.
        bt4.optimize.InfeasibleError: If the LOCAL constraints admit no feasible
            sequence.
    """
    config = config or OptimizeConfig()
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if temperature <= 0.0:
        raise ValueError(f"temperature must be > 0, got {temperature}")
    effective_seed = config.seed if seed is None else seed
    # Thread the effective seed through the config so the manifest, the audit, and
    # the sampler all agree on it (a stamp must reproduce the library it names).
    if effective_seed != config.seed:
        config = dataclasses.replace(config, seed=effective_seed)

    p = validate_protein(protein)
    table = load_table(config.organism, reference_set=config.reference_set)
    terms = [term for term, _ in _active_terms(table, config)]
    constraints = _build_constraints(config)
    global_constraints = _build_global_constraints(config)
    # Local + global together are the reporting set: validate() audits both, so a
    # residual non-local violation (a max-repeat / uORF the sampler did not avoid)
    # surfaces honestly in each member's violations and hard-violation count.
    report_constraints = [*constraints, *global_constraints]
    residues = [*p, STOP]

    dnas = sample_sequences(
        residues,
        codon_weights=dict(table.frequency),
        constraints=constraints,
        n=n,
        seed=effective_seed,
        temperature=temperature,
    )

    # The manifest records that this was a library run and on what terms, so two
    # libraries differing only in n or temperature stamp differently (invariant #9).
    manifest_extra: dict[str, object] = {
        "mode": "library",
        "library_n": n,
        "temperature": temperature,
    }
    certificate = OptimalityCertificate(
        status=OptimalityStatus.SAMPLED,
        solver="library_sampler",
        detail=(
            f"stochastic draw from the codon distribution (temperature={temperature}); "
            "no optimality or expression claim"
        ),
    )
    results = tuple(
        _make_result(
            protein=p,
            dna=dna,
            table=table,
            terms=terms,
            constraints=report_constraints,
            certificate=certificate,
            config=config,
            alpha=None,
            extra_audit=dict(manifest_extra),
            manifest_extra=manifest_extra,
        )
        for dna in dnas
    )
    distinct, mean_hamming = _diversity(list(dnas))
    # The top-level manifest is built by the same content-addressed path each
    # member's manifest used, so it is the single stamp the whole library
    # reproduces from (every member shares config, seed, n, and temperature).
    return LibraryResult(
        results=results,
        manifest=_manifest(config, manifest_extra),
        distinct=distinct,
        mean_pairwise_hamming=mean_hamming,
    )
