"""The run orchestration: compose biomodels + objectives + constraints + solver.

This layer turns a protein and an :class:`OptimizeConfig` into an honest
:class:`~bt4.domain.result.Result` (single solve) or a Pareto
:class:`FrontierResult` (a scalarization sweep of exact solves). It is where the
honesty invariants are enforced in practice:

* every metric on the returned result is **recomputed from the DNA** by its
  owning model (invariant #2), never read back from the solver;
* the appended stop codon is optimized through the same trellis, so it is
  re-validated through ``ok_suffix`` (invariant #8);
* the run carries a content-addressed :class:`~bt4.provenance.manifest.Manifest`
  (invariant #9) built from the codon table's own SHA-256.

Only :mod:`bt4.domain`, the pure/biomodel layers, and the solver are imported
here; nothing above ``pipeline`` is referenced.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass

from bt4 import __version__
from bt4.biomodels.codon.tables import CodonUsageTable, load_provenance, load_table
from bt4.biomodels.codon.tai import load_tai_provenance, load_tai_table
from bt4.constraints.kozak import InternalStartConstraint
from bt4.constraints.repeats import InvertedRepeatConstraint, TandemRepeatConstraint
from bt4.constraints.restriction import RestrictionSiteConstraint
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint
from bt4.domain import (
    CODON_TABLE,
    STOP,
    Constraint,
    Frontier,
    Metrics,
    ObjectiveTerm,
    ObjectiveVector,
    OptimalityCertificate,
    Result,
    Scope,
    Severity,
    Violation,
    dominates,
    gc_fraction,
    translate,
    validate_dna,
    validate_protein,
)
from bt4.objectives.dinucleotide import DinucleotideTerm
from bt4.objectives.minmax import MinMaxTerm
from bt4.objectives.ramp import RampTerm
from bt4.objectives.tai import TaiTerm
from bt4.objectives.terms import CaiTerm, GcProximityTerm
from bt4.optimize import SolveResult, solve_exact
from bt4.provenance import Manifest, build_manifest

__all__ = [
    "FrontierResult",
    "OptimizeConfig",
    "ValidationReport",
    "run_frontier",
    "run_optimize",
    "run_validate",
]

_NON_SCORED_AA = frozenset({"M", "W"})

# Cap on the number of scalarization points swept for a multi-objective frontier.
# Each point is a full exact solve, so this bounds the frontier's runtime; the
# grid resolution is reduced to stay under it (see _simplex_grid).
_MAX_FRONTIER_POINTS = 60


@dataclass(frozen=True, slots=True)
class OptimizeConfig:
    """User-facing knobs for a back-translation run.

    Attributes:
        organism: Codon-usage table key or alias (default human).
        gc_target: Desired GC fraction in ``[0, 1]`` for the GC-proximity term.
        cai_weight: Weight on the CAI (log-w) objective in a single solve.
        tai_weight: Weight on the tRNA-adaptation-index (log-w) objective (0
            disables it). Requires a bundled tAI table for ``organism`` (human
            ships one); other organisms raise until their tRNA data is added.
        gc_weight: Weight on the GC-proximity objective in a single solve.
        max_homopolymer: Longest allowed single-base run, or ``None`` to disable.
        forbidden_motifs: Substrings (and, when enabled, their reverse
            complements) that may not appear in the coding sequence.
        avoid_reverse_complement: Also forbid the reverse complement of each
            motif.
        restriction_enzymes: Names of restriction enzymes whose recognition
            sites (and their reverse complements) may not appear.
        ramp_weight: Weight on the 5' translation-ramp term (0 disables it).
        ramp_codons: Length of the 5' ramp window in codons.
        cpg_weight: Weight on the CpG-dinucleotide term (0 disables it).
        cpg_mode: ``"deplete"`` (fewer CpGs, stealth) or ``"elevate"`` (more,
            immunostimulatory) -- only used when ``cpg_weight`` is non-zero.
        minmax_weight: Weight on the %MinMax codon-commonness term (0 disables
            it). Positive values push the sequence toward the ``minmax_direction``
            end of the synonymous-usage range.
        minmax_direction: ``"max"`` (favour common codons) or ``"min"`` (favour
            rare codons) -- only used when ``minmax_weight`` is non-zero.
        tandem_unit: Repeated-unit length whose ``tandem_copies``-fold tandem
            repeat is banned, or ``None`` to disable the tandem-repeat constraint.
        tandem_copies: Number of back-to-back copies that constitutes a banned
            tandem repeat (only used when ``tandem_unit`` is set).
        inverted_stem: Arm length of a banned hairpin (stem-loop) inverted
            repeat, or ``None`` to disable the inverted-repeat constraint.
        inverted_loop: Maximum loop length between the hairpin arms (only used
            when ``inverted_stem`` is set).
        avoid_internal_start: When ``True``, forbid any internal ATG (past the
            first codon) sitting in a strong Kozak context (purine at -3 and G at
            +4) to suppress spurious re-initiation. Off by default.
        gc_min: Optional lower bound on the total GC nucleotide count. Setting a
            GC budget routes the solve through the OR-Tools CP-SAT backend.
        gc_max: Optional upper bound on the total GC nucleotide count.
        refine: When ``True``, run a simulated-annealing refinement pass over the
            exact-DP seed that also maximizes a 5' mRNA-folding score (opening up
            start-proximal structure). The result's certificate becomes
            ``HEURISTIC``. With the default (uncalibrated) folding backend the
            reported folding number is a labeled proxy, not real deltaG -- see
            ``folding_calibrated`` in the result audit. Incompatible with a GC
            budget (the budget is not re-enforced during refinement).
        refine_iterations: Number of SA proposals when ``refine`` is set.
        folding_weight: Weight on the folding score in the refinement objective.
        beam: ``None`` for an exact DP; an int caps the trellis beam width
            (certificate then reports ``beam_truncated``).
        seed: Master seed recorded in the manifest (the solver is deterministic).
    """

    organism: str = "homo_sapiens"
    gc_target: float = 0.55
    cai_weight: float = 1.0
    tai_weight: float = 0.0
    gc_weight: float = 0.0
    max_homopolymer: int | None = 6
    forbidden_motifs: tuple[str, ...] = ()
    avoid_reverse_complement: bool = True
    restriction_enzymes: tuple[str, ...] = ()
    ramp_weight: float = 0.0
    ramp_codons: int = 35
    cpg_weight: float = 0.0
    cpg_mode: str = "deplete"
    minmax_weight: float = 0.0
    minmax_direction: str = "max"
    tandem_unit: int | None = None
    tandem_copies: int = 3
    inverted_stem: int | None = None
    inverted_loop: int = 0
    avoid_internal_start: bool = False
    gc_min: int | None = None
    gc_max: int | None = None
    refine: bool = False
    refine_iterations: int = 2000
    folding_weight: float = 1.0
    beam: int | None = None
    seed: int = 0


@dataclass(frozen=True, slots=True)
class FrontierResult:
    """A Pareto sweep: the frontier plus the full result behind each point."""

    frontier: Frontier
    results: tuple[Result, ...]
    manifest: Manifest

    def delivered(self) -> Result | None:
        """The result at the frontier's chosen index, or ``None``."""
        i = self.frontier.chosen
        if 0 <= i < len(self.results):
            return self.results[i]
        return None


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """A whole-sequence audit of a caller-supplied DNA."""

    dna: str
    violations: tuple[Violation, ...]
    metrics: Metrics

    @property
    def is_feasible(self) -> bool:
        """True when no hard violation is present."""
        return not any(v.severity is Severity.HARD for v in self.violations)


def _active_terms(
    table: CodonUsageTable,
    config: OptimizeConfig,
) -> list[tuple[ObjectiveTerm, float]]:
    """Objective terms active for this config, paired with their weights.

    CAI and GC-proximity are always present; the 5' ramp, CpG, and %MinMax terms
    join only when their weight is non-zero.
    """
    w = table.relative_adaptiveness()
    active: list[tuple[ObjectiveTerm, float]] = [
        (CaiTerm(w), config.cai_weight),
        (GcProximityTerm(config.gc_target), config.gc_weight),
    ]
    if config.tai_weight != 0.0:
        # Raises a clear ValueError if the organism has no bundled tRNA data.
        tai_table = load_tai_table(config.organism)
        active.append((TaiTerm(tai_table.relative_adaptiveness()), config.tai_weight))
    if config.ramp_weight != 0.0:
        active.append((RampTerm(w, config.ramp_codons), config.ramp_weight))
    if config.cpg_weight != 0.0:
        active.append((DinucleotideTerm("CG", config.cpg_mode), config.cpg_weight))
    if config.minmax_weight != 0.0:
        active.append(
            (MinMaxTerm(table.frequency, config.minmax_direction), config.minmax_weight)
        )
    return active


def _frontier_axes(table: CodonUsageTable, config: OptimizeConfig) -> list[ObjectiveTerm]:
    """The objective terms the frontier sweep trades off, in a stable order.

    Exactly the terms :func:`_active_terms` activates -- CAI and GC-proximity
    always, plus any ramp/CpG/%MinMax term whose config weight is non-zero -- but
    stripped of their config weights, because the sweep assigns its own simplex
    weights to each axis.
    """
    return [term for term, _ in _active_terms(table, config)]


def _compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    """Yield every non-negative integer ``parts``-tuple summing to ``total``."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first, *rest)


def _simplex_grid(n_axes: int, steps: int) -> list[tuple[float, ...]]:
    """Weight vectors on the ``n_axes`` unit simplex, at resolution ``steps``.

    Each vector is a normalized integer composition of ``steps - 1`` into
    ``n_axes`` non-negative parts, so the grid includes every axis-corner (all
    weight on one objective) plus the interior mixtures. For two axes this is
    exactly the classic ``[0, 1]`` alpha sweep ``(alpha, 1 - alpha)``, so the
    CAI/GC frontier is unchanged. For three or more axes the raw grid grows
    combinatorially, so the resolution is reduced until the point count is at
    most :data:`_MAX_FRONTIER_POINTS`; each surviving point is still an exact,
    proven-optimal solve -- only the sampling of the trade-off surface is
    coarser, never the optimality of any point.
    """
    if n_axes < 1:
        return []
    if steps <= 1:
        return [tuple(1.0 / n_axes for _ in range(n_axes))]
    res = steps
    while res > 2 and math.comb(res - 1 + n_axes - 1, n_axes - 1) > _MAX_FRONTIER_POINTS:
        res -= 1
    total = res - 1
    return [tuple(c / total for c in comp) for comp in _compositions(total, n_axes)]


def _build_constraints(config: OptimizeConfig) -> list[Constraint]:
    constraints: list[Constraint] = []
    if config.max_homopolymer is not None and config.max_homopolymer > 0:
        constraints.append(HomopolymerConstraint(config.max_homopolymer))
    if config.forbidden_motifs:
        constraints.append(
            ForbiddenMotifConstraint(
                tuple(config.forbidden_motifs),
                reverse_complement=config.avoid_reverse_complement,
            )
        )
    if config.restriction_enzymes:
        constraints.append(
            RestrictionSiteConstraint(enzymes=tuple(config.restriction_enzymes))
        )
    if config.tandem_unit is not None:
        constraints.append(
            TandemRepeatConstraint(config.tandem_unit, config.tandem_copies)
        )
    if config.inverted_stem is not None:
        constraints.append(
            InvertedRepeatConstraint(config.inverted_stem, config.inverted_loop)
        )
    if config.avoid_internal_start:
        constraints.append(InternalStartConstraint())
    return constraints


def _scalar_delta(
    active: Sequence[tuple[ObjectiveTerm, float]],
) -> Callable[[str, str, int], float]:
    pairs = tuple(active)

    def delta(prefix: str, codon: str, pos: int) -> float:
        return sum((w * term.delta(prefix, codon, pos) for term, w in pairs), 0.0)

    return delta


def _violations(dna: str, constraints: Sequence[Constraint]) -> tuple[Violation, ...]:
    out: list[Violation] = []
    for constraint in constraints:
        out.extend(constraint.validate(dna))
    return tuple(out)


def _n_scored_codons(dna: str) -> int:
    return sum(
        1
        for i in range(0, len(dna), 3)
        if (aa := CODON_TABLE.get(dna[i : i + 3])) is not None
        and aa != STOP
        and aa not in _NON_SCORED_AA
    )


def _config_dict(config: OptimizeConfig) -> dict[str, object]:
    return {
        "organism": config.organism,
        "gc_target": config.gc_target,
        "cai_weight": config.cai_weight,
        "tai_weight": config.tai_weight,
        "gc_weight": config.gc_weight,
        "max_homopolymer": config.max_homopolymer,
        "forbidden_motifs": list(config.forbidden_motifs),
        "avoid_reverse_complement": config.avoid_reverse_complement,
        "restriction_enzymes": list(config.restriction_enzymes),
        "ramp_weight": config.ramp_weight,
        "ramp_codons": config.ramp_codons,
        "cpg_weight": config.cpg_weight,
        "cpg_mode": config.cpg_mode,
        "minmax_weight": config.minmax_weight,
        "minmax_direction": config.minmax_direction,
        "tandem_unit": config.tandem_unit,
        "tandem_copies": config.tandem_copies,
        "inverted_stem": config.inverted_stem,
        "inverted_loop": config.inverted_loop,
        "avoid_internal_start": config.avoid_internal_start,
        "gc_min": config.gc_min,
        "gc_max": config.gc_max,
        "refine": config.refine,
        "refine_iterations": config.refine_iterations,
        "folding_weight": config.folding_weight,
        "beam": config.beam,
        "seed": config.seed,
    }


def _manifest(config: OptimizeConfig, extra: dict[str, object]) -> Manifest:
    prov = load_provenance(config.organism)
    cfg = _config_dict(config)
    cfg.update(extra)
    inputs: dict[str, str] = {"codon_table_sha256": prov.sha256}
    if config.tai_weight != 0.0:
        # tAI in play => the tRNA-count table's content hash enters the stamp too.
        inputs["trna_table_sha256"] = load_tai_provenance(config.organism).sha256
    return build_manifest(
        bt4_version=__version__,
        config=cfg,
        inputs=inputs,
        seed=config.seed,
    )


def _metrics(
    dna: str, terms: Sequence[ObjectiveTerm], violations: Sequence[Violation]
) -> Metrics:
    objective = ObjectiveVector({term.name: term.score(dna) for term in terms})
    hard = sum(1 for v in violations if v.severity is Severity.HARD)
    return Metrics(
        objective=objective,
        gc=gc_fraction(dna),
        length_nt=len(dna),
        hard_violations=hard,
        soft_violations=len(violations) - hard,
    )


def _make_result(
    *,
    protein: str,
    dna: str,
    table: CodonUsageTable,
    terms: Sequence[ObjectiveTerm],
    constraints: Sequence[Constraint],
    certificate: OptimalityCertificate,
    config: OptimizeConfig,
    alpha: float | None,
    extra_audit: dict[str, object] | None = None,
) -> Result:
    # Enforce invariant #1 at the boundary: the returned DNA must translate back.
    if translate(dna) != protein + STOP:
        raise AssertionError("round-trip invariant violated: translate(dna) != protein + stop")
    violations = _violations(dna, constraints)
    metrics = _metrics(dna, terms, violations)
    audit: dict[str, object] = {
        "cai": table.cai(dna),
        "gc_percent": metrics.gc * 100.0,
        "n_scored_codons": _n_scored_codons(dna),
        "solver": certificate.solver,
        "seed": config.seed,
        "manifest": _manifest(config, {} if alpha is None else {"alpha": alpha}).to_dict(),
    }
    if config.tai_weight != 0.0:
        audit["tai"] = load_tai_table(config.organism).tai(dna)
    if alpha is not None:
        audit["alpha"] = alpha
    if extra_audit:
        audit.update(extra_audit)
    return Result(
        protein=protein,
        dna=dna,
        metrics=metrics,
        certificate=certificate,
        violations=violations,
        audit=audit,
    )


def _solve_with_gc_budget(
    residues: Sequence[str],
    active: Sequence[tuple[ObjectiveTerm, float]],
    constraints: Sequence[Constraint],
    config: OptimizeConfig,
) -> tuple[SolveResult, OptimalityCertificate]:
    """Solve with a global GC-count budget, choosing an honest backend.

    With neither a local sequence constraint nor a pairwise/positional objective
    term, the problem is a pure additive, context-free objective under a linear
    budget: the CP-SAT backend (see cpsat.py) solves it and proves optimality of
    the integer-scaled objective. When a local constraint or a pairwise term is in
    force -- which CP-SAT cannot encode -- the Lagrangian backend (see
    lagrangian.py) dualizes the GC budget into the exact DP so those constraints
    and terms stay honored, at the cost of a gap-bounded (not proven-optimal)
    certificate. Both backends recompute the delivered metrics from the DNA, so a
    budget is always honestly reported.
    """
    non_local = any(
        term.scope() is Scope.PAIRWISE or term.context_len() > 0 for term, _ in active
    )
    if constraints or non_local:
        # Local constraints and pairwise terms are exactly what CP-SAT drops; the
        # Lagrangian backend keeps them by dualizing the budget into the exact DP.
        from bt4._accel import gc_count  # lazy: matches the accel import site
        from bt4.optimize.lagrangian import solve_lagrangian

        solve = solve_lagrangian(
            residues,
            scalar_delta=_scalar_delta(active),
            constraints=constraints,
            amount=gc_count,
            budget_min=config.gc_min,
            budget_max=config.gc_max,
            beam=config.beam,
            objective_context=max((term.context_len() for term, _ in active), default=0),
            budget_name="gc_budget",
        )
        return solve, solve.certificate

    from bt4.optimize.cpsat import solve_cpsat  # lazy: keeps OR-Tools optional

    pairs = tuple(active)

    def codon_score(codon: str, pos: int) -> float:
        return sum((w * term.delta("", codon, pos) for term, w in pairs), 0.0)

    solve = solve_cpsat(
        residues, codon_score=codon_score, gc_min=config.gc_min, gc_max=config.gc_max
    )
    return solve, solve.certificate


def _refine(
    seed_dna: str,
    residues: Sequence[str],
    active: Sequence[tuple[ObjectiveTerm, float]],
    constraints: Sequence[Constraint],
    config: OptimizeConfig,
) -> tuple[str, OptimalityCertificate, dict[str, object]]:
    """Run a folding-aware SA refinement pass over the exact-DP seed.

    Maximizes the additive objective plus a 5' mRNA-folding score by
    synonymous-only simulated annealing (invariants #1/#5/#7 hold by
    construction). The folding backend is the best available
    (:func:`bt4.biomodels.folding.default`); when that is the uncalibrated
    baseline the returned ``folding_dg`` is a labeled proxy, not real deltaG, and
    is flagged ``folding_calibrated=False`` (CLAUDE.md sections 6 and 10.6). The
    returned certificate is ``HEURISTIC`` -- refinement never claims optimality.
    """
    from bt4.biomodels.folding import default as folding_default  # lazy: keep core light
    from bt4.optimize.anneal_refine import anneal_refine  # lazy

    model = folding_default()
    pairs = tuple(active)
    weight = config.folding_weight

    def score(dna: str) -> float:
        additive = sum(w * term.score(dna) for term, w in pairs)
        return additive + weight * model.score_sequence(dna)

    result = anneal_refine(
        seed_dna,
        residues,
        score,
        constraints,
        iterations=config.refine_iterations,
        seed=config.seed,
    )
    extra_audit: dict[str, object] = {
        "refined": True,
        "folding_model": model.name,
        "folding_calibrated": model.calibrated,
        "folding_dg": model.five_prime_dg(result.dna),
    }
    return result.dna, result.certificate, extra_audit


def run_optimize(protein: str, config: OptimizeConfig | None = None) -> Result:
    """Optimize ``protein`` into a coding sequence under ``config`` (single solve).

    A GC budget (``gc_min``/``gc_max``) routes the solve through a budget backend:
    the OR-Tools CP-SAT backend for a pure additive, context-free objective with
    no local constraints (proven optimal), or the Lagrangian backend when a local
    constraint or pairwise term is in force (gap-bounded, but those are honored).
    With no GC budget the exact codon-trellis DP is used.

    Raises:
        ValueError: On an invalid protein.
        bt4.optimize.InfeasibleError: If no feasible sequence exists (including a
            GC budget that admits no assignment).
    """
    config = config or OptimizeConfig()
    has_budget = config.gc_min is not None or config.gc_max is not None
    if config.refine and has_budget:
        raise ValueError(
            "refine is not supported together with a GC budget (gc_min/gc_max): "
            "the budget is a global constraint that refinement would not re-enforce"
        )
    p = validate_protein(protein)
    table = load_table(config.organism)
    active = _active_terms(table, config)
    terms = [term for term, _ in active]
    constraints = _build_constraints(config)
    residues = [*p, STOP]
    if has_budget:
        solve, certificate = _solve_with_gc_budget(residues, active, constraints, config)
    else:
        solve = solve_exact(
            residues,
            scalar_delta=_scalar_delta(active),
            constraints=constraints,
            beam=config.beam,
            objective_context=max((term.context_len() for term, _ in active), default=0),
        )
        certificate = solve.certificate
    dna = solve.dna
    extra_audit: dict[str, object] | None = None
    if config.refine:
        dna, certificate, extra_audit = _refine(dna, residues, active, constraints, config)
    return _make_result(
        protein=p,
        dna=dna,
        table=table,
        terms=terms,
        constraints=constraints,
        certificate=certificate,
        config=config,
        alpha=None,
        extra_audit=extra_audit,
    )


def _delivered_index(points: Sequence[ObjectiveVector]) -> int:
    """Index of the delivered frontier point (highest CAI), or ``-1`` if empty."""
    if not points:
        return -1
    return max(range(len(points)), key=lambda i: points[i].get("cai_logw"))


def _nondominated_indices(vectors: Sequence[ObjectiveVector]) -> list[int]:
    keep: list[int] = []
    for i, p in enumerate(vectors):
        if not any(dominates(vectors[j], p) for j in range(len(vectors)) if j != i):
            keep.append(i)
    return keep


def run_frontier(
    protein: str,
    config: OptimizeConfig | None = None,
    steps: int = 11,
    *,
    on_progress: Callable[[int, int], None] | None = None,
    should_cancel: Callable[[], bool] | None = None,
) -> FrontierResult:
    """Sweep the objective trade-offs and return the Pareto frontier of exact solves.

    The frontier trades off every *active* objective term: CAI and GC-proximity
    always, plus any 5' ramp, CpG, or %MinMax term whose config weight is
    non-zero. Scalarization weights are drawn from a uniform grid on the unit
    simplex over those axes (see :func:`_simplex_grid`); with only CAI and GC
    active this is exactly the classic ``(alpha, 1 - alpha)`` sweep, so the
    two-objective frontier is unchanged, while three or more active objectives
    now trace a genuine multi-dimensional trade-off surface rather than holding
    the extra terms at a fixed weight.

    Each grid point is an exact (or beam-capped) solve; duplicate sequences are
    collapsed and the non-dominated subset (over the *full* objective vector) is
    returned, with the delivered point set to the highest-CAI frontier member.
    The frontier uses the exact DP; a GC budget (``gc_min``/``gc_max``) applies
    to :func:`run_optimize` only, not here.

    Raises:
        ValueError: On an invalid protein or ``steps < 1``.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    config = config or OptimizeConfig()
    if steps < 1:
        raise ValueError(f"steps must be >= 1, got {steps}")
    p = validate_protein(protein)
    table = load_table(config.organism)
    constraints = _build_constraints(config)
    residues = [*p, STOP]

    axes = _frontier_axes(table, config)
    objective_context = max((term.context_len() for term in axes), default=0)
    grid = _simplex_grid(len(axes), steps)
    total = len(grid)
    by_dna: dict[str, Result] = {}
    for i, weights in enumerate(grid):
        if should_cancel is not None and should_cancel():
            # Cooperative cancellation: stop the sweep and return whatever points
            # completed. Each computed point is a real exact solve, so a partial
            # frontier is honest -- just a coarser sample of the trade-off surface.
            break
        active = list(zip(axes, weights, strict=True))
        solve = solve_exact(
            residues,
            scalar_delta=_scalar_delta(active),
            constraints=constraints,
            beam=config.beam,
            objective_context=objective_context,
        )
        if solve.dna not in by_dna:
            by_dna[solve.dna] = _make_result(
                protein=p,
                dna=solve.dna,
                table=table,
                terms=axes,
                constraints=constraints,
                certificate=solve.certificate,
                config=config,
                alpha=weights[0],
            )
        if on_progress is not None:
            on_progress(i + 1, total)

    if not by_dna:
        # Only reachable if cancelled before the first grid point completed.
        raise ValueError("optimization was cancelled before any result was produced")

    uniq = list(by_dna.values())
    vectors = [r.metrics.objective for r in uniq]
    keep = _nondominated_indices(vectors)
    kept = sorted((uniq[i] for i in keep), key=lambda r: r.metrics.objective.get("cai_logw"))
    points = tuple(r.metrics.objective for r in kept)
    frontier = Frontier(points=points, chosen=_delivered_index(points))
    return FrontierResult(
        frontier=frontier,
        results=tuple(kept),
        manifest=_manifest(
            config, {"frontier_steps": steps, "frontier_objectives": len(axes)}
        ),
    )


def run_validate(dna: str, config: OptimizeConfig | None = None) -> ValidationReport:
    """Audit a caller-supplied ``dna`` under ``config`` (no optimization).

    Raises:
        ValueError: On non-ACGT input.
    """
    config = config or OptimizeConfig()
    d = validate_dna(dna)
    table = load_table(config.organism)
    terms = [term for term, _ in _active_terms(table, config)]
    constraints = _build_constraints(config)
    violations = _violations(d, constraints)
    if len(d) % 3 == 0:
        metrics = _metrics(d, terms, violations)
    else:
        hard = sum(1 for v in violations if v.severity is Severity.HARD)
        metrics = Metrics(
            objective=ObjectiveVector({}),
            gc=gc_fraction(d),
            length_nt=len(d),
            hard_violations=hard,
            soft_violations=len(violations) - hard,
        )
    return ValidationReport(dna=d, violations=violations, metrics=metrics)
