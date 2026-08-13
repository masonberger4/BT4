"""Candidate-set assembly + expression rerank (design-flow step 3).

This composes the already-shipped pieces into the finalist set the expression
head ranks (``docs/DESIGN_expression_splice_flow.md`` step 3): the **Pareto
frontier** (objective trade-off diversity) plus, when a GLOBAL rule is active and
the exact-DP seed actually violates it, a small deterministic library of
**repeat-refined variants** (structural diversity). The combined set is
deduplicated, scored by an :class:`~bt4.biomodels.expression.ExpressionPredictor`
(batched in one call when the backend supports it), and delivered under the same
honesty rule as :func:`~bt4.pipeline.rerank.rerank_by_expression`.

The load-bearing honesty rules (CLAUDE.md §10.5/§10.6), hardened after a design
review:

* **An uncalibrated score never steers delivery, and an uncalibrated ordering is
  never presented as a ranking.** With the default placeholder (or any
  ``calibrated is False`` head, including the shipped RiboNN adapter) the set is
  returned in **discovery order** (the solver-delivered sequence first, then the
  rest of the frontier, then the repeat-refined variants), ``order_basis`` is
  ``"discovery"``, and ``chosen`` is the solver-delivered sequence -- expression
  scores are pure annotation. Only when ``predictor.calibrated`` is ``True`` is the
  set reordered by predicted expression (``order_basis == "expression_rank"``) and
  ``chosen`` moved to the top-scoring member.
* **The delivered (``chosen``) sequence is invariant to ``n``.** Uncalibrated, the
  sequence that ships is the solver-delivered one; it is pinned **first** in
  discovery order so the cap (keep first-``n``) can never drop it. Calibrated, the
  sequence that ships is the head's top pick; the cap keeps the top-``n`` by score,
  so that pick (the top of the keep) is always retained. Either way the cap is
  applied **after** scoring, so a calibrated reranker never has its best candidate
  truncated away before it is ranked.
* **No silent truncation, no unenforced claim.** How many members were dropped by
  de-duplication (``n_dedup_dropped``) and by the size cap (``n_dropped_cap``) are
  both reported. Variants are labelled ``"repeat_refined"`` (the *process*, not a
  guaranteed outcome): each carries its own :class:`~bt4.domain.Result` whose
  ``violations`` disclose any residual GLOBAL rule it did *not* clear, so a
  still-violating variant is never presented as clean.
* **Determinism (#7).** Frontier and variants are seeded; the calibrated reorder
  uses the total order ``(score desc, discovery index asc)``.

Scoping note: ranking by the protocol-generic ``ExpressionResult.score`` is valid
here precisely because **every candidate encodes the same protein** with the same
UTR context, so the absolute score and the CDS-attributable Δ are rank-equivalent.
Cross-protein ranking (which would need the UTR-fixed ``delta_logte`` framing) is
out of scope for this assembler.

This layer composes :mod:`bt4.domain`, :mod:`bt4.biomodels`, and its sibling
:mod:`bt4.pipeline.optimize`; it imports nothing above ``pipeline``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from bt4.biomodels.codon.tables import CodonUsageTable, load_table
from bt4.biomodels.expression import (
    BatchExpressionPredictor,
    ExpressionPredictor,
    ExpressionResult,
)
from bt4.biomodels.expression import default as expression_default
from bt4.domain import STOP, Constraint, ObjectiveTerm, Result, Severity, validate_protein
from bt4.pipeline.optimize import (
    OptimizeConfig,
    _active_terms,
    _build_constraints,
    _build_global_constraints,
    _make_result,
    _manifest,
    _refine,
    run_frontier,
)
from bt4.provenance import Manifest

__all__ = ["Candidate", "CandidateSet", "assemble_and_rank_candidates"]

# Default cap on the candidate set (each member is a heavy-CNN input downstream).
_DEFAULT_MAX_CANDIDATES = 24
# Default number of repeat-refined variants to draw when the seed violates a GLOBAL rule.
_DEFAULT_REPEAT_VARIANTS = 4


@dataclass(frozen=True, slots=True)
class Candidate:
    """One member of a candidate set: a full result plus its expression annotation.

    Attributes:
        result: The candidate sequence as a full :class:`~bt4.domain.Result` --
            it round-trips (invariant #1), carries metrics recomputed from its own
            DNA (invariant #2), a certificate, and ``violations`` that disclose any
            residual (including GLOBAL) rule it does not satisfy.
        source: ``"frontier"`` (a Pareto-frontier point) or ``"repeat_refined"``
            (a synonymous SA-refined variant of the delivered seed -- the label
            names the *process*, not a guarantee the repeat was fully removed;
            check ``result.violations``).
        expression_score: The predictor's score for this sequence (larger is
            better). **Only** a validated prediction when ``expression_calibrated``
            is ``True``; otherwise an annotation that must not be read as one.
        expression_model: The producing backend's name.
        expression_calibrated: Mirror of the predictor's ``calibrated`` flag.
        expression_units: The score's units label.
    """

    result: Result
    source: str
    expression_score: float
    expression_model: str
    expression_calibrated: bool
    expression_units: str


@dataclass(frozen=True, slots=True)
class CandidateSet:
    """An assembled, expression-annotated candidate set for one protein.

    Attributes:
        candidates: The members. Ordered by predicted expression (best first) when
            ``order_basis == "expression_rank"``, else in discovery order (the
            delivered sequence first). Never read discovery order as a ranking.
        chosen: Index of the delivered candidate -- the top-scoring member when the
            predictor is calibrated, else the solver-delivered sequence.
        calibrated: Whether the predictor was calibrated (and thus steered
            delivery/ordering). ``False`` for the default placeholder and the
            shipped RiboNN adapter.
        order_basis: ``"expression_rank"`` (calibrated) or ``"discovery"``.
        manifest: Content-addressed provenance stamp; folds in the predictor
            identity and the assembly knobs, so two sets differing only in
            predictor or ``n`` stamp differently (invariant #9).
        n_frontier: Frontier points assembled (before de-dup).
        n_repeat_refined: Repeat-refined variants assembled (before de-dup).
        n_dedup_dropped: Members dropped as duplicate sequences.
        n_dropped_cap: Members dropped by the ``n`` size cap (0 if none).
        scored_batched: Whether scoring used the backend's one-call batch path.
        repeat_note: Honest one-line reason for the repeat-refined count (e.g. the
            seed already satisfied the GLOBAL rules, or none were active).
    """

    candidates: tuple[Candidate, ...]
    chosen: int
    calibrated: bool
    order_basis: str
    manifest: Manifest
    n_frontier: int
    n_repeat_refined: int
    n_dedup_dropped: int
    n_dropped_cap: int
    scored_batched: bool
    repeat_note: str

    def delivered(self) -> Candidate | None:
        """The candidate at :attr:`chosen`, or ``None`` if the set is empty."""
        if 0 <= self.chosen < len(self.candidates):
            return self.candidates[self.chosen]
        return None


def _score_candidates(
    backend: ExpressionPredictor, dnas: list[str]
) -> tuple[list[ExpressionResult], bool]:
    """Score ``dnas`` (in order); use the backend's batch path when it has one.

    Returns ``(results, batched)``. A :class:`BatchExpressionPredictor` scores the
    whole set in one invocation (the reason step 1 exists); any other predictor is
    scored per sequence via ``score_sequence``. The two paths must agree per the
    batch contract, so ordering and values are identical either way.
    """
    if not dnas:
        return [], False
    if isinstance(backend, BatchExpressionPredictor):
        results = list(backend.score_many(dnas))
        if len(results) != len(dnas):
            raise RuntimeError(
                f"score_many returned {len(results)} results for {len(dnas)} inputs"
            )
        return results, True
    return [backend.score_sequence(dna) for dna in dnas], False


def _repeat_refined_variants(
    protein: str,
    seed_dna: str,
    config: OptimizeConfig,
    table: CodonUsageTable,
    active: list[tuple[ObjectiveTerm, float]],
    constraints: list[Constraint],
    global_constraints: list[Constraint],
    k: int,
) -> list[Result]:
    """Return up to ``k`` deterministic SA-refined variants of ``seed_dna``.

    Each variant reuses :func:`~bt4.pipeline.optimize._refine` (no duplicated
    score/penalty/folding construction) with a distinct derived seed
    (``config.seed + 1 .. + k``), so the trajectories differ yet every run is
    reproducible (#7). Variants are *not* guaranteed repeat-free: single-codon SA
    can hit a feasibility floor, so each variant's residual GLOBAL violations are
    reported through its own :class:`~bt4.domain.Result`.
    """
    residues = [*protein, STOP]
    variants: list[Result] = []
    for j in range(1, k + 1):
        dna, certificate, extra_audit, manifest_extra = _refine(
            seed_dna,
            residues,
            active,
            constraints,
            global_constraints,
            config,
            with_folding=config.refine,
            seed=config.seed + j,
        )
        variants.append(
            _make_result(
                protein=protein,
                dna=dna,
                table=table,
                terms=[term for term, _ in active],
                constraints=[*constraints, *global_constraints],
                certificate=certificate,
                config=config,
                alpha=None,
                extra_audit=extra_audit,
                manifest_extra=manifest_extra,
            )
        )
    return variants


def assemble_and_rank_candidates(
    protein: str,
    config: OptimizeConfig | None = None,
    *,
    steps: int = 11,
    n: int = _DEFAULT_MAX_CANDIDATES,
    repeat_variants: int = _DEFAULT_REPEAT_VARIANTS,
    predictor: ExpressionPredictor | None = None,
    on_progress: Callable[[int, int], None] | None = None,
    should_cancel: Callable[[], bool] | None = None,
) -> CandidateSet:
    """Assemble the frontier (+ repeat-refined variants) and rank it by expression.

    Runs the Pareto frontier, optionally augments it with deterministic
    repeat-refined variants of the delivered seed (only when a GLOBAL rule is
    active *and* the seed actually violates it -- mirroring
    :func:`~bt4.pipeline.optimize.run_optimize`), de-duplicates by sequence,
    scores every member with ``predictor`` (batched in one call when supported),
    and delivers under the calibrated-gating rule (see the module docstring): an
    uncalibrated head annotates only, a calibrated head reorders and re-picks.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`.
        steps: Frontier scalarization grid resolution (see
            :func:`~bt4.pipeline.optimize.run_frontier`).
        n: Maximum candidates to keep (``>= 1``). The delivered sequence is always
            retained; the cap is applied after scoring.
        repeat_variants: Number of repeat-refined variants to attempt when the
            delivered seed violates a GLOBAL rule (``>= 0``).
        predictor: Expression backend; defaults to
            :func:`bt4.biomodels.expression.default` (the neutral placeholder, so
            ranking is a pure reporting no-op).
        on_progress: Forwarded to the frontier sweep.
        should_cancel: Forwarded to the frontier sweep.

    Returns:
        A :class:`CandidateSet`.

    Raises:
        ValueError: On an invalid protein, ``n < 1``, ``repeat_variants < 0``, or
            ``steps < 1``.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if repeat_variants < 0:
        raise ValueError(f"repeat_variants must be >= 0, got {repeat_variants}")
    config = config or OptimizeConfig()
    backend = predictor or expression_default()
    p = validate_protein(protein)
    table = load_table(config.organism, reference_set=config.reference_set)

    fr = run_frontier(
        protein, config, steps, on_progress=on_progress, should_cancel=should_cancel
    )
    frontier_results = list(fr.results)
    delivered = fr.delivered()
    delivered_dna = delivered.dna if delivered is not None else None

    # Repeat-refined variants: only when a GLOBAL rule is active AND the delivered
    # exact-DP seed actually violates it (a clean seed needs no repair) -- exactly
    # run_optimize's gate. Otherwise record honestly why K was not used.
    global_constraints = _build_global_constraints(config)
    repeat_results: list[Result] = []
    if not global_constraints:
        repeat_note = "no GLOBAL rule active (max_repeat_length / avoid_uorf)"
    elif repeat_variants == 0:
        repeat_note = "repeat_variants=0"
    elif delivered_dna is None:
        repeat_note = "no delivered frontier seed"
    else:
        seed_hard = sum(
            1
            for c in global_constraints
            for v in c.validate(delivered_dna)
            if v.severity is Severity.HARD
        )
        if seed_hard == 0:
            repeat_note = "delivered seed already satisfies the GLOBAL rules"
        else:
            active = _active_terms(table, config)
            constraints = _build_constraints(config)
            repeat_results = _repeat_refined_variants(
                p, delivered_dna, config, table, active, constraints,
                global_constraints, repeat_variants,
            )
            repeat_note = f"{len(repeat_results)} variant(s) refined from the delivered seed"

    n_frontier = len(frontier_results)
    n_repeat_refined = len(repeat_results)

    # Discovery order: delivered sequence FIRST (pinned), then the rest of the
    # frontier, then the repeat-refined variants; de-duplicated by DNA.
    discovery: list[tuple[Result, str]] = []
    seen: set[str] = set()

    def _add(result: Result, source: str) -> None:
        if result.dna not in seen:
            seen.add(result.dna)
            discovery.append((result, source))

    if delivered is not None:
        _add(delivered, "frontier")
    for r in frontier_results:
        _add(r, "frontier")
    for r in repeat_results:
        _add(r, "repeat_refined")

    n_dedup_dropped = (n_frontier + n_repeat_refined) - len(discovery)

    # Score the FULL de-duplicated set before any cap (so a calibrated reranker
    # never loses its best candidate to truncation).
    ers, batched = _score_candidates(backend, [r.dna for r, _ in discovery])
    scored = [
        Candidate(
            result=r,
            source=source,
            expression_score=er.score,
            expression_model=er.model_name,
            expression_calibrated=er.calibrated,
            expression_units=er.units,
        )
        for (r, source), er in zip(discovery, ers, strict=True)
    ]

    calibrated = backend.calibrated and bool(scored)
    if calibrated:
        # A calibrated head delivers its top-scoring pick, so THAT is the sequence
        # that must survive the cap -- and it trivially does, being the top of a
        # top-n-by-score keep. The solver-delivered original is then just one
        # candidate, retained iff it ranks within n (no special pin needed here).
        order_basis = "expression_rank"
        # Total order: score desc, then discovery index asc (deterministic ties, #7).
        order = sorted(range(len(scored)), key=lambda i: (-scored[i].expression_score, i))
        candidates = tuple(scored[i] for i in order[:n])
        chosen = 0  # the top predicted expression
    else:
        # An uncalibrated head only annotates: the sequence that ships is the
        # solver-delivered one, which is pinned FIRST in discovery order, so the cap
        # (keep first-n) can never drop it -- it is invariant to n.
        order_basis = "discovery"
        candidates = tuple(scored[:n])
        chosen = 0  # the solver-delivered sequence

    n_dropped_cap = len(scored) - len(candidates)

    manifest = _manifest(
        config,
        {
            "mode": "candidates",
            "predictor": backend.name,
            "predictor_calibrated": backend.calibrated,
            "candidate_n": n,
            "repeat_variants": repeat_variants,
            "frontier_steps": steps,
        },
    )
    return CandidateSet(
        candidates=candidates,
        chosen=chosen,
        calibrated=calibrated,
        order_basis=order_basis,
        manifest=manifest,
        n_frontier=n_frontier,
        n_repeat_refined=n_repeat_refined,
        n_dedup_dropped=n_dedup_dropped,
        n_dropped_cap=n_dropped_cap,
        scored_batched=batched,
        repeat_note=repeat_note,
    )
