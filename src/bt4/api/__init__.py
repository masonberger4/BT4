"""BT4 api: the stable, print-free entry point every frontend calls.

The CLI, the BT4 Studio desktop app, and the optional HTTP service all go
through this module and nothing below it. Results carry their own recomputed
metrics, optimality certificate, and provenance manifest, so a caller never has
to reach into the optimizer or the biomodels to understand what it got back.

This layer never prints; it raises on error and returns immutable results.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.tables import (
    GENOME_WIDE,
    HIGHLY_EXPRESSED,
    REFERENCE_SETS,
    TableProvenance,
    available_organisms,
    available_reference_sets,
    default_reference_set,
    load_provenance,
)
from bt4.biomodels.codon.tai import available_tai_organisms
from bt4.biomodels.expression import (
    ExpressionPanel,
    ExpressionPredictor,
    ExpressionResult,
    PanelRow,
    read_panel,
)
from bt4.biomodels.expression import available_backends as available_expression_backends
from bt4.biomodels.expression import default as expression_model
from bt4.biomodels.expression import resolve_backend as resolve_expression_backend
from bt4.biomodels.splice import SpliceAuditReport, SpliceFlag, SplicePredictor
from bt4.constraints import (
    ForbiddenPreset,
    available_enzymes,
    available_forbidden_presets,
    enzyme_provenance,
    enzyme_suggestions,
    resolve_enzyme,
    unknown_enzyme_message,
)
from bt4.domain import (
    AMINO_ACIDS,
    ConstructContext,
    Result,
    Severity,
    Violation,
    validate_protein,
)
from bt4.io import (
    GenBankFeature,
    GenBankRecord,
    context_from_genbank,
    parse_fasta,
    parse_genbank,
    read_fasta,
    read_genbank,
    result_to_dict,
    result_to_json,
    to_fasta,
    write_genbank,
)
from bt4.pipeline import (
    APPLICATION_PRESETS,
    PACKAGING_LIMITS,
    ApplicationPreset,
    Candidate,
    CandidateSet,
    ConstructAudit,
    CrossCheckSite,
    EnzymeOccurrence,
    FrontierResult,
    GateComparison,
    GateSettings,
    InfeasibleError,
    LibraryResult,
    OptimizeConfig,
    PackagingReport,
    SpliceCrossCheck,
    Track,
    TracksResult,
    ValidationReport,
    apply_preset,
    assemble_and_rank_candidates,
    audit_candidate_set,
    audit_construct,
    available_presets,
    available_splice_backends,
    packaging_report,
    rerank_by_expression,
    resolve_preset,
    resolve_splice_backend,
    run_frontier,
    run_library,
    run_optimize,
    run_splice_crosscheck,
    run_tracks,
    run_validate,
    summarize,
)
from bt4.pipeline import (
    run_panel_gate as expression_gate,
)

__all__ = [
    "AMINO_ACIDS",
    "APPLICATION_PRESETS",
    "GENOME_WIDE",
    "HIGHLY_EXPRESSED",
    "PACKAGING_LIMITS",
    "REFERENCE_SETS",
    "ApplicationPreset",
    "Candidate",
    "CandidateSet",
    "ConstructAudit",
    "ConstructContext",
    "CrossCheckSite",
    "EnzymeOccurrence",
    "ExpressionPanel",
    "ExpressionPredictor",
    "ExpressionResult",
    "ForbiddenPreset",
    "FrontierResult",
    "GateComparison",
    "GateSettings",
    "GenBankFeature",
    "GenBankRecord",
    "InfeasibleError",
    "LibraryResult",
    "OptimizeConfig",
    "PackagingReport",
    "PanelRow",
    "Result",
    "Severity",
    "SpliceAuditReport",
    "SpliceCrossCheck",
    "SpliceFlag",
    "TableProvenance",
    "Track",
    "TracksResult",
    "ValidationReport",
    "Violation",
    "apply_preset",
    "assemble_and_rank_candidates",
    "audit_candidate_set",
    "audit_construct",
    "available_enzymes",
    "available_expression_backends",
    "available_forbidden_presets",
    "available_organisms",
    "available_presets",
    "available_reference_sets",
    "available_splice_backends",
    "available_tai_organisms",
    "build_table",
    "candidates",
    "context_from_genbank",
    "count_codons",
    "default_reference_set",
    "enzyme_provenance",
    "enzyme_suggestions",
    "expression_gate",
    "expression_model",
    "frontier",
    "library",
    "load_provenance",
    "optimize",
    "packaging_report",
    "parse_fasta",
    "parse_genbank",
    "read_fasta",
    "read_genbank",
    "read_panel",
    "rerank_by_expression",
    "resolve_enzyme",
    "resolve_expression_backend",
    "resolve_preset",
    "resolve_splice_backend",
    "result_to_dict",
    "result_to_json",
    "splice_audit",
    "splice_crosscheck",
    "summarize",
    "to_fasta",
    "tracks",
    "unknown_enzyme_message",
    "validate",
    "validate_protein",
    "write_genbank",
    "write_table",
]


def optimize(protein: str, config: OptimizeConfig | None = None) -> Result:
    """Back-translate ``protein`` into an optimized coding sequence (single solve).

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`.

    Returns:
        A :class:`~bt4.domain.result.Result` with recomputed metrics, an
        optimality certificate, any violations, and a provenance manifest.

    Raises:
        ValueError: On an invalid protein or unknown organism.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    return run_optimize(protein, config)


def frontier(
    protein: str,
    config: OptimizeConfig | None = None,
    steps: int = 11,
    *,
    on_progress: Callable[[int, int], None] | None = None,
    should_cancel: Callable[[], bool] | None = None,
) -> FrontierResult:
    """Compute the multi-objective Pareto frontier for ``protein``.

    Trades off every active objective axis (CAI and GC always, plus any
    ramp/CpG/%MinMax/tAI term whose config weight is non-zero); with only CAI and
    GC active this is the classic CAI-vs-GC frontier.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`.
        steps: Resolution of the scalarization weight grid on the objective
            simplex (for two objectives, the ``[0, 1]`` alpha sweep).
        on_progress: Optional callback invoked as ``on_progress(done, total)``
            after each frontier grid point completes. Lets a UI show progress
            without blocking the engine.
        should_cancel: Optional predicate polled before each grid point; when it
            returns ``True`` the sweep stops early and returns the frontier of the
            points computed so far (raising only if none completed).

    Returns:
        A :class:`FrontierResult`: the non-dominated frontier plus the full
        result behind each point and a top-level manifest.

    Raises:
        ValueError: On an invalid protein, unknown organism, ``steps < 1``, or a
            cancellation before any point completed.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    return run_frontier(
        protein, config, steps, on_progress=on_progress, should_cancel=should_cancel
    )


def validate(dna: str, config: OptimizeConfig | None = None) -> ValidationReport:
    """Audit a caller-supplied coding sequence (no optimization).

    Args:
        dna: An ACGT coding sequence.
        config: Run configuration whose constraints define the audit.

    Returns:
        A :class:`ValidationReport` with the whole-sequence violations and
        recomputed metrics.

    Raises:
        ValueError: On non-ACGT input or unknown organism.
    """
    return run_validate(dna, config)


def tracks(
    dna: str,
    organism: str = "homo_sapiens",
    *,
    reference_set: str | None = None,
    nt_window: int = 50,
    codon_window: int = 18,
    splice: bool = False,
    splice_predictor: object | None = None,
    context: ConstructContext | None = None,
) -> TracksResult:
    """Compute per-site tracks for a coding sequence.

    Sliding-window **reporting** profiles (GC fraction, CpG density, and -- when
    codon-aligned -- %MinMax) so a delivered sequence's composition can be
    audited or plotted position-by-position, plus an opt-in per-nucleotide
    cryptic-splice-site track. Nothing here feeds the optimizer.

    Args:
        dna: An ACGT coding sequence (case-insensitive).
        organism: Codon-usage table key/alias for the %MinMax reference.
        reference_set: That organism's reference set for the %MinMax
            frequencies (``None`` = its default).
        nt_window: Window (nucleotides) for the GC and CpG tracks.
        codon_window: Window (codons) for the %MinMax track.
        splice: Add the per-nucleotide ``splice_site`` track. Opt-in because it
            runs a model, and its numbers inherit that model's calibration
            status -- with the default PWM baseline they are UNCALIBRATED
            pseudo-scores showing where consensus-like positions sit, not
            probabilities.
        splice_predictor: Backend to score with (``None`` = the baseline).
        context: Known flanking sequence. When given, the splice track is scored
            with the CDS in its real flanks rather than alone; values stay aligned
            to the CDS.

    Returns:
        A :class:`~bt4.pipeline.tracks.TracksResult` bundling the named tracks,
        recording which splice model ran and whether it is calibrated.

    Raises:
        ValueError: On non-ACGT input, a non-positive window, or unknown organism.
    """
    return run_tracks(
        dna,
        organism,
        reference_set=reference_set,
        nt_window=nt_window,
        codon_window=codon_window,
        splice=splice,
        splice_predictor=splice_predictor,
        context=context,
    )


def library(
    protein: str,
    config: OptimizeConfig | None = None,
    n: int = 8,
    *,
    seed: int | None = None,
    temperature: float = 1.0,
) -> LibraryResult:
    """Sample a library of coding sequences for ``protein`` (stochastic, not optimal).

    Library / degenerate-design mode (CLAUDE.md §9, Phase 5). Instead of a single
    most-favored-codon optimum, this draws ``n`` sequences by **sampling** each
    residue's synonymous-codon distribution (organism usage frequencies, raised to
    ``1 / temperature``), keeping only codons that satisfy every LOCAL constraint.
    The delivered sequences are **sampled, not optimized**: each carries the
    ``SAMPLED`` certificate and makes no optimality or expression claim. GLOBAL
    constraints (``max_repeat_length``, ``avoid_uorf``) are not enforced during
    sampling but are validated and reported honestly on every member.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`. Objective
            weights do not steer the draw (this is a sampler, not a solver); the
            codon table and the LOCAL constraints shape it.
        n: Number of sequences to sample (``>= 1``).
        seed: Master sampling seed; when ``None`` the run uses ``config.seed``.
            The effective seed enters the manifest, so the library reproduces from
            its stamp.
        temperature: Sampling temperature (``> 0``). ``-> 0`` approaches the
            per-residue argmax, ``1.0`` is the natural distribution, large values
            approach uniform.

    Returns:
        A :class:`LibraryResult`: the sampled members (each a full
        :class:`~bt4.domain.Result` with recomputed metrics, a ``SAMPLED``
        certificate, and any residual violations), a shared provenance manifest,
        and honest diversity statistics.

    Raises:
        ValueError: On an invalid protein, unknown organism, ``n < 1``, or
            ``temperature <= 0``.
        bt4.optimize.InfeasibleError: If the LOCAL constraints admit no feasible
            sequence.
    """
    return run_library(protein, config, n, seed=seed, temperature=temperature)


def candidates(
    protein: str,
    config: OptimizeConfig | None = None,
    *,
    steps: int = 11,
    n: int = 24,
    repeat_variants: int = 4,
    predictor: ExpressionPredictor | None = None,
) -> CandidateSet:
    """Assemble the frontier (+ repeat-refined variants) and rank it by expression.

    Design-flow step 3 (``docs/DESIGN_expression_splice_flow.md``): builds the
    finalist set an expression head ranks -- the Pareto frontier plus, when a
    GLOBAL rule is active and the delivered exact-DP seed violates it, a small
    deterministic library of repeat-refined variants -- de-duplicates it, scores
    every member with ``predictor`` (in one batched call when the backend supports
    it), and delivers under the **calibrated-gating** rule: an uncalibrated head
    (the default placeholder, and the shipped RiboNN adapter) only *annotates* --
    the set stays in discovery order and the solver-delivered sequence is
    ``chosen`` -- while a calibrated head reorders by predicted expression and
    re-picks the top. The delivered (``chosen``) sequence is invariant to ``n``
    (the cap is applied after scoring and never drops it) (CLAUDE.md §10.5/§10.6).

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to :class:`OptimizeConfig`.
        steps: Frontier scalarization grid resolution.
        n: Maximum candidates to keep (``>= 1``); the delivered sequence is always
            retained and the cap is applied after scoring.
        repeat_variants: Repeat-refined variants to attempt when the delivered seed
            violates a GLOBAL rule (``>= 0``).
        predictor: Expression backend; defaults to the neutral placeholder, so
            ranking is a pure reporting no-op.

    Returns:
        A :class:`CandidateSet`: the (possibly reranked) candidates, the delivered
        index, the calibration/order-basis flags, honest de-dup/cap counts, and a
        provenance manifest that folds in the predictor identity.

    Raises:
        ValueError: On an invalid protein, ``n < 1``, ``repeat_variants < 0``, or
            ``steps < 1``.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    return assemble_and_rank_candidates(
        protein,
        config,
        steps=steps,
        n=n,
        repeat_variants=repeat_variants,
        predictor=predictor,
    )


def splice_audit(
    candidate_set: CandidateSet,
    *,
    reference: str | None = None,
    predictors: Sequence[SplicePredictor] | None = None,
    threshold: float = 0.5,
    match_window: int = 3,
    context: ConstructContext | None = None,
) -> SpliceAuditReport:
    """Localize-and-flag cryptic splice sites across a candidate set (no editing).

    Design-flow step 4 (``docs/DESIGN_expression_splice_flow.md`` Stage C): runs the
    available splice backends over the step-3 candidate set to **localize** residual
    cryptic sites and attach whole-panel **backend agreement** -- an advisory
    annotation pass that never edits the sequences. Every shipped backend is
    ``calibrated is False`` today, so ``report.all_calibrated`` is ``False`` and
    every flag is advisory (a targeted auto-edit is a future, calibrated-gated step;
    CLAUDE.md §6/§10.6).

    Args:
        candidate_set: A step-3 :class:`~bt4.pipeline.candidates.CandidateSet` (from
            :func:`candidates`).
        reference: Sequence each candidate's added risk is measured against;
            defaults to the delivered (``chosen``) candidate.
        predictors: Splice backends to run; defaults to the honest baseline. Pass
            :func:`available_splice_backends` to include the wrapped SpliceAI /
            Pangolin CNNs when installed.
        threshold: Site-localization threshold (a heuristic display knob, not a
            calibrated cutoff).
        match_window: +/- nt window for the approximate cross-backend co-occurrence.
        context: Known flanking sequence. When given, every candidate (and the
            reference) is scored with the CDS in its real flanks instead of alone --
            the wrapped CNNs otherwise pad their wide context with literal ``N``.
            Reported positions stay in CDS coordinates. Better input is not
            calibration: an uncalibrated backend stays uncalibrated.

    Returns:
        A :class:`~bt4.biomodels.splice.audit.SpliceAuditReport`.

    Raises:
        ValueError: If the candidate set is empty (or per ``audit_candidate_set``).
    """
    return audit_candidate_set(
        candidate_set,
        reference=reference,
        predictors=predictors,
        threshold=threshold,
        match_window=match_window,
        context=context,
    )


def splice_crosscheck(
    dna: str,
    *,
    backend: str = "assp",
    predictor: SplicePredictor | None = None,
    threshold: float = 0.5,
    top_k: int = 3,
) -> SpliceCrossCheck:
    """Cross-check one delivered sequence with an opt-in, out-of-loop splice backend.

    The *final audit / validation pass* of BT4's splice story (CLAUDE.md §6): runs
    the named backend over ``dna`` and reports its predicted cryptic-splice sites and
    pooled risk. It is what ``bt4 validate --splice-backend`` and ``bt4 optimize
    --check-splice`` call. **Never blocking:** a backend failure -- an **ASSP** outage
    or garbled response, or a wrapped CNN's missing deps / weights -- is caught and
    reported as ``available is False`` with the reason, so a cross-check can never
    fail a run; only a genuinely invalid sequence (non-ACGT) raises. For ASSP the
    result is stamped ``network_derived is True`` and its numbers are **excluded from
    the reproducible-from-manifest guarantee** (report them as a separate advisory
    section, never in a :class:`Result` audit or manifest).

    Args:
        dna: The delivered coding sequence to audit.
        backend: Backend to construct when ``predictor`` is not given -- ``"assp"``
            (the opt-in online cross-check, requires the ``bt4[assp]`` extra),
            ``"pwm"`` (the offline baseline), or ``"pangolin"`` / ``"spliceai"``
            (the wrapped CNNs, when installed).
        predictor: An explicit backend overriding ``backend`` (e.g. a
            fixture-backed ASSP predictor).
        threshold: Localization threshold for dense-track backends (heuristic).
        top_k: Pooling depth for the reported pooled risk.

    Returns:
        A :class:`~bt4.pipeline.splice_crosscheck.SpliceCrossCheck`.

    Raises:
        ValueError: If ``dna`` is empty / non-ACGT, or ``backend`` is unknown.
    """
    return run_splice_crosscheck(
        dna, backend=backend, predictor=predictor, threshold=threshold, top_k=top_k
    )
