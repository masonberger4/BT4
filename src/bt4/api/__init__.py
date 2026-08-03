"""BT4 api: the stable, print-free entry point every frontend calls.

The CLI, the BT4 Studio desktop app, and the optional HTTP service all go
through this module and nothing below it. Results carry their own recomputed
metrics, optimality certificate, and provenance manifest, so a caller never has
to reach into the optimizer or the biomodels to understand what it got back.

This layer never prints; it raises on error and returns immutable results.
"""

from __future__ import annotations

from collections.abc import Callable

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.tables import available_organisms
from bt4.biomodels.codon.tai import available_tai_organisms
from bt4.biomodels.expression import ExpressionPredictor, ExpressionResult
from bt4.biomodels.expression import default as expression_model
from bt4.constraints import (
    ForbiddenPreset,
    available_enzymes,
    available_forbidden_presets,
)
from bt4.domain import AMINO_ACIDS, Result, Severity, Violation, validate_protein
from bt4.io import parse_fasta, read_fasta, result_to_dict, result_to_json, to_fasta
from bt4.pipeline import (
    FrontierResult,
    InfeasibleError,
    LibraryResult,
    OptimizeConfig,
    Track,
    TracksResult,
    ValidationReport,
    rerank_by_expression,
    run_frontier,
    run_library,
    run_optimize,
    run_tracks,
    run_validate,
    summarize,
)

__all__ = [
    "AMINO_ACIDS",
    "ExpressionPredictor",
    "ExpressionResult",
    "ForbiddenPreset",
    "FrontierResult",
    "InfeasibleError",
    "LibraryResult",
    "OptimizeConfig",
    "Result",
    "Severity",
    "Track",
    "TracksResult",
    "ValidationReport",
    "Violation",
    "available_enzymes",
    "available_forbidden_presets",
    "available_organisms",
    "available_tai_organisms",
    "build_table",
    "count_codons",
    "expression_model",
    "frontier",
    "library",
    "optimize",
    "parse_fasta",
    "read_fasta",
    "rerank_by_expression",
    "result_to_dict",
    "result_to_json",
    "summarize",
    "to_fasta",
    "tracks",
    "validate",
    "validate_protein",
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
    nt_window: int = 50,
    codon_window: int = 18,
) -> TracksResult:
    """Compute per-site composition tracks for a coding sequence.

    Sliding-window **reporting** profiles (GC fraction, CpG density, and -- when
    codon-aligned -- %MinMax) so a delivered sequence's composition can be
    audited or plotted position-by-position. Nothing here feeds the optimizer.

    Args:
        dna: An ACGT coding sequence (case-insensitive).
        organism: Codon-usage table key/alias for the %MinMax reference.
        nt_window: Window (nucleotides) for the GC and CpG tracks.
        codon_window: Window (codons) for the %MinMax track.

    Returns:
        A :class:`~bt4.pipeline.tracks.TracksResult` bundling the named tracks.

    Raises:
        ValueError: On non-ACGT input, a non-positive window, or unknown organism.
    """
    return run_tracks(dna, organism, nt_window=nt_window, codon_window=codon_window)


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
