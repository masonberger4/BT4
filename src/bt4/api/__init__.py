"""BT4 api: the stable, print-free entry point every frontend calls.

The CLI, the BT4 Studio desktop app, and the optional HTTP service all go
through this module and nothing below it. Results carry their own recomputed
metrics, optimality certificate, and provenance manifest, so a caller never has
to reach into the optimizer or the biomodels to understand what it got back.

This layer never prints; it raises on error and returns immutable results.
"""

from __future__ import annotations

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.tables import available_organisms
from bt4.constraints import available_enzymes
from bt4.domain import Result
from bt4.io import parse_fasta, read_fasta, result_to_dict, result_to_json, to_fasta
from bt4.pipeline import (
    FrontierResult,
    InfeasibleError,
    OptimizeConfig,
    ValidationReport,
    run_frontier,
    run_optimize,
    run_validate,
)

__all__ = [
    "FrontierResult",
    "InfeasibleError",
    "OptimizeConfig",
    "Result",
    "ValidationReport",
    "available_enzymes",
    "available_organisms",
    "build_table",
    "count_codons",
    "frontier",
    "optimize",
    "parse_fasta",
    "read_fasta",
    "result_to_dict",
    "result_to_json",
    "to_fasta",
    "validate",
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
    protein: str, config: OptimizeConfig | None = None, steps: int = 11
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

    Returns:
        A :class:`FrontierResult`: the non-dominated frontier plus the full
        result behind each point and a top-level manifest.

    Raises:
        ValueError: On an invalid protein, unknown organism, or ``steps < 1``.
        bt4.optimize.InfeasibleError: If the constraints admit no feasible codon.
    """
    return run_frontier(protein, config, steps)


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
