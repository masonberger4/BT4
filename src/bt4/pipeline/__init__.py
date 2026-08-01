"""BT4 pipeline: composes biomodels + objectives + constraints + solver into runs.

This layer owns the two-stage orchestration and the Pareto sweep, and is where
metrics are recomputed from the delivered sequence and the provenance manifest is
stamped. It imports only :mod:`bt4.domain`, the pure/biomodel layers, and the
solver -- never anything above it.
"""

from __future__ import annotations

from bt4.optimize import InfeasibleError
from bt4.pipeline.optimize import (
    FrontierResult,
    OptimizeConfig,
    ValidationReport,
    run_frontier,
    run_optimize,
    run_validate,
)
from bt4.pipeline.tracks import Track, TracksResult, run_tracks, summarize

__all__ = [
    "FrontierResult",
    "InfeasibleError",
    "OptimizeConfig",
    "Track",
    "TracksResult",
    "ValidationReport",
    "run_frontier",
    "run_optimize",
    "run_tracks",
    "run_validate",
    "summarize",
]
