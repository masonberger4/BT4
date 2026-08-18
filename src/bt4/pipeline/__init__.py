"""BT4 pipeline: composes biomodels + objectives + constraints + solver into runs.

This layer owns the two-stage orchestration and the Pareto sweep, and is where
metrics are recomputed from the delivered sequence and the provenance manifest is
stamped. It imports only :mod:`bt4.domain`, the pure/biomodel layers, and the
solver -- never anything above it.
"""

from __future__ import annotations

from bt4.optimize import InfeasibleError
from bt4.pipeline.candidates import (
    Candidate,
    CandidateSet,
    assemble_and_rank_candidates,
)
from bt4.pipeline.construct import ConstructAudit, EnzymeOccurrence, audit_construct
from bt4.pipeline.expression_gate import (
    BASELINES,
    GateComparison,
    GateSettings,
    run_panel_gate,
    score_panel,
)
from bt4.pipeline.library import LibraryResult, run_library
from bt4.pipeline.optimize import (
    FrontierResult,
    OptimizeConfig,
    ValidationReport,
    run_frontier,
    run_optimize,
    run_validate,
)
from bt4.pipeline.packaging import (
    PACKAGING_LIMITS,
    PackagingReport,
    packaging_report,
)
from bt4.pipeline.presets import (
    APPLICATION_PRESETS,
    ApplicationPreset,
    apply_preset,
    available_presets,
    resolve_preset,
)
from bt4.pipeline.rerank import rerank_by_expression
from bt4.pipeline.splice_audit import audit_candidate_set, available_splice_backends
from bt4.pipeline.splice_crosscheck import (
    CrossCheckSite,
    SpliceCrossCheck,
    resolve_splice_backend,
    run_splice_crosscheck,
)
from bt4.pipeline.splice_gate import (
    SPLICE_BASELINES,
    AlignmentDiagnostic,
    SpliceGateComparison,
    SpliceGateSettings,
    run_splice_panel_gate,
    score_splice_panel,
)
from bt4.pipeline.tracks import Track, TracksResult, run_tracks, summarize

__all__ = [
    "APPLICATION_PRESETS",
    "BASELINES",
    "PACKAGING_LIMITS",
    "SPLICE_BASELINES",
    "AlignmentDiagnostic",
    "ApplicationPreset",
    "Candidate",
    "CandidateSet",
    "ConstructAudit",
    "CrossCheckSite",
    "EnzymeOccurrence",
    "FrontierResult",
    "GateComparison",
    "GateSettings",
    "InfeasibleError",
    "LibraryResult",
    "OptimizeConfig",
    "PackagingReport",
    "SpliceCrossCheck",
    "SpliceGateComparison",
    "SpliceGateSettings",
    "Track",
    "TracksResult",
    "ValidationReport",
    "apply_preset",
    "assemble_and_rank_candidates",
    "audit_candidate_set",
    "audit_construct",
    "available_presets",
    "available_splice_backends",
    "packaging_report",
    "rerank_by_expression",
    "resolve_preset",
    "resolve_splice_backend",
    "run_frontier",
    "run_library",
    "run_optimize",
    "run_panel_gate",
    "run_splice_crosscheck",
    "run_splice_panel_gate",
    "run_tracks",
    "run_validate",
    "score_panel",
    "score_splice_panel",
    "summarize",
]
