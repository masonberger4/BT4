"""BT4 constraints layer: the feasibility rules (CLAUDE.md §4.2, §6).

Exposes the :class:`~bt4.constraints.base.Constraint` contract together with the
concrete local constraints.
"""

from __future__ import annotations

from bt4.constraints.base import Constraint
from bt4.constraints.forbidden import (
    FORBIDDEN_PRESETS,
    ForbiddenPreset,
    available_forbidden_presets,
    resolve_forbidden_motifs,
)
from bt4.constraints.gc_run import GcRunConstraint
from bt4.constraints.kozak import InternalStartConstraint
from bt4.constraints.max_repeat import MaxRepeatConstraint
from bt4.constraints.repeats import InvertedRepeatConstraint, TandemRepeatConstraint
from bt4.constraints.restriction import (
    RestrictionSiteConstraint,
    available_enzymes,
    enzyme_provenance,
    enzyme_suggestions,
    resolve_enzyme,
    unknown_enzyme_message,
)
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint
from bt4.constraints.splice_motif import (
    DEFAULT_ACCEPTOR_MOTIFS,
    DEFAULT_DONOR_MOTIFS,
    SpliceSiteMotifConstraint,
)
from bt4.constraints.uorf import UorfConstraint

__all__ = [
    "DEFAULT_ACCEPTOR_MOTIFS",
    "DEFAULT_DONOR_MOTIFS",
    "FORBIDDEN_PRESETS",
    "Constraint",
    "ForbiddenMotifConstraint",
    "ForbiddenPreset",
    "GcRunConstraint",
    "HomopolymerConstraint",
    "InternalStartConstraint",
    "InvertedRepeatConstraint",
    "MaxRepeatConstraint",
    "RestrictionSiteConstraint",
    "SpliceSiteMotifConstraint",
    "TandemRepeatConstraint",
    "UorfConstraint",
    "available_enzymes",
    "available_forbidden_presets",
    "enzyme_provenance",
    "enzyme_suggestions",
    "resolve_enzyme",
    "resolve_forbidden_motifs",
    "unknown_enzyme_message",
]
