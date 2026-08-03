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
from bt4.constraints.restriction import RestrictionSiteConstraint, available_enzymes
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint
from bt4.constraints.uorf import UorfConstraint

__all__ = [
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
    "TandemRepeatConstraint",
    "UorfConstraint",
    "available_enzymes",
    "available_forbidden_presets",
    "resolve_forbidden_motifs",
]
