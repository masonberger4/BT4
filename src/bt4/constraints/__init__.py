"""BT4 constraints layer: the feasibility rules (CLAUDE.md §4.2, §6).

Exposes the :class:`~bt4.constraints.base.Constraint` contract together with the
concrete local constraints.
"""

from __future__ import annotations

from bt4.constraints.base import Constraint
from bt4.constraints.rules import ForbiddenMotifConstraint, HomopolymerConstraint

__all__ = ["Constraint", "ForbiddenMotifConstraint", "HomopolymerConstraint"]
