"""BT4 objectives layer: the objective-vector terms (CLAUDE.md §4.1, §6).

Exposes the :class:`~bt4.objectives.base.ObjectiveTerm` contract together with
the concrete additive terms and the :func:`~bt4.objectives.base.iter_codons`
helper used by every term's ``score``.
"""

from __future__ import annotations

from bt4.objectives.base import ObjectiveTerm, iter_codons
from bt4.objectives.terms import CaiTerm, GcProximityTerm

__all__ = ["CaiTerm", "GcProximityTerm", "ObjectiveTerm", "iter_codons"]
