"""BT4 objectives layer: the objective-vector terms (CLAUDE.md §4.1, §6).

Exposes the :class:`~bt4.objectives.base.ObjectiveTerm` contract together with
the concrete additive terms and the :func:`~bt4.objectives.base.iter_codons`
helper used by every term's ``score``.
"""

from __future__ import annotations

from bt4.objectives.base import ObjectiveTerm, iter_codons
from bt4.objectives.codon_pair import CpbTerm
from bt4.objectives.dinucleotide import DinucleotideTerm
from bt4.objectives.minmax import MinMaxTerm, min_max_profile
from bt4.objectives.ramp import RampTerm
from bt4.objectives.terms import CaiTerm, GcProximityTerm

__all__ = [
    "CaiTerm",
    "CpbTerm",
    "DinucleotideTerm",
    "GcProximityTerm",
    "MinMaxTerm",
    "ObjectiveTerm",
    "RampTerm",
    "iter_codons",
    "min_max_profile",
]
