"""Concrete additive objective terms (CLAUDE.md §6, invariant #4).

Both terms here are strictly ``LOCAL`` with ``context_len == 0``: each codon's
contribution depends only on the codon itself, so ``delta`` and the whole-
sequence ``score`` agree by construction (invariant #4, "delta == score"). Every
term is oriented so that *larger is better*.

* :class:`CaiTerm` scores log relative-adaptiveness (``log w``). Summed over a
  sequence this is ``log(CAI) * n_scored``; the Codon Adaptation Index itself is
  ``exp(score / n_scored)``. Met/Trp/stop codons are non-degenerate and score
  ``0.0`` (they carry no coding-choice information), matching
  :meth:`bt4.biomodels.codon.tables.CodonUsageTable.cai`.
* :class:`GcProximityTerm` is an honest *per-codon* proxy for a whole-sequence
  GC target: it penalizes the squared distance between a codon's GC fraction and
  the target. It is additive by construction, so - unlike BT3's window-rescoring
  GC penalty - its ``delta`` sum equals its ``score`` exactly.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

from bt4._accel import gc_count
from bt4.domain.genetic_code import CODON_TABLE, STOP
from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons

__all__ = ["CaiTerm", "GcProximityTerm"]

# Amino acids with a single codon: no coding choice, so they carry no CAI signal.
_NON_DEGENERATE: frozenset[str] = frozenset({"M", "W"})


@dataclass(frozen=True, slots=True)
class CaiTerm:
    """Log relative-adaptiveness objective term (a weak, cheap CAI prior).

    Holds the relative-adaptiveness mapping directly (``codon -> w``) rather than
    a whole codon-usage table, so this pure objective term depends on nothing
    below ``domain``. Build it from a table via
    ``CaiTerm(table.relative_adaptiveness())``.

    Attributes:
        weights: Mapping ``codon -> w(codon)`` (relative adaptiveness, ``w`` in
            ``(0, 1]`` within each amino acid).
    """

    weights: Mapping[str, float]
    name: str = field(default="cai_logw", init=False)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``0`` - each codon is scored independently."""
        return 0

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return ``log w(codon)`` (``0.0`` for Met, Trp, and stop codons).

        Args:
            prefix: Unused (this term needs no context).
            codon: The 3-nt codon being placed.
            pos: Unused 0-based codon index.
        """
        aa = CODON_TABLE[codon.upper()]
        if aa == STOP or aa in _NON_DEGENERATE:
            return 0.0
        return math.log(self.weights[codon.upper()])

    def score(self, dna: str) -> float:
        """Return the sum of :meth:`delta` over the codons of ``dna``."""
        return sum(self.delta("", codon, pos) for pos, codon in iter_codons(dna))


@dataclass(frozen=True, slots=True)
class GcProximityTerm:
    """Per-codon GC-target proximity term (larger is closer to ``target``).

    Attributes:
        target: Desired GC fraction in ``[0, 1]``.
    """

    target: float
    name: str = field(default="gc_proximity", init=False)

    def __post_init__(self) -> None:
        """Validate the GC target.

        Raises:
            ValueError: If ``target`` is not in ``[0, 1]``.
        """
        if not 0.0 <= self.target <= 1.0:
            raise ValueError(f"gc target must be in [0, 1], got {self.target!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``0`` - each codon is scored independently."""
        return 0

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return ``-(gc_fraction(codon) - target) ** 2`` for the codon.

        Args:
            prefix: Unused (this term needs no context).
            codon: The 3-nt codon being placed.
            pos: Unused 0-based codon index.
        """
        gc = gc_count(codon) / 3.0
        return -((gc - self.target) ** 2)

    def score(self, dna: str) -> float:
        """Return the sum of :meth:`delta` over the codons of ``dna``."""
        return sum(self.delta("", codon, pos) for pos, codon in iter_codons(dna))
