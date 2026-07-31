"""Codon-pair bias objective term (CLAUDE.md §6, Scope.PAIRWISE).

:class:`CpbTerm` scores a coding sequence by the codon-pair scores of its
in-frame adjacent codon pairs (see :mod:`bt4.biomodels.codon.pairs`). It is the
canonical ``PAIRWISE`` objective term: each codon's contribution depends on the
single immediately-preceding codon, so the codon trellis can solve it exactly by
extending its state with the previous codon (``context_len == 3``).

The term is oriented **larger is better**: a higher total codon-pair score
prefers over-represented (positive-CPS) pairs. A pipeline that instead wants to
*attenuate* codon-pair bias -- e.g. codon-pair deoptimization for attenuated
vaccine design (Coleman et al. 2008) -- flips this term's weight sign at the
objective-vector level, which turns "prefer over-represented pairs" into "prefer
under-represented pairs"; this term itself stays honestly larger-is-better.

Invariant #4 (``delta == score``): :meth:`CpbTerm.score` sums exactly the same
per-pair contributions that :meth:`CpbTerm.delta` yields as the DP grows the
prefix codon by codon, so the two agree by construction.

This module imports only :mod:`bt4.domain` and its own objectives layer. Like
:class:`~bt4.objectives.terms.CaiTerm` (which holds a weights mapping, not a
whole codon-usage table), ``CpbTerm`` holds the codon-pair *scores mapping*
directly rather than a :class:`~bt4.biomodels.codon.pairs.CodonPairTable`, so
``bt4.objectives`` carries no dependency on ``bt4.biomodels`` (the strict-
layering rule of CLAUDE.md §3). Build it from a table via ``CpbTerm(table.scores)``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons

__all__ = ["CpbTerm"]


@dataclass(frozen=True, slots=True)
class CpbTerm:
    """Codon-pair bias objective term (larger total CPS is better).

    Attributes:
        scores: Mapping ``(codon_a, codon_b) -> codon-pair score`` with upper-case
            codon keys (e.g. from ``CodonPairTable.scores``). Unknown pairs score
            ``0.0`` (neutral), so a partial table degrades gracefully.
    """

    scores: Mapping[tuple[str, str], float]
    name: str = field(default="codon_pair", init=False)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.PAIRWISE`."""
        return Scope.PAIRWISE

    def context_len(self) -> int:
        """Return ``3`` -- the single immediately-preceding codon."""
        return 3

    def _pair(self, prev: str, codon: str) -> float:
        return self.scores.get((prev.upper(), codon.upper()), 0.0)

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return the codon-pair score of the previous codon then ``codon``.

        Args:
            prefix: The DNA placed so far; its last three characters are the
                previous codon.
            codon: The 3-nt codon being placed at index ``pos``.
            pos: The 0-based codon index. The first codon (``pos == 0``) has no
                predecessor, so it contributes ``0.0``.

        Returns:
            The pair score of ``(previous_codon, codon)``, or ``0.0`` at ``pos == 0``.
        """
        if pos == 0:
            return 0.0
        return self._pair(prefix[-3:], codon)

    def score(self, dna: str) -> float:
        """Return the summed codon-pair score over adjacent codon pairs of ``dna``.

        Equivalent to accumulating :meth:`delta` over the codons of ``dna`` with a
        growing prefix, so ``score`` equals the sum of ``delta`` (invariant #4).

        Args:
            dna: Coding DNA whose length is a multiple of three.

        Returns:
            The total codon-pair score.

        Raises:
            ValueError: If ``len(dna)`` is not a multiple of three.
        """
        total = 0.0
        prefix = ""
        for pos, codon in iter_codons(dna):
            if pos != 0:
                total += self._pair(prefix[-3:], codon)
            prefix += codon
        return total
