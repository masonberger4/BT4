"""The tRNA-adaptation-index objective term (CLAUDE.md sections 6, 4.1).

:class:`TaiTerm` is the tAI analogue of :class:`~bt4.objectives.terms.CaiTerm`: a
strictly ``LOCAL``, ``context_len == 0`` additive term whose per-codon
contribution is ``log w(codon)``, where ``w`` is the dos Reis relative
adaptiveness derived from real tRNA gene copy numbers (built by
:mod:`bt4.biomodels.codon.tai`). Summed over a sequence this is
``log(tAI) * n_scored``; the index itself is ``exp(score / n_scored)``.

Like :class:`CaiTerm`, it holds the ``codon -> w`` mapping directly rather than a
table, so this pure objective term depends on nothing below ``domain``. Build it
with ``TaiTerm(tai_table.relative_adaptiveness())``. Codons absent from the
mapping (Met and stops, which carry no synonymous choice and are excluded from
tAI exactly as in the reference) score ``0.0``. Oriented so larger is better.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons

__all__ = ["TaiTerm"]


@dataclass(frozen=True, slots=True)
class TaiTerm:
    """Log tRNA-adaptation objective term (a mechanistic, tRNA-grounded prior).

    Attributes:
        weights: Mapping ``codon -> w(codon)`` (dos Reis relative adaptiveness,
            ``w`` in ``(0, 1]``). Typically the 60 scoreable codons from
            :meth:`bt4.biomodels.codon.tai.TaiTable.relative_adaptiveness`;
            codons not present (Met, stops) contribute ``0.0``.
    """

    weights: Mapping[str, float]
    name: str = field(default="tai_logw", init=False)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``0`` - each codon is scored independently."""
        return 0

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return ``log w(codon)``, or ``0.0`` for codons absent from ``weights``.

        Args:
            prefix: Unused (this term needs no context).
            codon: The 3-nt codon being placed.
            pos: Unused 0-based codon index.
        """
        w = self.weights.get(codon.upper())
        return 0.0 if w is None else math.log(w)

    def score(self, dna: str) -> float:
        """Return the sum of :meth:`delta` over the codons of ``dna``."""
        return sum(self.delta("", codon, pos) for pos, codon in iter_codons(dna))
