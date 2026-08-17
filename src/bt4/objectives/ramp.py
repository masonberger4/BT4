"""5' shaping term over the first few dozen codons (a prior, NOT a mechanism).

Natural highly-expressed genes often show a "ramp" of slower (less codon-adapted)
codons across roughly the first 30-50 codons. BT4 models that *observation* as a
deliberately simple shaping term -- and is explicit that the mechanism originally
proposed for it does not survive the experiments.

**The causal claim is falsified; the 5' effect is real.** Goodman, Church & Kosuri
(2013) varied 5' codon usage and mRNA structure independently across a large
synthetic library and concluded that **reduced RNA secondary structure, not codon
rarity itself**, is responsible for the expression benefit of a "slow" 5' region.
So this term must never be described as modelling ribosome-loading dynamics: rare
5' codons are *correlated* with the effect largely because they tend to be
AT-richer and therefore less structured. The lever BT4 offers for the *validated*
mechanism is the 5' folding objective (:mod:`bt4.biomodels.folding`, reached via
``refine``/``folding_weight``), which acts on start-proximal structure directly.

What :class:`RampTerm` is, then: a **shaping prior**. It rewards lower codon
adaptiveness early in the sequence and fades linearly to no effect at and beyond
``ramp_codons``. It is oriented so that larger is better (CLAUDE.md invariant #4,
"delta == score"), which -- combined with a CAI/tAI term pulling the other way --
creates a genuine 5' trade-off the frontier can expose. It is kept because that
trade-off is real and useful to see, not because the ramp hypothesis is settled.

Because its ``delta`` reads only ``pos`` and the codon itself (never the growing
``prefix``), it is a ``POSITIONAL`` term with ``context_len == 0``, and its
per-codon deltas sum exactly to its whole-sequence ``score``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from bt4.domain.genetic_code import CODON_TABLE, STOP
from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons

__all__ = ["RampTerm"]

# Amino acids with a single codon: no coding choice, so no ramp signal to shape.
_NON_DEGENERATE: frozenset[str] = frozenset({"M", "W"})


@dataclass(frozen=True, slots=True)
class RampTerm:
    """Heuristic 5' translation-ramp shaping term (larger is better).

    Rewards *lower* codon adaptiveness in the first ``ramp_codons`` codons and
    fades linearly to zero effect at and beyond that index. Holds the relative-
    adaptiveness mapping directly (``codon -> w``) rather than a whole codon-
    usage table, so this pure objective term depends on nothing below ``domain``.
    Build it from a table via ``RampTerm(table.relative_adaptiveness())``.

    This is a heuristic prior, not a validated model of ribosome loading; see the
    module docstring.

    Attributes:
        weights: Mapping ``codon -> w(codon)`` (relative adaptiveness, ``w`` in
            ``(0, 1]`` within each amino acid).
        ramp_codons: Length of the ramp in codons; the shaping weight decays
            linearly from ``1.0`` at ``pos == 0`` to ``0.0`` at ``pos ==
            ramp_codons`` (and stays ``0.0`` beyond).
    """

    weights: Mapping[str, float]
    ramp_codons: int = 35
    name: str = field(default="ramp", init=False)

    def __post_init__(self) -> None:
        """Validate the ramp length.

        Raises:
            ValueError: If ``ramp_codons`` is less than one.
        """
        if self.ramp_codons < 1:
            raise ValueError(f"ramp_codons must be >= 1, got {self.ramp_codons!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.POSITIONAL`."""
        return Scope.POSITIONAL

    def context_len(self) -> int:
        """Return ``0`` - the delta reads only ``pos`` and the codon."""
        return 0

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return the ramp contribution of ``codon`` at codon index ``pos``.

        The value is ``-w(codon) * ramp`` where ``ramp`` decays linearly from
        ``1.0`` at ``pos == 0`` to ``0.0`` at ``pos == ramp_codons``. It is
        negative, strongest at ``pos == 0``, so maximizing it prefers *slower*
        (lower-``w``) codons early. Met, Trp, and stop codons carry no coding
        choice and score ``0.0`` at any position.

        Args:
            prefix: Unused (this term needs no context).
            codon: The 3-nt codon being placed.
            pos: 0-based codon index (drives the positional decay).
        """
        aa = CODON_TABLE[codon.upper()]
        if aa == STOP or aa in _NON_DEGENERATE:
            return 0.0
        ramp = max(0.0, 1.0 - pos / self.ramp_codons)
        return -self.weights[codon.upper()] * ramp

    def score(self, dna: str) -> float:
        """Return the sum of :meth:`delta` over the codons of ``dna``."""
        return sum(self.delta("", codon, pos) for pos, codon in iter_codons(dna))
