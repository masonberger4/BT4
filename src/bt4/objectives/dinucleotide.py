"""Dinucleotide-content objective term (CpG / UpA control) -- exact counting.

Dinucleotide composition is a real, controllable knob in mRNA design. The CpG
dinucleotide (the ``CG`` step) is sensed by innate-immunity pathways, so
*depleting* CpG makes a transcript stealthier while *elevating* it can raise
immunogenicity for a vaccine; the UpA dinucleotide (the ``TA`` step) is likewise
a design lever. :class:`DinucleotideTerm` lets the optimizer push a chosen 2-mer
up or down.

This term is deliberately dull: it is exact counting, not a learned or fitted
model. There is no biology to calibrate and nothing to validate beyond
arithmetic -- ``score`` is just ``sign x (number of occurrences of the 2-mer)``,
and every reported number is recomputed from the delivered sequence.

Because a dinucleotide can straddle a codon boundary (the ``C`` ending one codon
and the ``G`` beginning the next), the term is ``PAIRWISE`` with
``context_len == 1``: :meth:`DinucleotideTerm.delta` needs the final base of the
prefix to see a boundary-straddling occurrence. Following BT4's orientation
convention every term is written so that **larger is better**, so the sign is
``-1`` for depletion (fewer occurrences scores higher) and ``+1`` for elevation.
The per-codon deltas sum exactly to the whole-sequence ``score`` (CLAUDE.md
invariant #4, "delta == score"): each occurrence is attributed to the single
codon that contains its second (end) base, so it is counted once and only once.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt4.domain.scope import Scope
from bt4.domain.sequence import validate_dna

__all__ = ["DinucleotideTerm"]

_DIRECTIONS: frozenset[str] = frozenset({"deplete", "elevate"})


@dataclass(frozen=True, slots=True)
class DinucleotideTerm:
    """Push a chosen dinucleotide's count up or down (larger is better).

    The term counts overlapping occurrences of a fixed 2-mer (e.g. ``CG`` for
    CpG, ``TA`` for UpA) across the coding sequence and orients the score so
    that maximizing it moves the count in the requested ``direction``. It is
    exact counting only -- no learned model, no fitted parameters.

    Attributes:
        dinucleotide: The 2-mer to control, exactly two ``ACGT`` characters
            (case-insensitive; stored upper-cased).
        direction: ``"deplete"`` to reward *fewer* occurrences (sign ``-1``) or
            ``"elevate"`` to reward *more* (sign ``+1``).
    """

    dinucleotide: str
    direction: str = "deplete"

    def __post_init__(self) -> None:
        """Normalize and validate the dinucleotide and direction.

        Raises:
            ValueError: If ``dinucleotide`` is not exactly two ``ACGT``
                characters, or ``direction`` is not ``"deplete"`` or
                ``"elevate"``.
        """
        dinuc = validate_dna(self.dinucleotide)
        if len(dinuc) != 2:
            raise ValueError(
                f"dinucleotide must be exactly 2 ACGT characters, got {self.dinucleotide!r}"
            )
        object.__setattr__(self, "dinucleotide", dinuc)
        if self.direction not in _DIRECTIONS:
            raise ValueError(
                f"direction must be one of {sorted(_DIRECTIONS)}, got {self.direction!r}"
            )

    @property
    def name(self) -> str:
        """Return the stable identifier, e.g. ``"dinuc_cg_deplete"``."""
        return f"dinuc_{self.dinucleotide.lower()}_{self.direction}"

    @property
    def sign(self) -> float:
        """Return ``-1.0`` for depletion or ``+1.0`` for elevation."""
        return -1.0 if self.direction == "deplete" else 1.0

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.PAIRWISE`."""
        return Scope.PAIRWISE

    def context_len(self) -> int:
        """Return ``1`` - the delta needs the last base of the prefix."""
        return 1

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return this codon's signed contribution to the dinucleotide count.

        Counts occurrences of the 2-mer in ``prefix[-1:] + codon`` whose end
        base lies inside ``codon`` -- i.e. occurrences the new codon is
        responsible for, including one straddling the codon boundary (the last
        base of the prefix paired with the first base of ``codon``). Attributing
        each occurrence to the codon holding its end base counts it exactly once
        across the sequence, so these deltas sum to :meth:`score`.

        Args:
            prefix: The DNA placed so far; only its final base is read.
            codon: The 3-nt codon being placed.
            pos: Unused 0-based codon index.
        """
        lead = prefix[-1:]
        window = (lead + codon).upper()
        offset = len(lead)  # index in ``window`` where the codon portion begins
        count = sum(
            1
            for start in range(len(window) - 1)
            if window[start : start + 2] == self.dinucleotide and start + 1 >= offset
        )
        return self.sign * count

    def score(self, dna: str) -> float:
        """Return ``sign x`` the total overlapping count of the 2-mer in ``dna``.

        Equals the running sum of :meth:`delta` over the codons of ``dna``
        (CLAUDE.md invariant #4).

        Args:
            dna: The coding sequence to score.
        """
        seq = validate_dna(dna)
        count = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == self.dinucleotide)
        return self.sign * count
