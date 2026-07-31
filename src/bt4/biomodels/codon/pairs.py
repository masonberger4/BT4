"""Codon-pair bias tables built from a real reference CDS set (Coleman 2008).

A :class:`CodonPairTable` holds a *codon-pair score* (CPS) for each observed
in-frame adjacent codon pair. The score follows Coleman et al. (2008,
*Science* 320:1784), which measures whether a codon pair is over- or
under-represented relative to what its two amino acids and their two codons
would predict under independence:

    CPS(AB) = ln( N(AB) / expected(AB) )

    expected(AB) = ( N(A) * N(B) ) / ( N(X) * N(Y) ) * N(XY)

where, over the reference CDS set:

* ``N(AB)`` is the observed count of the in-frame codon pair ``A`` then ``B``;
* ``N(A)`` / ``N(B)`` are the single-codon counts of ``A`` / ``B``;
* ``X = aa(A)`` and ``Y = aa(B)`` are the amino acids the codons encode;
* ``N(X)`` / ``N(Y)`` are the amino-acid counts of ``X`` / ``Y``;
* ``N(XY)`` is the count of the adjacent amino-acid pair ``X`` then ``Y``.

A positive CPS marks an over-represented (preferred) pair; a negative CPS marks
an under-represented (avoided) pair. Larger is "more preferred", so the derived
objective term is oriented larger-is-better.

**Honesty (CLAUDE.md).** There is no bundled default codon-pair table: the
scores are only meaningful for the organism whose coding sequences they were
counted over, so callers build one from their own reference CDS set via
:func:`build_codon_pair_table`. Nothing here fabricates a default.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from bt4.domain.genetic_code import CODON_TABLE
from bt4.domain.sequence import validate_dna

__all__ = ["CodonPairTable", "build_codon_pair_table"]


@dataclass(frozen=True, slots=True)
class CodonPairTable:
    """Read-only codon-pair scores for one reference CDS set.

    Attributes:
        scores: Read-only mapping ``(codon_a, codon_b) -> cps`` for every codon
            pair observed in the reference set. Keys are upper-cased on
            construction; pairs absent from the mapping score ``0.0`` (treated as
            neutral / unobserved).
    """

    scores: Mapping[tuple[str, str], float]
    name: str = field(default="codon_pair", init=False)

    def __post_init__(self) -> None:
        """Normalize keys to upper case and freeze the mapping."""
        normalized: dict[tuple[str, str], float] = {
            (codon_a.upper(), codon_b.upper()): float(cps)
            for (codon_a, codon_b), cps in self.scores.items()
        }
        object.__setattr__(self, "scores", MappingProxyType(normalized))

    def score(self, codon_a: str, codon_b: str) -> float:
        """Return the codon-pair score of ``codon_a`` then ``codon_b``.

        Args:
            codon_a: The 5' (first) codon of the adjacent pair (case-insensitive).
            codon_b: The 3' (second) codon of the adjacent pair (case-insensitive).

        Returns:
            The stored CPS, or ``0.0`` for a pair not present in the table.
        """
        return self.scores.get((codon_a.upper(), codon_b.upper()), 0.0)


def build_codon_pair_table(
    cds_sequences: Iterable[str], *, pseudocount: float = 1.0
) -> CodonPairTable:
    """Build a :class:`CodonPairTable` from a reference set of coding sequences.

    Counts single codons, amino acids, in-frame adjacent codon pairs, and
    in-frame adjacent amino-acid pairs across ``cds_sequences``, then assigns
    each *observed* codon pair its Coleman codon-pair score::

        CPS(AB) = ln( N(AB) / ( (N(A) * N(B)) / (N(X) * N(Y)) * N(XY) ) )

    with ``X = aa(A)`` and ``Y = aa(B)`` (see the module docstring). For a CDS of
    ``k`` codons the adjacent pairs counted are ``(codon[0], codon[1])``,
    ``(codon[1], codon[2])``, ..., ``(codon[k-2], codon[k-1])`` -- successive
    codons in the same reading frame only.

    Args:
        cds_sequences: An iterable of coding DNA strings. Each is validated with
            :func:`bt4.domain.sequence.validate_dna` and must have a length that
            is a multiple of three.
        pseudocount: Additive (Laplace) smoothing added to every count that
            enters the formula -- ``N(AB)``, ``N(A)``, ``N(B)``, ``N(X)``,
            ``N(Y)``, and ``N(XY)`` -- to avoid division-by-zero and ``log(0)``.
            Defaults to ``1.0``; ``0.0`` disables smoothing (safe because only
            observed pairs, whose marginal counts are all positive, are scored).

    Returns:
        A :class:`CodonPairTable` with one CPS per observed codon pair.

    Raises:
        ValueError: If ``pseudocount`` is negative, or if any CDS fails
            :func:`~bt4.domain.sequence.validate_dna` or has a length that is not
            a multiple of three.
    """
    if pseudocount < 0:
        raise ValueError(f"pseudocount must be non-negative, got {pseudocount!r}")

    codon_counts: Counter[str] = Counter()
    aa_counts: Counter[str] = Counter()
    pair_counts: Counter[tuple[str, str]] = Counter()
    aa_pair_counts: Counter[tuple[str, str]] = Counter()

    for cds in cds_sequences:
        seq = validate_dna(cds)
        if len(seq) % 3 != 0:
            raise ValueError(f"CDS length {len(seq)} is not a multiple of three")
        codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
        amino_acids = [CODON_TABLE[codon] for codon in codons]
        codon_counts.update(codons)
        aa_counts.update(amino_acids)
        for i in range(len(codons) - 1):
            pair_counts[(codons[i], codons[i + 1])] += 1
            aa_pair_counts[(amino_acids[i], amino_acids[i + 1])] += 1

    scores: dict[tuple[str, str], float] = {}
    for (codon_a, codon_b), n_ab in pair_counts.items():
        x = CODON_TABLE[codon_a]
        y = CODON_TABLE[codon_b]
        observed = n_ab + pseudocount
        expected = (
            (codon_counts[codon_a] + pseudocount)
            * (codon_counts[codon_b] + pseudocount)
            / ((aa_counts[x] + pseudocount) * (aa_counts[y] + pseudocount))
            * (aa_pair_counts[(x, y)] + pseudocount)
        )
        scores[(codon_a, codon_b)] = math.log(observed / expected)

    return CodonPairTable(scores)
