"""The standard genetic code — the pure, authoritative source of truth.

This module knows how to translate coding DNA to protein and how to enumerate
the synonymous codons of each amino acid. It depends on nothing (stdlib only)
and holds no optimization logic.

The one invariant everything downstream leans on:
``translate(dna)`` is the exact inverse relation of a valid back-translation —
for any protein ``p`` and any synonymous DNA ``d`` that BT4 emits,
``translate(d) == p (+ stop)``.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Final

__all__ = [
    "AA_TO_CODONS",
    "AMINO_ACIDS",
    "CODON_TABLE",
    "STOP",
    "is_stop",
    "synonymous_codons",
    "translate",
]

STOP: Final = "*"
"""Sentinel used for stop codons in translated protein strings."""

# The standard genetic code (DNA codons -> single-letter amino acid; '*' = stop).
_CODON_TABLE: Final[dict[str, str]] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODON_TABLE: Final = MappingProxyType(_CODON_TABLE)
"""Read-only view of the standard genetic code (codon -> amino acid)."""


def _build_aa_to_codons() -> MappingProxyType[str, tuple[str, ...]]:
    grouped: dict[str, list[str]] = {}
    for codon, aa in _CODON_TABLE.items():
        grouped.setdefault(aa, []).append(codon)
    # Sort codons within each amino acid for determinism (BT3 invariant kept).
    frozen: dict[str, tuple[str, ...]] = {
        aa: tuple(sorted(codons)) for aa, codons in grouped.items()
    }
    return MappingProxyType(frozen)


AA_TO_CODONS: Final = _build_aa_to_codons()
"""Read-only map amino acid -> sorted tuple of synonymous codons (includes '*')."""

AMINO_ACIDS: Final = frozenset(aa for aa in AA_TO_CODONS if aa != STOP)
"""The 20 standard amino acids (stop excluded)."""


def is_stop(codon: str) -> bool:
    """Return True if ``codon`` is a stop codon."""
    return _CODON_TABLE.get(codon.upper()) == STOP


def synonymous_codons(amino_acid: str) -> tuple[str, ...]:
    """Return the sorted synonymous codons for ``amino_acid`` (or stop ``'*'``).

    Args:
        amino_acid: A single-letter amino acid code, or ``'*'`` for stop.

    Returns:
        A deterministic (sorted) tuple of codons.

    Raises:
        KeyError: If ``amino_acid`` is not a valid single-letter code.
    """
    return AA_TO_CODONS[amino_acid.upper()]


def translate(dna: str) -> str:
    """Translate coding DNA (5'->3', frame 0) to a protein string.

    Stop codons translate to ``'*'``. This is the reference implementation used
    to enforce the round-trip invariant everywhere in BT4.

    Args:
        dna: A DNA string whose length is a multiple of three.

    Returns:
        The translated protein, with ``'*'`` marking any stop codon.

    Raises:
        ValueError: If ``dna`` length is not a multiple of three or contains an
            unknown codon.
    """
    seq = dna.upper()
    if len(seq) % 3 != 0:
        raise ValueError(f"DNA length {len(seq)} is not a multiple of three")
    out: list[str] = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        aa = _CODON_TABLE.get(codon)
        if aa is None:
            raise ValueError(f"unknown codon {codon!r} at position {i}")
        out.append(aa)
    return "".join(out)
