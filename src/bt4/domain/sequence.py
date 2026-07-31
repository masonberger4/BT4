"""Sequence validation helpers for the pure domain layer."""

from __future__ import annotations

from .genetic_code import AMINO_ACIDS

__all__ = ["DNA_BASES", "gc_fraction", "validate_dna", "validate_protein"]

DNA_BASES: frozenset[str] = frozenset("ACGT")


def validate_protein(protein: str) -> str:
    """Validate and normalize a protein string.

    Args:
        protein: A single-letter amino-acid string. May carry a trailing ``'*'``
            stop marker, which is rejected here (proteins are stop-free; the
            optimizer appends the stop).

    Returns:
        The upper-cased protein.

    Raises:
        ValueError: If empty or containing any non-amino-acid character.
    """
    p = protein.strip().upper()
    if not p:
        raise ValueError("protein sequence is empty")
    bad = sorted({ch for ch in p if ch not in AMINO_ACIDS})
    if bad:
        raise ValueError(f"protein contains non-amino-acid characters: {bad}")
    return p


def validate_dna(dna: str) -> str:
    """Validate and normalize a DNA string.

    Args:
        dna: A string over the alphabet ``{A,C,G,T}`` (case-insensitive).

    Returns:
        The upper-cased DNA.

    Raises:
        ValueError: If empty or containing any non-ACGT character.
    """
    d = dna.strip().upper()
    if not d:
        raise ValueError("DNA sequence is empty")
    bad = sorted({ch for ch in d if ch not in DNA_BASES})
    if bad:
        raise ValueError(f"DNA contains non-ACGT characters: {bad}")
    return d


def gc_fraction(dna: str) -> float:
    """Return the GC fraction of ``dna`` in ``[0, 1]`` (0.0 for empty input)."""
    if not dna:
        return 0.0
    gc = sum(1 for ch in dna.upper() if ch in ("G", "C"))
    return gc / len(dna)
