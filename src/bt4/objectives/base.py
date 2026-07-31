"""Objective-term vocabulary for the objectives package.

The :class:`~bt4.domain.contracts.ObjectiveTerm` protocol itself lives in
``domain`` (so the optimizer and pipeline can speak it while every pure layer
imports only ``domain``); it is re-exported here so ``from bt4.objectives import
ObjectiveTerm`` keeps working. This module also owns the small ``iter_codons``
helper the concrete terms share.
"""

from __future__ import annotations

from bt4.domain.contracts import ObjectiveTerm

__all__ = ["ObjectiveTerm", "iter_codons"]


def iter_codons(dna: str) -> list[tuple[int, str]]:
    """Return ``(codon_index, codon)`` pairs for a length-multiple-of-3 ``dna``.

    Raises:
        ValueError: If ``len(dna)`` is not a multiple of three.
    """
    if len(dna) % 3 != 0:
        raise ValueError(f"DNA length {len(dna)} is not a multiple of three")
    return [(i // 3, dna[i : i + 3]) for i in range(0, len(dna), 3)]
