"""IUPAC nucleotide code matching for degenerate DNA patterns.

Restriction recognition sites (and other motif families) are often written with
IUPAC ambiguity codes -- ``N`` for any base, ``R`` for a purine, ``S`` for a
strong pair, and so on. This module is the small, pure engine that decides when
an ACGT sequence matches such a degenerate pattern and how to take an
IUPAC-aware reverse complement (needed because restriction sites are
double-stranded).

It is stdlib-only and depends on nothing else in BT4, so every layer that speaks
motifs can lean on it. All patterns are validated: a pattern must be a non-empty
string of IUPAC codes, otherwise a :class:`ValueError` is raised.
"""

from __future__ import annotations

__all__ = [
    "IUPAC",
    "IUPAC_COMPLEMENT",
    "find_iupac",
    "is_iupac",
    "matches_at",
    "reverse_complement_iupac",
]

IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
"""Each IUPAC code mapped to the set of concrete ACGT bases it stands for."""

# Watson-Crick complement of the four concrete bases.
_BASE_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Inverse of IUPAC: a base-set maps back to its (unique) IUPAC code. Used to name
# the complement of an ambiguity code without hand-maintaining a second table.
_SET_TO_CODE: dict[frozenset[str], str] = {bases: code for code, bases in IUPAC.items()}


def _complement_code(code: str) -> str:
    """Return the IUPAC code standing for the complement of ``code``'s base set."""
    complemented = frozenset(_BASE_COMPLEMENT[base] for base in IUPAC[code])
    return _SET_TO_CODE[complemented]


IUPAC_COMPLEMENT: dict[str, str] = {code: _complement_code(code) for code in IUPAC}
"""Each IUPAC code mapped to the code for its complementary base set."""


def is_iupac(pattern: str) -> bool:
    """Return ``True`` iff ``pattern`` is a non-empty string of IUPAC codes.

    Args:
        pattern: The candidate pattern (case-sensitive; codes are upper-case).

    Returns:
        ``True`` when every character is a recognized IUPAC code and the pattern
        is non-empty; ``False`` otherwise.
    """
    return bool(pattern) and all(ch in IUPAC for ch in pattern)


def matches_at(window: str, pattern: str, start: int) -> bool:
    """Return ``True`` iff ``pattern`` matches ``window`` beginning at ``start``.

    Each ACGT base of ``window`` in the aligned span must belong to the base set
    of the corresponding IUPAC code in ``pattern``. A span that runs past either
    end of ``window`` simply does not match (no exception).

    Args:
        window: An ACGT sequence to test against the pattern.
        pattern: A non-empty IUPAC pattern.
        start: The 0-based index in ``window`` where the pattern is aligned.

    Returns:
        ``True`` when the pattern matches at ``start``, else ``False``.

    Raises:
        ValueError: If ``pattern`` is not a valid IUPAC string.
    """
    if not is_iupac(pattern):
        raise ValueError(f"not a valid IUPAC pattern: {pattern!r}")
    end = start + len(pattern)
    if start < 0 or end > len(window):
        return False
    return all(window[start + i] in IUPAC[pattern[i]] for i in range(len(pattern)))


def find_iupac(text: str, pattern: str) -> list[int]:
    """Return every start index where ``pattern`` occurs in ACGT ``text``.

    Args:
        text: The ACGT sequence to search.
        pattern: A non-empty IUPAC pattern.

    Returns:
        Ascending list of 0-based start indices of all (possibly overlapping)
        occurrences.

    Raises:
        ValueError: If ``pattern`` is not a valid IUPAC string.
    """
    if not is_iupac(pattern):
        raise ValueError(f"not a valid IUPAC pattern: {pattern!r}")
    span = len(pattern)
    return [i for i in range(len(text) - span + 1) if matches_at(text, pattern, i)]


def reverse_complement_iupac(pattern: str) -> str:
    """Return the IUPAC-aware reverse complement of ``pattern``.

    Each code is replaced by the code for its complementary base set, then the
    whole string is reversed -- so a palindromic site such as ``GAATTC`` (or
    ``GANTC``, palindromic through the ambiguous ``N``) maps back to itself.

    Args:
        pattern: A non-empty IUPAC pattern (upper-cased before processing).

    Returns:
        The reverse-complemented IUPAC pattern.

    Raises:
        ValueError: If ``pattern`` is not a valid IUPAC string.
    """
    upper = pattern.upper()
    if not is_iupac(upper):
        raise ValueError(f"not a valid IUPAC pattern: {pattern!r}")
    return "".join(IUPAC_COMPLEMENT[ch] for ch in reversed(upper))
