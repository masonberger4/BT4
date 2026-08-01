"""Per-window dinucleotide (CpG / general 2-mer) reporting profiles.

Dinucleotide composition is a real mRNA-design knob: the CpG step (``CG``) is
sensed by innate immunity, and its *local* density -- not just its whole-sequence
total -- is what matters for stealth vs. immunogenicity. This module reports how
that density varies along the sequence, as a sliding-window profile.

Like :mod:`bt4.objectives.minmax`, this module keeps the honest split between an
optimizer objective and a reporting metric:

* :class:`~bt4.objectives.dinucleotide.DinucleotideTerm` (in the sibling
  ``dinucleotide`` module) is the **DP-participating objective term**: a
  genuinely additive count that the trellis optimizes exactly (``delta ==
  score``). It controls the *whole-sequence total* of a 2-mer.
* :func:`dinucleotide_profile` (here) is a **reporting metric only**: the classic
  per-window 2-mer count/density used for display and auditing. A window count is
  a non-additive-per-codon window statistic (a single occurrence sits in up to
  ``window - 1`` overlapping windows), so -- exactly as with
  :func:`~bt4.objectives.minmax.min_max_profile` -- it is deliberately **not** an
  ``ObjectiveTerm`` and never enters the DP. Folding it in per codon would repeat
  the BT3 sin CLAUDE.md invariant #4 forbids. The optimizer's additive CpG
  handling stays in :class:`~bt4.objectives.dinucleotide.DinucleotideTerm`.

Both functions are pure, deterministic, and recompute every number from the input
sequence (CLAUDE.md invariant #2).
"""

from __future__ import annotations

from bt4.domain.sequence import validate_dna

__all__ = ["cpg_profile", "dinucleotide_profile"]


def dinucleotide_profile(
    dna: str,
    dinucleotide: str = "CG",
    window: int = 50,
    *,
    density: bool = False,
) -> list[float]:
    """Compute the per-window dinucleotide count/density profile (a metric).

    A window of ``window`` consecutive nucleotides contains ``window - 1``
    adjacent base-pair *slots*, and each window's value is the number of those
    slots that equal ``dinucleotide`` (overlapping occurrences fully inside the
    window). The window slides one nucleotide at a time, so consecutive profile
    entries share ``window - 1`` nucleotides; the value at index ``i`` describes
    ``dna[i : i + window]``.

    Occurrences are counted exactly, over the upper-cased sequence. With
    ``density=False`` (the default) each value is the raw count, an integer stored
    as a ``float`` in ``[0, window - 1]``. With ``density=True`` each value is that
    count divided by the number of slots (``window - 1``), a fraction in
    ``[0, 1]``; a degenerate ``window == 1`` (no slots) yields ``0.0`` for every
    window rather than dividing by zero.

    This is a **reporting/profile** helper only: it is deterministic, read-only,
    and does **not** participate in the optimizer (see the module docstring and
    :class:`~bt4.objectives.dinucleotide.DinucleotideTerm` for the additive DP
    term). Unlike codon-window metrics it operates purely on nucleotides, so
    ``len(dna)`` need not be a multiple of three.

    Args:
        dna: The sequence to profile, over ``{A,C,G,T}`` (case-insensitive).
        dinucleotide: The 2-mer to count, exactly two ``ACGT`` characters
            (case-insensitive; e.g. ``"CG"`` for CpG, ``"TA"`` for UpA).
        window: Number of consecutive nucleotides per window (must be ``>= 1``).
        density: When ``True`` return per-slot fractions in ``[0, 1]``; when
            ``False`` (default) return raw counts.

    Returns:
        The per-window values, of length ``max(0, len(dna) - window + 1)`` (empty
        when the sequence is shorter than one window).

    Raises:
        ValueError: If ``dna`` is empty or non-ACGT, if ``dinucleotide`` is not
            exactly two ``ACGT`` characters, or if ``window`` is less than one.
    """
    seq = validate_dna(dna)
    dinuc = validate_dna(dinucleotide)
    if len(dinuc) != 2:
        raise ValueError(
            f"dinucleotide must be exactly 2 ACGT characters, got {dinucleotide!r}"
        )
    if window < 1:
        raise ValueError(f"window must be >= 1, got {window!r}")

    n_windows = len(seq) - window + 1
    if n_windows <= 0:
        return []

    slots = window - 1  # adjacent base-pair positions inside one window
    hits = [1 if seq[i : i + 2] == dinuc else 0 for i in range(len(seq) - 1)]

    # Slide a running sum over `slots` slots; when slots == 0 every value is 0.
    current = sum(hits[0:slots])
    counts = [current]
    for start in range(1, n_windows):
        current += hits[start - 1 + slots] - hits[start - 1]
        counts.append(current)

    if density:
        if slots == 0:
            return [0.0] * n_windows
        return [count / slots for count in counts]
    return [float(count) for count in counts]


def cpg_profile(dna: str, window: int = 50, *, density: bool = False) -> list[float]:
    """Compute the per-window CpG (``CG``) profile (a reporting metric).

    Convenience wrapper for :func:`dinucleotide_profile` with
    ``dinucleotide="CG"``; see that function for the exact window convention,
    return length, and the ``density`` option.

    Args:
        dna: The sequence to profile, over ``{A,C,G,T}`` (case-insensitive).
        window: Number of consecutive nucleotides per window (must be ``>= 1``).
        density: When ``True`` return per-slot fractions in ``[0, 1]``; when
            ``False`` (default) return raw CpG counts.

    Returns:
        The per-window CpG values, of length ``max(0, len(dna) - window + 1)``.

    Raises:
        ValueError: If ``dna`` is empty or non-ACGT, or if ``window`` is less
            than one.
    """
    return dinucleotide_profile(dna, "CG", window, density=density)
