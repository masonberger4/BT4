"""Tests for the per-window dinucleotide reporting profile.

These pin the reporting metric :func:`~bt4.objectives.dinuc_profile.dinucleotide_profile`
(and its :func:`~bt4.objectives.dinuc_profile.cpg_profile` wrapper): the window-count
length convention (including the empty result when the sequence is shorter than one
window), an exact hand-counted profile over a small sequence, the density option, the
ordering property (a CpG-dense sequence profiles higher than a CpG-poor one), and the
input validation for the dinucleotide and the window. It is a metric only -- there is
no ``ObjectiveTerm`` surface to test here (that lives in ``dinucleotide.py``).
"""

from __future__ import annotations

import pytest

from bt4.objectives.dinuc_profile import cpg_profile, dinucleotide_profile


def _naive_profile(dna: str, dinuc: str, window: int) -> list[float]:
    """Recount each window independently (the O(n*window) reference)."""
    seq = dna.upper()
    n_windows = len(seq) - window + 1
    if n_windows <= 0:
        return []
    out: list[float] = []
    for start in range(n_windows):
        chunk = seq[start : start + window]
        out.append(float(sum(1 for i in range(len(chunk) - 1) if chunk[i : i + 2] == dinuc)))
    return out


# ---------------------------------------------------------------------------
# Length convention
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("window", [1, 2, 3, 10])
def test_profile_length_convention(window: int) -> None:
    dna = "ACGTACGTACGT"  # length 12
    profile = dinucleotide_profile(dna, "CG", window=window)
    assert len(profile) == len(dna) - window + 1


def test_profile_empty_when_shorter_than_window() -> None:
    assert dinucleotide_profile("ACGTA", "CG", window=50) == []
    # Exactly one nucleotide short of a full window is still empty.
    assert dinucleotide_profile("ACGT", "CG", window=5) == []


def test_profile_single_window_when_exactly_window_long() -> None:
    dna = "CGCG"  # length 4 == window
    profile = dinucleotide_profile(dna, "CG", window=4)
    assert len(profile) == 1


# ---------------------------------------------------------------------------
# Exact hand-counted profile
# ---------------------------------------------------------------------------


def test_hand_counted_profile() -> None:
    # seq = C G T T G C G A T  (length 9), window = 4 -> 3 CG slots per window,
    # 6 sliding windows. CG occurs at index 0 (CG..) and index 5 (..CG..).
    #   window 0 "CGTT": slots at 0,1,2 -> 1 (the CG at index 0)
    #   window 1 "GTTG": slots at 1,2,3 -> 0
    #   window 2 "TTGC": slots at 2,3,4 -> 0
    #   window 3 "TGCG": slots at 3,4,5 -> 1 (the CG at index 5)
    #   window 4 "GCGA": slots at 4,5,6 -> 1
    #   window 5 "CGAT": slots at 5,6,7 -> 1
    dna = "CGTTGCGAT"
    profile = dinucleotide_profile(dna, "CG", window=4)
    assert profile == [1.0, 0.0, 0.0, 1.0, 1.0, 1.0]


def test_matches_naive_reference() -> None:
    dna = "CGATCGCGTTACGGCATGCGCGAT"
    for window in (2, 3, 7, 12):
        assert dinucleotide_profile(dna, "CG", window=window) == _naive_profile(
            dna, "CG", window
        )
        assert dinucleotide_profile(dna, "TA", window=window) == _naive_profile(
            dna, "TA", window
        )


def test_case_insensitive_and_wrapper_match() -> None:
    dna = "cgttgcgat"
    assert dinucleotide_profile(dna, "cg", window=4) == dinucleotide_profile(
        dna.upper(), "CG", window=4
    )
    # cpg_profile is exactly the CG specialization.
    assert cpg_profile(dna, window=4) == dinucleotide_profile(dna, "CG", window=4)


def test_window_one_is_all_zero() -> None:
    # A one-nucleotide window has no adjacent-base slots, so every count is 0.
    profile = dinucleotide_profile("CGCGCG", "CG", window=1)
    assert profile == [0.0] * 6


# ---------------------------------------------------------------------------
# Density option
# ---------------------------------------------------------------------------


def test_density_is_count_over_slots() -> None:
    dna = "CGCGCGCGCG"  # length 10; every adjacent even slot is CG
    counts = dinucleotide_profile(dna, "CG", window=5)
    densities = dinucleotide_profile(dna, "CG", window=5, density=True)
    slots = 5 - 1
    assert densities == pytest.approx([c / slots for c in counts])
    for value in densities:
        assert 0.0 <= value <= 1.0


def test_density_of_cg_repeat_alternates() -> None:
    # Over a pure "CGCG..." repeat the 4-nt windows alternate: an even-start
    # "CGCG" has 2 of its 3 slots CG (density 2/3), an odd-start "GCGC" has 1
    # (density 1/3). Overlapping CG can never fill every slot, so 2/3 is the max.
    dna = "CG" * 25
    densities = dinucleotide_profile(dna, "CG", window=4, density=True)
    assert max(densities) == pytest.approx(2.0 / 3.0)
    assert min(densities) == pytest.approx(1.0 / 3.0)


def test_density_window_one_is_zero_not_division_error() -> None:
    assert dinucleotide_profile("ACGT", "CG", window=1, density=True) == [0.0] * 4


# ---------------------------------------------------------------------------
# Ordering: CpG-dense profiles higher than CpG-poor
# ---------------------------------------------------------------------------


def test_dense_profiles_higher_than_poor() -> None:
    window = 20
    dense = "CG" * 40  # CpG at every other position
    poor = "AT" * 40  # no CpG at all
    dense_profile = cpg_profile(dense, window=window)
    poor_profile = cpg_profile(poor, window=window)
    assert dense_profile and poor_profile
    assert min(dense_profile) > max(poor_profile)
    assert all(value == 0.0 for value in poor_profile)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["C", "CGT", "", "CX", "XY", "N"])
def test_invalid_dinucleotide_rejected(bad: str) -> None:
    with pytest.raises(ValueError):
        dinucleotide_profile("ACGTACGT", bad, window=4)


@pytest.mark.parametrize("bad", [0, -1, -50])
def test_non_positive_window_rejected(bad: int) -> None:
    with pytest.raises(ValueError, match="window"):
        dinucleotide_profile("ACGTACGT", "CG", window=bad)


def test_empty_dna_rejected() -> None:
    with pytest.raises(ValueError):
        dinucleotide_profile("", "CG", window=4)


def test_non_acgt_dna_rejected() -> None:
    with pytest.raises(ValueError):
        dinucleotide_profile("ACGTN", "CG", window=4)


def test_cpg_profile_rejects_bad_window() -> None:
    with pytest.raises(ValueError, match="window"):
        cpg_profile("ACGTACGT", window=0)
