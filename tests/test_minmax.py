"""Tests for %MinMax: the additive DP term and the windowed reporting metric.

The load-bearing property is invariant #4, ``delta == score``: the running sum
of :class:`~bt4.objectives.minmax.MinMaxTerm`'s per-codon ``delta`` (as the DP
would accumulate it over real growing prefixes) must equal its whole-sequence
``score``. We also pin the contract surface (``name``/``scope``/``context_len``),
the zeroing of non-degenerate and stop codons, the ``direction`` orientation, the
``direction`` validation, and the reporting-metric bounds and sign.
"""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4.biomodels.codon.tables import load_table
from bt4.domain.genetic_code import AMINO_ACIDS, STOP, synonymous_codons
from bt4.domain.scope import Scope
from bt4.objectives import iter_codons
from bt4.objectives.minmax import MinMaxTerm, min_max_profile

_TABLE = load_table("homo_sapiens")
_FREQ = _TABLE.frequency
_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=120)

# Degenerate residues (each has a clear most/least common synonymous codon).
_DEGENERATE: tuple[str, ...] = tuple(
    sorted(aa for aa in AMINO_ACIDS if aa not in {"M", "W"})
)


def _backtranslate(protein: str) -> str:
    """Deterministically pick the first synonymous codon per residue."""
    return "".join(synonymous_codons(aa)[0] for aa in protein)


def _most_common_codon(aa: str) -> str:
    """Return the synonymous codon of ``aa`` with the largest table frequency."""
    return max(synonymous_codons(aa), key=lambda c: _FREQ[c])


def _least_common_codon(aa: str) -> str:
    """Return the synonymous codon of ``aa`` with the smallest table frequency."""
    return min(synonymous_codons(aa), key=lambda c: _FREQ[c])


def _accumulated_delta(term: MinMaxTerm, dna: str) -> float:
    """Sum the term's deltas over growing real prefixes (as the DP would)."""
    acc = 0.0
    prefix = ""
    for pos, codon in iter_codons(dna):
        acc += term.delta(prefix, codon, pos)
        prefix += codon
    return acc


# ---------------------------------------------------------------------------
# MinMaxTerm: contract surface
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("direction", "expected"),
    [("max", "minmax_max"), ("min", "minmax_min")],
)
def test_name_scope_and_context(direction: str, expected: str) -> None:
    term = MinMaxTerm(_FREQ, direction=direction)
    assert term.name == expected
    assert term.scope() is Scope.LOCAL
    assert term.context_len() == 0


def test_bad_direction_rejected() -> None:
    with pytest.raises(ValueError, match="direction"):
        MinMaxTerm(_FREQ, direction="both")


def test_default_direction_is_max() -> None:
    assert MinMaxTerm(_FREQ).direction == "max"


def test_frozen_term_is_immutable() -> None:
    term = MinMaxTerm(_FREQ)
    with pytest.raises((AttributeError, TypeError)):
        term.direction = "min"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# MinMaxTerm: invariant #4 (delta == score) -- the load-bearing test
# ---------------------------------------------------------------------------


@given(protein=_PROTEIN, direction=st.sampled_from(["max", "min"]))
def test_delta_equals_score(protein: str, direction: str) -> None:
    term = MinMaxTerm(_FREQ, direction=direction)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


# ---------------------------------------------------------------------------
# MinMaxTerm: orientation and non-degenerate handling
# ---------------------------------------------------------------------------


def test_direction_flips_sign_for_common_codon() -> None:
    # Leucine's most-common codon is (well) above its synonymous average.
    codon = _most_common_codon("L")
    up = MinMaxTerm(_FREQ, direction="max").delta("", codon, 0)
    down = MinMaxTerm(_FREQ, direction="min").delta("", codon, 0)
    assert up > 0.0
    assert down < 0.0
    assert up == pytest.approx(-down)


def test_rare_codon_is_negative_under_max() -> None:
    codon = _least_common_codon("L")
    assert MinMaxTerm(_FREQ, direction="max").delta("", codon, 0) < 0.0


@pytest.mark.parametrize("direction", ["max", "min"])
def test_non_degenerate_and_stop_codons_score_zero(direction: str) -> None:
    term = MinMaxTerm(_FREQ, direction=direction)
    # Met (ATG) and Trp (TGG) carry no coding choice.
    assert term.delta("", "ATG", 0) == 0.0
    assert term.delta("", "TGG", 0) == 0.0
    for stop in synonymous_codons(STOP):
        assert term.delta("", stop, 0) == 0.0


def test_delta_is_case_insensitive() -> None:
    term = MinMaxTerm(_FREQ)
    codon = _most_common_codon("R")
    assert term.delta("", codon.lower(), 3) == pytest.approx(term.delta("", codon, 3))


# ---------------------------------------------------------------------------
# min_max_profile: reporting metric
# ---------------------------------------------------------------------------


def test_profile_length_convention() -> None:
    protein = ("LRSVAGTP" * 4)[:30]  # 30 degenerate residues
    dna = "".join(_most_common_codon(aa) for aa in protein)
    window = 18
    profile = min_max_profile(dna, _FREQ, window=window)
    assert len(profile) == 30 - window + 1


def test_profile_empty_when_shorter_than_window() -> None:
    dna = "".join(_most_common_codon(aa) for aa in "LRS")  # 3 codons
    assert min_max_profile(dna, _FREQ, window=18) == []


def test_all_most_common_is_maximally_positive() -> None:
    protein = "LRSVAGTP" * 3  # 24 degenerate residues, > one window
    dna = "".join(_most_common_codon(aa) for aa in protein)
    profile = min_max_profile(dna, _FREQ, window=18)
    assert profile  # non-empty
    for value in profile:
        assert value == pytest.approx(100.0)


def test_all_least_common_is_maximally_negative() -> None:
    protein = "LRSVAGTP" * 3
    dna = "".join(_least_common_codon(aa) for aa in protein)
    profile = min_max_profile(dna, _FREQ, window=18)
    assert profile
    for value in profile:
        assert value == pytest.approx(-100.0)


@given(protein=st.text(alphabet=sorted(AMINO_ACIDS), min_size=18, max_size=90))
def test_profile_values_within_bounds(protein: str) -> None:
    dna = _backtranslate(protein)
    profile = min_max_profile(dna, _FREQ, window=18)
    assert len(profile) == len(protein) - 18 + 1
    for value in profile:
        assert -100.0 <= value <= 100.0


def test_single_codon_window_reports_zero() -> None:
    # A window of only Met/Trp has no synonymous choice: Xmax == Xavg == Xmin.
    dna = "ATG" * 10 + "TGG" * 10  # 20 single-codon residues
    profile = min_max_profile(dna, _FREQ, window=18)
    assert profile
    for value in profile:
        assert value == 0.0


def test_profile_rejects_bad_window() -> None:
    dna = _backtranslate("LRSVAGTP")
    with pytest.raises(ValueError, match="window"):
        min_max_profile(dna, _FREQ, window=0)


def test_profile_rejects_missing_codon() -> None:
    sparse = {codon: _FREQ[codon] for codon in synonymous_codons("A")}
    with pytest.raises(ValueError, match="missing"):
        min_max_profile(_backtranslate("LR"), sparse, window=1)
