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
from bt4.domain.genetic_code import AMINO_ACIDS, CODON_TABLE, STOP, synonymous_codons
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


def test_minmax_term_is_scale_invariant() -> None:
    """The same table on a different scale must produce identical deltas.

    Load-bearing, not cosmetic: BT4's bundled tables are NOT all on one scale --
    the recounted organisms hold raw codon counts (~1e5) while the older
    hand-curated ones hold per-thousand values (~1e1) -- and ``bt4 build-table``
    hands users raw counts too. Without normalizing within each synonymous
    family, the same ``minmax_weight`` would mean ~10,000x more on one organism
    than another and silently swamp every other objective (CLAUDE.md §10.5).
    """
    base = load_table("homo_sapiens").frequency
    for factor in (1e-3, 7.0, 1e4):
        scaled = {codon: value * factor for codon, value in base.items()}
        for direction in ("max", "min"):
            plain = MinMaxTerm(base, direction)
            rescaled = MinMaxTerm(scaled, direction)
            for codon in base:
                assert plain.delta("", codon, 0) == pytest.approx(
                    rescaled.delta("", codon, 0), abs=1e-12
                )


def test_minmax_delta_magnitudes_are_comparable_across_organisms() -> None:
    """A raw-count organism and a per-thousand organism must weigh alike.

    This is the property a user relies on when they change the organism dropdown
    and expect ``minmax_weight`` to keep meaning the same thing.
    """
    scales = []
    for organism in ("homo_sapiens", "mus_musculus", "danio_rerio", "escherichia_coli"):
        term = MinMaxTerm(load_table(organism).frequency, "max")
        deltas = [abs(term.delta("", c, 0)) for c in load_table(organism).frequency]
        nonzero = [d for d in deltas if d > 0.0]
        scales.append(sum(nonzero) / len(nonzero))
    # Every organism's mean |delta| is a within-family fraction, so they sit in
    # the same order of magnitude rather than four apart.
    assert max(scales) / min(scales) < 3.0


def test_minmax_preserves_within_family_preference_order() -> None:
    """Normalizing must not reorder codons inside a synonymous family.

    The term still answers "which synonymous codon is more common"; only the
    units change. A reordering here would silently alter what the optimizer picks.
    """
    freq = load_table("homo_sapiens").frequency
    term = MinMaxTerm(freq, "max")
    for amino_acid in ("L", "R", "S", "A", "V"):
        synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
        by_raw = sorted(synonyms, key=lambda c: freq[c])
        by_delta = sorted(synonyms, key=lambda c: term.delta("", c, 0))
        assert by_raw == by_delta
