"""Tests for the 5' translation-ramp shaping term (CLAUDE.md invariant #4).

The load-bearing property is ``delta == score``: the running sum of the term's
per-codon ``delta`` (as the DP would accumulate it, over real growing prefixes)
must equal the whole-sequence ``score``. We also pin the contract surface
(``name``/``scope``/``context_len``), the zeroing of non-degenerate and stop
codons, the positional decay shape (stronger early, exactly zero at and beyond
``ramp_codons``), and the ``ramp_codons`` validation.
"""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4.biomodels.codon.tables import load_table
from bt4.domain.genetic_code import AMINO_ACIDS, STOP, synonymous_codons
from bt4.domain.scope import Scope
from bt4.objectives import iter_codons
from bt4.objectives.ramp import RampTerm

_TABLE = load_table("homo_sapiens")
_W = _TABLE.relative_adaptiveness()
_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=120)


def _backtranslate(protein: str) -> str:
    """Deterministically pick the first synonymous codon per residue."""
    return "".join(synonymous_codons(aa)[0] for aa in protein)


def _accumulated_delta(term: RampTerm, dna: str) -> float:
    """Sum the term's deltas over growing real prefixes (as the DP would)."""
    acc = 0.0
    prefix = ""
    for pos, codon in iter_codons(dna):
        acc += term.delta(prefix, codon, pos)
        prefix += codon
    return acc


def test_name_scope_and_context() -> None:
    term = RampTerm(_W)
    assert term.name == "ramp"
    assert term.scope() is Scope.POSITIONAL
    assert term.context_len() == 0


@given(protein=_PROTEIN, ramp_codons=st.integers(min_value=1, max_value=60))
def test_delta_equals_score(protein: str, ramp_codons: int) -> None:
    term = RampTerm(_W, ramp_codons=ramp_codons)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


@pytest.mark.parametrize("pos", [0, 1, 10, 34, 35, 100])
def test_non_degenerate_and_stop_codons_score_zero(pos: int) -> None:
    term = RampTerm(_W)
    # Met (ATG), Trp (TGG) and every stop codon carry no coding choice.
    assert term.delta("", "ATG", pos) == 0.0
    assert term.delta("", "TGG", pos) == 0.0
    for stop in synonymous_codons(STOP):
        assert term.delta("", stop, pos) == 0.0


def test_ramp_decays_and_is_zero_past_ramp() -> None:
    ramp_codons = 35
    term = RampTerm(_W, ramp_codons=ramp_codons)
    # CTG (Leu) is degenerate and non-degenerate-excluded, so it carries a signal.
    codon = "CTG"
    early = abs(term.delta("", codon, 0))
    later = abs(term.delta("", codon, 20))
    assert early > later > 0.0
    # The ramp weight hits zero exactly at ramp_codons and stays there beyond.
    assert term.delta("", codon, ramp_codons) == 0.0
    assert term.delta("", codon, ramp_codons + 5) == 0.0


def test_delta_at_pos_zero_equals_negative_weight() -> None:
    term = RampTerm(_W)
    codon = "CTG"
    assert term.delta("", codon, 0) == pytest.approx(-_TABLE.weight(codon))


def test_delta_is_case_insensitive() -> None:
    term = RampTerm(_W)
    assert term.delta("", "ctg", 3) == pytest.approx(term.delta("", "CTG", 3))


@pytest.mark.parametrize("bad", [0, -1, -35])
def test_ramp_codons_below_one_rejected(bad: int) -> None:
    with pytest.raises(ValueError, match="ramp_codons"):
        RampTerm(_W, ramp_codons=bad)


def test_frozen_term_is_immutable() -> None:
    term = RampTerm(_W)
    with pytest.raises((AttributeError, TypeError)):
        term.name = "other"  # type: ignore[misc]
    with pytest.raises((AttributeError, TypeError)):
        term.ramp_codons = 10  # type: ignore[misc]
