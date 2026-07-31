"""Tests for the dinucleotide-content objective term (CLAUDE.md invariant #4).

The load-bearing property is ``delta == score``: the running sum of the term's
per-codon ``delta`` (as the DP would accumulate it, over real growing prefixes)
must equal the whole-sequence ``score``, for both depletion and elevation and
across dinucleotides that do and do not straddle codon boundaries. We also pin
exact counting against a hand-built sequence (including a boundary-straddling
occurrence), the depletion ordering (CpG-rich scores lower than CpG-poor), the
contract surface (``name``/``scope``/``context_len``/``sign``), and the input
validation for the dinucleotide and the direction.
"""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons
from bt4.domain.scope import Scope
from bt4.objectives import iter_codons
from bt4.objectives.dinucleotide import DinucleotideTerm

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=120)


def _backtranslate(protein: str) -> str:
    """Deterministically pick the first synonymous codon per residue."""
    return "".join(synonymous_codons(aa)[0] for aa in protein)


def _accumulated_delta(term: DinucleotideTerm, dna: str) -> float:
    """Sum the term's deltas over growing real prefixes (as the DP would)."""
    acc = 0.0
    prefix = ""
    for pos, codon in iter_codons(dna):
        acc += term.delta(prefix, codon, pos)
        prefix += codon
    return acc


def test_name_scope_and_context() -> None:
    term = DinucleotideTerm("CG")
    assert term.name == "dinuc_cg_deplete"
    assert term.scope() is Scope.PAIRWISE
    assert term.context_len() == 1
    assert term.sign == -1.0

    elevate = DinucleotideTerm("TA", direction="elevate")
    assert elevate.name == "dinuc_ta_elevate"
    assert elevate.sign == 1.0


def test_dinucleotide_is_normalized_to_upper() -> None:
    term = DinucleotideTerm("cg")
    assert term.dinucleotide == "CG"
    assert term.name == "dinuc_cg_deplete"


@pytest.mark.parametrize("direction", ["deplete", "elevate"])
@pytest.mark.parametrize("dinuc", ["CG", "TA"])
@given(protein=_PROTEIN)
def test_delta_equals_score(dinuc: str, direction: str, protein: str) -> None:
    term = DinucleotideTerm(dinuc, direction=direction)
    dna = _backtranslate(protein)
    assert term.score(dna) == pytest.approx(_accumulated_delta(term, dna))


def test_counting_including_straddle() -> None:
    # CGTTGCGAT: one CG inside codon 0 (CGT) and one straddling codon 1|2 (TGC|GAT).
    #   index:  0 1 2 3 4 5 6 7 8
    #   base:   C G T T G C G A T
    #   CG at:  ^         ^
    dna = "CGTTGCGAT"
    deplete = DinucleotideTerm("CG", direction="deplete")
    elevate = DinucleotideTerm("CG", direction="elevate")
    assert deplete.score(dna) == -2.0
    assert elevate.score(dna) == 2.0
    # The straddling occurrence is only visible with the 1-base prefix context.
    assert deplete.score(dna) == pytest.approx(_accumulated_delta(deplete, dna))


def test_single_straddling_occurrence() -> None:
    # TGCGAT: the only CG spans the codon boundary between TGC and GAT.
    dna = "TGCGAT"
    term = DinucleotideTerm("CG", direction="deplete")
    assert term.score(dna) == -1.0
    # It must be attributed to codon 1 (the codon holding the CG's end base G).
    assert term.delta("TGC", "GAT", 1) == -1.0
    assert term.delta("", "TGC", 0) == 0.0


def test_deplete_scores_lower_on_rich_than_poor() -> None:
    term = DinucleotideTerm("CG", direction="deplete")
    rich = "CGCGCGCGCG"  # five overlapping CG occurrences
    poor = "ATATATATAT"  # no CG at all
    assert term.score(rich) == -5.0
    assert term.score(poor) == 0.0
    assert term.score(rich) < term.score(poor)


def test_elevate_scores_higher_on_rich_than_poor() -> None:
    term = DinucleotideTerm("CG", direction="elevate")
    rich = "CGCGCGCGCG"
    poor = "ATATATATAT"
    assert term.score(rich) > term.score(poor)


@pytest.mark.parametrize("bad", ["C", "CGT", "", "CX", "XY", "CGA"])
def test_invalid_dinucleotide_rejected(bad: str) -> None:
    with pytest.raises(ValueError):
        DinucleotideTerm(bad)


@pytest.mark.parametrize("bad", ["both", "up", "DEPLETE", "", "raise"])
def test_invalid_direction_rejected(bad: str) -> None:
    with pytest.raises(ValueError, match="direction"):
        DinucleotideTerm("CG", direction=bad)


def test_frozen_term_is_immutable() -> None:
    term = DinucleotideTerm("CG")
    with pytest.raises((AttributeError, TypeError)):
        term.dinucleotide = "TA"  # type: ignore[misc]
    with pytest.raises((AttributeError, TypeError)):
        term.direction = "elevate"  # type: ignore[misc]
