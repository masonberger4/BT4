"""Tests for the strong splice-consensus motif constraint (CLAUDE.md invariant #3).

The load-bearing property, as for every LOCAL constraint: a sequence built
respecting ``ok_suffix`` has zero hard violations under ``validate``
(``ok_suffix <=> validate``), and ``context_len`` must actually suffice for a
motif that straddles a codon boundary. These also pin the honest design choices:
sense strand only (no reverse complement -- splicing is strand-specific), IUPAC
degeneracy, and that the bare ``GT``/``AG`` are never banned.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.base import Constraint
from bt4.constraints.splice_motif import (
    DEFAULT_ACCEPTOR_MOTIFS,
    DEFAULT_DONOR_MOTIFS,
    SpliceSiteMotifConstraint,
)
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)

_DONOR_ONLY = SpliceSiteMotifConstraint(donor_motifs=("GTRAGT",), acceptor_motifs=())
_ACCEPTOR_ONLY = SpliceSiteMotifConstraint(donor_motifs=(), acceptor_motifs=("YYYYYYNYAGG",))


def _build_feasible(protein: str, constraints: list[Constraint]) -> str:
    """Greedily first-fit a codon per residue that passes every ``ok_suffix``."""
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if all(c.ok_suffix(dna, codon) for c in constraints):
                dna += codon
                break
        else:  # dead end - stop with the valid prefix.
            break
    return dna


def _hard(constraint: Constraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN)
def test_splice_ok_suffix_implies_validate_clean(protein: str) -> None:
    constraint = SpliceSiteMotifConstraint()  # both donor + acceptor defaults
    dna = _build_feasible(protein, [constraint])
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # A genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


@given(protein=_PROTEIN)
def test_splice_donor_only_ok_suffix_implies_validate_clean(protein: str) -> None:
    dna = _build_feasible(protein, [_DONOR_ONLY])
    assume(len(dna) >= 3)
    assert _hard(_DONOR_ONLY, dna) == []


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing motifs are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_donor_context_len_is_maxlen_minus_one() -> None:
    assert _DONOR_ONLY.context_len() == 5  # len("GTRAGT") - 1
    assert _ACCEPTOR_ONLY.context_len() == 10  # len("YYYYYYNYAGG") - 1
    # With both arms active the context is the larger of the two.
    assert SpliceSiteMotifConstraint().context_len() == 10


def test_donor_ok_suffix_vetoes_boundary_crossing_motif() -> None:
    # prefix ends in GTAAG (matches GTRAG with R=A); the codon 'TCC' completes
    # GTAAGT across the seam.
    assert _DONOR_ONLY.ok_suffix("AAAGTAAG", "TCC") is False
    # A codon that does not complete the donor core is fine.
    assert _DONOR_ONLY.ok_suffix("AAAGTAAG", "CCC") is True


def test_ok_suffix_allows_motif_fully_in_prefix() -> None:
    # A motif entirely inside the already-feasible prefix, not touching the new
    # codon, is not this extension's fault - ok_suffix only judges next_codon.
    assert _DONOR_ONLY.ok_suffix("GTGAGT", "AAA") is True
    # But a motif overlapping the incoming codon is vetoed.
    assert _DONOR_ONLY.ok_suffix("AAAGTGAG", "TCC") is False


# --------------------------------------------------------------------------- #
# Positive detection + IUPAC degeneracy.
# --------------------------------------------------------------------------- #


def test_validate_flags_donor_both_purines_at_plus3() -> None:
    # R = A or G at donor +3: GTAAGT and GTGAGT are both strong donors.
    for core in ("GTAAGT", "GTGAGT"):
        vios = list(_DONOR_ONLY.validate("CCC" + core + "CCC"))
        assert [v.constraint for v in vios] == ["splice_site"]
        assert vios[0].severity is Severity.HARD
        assert (vios[0].start, vios[0].end) == (3, 9)
        assert "donor" in vios[0].detail
    # GTCAGT (pyrimidine at +3) is NOT the strong core and is not flagged.
    assert list(_DONOR_ONLY.validate("CCCGTCAGTCCC")) == []


def test_validate_flags_strong_acceptor() -> None:
    # YYYYYY N Y A G | G with a T/C tract, N=A, Y=C: TTTTTT A C AG G.
    vios = list(_ACCEPTOR_ONLY.validate("AA" + "TTTTTTACAGG" + "AA"))
    assert [v.constraint for v in vios] == ["splice_site"]
    assert "acceptor" in vios[0].detail
    assert (vios[0].start, vios[0].end) == (2, 13)


def test_bare_dinucleotides_are_never_banned() -> None:
    # The ubiquitous bare donor GT and acceptor AG must never be flagged.
    constraint = SpliceSiteMotifConstraint()
    assert list(constraint.validate("ATGGTAGCAGATGGTTAA")) == []


# --------------------------------------------------------------------------- #
# Sense strand only: no reverse-complement banning.
# --------------------------------------------------------------------------- #


def test_no_reverse_complement_banning() -> None:
    # reverse_complement("GTAAGT") == "ACTTAC"; splicing is strand-specific, so
    # the RC of a donor must NOT be flagged (unlike restriction sites).
    assert list(_DONOR_ONLY.validate("CCCACTTACCCC")) == []
    # The forward motif still is.
    assert len(list(_DONOR_ONLY.validate("CCCGTAAGTCCC"))) == 1


# --------------------------------------------------------------------------- #
# Degenerate / configuration cases.
# --------------------------------------------------------------------------- #


def test_empty_motifs_is_inert() -> None:
    constraint = SpliceSiteMotifConstraint(donor_motifs=(), acceptor_motifs=())
    assert constraint.context_len() == 0
    assert constraint.ok_suffix("", "ATG") is True
    assert constraint.ok_suffix("CCCGTAAGT", "CCC") is True
    assert list(constraint.validate("GTAAGTTTTTTTACAGG")) == []


def test_rejects_invalid_iupac_motif() -> None:
    with pytest.raises(ValueError, match=r"donor"):
        SpliceSiteMotifConstraint(donor_motifs=("GTXAGT",))
    with pytest.raises(ValueError, match=r"acceptor"):
        SpliceSiteMotifConstraint(acceptor_motifs=("",))


def test_name_scope_penalty_and_defaults() -> None:
    constraint = SpliceSiteMotifConstraint()
    assert constraint.name == "splice_site"
    assert constraint.scope() is Scope.LOCAL
    assert constraint.penalty("AAA", "CCC") == 0.0
    assert DEFAULT_DONOR_MOTIFS == ("GTRAGT",)
    assert DEFAULT_ACCEPTOR_MOTIFS == ("YYYYYYNYAGG",)
    assert isinstance(constraint, Constraint)


def test_duplicate_pattern_across_arms_not_double_counted() -> None:
    # If the same pattern is given to both arms, validate reports it once (donor
    # wins the label since it is collected first).
    constraint = SpliceSiteMotifConstraint(donor_motifs=("GTRAGT",), acceptor_motifs=("GTRAGT",))
    vios = list(constraint.validate("CCCGTAAGTCCC"))
    assert len(vios) == 1
    assert "donor" in vios[0].detail
