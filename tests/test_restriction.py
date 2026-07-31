"""Tests for the restriction-site constraint and the IUPAC matcher.

Two things are proven here. First, the IUPAC engine matches degenerate patterns
correctly and takes a strand-correct reverse complement (palindromic sites map
to themselves). Second, the restriction constraint obeys invariant #3: a
sequence built respecting ``ok_suffix`` has zero hard violations under
``validate``, including sites that straddle a codon boundary and sites on the
reverse strand.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.iupac import find_iupac, matches_at, reverse_complement_iupac
from bt4.constraints.restriction import (
    ENZYMES,
    RestrictionSiteConstraint,
    available_enzymes,
)
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)


def _build_feasible(protein: str, constraint: RestrictionSiteConstraint) -> str:
    """Greedily first-fit a codon per residue that passes ``ok_suffix``.

    Returns the DNA built so far; stops early (returns the valid prefix) at the
    first residue where no synonymous codon is feasible.
    """
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if constraint.ok_suffix(dna, codon):
                dna += codon
                break
        else:  # dead end - stop with the feasible prefix built so far.
            break
    return dna


def _hard(constraint: RestrictionSiteConstraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# IUPAC matcher: concrete and degenerate matching.
# --------------------------------------------------------------------------- #


def test_matches_at_concrete() -> None:
    assert matches_at("GAATTC", "GAATTC", 0) is True
    assert matches_at("GAATTC", "GAATTC", 1) is False
    assert matches_at("AGAATTCA", "GAATTC", 1) is True


def test_matches_at_degenerate_gantc() -> None:
    # GANTC (N = any base) matches each concrete filling of the third position.
    for third in "ACGT":
        assert matches_at(f"GA{third}TC", "GANTC", 0) is True
    # A base that violates a fixed position does not match.
    assert matches_at("GGATC", "GANTC", 0) is False  # pos 1 must be A


def test_matches_at_out_of_bounds_is_false() -> None:
    assert matches_at("GAAT", "GANTC", 0) is False  # pattern runs past the window
    assert matches_at("GAATTC", "GAATTC", -1) is False


def test_find_iupac_concrete_and_overlapping() -> None:
    assert find_iupac("AAGAATTCAA", "GAATTC") == [2]
    # Overlapping occurrences are all reported.
    assert find_iupac("AAAA", "AA") == [0, 1, 2]
    assert find_iupac("CCCC", "GAATTC") == []


def test_find_iupac_degenerate() -> None:
    # GANTC occurs once for each of the concrete third-position bases.
    assert find_iupac("TTGACTCTT", "GANTC") == [2]
    assert find_iupac("GAATCGAGTC", "GANTC") == [0, 5]


def test_find_iupac_invalid_pattern_raises() -> None:
    with pytest.raises(ValueError, match="IUPAC"):
        find_iupac("ACGT", "GAXTC")
    with pytest.raises(ValueError, match="IUPAC"):
        find_iupac("ACGT", "")


def test_reverse_complement_iupac_palindromes() -> None:
    # EcoRI is a classic palindrome; it maps back onto itself.
    assert reverse_complement_iupac("GAATTC") == "GAATTC"
    # GANTC is palindromic through the ambiguous N.
    assert reverse_complement_iupac("GANTC") == "GANTC"
    # DraIII CACNNNGTG is palindromic across the NNN core.
    assert reverse_complement_iupac("CACNNNGTG") == "CACNNNGTG"


def test_reverse_complement_iupac_asymmetric() -> None:
    assert reverse_complement_iupac("GGGAAA") == "TTTCCC"
    assert reverse_complement_iupac("TTTCCC") == "GGGAAA"
    # K = {G,T} complements to M = {A,C}.
    assert reverse_complement_iupac("K") == "M"


# --------------------------------------------------------------------------- #
# Catalog surface.
# --------------------------------------------------------------------------- #


def test_available_enzymes_nonempty_and_sorted() -> None:
    names = available_enzymes()
    assert names
    assert list(names) == sorted(names)
    assert "EcoRI" in names
    assert "HinfI" in names  # degenerate site is catalogued
    assert set(names) == set(ENZYMES)


def test_unknown_enzyme_raises() -> None:
    with pytest.raises(ValueError, match="unknown enzyme"):
        RestrictionSiteConstraint(enzymes=("NotAnEnzyme",))


def test_invalid_extra_site_raises() -> None:
    with pytest.raises(ValueError, match="IUPAC"):
        RestrictionSiteConstraint(extra_sites=("GAXTC",))


def test_name_scope_penalty() -> None:
    constraint = RestrictionSiteConstraint(enzymes=("EcoRI",))
    assert constraint.name == "restriction_site"
    assert constraint.scope() is Scope.LOCAL
    assert constraint.penalty("AAA", "CCC") == 0.0


# --------------------------------------------------------------------------- #
# Positive detection: forward strand and reverse strand.
# --------------------------------------------------------------------------- #


def test_validate_flags_ecori_site() -> None:
    constraint = RestrictionSiteConstraint(enzymes=("EcoRI",))
    violations = list(constraint.validate("AAAGAATTCAAA"))
    assert len(violations) == 1
    v = violations[0]
    assert v.constraint == "restriction_site"
    assert v.severity is Severity.HARD
    assert (v.start, v.end) == (3, 9)


def test_validate_catches_reverse_complement_strand() -> None:
    # An asymmetric site: GGGAAA (RC = TTTCCC) must be banned on both strands.
    constraint = RestrictionSiteConstraint(extra_sites=("GGGAAA",))
    fwd = list(constraint.validate("AAAGGGAAACCC"))
    assert [v.constraint for v in fwd] == ["restriction_site"]
    assert (fwd[0].start, fwd[0].end) == (3, 9)
    # The reverse-complement occurrence is caught even though only GGGAAA was given.
    rc = list(constraint.validate("AAATTTCCCAAA"))
    assert [v.constraint for v in rc] == ["restriction_site"]
    assert (rc[0].start, rc[0].end) == (3, 9)


def test_validate_flags_degenerate_hinfi_site() -> None:
    constraint = RestrictionSiteConstraint(enzymes=("HinfI",))
    # GACTC is a concrete filling of HinfI's GANTC.
    starts = sorted(v.start for v in constraint.validate("AAGACTCAA"))
    assert starts == [2]


def test_validate_empty_sequence_is_clean() -> None:
    constraint = RestrictionSiteConstraint(enzymes=("EcoRI",))
    assert list(constraint.validate("")) == []


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing sites are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_ok_suffix_vetoes_boundary_crossing_site() -> None:
    constraint = RestrictionSiteConstraint(enzymes=("EcoRI",))
    # context_len is exactly maxlen - 1 == 5 for a single 6-nt site.
    assert constraint.context_len() == 5
    # prefix ends in GAATT; the codon 'CAA' completes GAATTC across the seam.
    assert constraint.ok_suffix("GGGGAATT", "CAA") is False
    # A codon that does not complete the site is allowed.
    assert constraint.ok_suffix("GGGGAATT", "GAA") is True


def test_ok_suffix_allows_site_fully_in_prefix() -> None:
    # A site entirely inside the already-feasible prefix, not touching the new
    # codon, is not this extension's fault - ok_suffix only judges next_codon.
    constraint = RestrictionSiteConstraint(enzymes=("EcoRI",))
    assert constraint.ok_suffix("GAATTC", "AAA") is True


def test_ok_suffix_vetoes_reverse_strand_boundary_crossing() -> None:
    constraint = RestrictionSiteConstraint(extra_sites=("GGGAAA",))
    # RC site TTTCCC straddles the seam: prefix tail 'TTT' + codon 'CCC'.
    assert constraint.ok_suffix("AAATTT", "CCC") is False
    assert constraint.ok_suffix("AAATTT", "GGG") is True


def test_empty_constraint_is_inert() -> None:
    constraint = RestrictionSiteConstraint()
    assert constraint.context_len() == 0
    assert constraint.ok_suffix("", "ATG") is True
    assert constraint.ok_suffix("GAATTC", "GAA") is True
    assert list(constraint.validate("GAATTCGGATCC")) == []


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN)
def test_ok_suffix_implies_validate_clean(protein: str) -> None:
    constraint = RestrictionSiteConstraint(
        enzymes=("EcoRI", "BamHI", "HindIII", "HinfI", "XhoI"),
    )
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []


@given(protein=_PROTEIN)
def test_ok_suffix_implies_validate_clean_with_extra_sites(protein: str) -> None:
    constraint = RestrictionSiteConstraint(
        enzymes=("EcoRV", "DraIII"),
        extra_sites=("GGGAAA",),
    )
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
