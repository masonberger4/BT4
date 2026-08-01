"""Tests for the internal strong-Kozak ATG constraint (CLAUDE.md invariant #3).

The load-bearing property: a sequence built respecting ``ok_suffix`` contains
zero hard violations under ``validate`` (``ok_suffix <=> validate``), and the
declared ``context_len`` actually suffices for the veto -- so a strong-Kozak ATG
whose +4 base is completed across a codon seam is still caught. This mirrors the
methodology of ``tests/test_constraints.py``.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.kozak import InternalStartConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=80)


def _build_feasible(protein: str, constraint: InternalStartConstraint) -> str:
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
        else:  # no synonymous codon worked - dead end, stop with the valid prefix.
            break
    return dna


def _hard(constraint: InternalStartConstraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN)
def test_ok_suffix_implies_validate_clean_default(protein: str) -> None:
    constraint = InternalStartConstraint()  # strong = purine@-3 AND G@+4
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


@given(
    protein=_PROTEIN,
    require_purine=st.booleans(),
    require_g=st.booleans(),
    min_start=st.integers(min_value=0, max_value=6),
)
def test_ok_suffix_implies_validate_clean_all_configs(
    protein: str, require_purine: bool, require_g: bool, min_start: int
) -> None:
    constraint = InternalStartConstraint(
        require_purine_m3=require_purine,
        require_g_p4=require_g,
        min_start=min_start,
    )
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []


# --------------------------------------------------------------------------- #
# Positive detection: validate flags strong internal ATGs; weak ones are clean.
# --------------------------------------------------------------------------- #


def test_validate_flags_strong_internal_atg() -> None:
    constraint = InternalStartConstraint()
    # AAA ATG G: ATG at index 3, purine A at -3 (index 0), G at +4 (index 6).
    violations = list(constraint.validate("AAAATGG"))
    assert len(violations) == 1
    v = violations[0]
    assert v.constraint == "internal_start"
    assert v.severity is Severity.HARD
    assert (v.start, v.end) == (3, 6)


def test_validate_ignores_weak_context() -> None:
    constraint = InternalStartConstraint()
    # Pyrimidine C at -3 -> not a strong context under the default.
    assert list(constraint.validate("CCCATGG")) == []
    # Purine at -3 but +4 is C (not G) -> also weak.
    assert list(constraint.validate("AAAATGC")) == []


# --------------------------------------------------------------------------- #
# context_len suffices: seam-crossing +4 base is vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_context_len_is_six() -> None:
    assert InternalStartConstraint().context_len() == 6


def test_ok_suffix_vetoes_seam_crossing_plus4() -> None:
    constraint = InternalStartConstraint()
    # prefix ends AAG ATG (purine at -3, ATG at index 3); a codon starting with G
    # completes the strong +4 across the seam -> vetoed.
    assert constraint.ok_suffix("AAGATG", "GCC") is False
    # A codon that does not put G at +4 leaves the context weak -> allowed.
    assert constraint.ok_suffix("AAGATG", "CCC") is True


def test_ok_suffix_allows_occurrence_fully_in_prefix() -> None:
    # The strong ATG (index 3, +4 G at index 6) lies wholly inside the feasible
    # prefix; completing an unrelated codon is not this extension's fault.
    constraint = InternalStartConstraint()
    assert constraint.ok_suffix("AAGATGG", "CCC") is True


def test_ok_suffix_vetoes_atg_formed_across_seam() -> None:
    # prefix 'AGAAT' ends ...AT; codon 'GGC' supplies the ATG's G (index 5) and
    # the +4 G (index 6). The ATG (A@3,T@4,G@5) straddles the seam, purine A sits
    # at -3 (index 0) -> a new strong internal ATG is vetoed.
    constraint = InternalStartConstraint()
    assert constraint.ok_suffix("AGAAT", "GGC") is False


# --------------------------------------------------------------------------- #
# min_start: the real start codon is exempt; the same pattern internal is not.
# --------------------------------------------------------------------------- #


def test_min_start_exempts_index_zero() -> None:
    # Isolate min_start: -3 is not required, so only the A-index distinguishes.
    constraint = InternalStartConstraint(require_purine_m3=False, require_g_p4=True)
    assert list(constraint.validate("ATGG")) == []  # ATG at index 0: real start.
    flagged = list(constraint.validate("CCCATGG"))  # same +4=G pattern, index 3.
    assert [v.start for v in flagged] == [3]


def test_ok_suffix_never_vetoes_first_codon() -> None:
    # An empty prefix means next_codon is codon 0; no internal ATG is possible.
    constraint = InternalStartConstraint(require_purine_m3=False, require_g_p4=False)
    assert constraint.ok_suffix("", "ATG") is True


# --------------------------------------------------------------------------- #
# Configuration: toggling each condition behaves as documented.
# --------------------------------------------------------------------------- #


def test_both_conditions_off_forbids_every_internal_atg() -> None:
    constraint = InternalStartConstraint(require_purine_m3=False, require_g_p4=False)
    # Weak context (C at -3, C at +4) is still forbidden when context is ignored.
    assert [v.start for v in constraint.validate("CCCATGCCC")] == [3]
    # min_start still exempts the start codon.
    assert list(constraint.validate("ATGCCC")) == []
    # ok_suffix vetoes any new internal ATG regardless of flanks.
    assert constraint.ok_suffix("CCC", "ATG") is False


def test_require_purine_only() -> None:
    constraint = InternalStartConstraint(require_purine_m3=True, require_g_p4=False)
    # Purine at -3, +4 irrelevant -> forbidden.
    assert [v.start for v in constraint.validate("AAAATGC")] == [3]
    # Pyrimidine at -3 -> allowed.
    assert list(constraint.validate("CCCATGC")) == []


def test_require_g_p4_only() -> None:
    constraint = InternalStartConstraint(require_purine_m3=False, require_g_p4=True)
    # G at +4, -3 irrelevant -> forbidden.
    assert [v.start for v in constraint.validate("CCCATGG")] == [3]
    # Non-G at +4 -> allowed.
    assert list(constraint.validate("CCCATGC")) == []


# --------------------------------------------------------------------------- #
# Edge cases: absent -3 / +4 flanks fail their condition (not-strong => allowed).
# --------------------------------------------------------------------------- #


def test_missing_minus3_base_is_not_a_purine() -> None:
    # min_start=0 makes the very-start ATG internal, but with no -3 base the
    # purine condition is unmet -> not forbidden.
    constraint = InternalStartConstraint(
        require_purine_m3=True, require_g_p4=False, min_start=0
    )
    assert list(constraint.validate("ATGCCC")) == []
    # An ATG far enough in to have a purine -3 base IS forbidden.
    assert [v.start for v in constraint.validate("AAAATG")] == [3]


def test_missing_plus4_base_is_not_a_g() -> None:
    # ATG at the very 3' end has no +4 base -> the G@+4 condition is unmet.
    constraint = InternalStartConstraint(require_purine_m3=False, require_g_p4=True)
    assert list(constraint.validate("CCCATG")) == []
    # The same ATG with a real G at +4 IS forbidden.
    assert [v.start for v in constraint.validate("CCCATGG")] == [3]


# --------------------------------------------------------------------------- #
# Identity / configuration validation.
# --------------------------------------------------------------------------- #


def test_name_scope_and_penalty() -> None:
    constraint = InternalStartConstraint()
    assert constraint.name == "internal_start"
    assert constraint.scope() is Scope.LOCAL
    assert constraint.penalty("AAA", "ATG") == 0.0


def test_rejects_negative_min_start() -> None:
    with pytest.raises(ValueError, match="min_start"):
        InternalStartConstraint(min_start=-1)
