"""Tests for :class:`~bt4.constraints.uorf.UorfConstraint`.

The load-bearing property (CLAUDE.md invariant #3, "ok_suffix <=> validate"): a
sequence built respecting ``ok_suffix`` contains zero hard violations under
``validate``. This holds because ``UorfConstraint`` is honest about being GLOBAL
-- its ``ok_suffix`` reads the *whole* prefix (a uORF's ATG and its closing stop
can lie far apart), so a uORF whose stop completes on the incoming codon is still
vetoed. The remaining tests pin the out-of-frame / in-frame distinction, the
downstream-in-frame-stop pairing, and the 5'-proximal scan window.
"""

from __future__ import annotations

import sys

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4.constraints.base import Constraint
from bt4.constraints.uorf import UorfConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=60)


def _build_feasible(protein: str, constraint: Constraint) -> str:
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


def _hard(constraint: Constraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN, region=st.sampled_from([30, 60, None]))
def test_ok_suffix_implies_validate_clean(protein: str, region: int | None) -> None:
    constraint = UorfConstraint(region_nt=region)
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


# --------------------------------------------------------------------------- #
# Positive detection: validate must catch a real out-of-frame uORF.
# --------------------------------------------------------------------------- #


def test_validate_flags_out_of_frame_uorf() -> None:
    # ATGAATGTAA: an out-of-frame ATG at index 4 (frame +1) whose first in-frame
    # stop is TAA at index 7 -> a uORF spanning [4, 10).
    constraint = UorfConstraint(min_start=3, region_nt=None)
    violations = list(constraint.validate("ATGAATGTAA"))
    assert len(violations) == 1
    v = violations[0]
    assert v.constraint == "uorf"
    assert v.severity is Severity.HARD
    assert (v.start, v.end) == (4, 10)


def test_in_frame_internal_atg_is_not_a_uorf() -> None:
    # ATGATGTAA: the internal ATG at index 3 is IN the main frame (3 % 3 == 0),
    # so it is a Met codon, not an out-of-frame uORF.
    assert list(UorfConstraint(min_start=3, region_nt=None).validate("ATGATGTAA")) == []


def test_out_of_frame_atg_without_downstream_stop_is_not_flagged() -> None:
    # ATGAATGAAA: out-of-frame ATG at index 4 but no in-frame (+1) stop downstream
    # (positions 7 = AAA), so there is no bounded uORF to ban.
    assert list(UorfConstraint(min_start=3, region_nt=None).validate("ATGAATGAAA")) == []


def test_region_window_bounds_the_atg() -> None:
    seq = "ATGAATGTAA"  # uORF ATG at index 4.
    # A window that excludes index 4 finds nothing...
    assert list(UorfConstraint(min_start=3, region_nt=4).validate(seq)) == []
    # ...one that includes it flags the uORF.
    assert len(list(UorfConstraint(min_start=3, region_nt=5).validate(seq))) == 1


def test_min_start_skips_the_real_start() -> None:
    # With min_start=0 an out-of-frame ATG at index 1 is a uORF; the default
    # min_start=3 would skip it.
    seq = "CATGTAA"  # ATG at index 1 (frame +1), stop TAA at index 4.
    assert len(list(UorfConstraint(min_start=0, region_nt=None).validate(seq))) == 1
    assert list(UorfConstraint(min_start=3, region_nt=None).validate(seq)) == []


# --------------------------------------------------------------------------- #
# ok_suffix reads the whole prefix (global): the closing stop is vetoed.
# --------------------------------------------------------------------------- #


def test_ok_suffix_vetoes_stop_that_closes_a_uorf() -> None:
    constraint = UorfConstraint(min_start=0, region_nt=None)
    # prefix "CATG" carries an out-of-frame ATG at index 1; the codon "TAA"
    # supplies its in-frame stop, completing the uORF -> vetoed.
    assert constraint.ok_suffix("CATG", "TAA") is False
    # A codon that does not close the uORF is accepted.
    assert constraint.ok_suffix("CATG", "AAA") is True


# --------------------------------------------------------------------------- #
# Contract surface: scope, context, name, penalty, configuration.
# --------------------------------------------------------------------------- #


def test_scope_is_global() -> None:
    assert UorfConstraint().scope() is Scope.GLOBAL


def test_context_len_is_unbounded_sentinel() -> None:
    assert UorfConstraint().context_len() == sys.maxsize


def test_name_and_penalty() -> None:
    constraint = UorfConstraint()
    assert constraint.name == "uorf"
    assert constraint.penalty("AAA", "CCC") == 0.0


def test_rejects_bad_params() -> None:
    with pytest.raises(ValueError, match="min_start"):
        UorfConstraint(min_start=-1)
    with pytest.raises(ValueError, match="region_nt"):
        UorfConstraint(region_nt=0)


# --------------------------------------------------------------------------- #
# Refinement enforcement: the global-constraint gate drives uORFs to zero
# without ever raising the count (invariant #5, global edition).
# --------------------------------------------------------------------------- #


def test_refinement_removes_a_uorf_from_a_dirty_seed() -> None:
    from bt4.optimize.anneal_refine import anneal_refine

    # ATG GAT GTA AAA TAA (= M D V K stop) carries an out-of-frame uORF: the ATG
    # at index 4 (D's "AT" + V's "G") pairs with the in-frame stop "TAA" at index
    # 7. A synonymous swap of D (GAT->GAC) or V (GTA->GTC/GTT) removes it.
    residues = ["M", "D", "V", "K", "*"]
    seed = "ATGGATGTAAAATAA"
    assert translate(seed) == "MDVK*"
    u = UorfConstraint(min_start=3, region_nt=None)
    assert _hard(u, seed)  # the seed genuinely has a uORF.

    def hard(dna: str) -> int:
        return sum(1 for v in u.validate(dna) if v.severity is Severity.HARD)

    def score(dna: str) -> float:
        return -1e9 * hard(dna)

    result = anneal_refine(
        seed, residues, score, (), global_constraints=(u,), iterations=500, seed=0
    )
    assert translate(result.dna) == "MDVK*"
    assert hard(result.dna) == 0  # refinement removed the uORF...
    assert hard(result.dna) <= hard(seed)  # ...and never raised the count (#5).
