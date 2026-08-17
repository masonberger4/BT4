"""Tests for per-system biology: functional poly(A) signals and packaging headroom.

Two rules with opposite postures, and the tests pin both:

* :class:`~bt4.constraints.polya.FunctionalPolyASignalConstraint` is a *constraint*
  -- it must be strictly more permissive than the blunt hexamer ban while still
  catching the bipartite signal the cleavage machinery actually recognises.
* :mod:`bt4.pipeline.packaging` is a *report* -- BT4 controls no lever over
  cassette size, so the test that matters is that it never overstates what it
  measured.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4 import api
from bt4.constraints.polya import FunctionalPolyASignalConstraint
from bt4.domain.context import ConstructContext
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope
from bt4.pipeline.packaging import PACKAGING_LIMITS, packaging_report

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=60)
_GAP = "GC" * 5  # 10 nt, spans the cleavage-site gap without being U/GU-rich


def _build_feasible(protein: str, constraint: FunctionalPolyASignalConstraint) -> str:
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if constraint.ok_suffix(dna, codon):
                dna += codon
                break
        else:
            break
    return dna


def _hard(constraint: FunctionalPolyASignalConstraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


# --------------------------------------------------------------------------- #
# The bipartite rule: hexamer alone is not a signal.
# --------------------------------------------------------------------------- #


def test_bare_hexamer_is_not_flagged() -> None:
    """The whole point: AATAAA with no downstream element is not a poly(A) site."""
    constraint = FunctionalPolyASignalConstraint()
    assert _hard(constraint, "AATAAA" + "GC" * 25) == []


def test_hexamer_with_u_rich_element_is_flagged() -> None:
    constraint = FunctionalPolyASignalConstraint()
    dna = "AATAAA" + _GAP + "T" * 10 + "GC" * 10
    violations = _hard(constraint, dna)
    assert violations
    assert violations[0].constraint == "polya_signal"  # type: ignore[attr-defined]
    assert "downstream element" in violations[0].detail  # type: ignore[attr-defined]


def test_hexamer_with_gu_rich_element_is_flagged() -> None:
    constraint = FunctionalPolyASignalConstraint()
    assert _hard(constraint, "AATAAA" + _GAP + "GT" * 5 + "GC" * 10)


def test_pure_g_run_is_not_a_gu_rich_element() -> None:
    """A GU-rich element needs its U: a poly-G stretch must not count."""
    constraint = FunctionalPolyASignalConstraint()
    assert _hard(constraint, "AATAAA" + _GAP + "G" * 10 + "GC" * 10) == []


def test_element_inside_the_cleavage_gap_is_not_counted() -> None:
    """The DSE lies beyond the cleavage site, so a U-run flush against the hexamer
    is not in the position a downstream element occupies."""
    constraint = FunctionalPolyASignalConstraint(dse_start=10, dse_end=40)
    assert _hard(constraint, "AATAAA" + "T" * 10 + "GC" * 20) == []


def test_variant_hexamer_is_recognised() -> None:
    constraint = FunctionalPolyASignalConstraint()
    assert _hard(constraint, "ATTAAA" + _GAP + "T" * 10 + "GC" * 10)


def test_is_more_permissive_than_the_blunt_hexamer_ban() -> None:
    """The claim that justifies this rule existing beside the preset."""
    functional = FunctionalPolyASignalConstraint()
    bare = "GGGCCC" + "AATAAA" + "GC" * 25
    # The blunt preset bans this outright; the functional rule leaves it alone.
    blunt_report = api.validate(
        bare, api.OptimizeConfig(forbidden_presets=("poly_a_signal",), max_homopolymer=None)
    )
    assert not blunt_report.is_feasible
    assert _hard(functional, bare) == []


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds carry no hard violations.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN)
def test_ok_suffix_implies_validate_clean(protein: str) -> None:
    constraint = FunctionalPolyASignalConstraint()
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    assert translate(dna) == protein[: len(dna) // 3]


def test_scope_context_and_penalty() -> None:
    constraint = FunctionalPolyASignalConstraint()
    assert constraint.scope() is Scope.LOCAL
    # Hexamer + the whole downstream search window -- far too wide for the trellis,
    # which is why the pipeline routes this rule to refinement.
    assert constraint.context_len() == 6 + 40 - 1
    assert constraint.penalty("AATAAA", "TTT") == 0.0


def test_rejects_invalid_geometry() -> None:
    with pytest.raises(ValueError, match="dse_window"):
        FunctionalPolyASignalConstraint(dse_window=0)
    with pytest.raises(ValueError, match="dse_end"):
        FunctionalPolyASignalConstraint(dse_start=40, dse_end=10)
    with pytest.raises(ValueError, match="min_u"):
        FunctionalPolyASignalConstraint(dse_window=4, min_u=9)
    with pytest.raises(ValueError, match="hexamers"):
        FunctionalPolyASignalConstraint(hexamers=())


# --------------------------------------------------------------------------- #
# End-to-end through the pipeline (refinement-enforced, honestly reported).
# --------------------------------------------------------------------------- #


def test_avoid_polya_runs_and_reports_enforcement() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(avoid_polya=True))
    assert "polya_signal_enforced" in result.audit
    assert result.audit["polya_signal_enforced"] in {"clean", "partial"}


def test_avoid_polya_is_audited_by_validate() -> None:
    """A GLOBAL-routed rule must be checked by validate too (the A.2 lesson)."""
    dna = "ATG" + "AATAAA" + _GAP + "T" * 10 + "GCGCGC"
    report = api.validate(dna, api.OptimizeConfig(avoid_polya=True, max_homopolymer=None))
    assert not report.is_feasible
    assert any(v.constraint == "polya_signal" for v in report.violations)


def test_avoid_polya_rejects_a_budget_combination() -> None:
    with pytest.raises(ValueError, match="avoid_polya"):
        api.optimize("MAALKHETQWY", api.OptimizeConfig(avoid_polya=True, gc_min=30))


def test_avoid_polya_changes_the_manifest() -> None:
    plain = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    guarded = api.optimize(
        "MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None, avoid_polya=True)
    )
    assert plain.audit["manifest"] != guarded.audit["manifest"]


# --------------------------------------------------------------------------- #
# Packaging: a report that must never overstate what it measured.
# --------------------------------------------------------------------------- #


def test_packaging_counts_only_the_cds_without_context() -> None:
    report = packaging_report("ATG" * 100, system="aav")
    assert report.counted_nt == 300
    assert report.counted == ("CDS",)
    assert not report.complete
    # The caveat is mandatory when the construct is only partly known.
    assert "not included" in report.summary()


def test_packaging_includes_supplied_context() -> None:
    context = ConstructContext(upstream="A" * 500, downstream="T" * 400)
    report = packaging_report("ATG" * 100, system="aav", context=context)
    assert report.context_nt == 900
    assert report.counted_nt == 1200
    assert "supplied context" in report.counted


def test_packaging_flags_an_oversized_cassette() -> None:
    context = ConstructContext(upstream="A" * 2500, downstream="T" * 2500)
    report = packaging_report("ATG" * 100, system="aav", context=context)
    assert report.over_limit
    assert report.headroom_nt < 0
    assert "OVER by" in report.summary()


def test_packaging_drops_the_caveat_only_when_completeness_is_asserted() -> None:
    context = ConstructContext(upstream="A" * 500, downstream="T" * 400)
    partial = packaging_report("ATG" * 100, context=context)
    whole = packaging_report("ATG" * 100, context=context, complete=True)
    assert "not included" in partial.summary()
    assert "not included" not in whole.summary()


def test_packaging_systems_and_overrides() -> None:
    assert packaging_report("ATG", system="lvv").limit_nt == PACKAGING_LIMITS["lvv"]
    assert packaging_report("ATG", limit_nt=123).limit_nt == 123
    with pytest.raises(ValueError, match="unknown vector system"):
        packaging_report("ATG", system="nonesuch")
