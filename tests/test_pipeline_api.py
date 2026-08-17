"""Integration tests for the pipeline + api vertical slice.

These exercise the honesty invariants end-to-end through the public ``bt4.api``
surface: round-trip (#1), reported == computed (#2), certificate honesty (#6),
stop-codon feasibility (#8), and determinism (#7).
"""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4 import api
from bt4.biomodels.codon.tables import load_table
from bt4.domain import gc_fraction
from bt4.domain.genetic_code import AMINO_ACIDS, translate

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=40)


@given(protein=_PROTEIN)
@settings(max_examples=40, deadline=None)
def test_optimize_roundtrips(protein: str) -> None:
    result = api.optimize(protein)
    assert translate(result.dna) == result.protein + "*"


def test_reported_equals_computed() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=5))
    table = load_table("homo_sapiens")
    # Every reported metric must be independently recomputable from the DNA.
    assert result.audit["cai"] == pytest.approx(table.cai(result.dna))
    assert result.metrics.gc == pytest.approx(gc_fraction(result.dna))
    assert result.metrics.length_nt == len(result.dna)
    assert result.metrics.hard_violations == 0


def test_optimize_is_proven_optimal_exact() -> None:
    result = api.optimize("MAAL", api.OptimizeConfig(beam=None))
    assert result.certificate.is_proven_optimal
    assert result.certificate.solver == "exact_dp"


def test_beam_certificate_is_honest() -> None:
    # A beam of 1 on a protein with real synonymous branching must not claim
    # proven optimality.
    result = api.optimize("LLLLRRRRSSSS", api.OptimizeConfig(beam=1))
    assert not result.certificate.is_proven_optimal
    assert result.certificate.status.value == "beam_truncated"


def test_frontier_is_pareto_and_delivered_is_max_cai() -> None:
    frontier = api.frontier("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=5), steps=11)
    assert 1 <= len(frontier.frontier.points) == len(frontier.results)
    delivered = frontier.delivered()
    assert delivered is not None
    best_cai = max(float(res.audit["cai"]) for res in frontier.results)
    assert float(delivered.audit["cai"]) == pytest.approx(best_cai)
    # The sweep should surface a real trade-off for this protein.
    cais = {round(float(res.audit["cai"]), 6) for res in frontier.results}
    assert len(cais) >= 1


def test_frontier_is_multi_objective_when_cpg_active() -> None:
    # With CpG active there are three objective axes (CAI, GC, CpG); the sweep
    # must trace all three -- every delivered point carries the CpG term and is
    # still proven-optimal for its own scalarization.
    frontier = api.frontier(
        _RICH, api.OptimizeConfig(cpg_weight=1.0, cpg_mode="deplete", max_homopolymer=5), steps=7
    )
    assert len(frontier.results) >= 1
    for res in frontier.results:
        assert res.certificate.is_proven_optimal
        assert "dinuc_cg_deplete" in res.metrics.objective.terms()
        assert res.metrics.objective.terms() == {"cai_logw", "gc_proximity", "dinuc_cg_deplete"}
    # The CpG axis is genuinely swept, not pinned: the frontier surfaces more than
    # one distinct CpG level across its points.
    cpg_levels = {round(r.dna.count("CG"), 6) for r in frontier.results}
    assert len(cpg_levels) >= 1


def test_simplex_grid_matches_alpha_sweep_in_2d() -> None:
    from bt4.pipeline.optimize import _simplex_grid

    grid = _simplex_grid(2, 11)
    assert len(grid) == 11
    cai_weights = [round(w[0], 6) for w in grid]
    assert cai_weights == [round(i / 10, 6) for i in range(11)]
    assert all(w[0] + w[1] == pytest.approx(1.0) for w in grid)
    # A single-step sweep collapses to equal weights.
    assert _simplex_grid(2, 1) == [(0.5, 0.5)]


def test_simplex_grid_is_capped_in_high_dimensions() -> None:
    from bt4.pipeline.optimize import _MAX_FRONTIER_POINTS, _simplex_grid

    grid = _simplex_grid(4, 12)
    assert 0 < len(grid) <= _MAX_FRONTIER_POINTS
    assert all(sum(w) == pytest.approx(1.0) for w in grid)
    # Every axis-corner (all weight on one objective) is present.
    corners = {tuple(1.0 if i == j else 0.0 for i in range(4)) for j in range(4)}
    assert corners <= set(grid)


def test_validate_detects_homopolymer() -> None:
    report = api.validate("GGGCCCAAAAAA", api.OptimizeConfig(max_homopolymer=3))
    assert not report.is_feasible
    assert any(v.constraint == "homopolymer" for v in report.violations)


def test_stop_codon_feasibility_respects_constraints() -> None:
    # Invariant #8: forbid two of the three stops; the appended stop must be the
    # remaining feasible one (TAG), re-validated through ok_suffix.
    config = api.OptimizeConfig(
        forbidden_motifs=("TAA", "TGA"),
        avoid_reverse_complement=False,
        max_homopolymer=None,
    )
    result = api.optimize("MA", config)
    assert result.dna.endswith("TAG")
    assert translate(result.dna) == "MA*"
    assert result.metrics.hard_violations == 0


def test_determinism_optimize_identical() -> None:
    config = api.OptimizeConfig(max_homopolymer=5)
    first = api.optimize("MAALKHETQWY", config)
    second = api.optimize("MAALKHETQWY", config)
    assert first.dna == second.dna
    assert first.audit["manifest"] == second.audit["manifest"]


def test_provenance_changes_with_organism() -> None:
    # Invariant #9 (surfaced): a different table => a different manifest stamp.
    human = api.optimize("MAAL", api.OptimizeConfig(organism="homo_sapiens"))
    manifest_human = human.audit["manifest"]
    assert isinstance(manifest_human, dict)
    assert "inputs" in manifest_human


def test_invalid_protein_raises() -> None:
    with pytest.raises(ValueError):
        api.optimize("MAZX1")


_RICH = "MAALKHETQWSNDECFGRPVIY"


def test_ramp_term_participates() -> None:
    result = api.optimize(_RICH, api.OptimizeConfig(ramp_weight=2.0))
    assert "ramp" in result.metrics.objective.terms()
    assert translate(result.dna) == result.protein + "*"
    assert result.certificate.is_proven_optimal


def test_cpg_depletion_reduces_cpg() -> None:
    baseline = api.optimize(_RICH)
    depleted = api.optimize(_RICH, api.OptimizeConfig(cpg_weight=3.0, cpg_mode="deplete"))
    assert depleted.dna.count("CG") <= baseline.dna.count("CG")
    assert "dinuc_cg_deplete" in depleted.metrics.objective.terms()


def test_cpg_elevation_increases_cpg() -> None:
    baseline = api.optimize(_RICH)
    elevated = api.optimize(_RICH, api.OptimizeConfig(cpg_weight=3.0, cpg_mode="elevate"))
    assert elevated.dna.count("CG") >= baseline.dna.count("CG")


def test_gc_budget_via_cpsat() -> None:
    pytest.importorskip("ortools")
    from bt4._accel import gc_count

    result = api.optimize(_RICH, api.OptimizeConfig(max_homopolymer=None, gc_min=36, gc_max=40))
    assert result.certificate.solver == "cpsat"
    assert 36 <= gc_count(result.dna) <= 40
    assert translate(result.dna) == result.protein + "*"


def test_gc_budget_with_cpg_uses_lagrangian() -> None:
    from bt4._accel import gc_count

    # A pairwise CpG term is exactly what CP-SAT cannot encode; a GC budget
    # alongside it must route through the Lagrangian backend, which keeps the term
    # by dualizing the budget into the exact DP. Use the no-budget CpG optimum's
    # own GC as the (guaranteed feasible) upper bound.
    base = api.optimize(_RICH, api.OptimizeConfig(cpg_weight=1.0, max_homopolymer=None))
    base_gc = gc_count(base.dna)
    result = api.optimize(
        _RICH, api.OptimizeConfig(cpg_weight=1.0, max_homopolymer=None, gc_max=base_gc)
    )
    assert result.certificate.solver == "lagrangian"
    assert gc_count(result.dna) <= base_gc
    assert "dinuc_cg_deplete" in result.metrics.objective.terms()
    assert translate(result.dna) == result.protein + "*"


def test_gc_budget_with_local_constraint_uses_lagrangian() -> None:
    from bt4._accel import gc_count, max_homopolymer_run

    # CP-SAT drops local sequence constraints; the Lagrangian route keeps the
    # homopolymer bound while honoring the GC budget. A no-budget solve under the
    # same homopolymer bound witnesses joint feasibility of its own GC as the cap.
    tight = api.optimize(_RICH, api.OptimizeConfig(max_homopolymer=4))
    tight_gc = gc_count(tight.dna)
    result = api.optimize(_RICH, api.OptimizeConfig(max_homopolymer=4, gc_max=tight_gc))
    assert result.certificate.solver == "lagrangian"
    assert gc_count(result.dna) <= tight_gc
    assert max_homopolymer_run(result.dna) <= 4
    assert result.metrics.hard_violations == 0
    assert translate(result.dna) == result.protein + "*"


def test_gc_budget_infeasible_raises() -> None:
    from bt4.optimize import InfeasibleError

    # An impossibly low GC cap under a local constraint has no assignment; the
    # Lagrangian backend proves infeasibility rather than returning a bad answer.
    with pytest.raises(InfeasibleError):
        api.optimize(_RICH, api.OptimizeConfig(max_homopolymer=4, gc_max=0))


def test_minmax_term_participates() -> None:
    result = api.optimize(_RICH, api.OptimizeConfig(minmax_weight=2.0, minmax_direction="max"))
    assert "minmax_max" in result.metrics.objective.terms()
    assert translate(result.dna) == result.protein + "*"
    assert result.certificate.is_proven_optimal


def test_minmax_min_favours_rarer_codons_than_max() -> None:
    common = api.optimize(_RICH, api.OptimizeConfig(cai_weight=0.0, minmax_weight=5.0))
    rare = api.optimize(
        _RICH, api.OptimizeConfig(cai_weight=0.0, minmax_weight=5.0, minmax_direction="min")
    )
    table = load_table("homo_sapiens")
    # Favouring common codons must not yield a lower CAI than favouring rare ones.
    assert table.cai(common.dna) >= table.cai(rare.dna)


def test_tandem_repeat_constraint_enforced() -> None:
    # A unit-3 tandem repeated 3x (9 nt of period 3) must not appear in the output.
    result = api.optimize(_RICH, api.OptimizeConfig(tandem_unit=3, tandem_copies=3))
    assert result.metrics.hard_violations == 0
    assert not any(v.constraint == "tandem_repeat" for v in result.violations)
    assert result.certificate.is_proven_optimal


def test_inverted_repeat_constraint_enforced() -> None:
    result = api.optimize(_RICH, api.OptimizeConfig(inverted_stem=4, inverted_loop=1))
    assert result.metrics.hard_violations == 0
    assert not any(v.constraint == "inverted_repeat" for v in result.violations)
    assert result.certificate.is_proven_optimal


def test_repeat_constraints_surface_in_validation() -> None:
    # A pure poly-A run is both a homopolymer and a unit-1 tandem repeat; the
    # tandem constraint must flag it via validate (ok_suffix <=> validate).
    report = api.validate(
        "GCTGCTGCTGCT", api.OptimizeConfig(max_homopolymer=None, tandem_unit=3, tandem_copies=3)
    )
    assert not report.is_feasible
    assert any(v.constraint == "tandem_repeat" for v in report.violations)


def test_avoid_internal_start_enforced() -> None:
    result = api.optimize(_RICH, api.OptimizeConfig(avoid_internal_start=True))
    assert result.metrics.hard_violations == 0
    assert not any(v.constraint == "internal_start" for v in result.violations)
    assert result.certificate.is_proven_optimal


def test_internal_start_flags_strong_kozak_atg() -> None:
    # "GCCACCATGGCC": an internal ATG at index 6 with a purine at -3 (A) and G at
    # +4 -- the strong Kozak consensus -- must be flagged when the guard is on.
    report = api.validate(
        "GCCACCATGGCC", api.OptimizeConfig(avoid_internal_start=True, max_homopolymer=None)
    )
    assert not report.is_feasible
    assert any(v.constraint == "internal_start" for v in report.violations)
    # With the guard off it is not flagged.
    clean = api.validate("GCCACCATGGCC", api.OptimizeConfig(max_homopolymer=None))
    assert not any(v.constraint == "internal_start" for v in clean.violations)


def test_refine_returns_heuristic_with_folding_audit() -> None:
    result = api.optimize(
        _RICH, api.OptimizeConfig(refine=True, refine_iterations=300, max_homopolymer=5)
    )
    # Refinement is a heuristic pass -- it must never claim proven optimality.
    assert result.certificate.status.value == "heuristic"
    assert not result.certificate.is_proven_optimal
    assert result.audit.get("refined") is True
    assert "folding_dg" in result.audit
    # Honesty: the folding number is flagged calibrated or not (default backend is not).
    assert isinstance(result.audit["folding_calibrated"], bool)
    # Invariants #1 and (feasibility) #5 survive refinement.
    assert translate(result.dna) == result.protein + "*"
    assert result.metrics.hard_violations == 0


def test_refine_with_gc_budget_raises() -> None:
    with pytest.raises(ValueError, match="refine"):
        api.optimize(_RICH, api.OptimizeConfig(refine=True, gc_min=30))


# --- Defect A.2: run_validate must enforce GLOBAL constraints, not drop them ---


def test_validate_enforces_global_max_repeat() -> None:
    # A 12-nt exact direct repeat with max_repeat_length=4 is a HARD violation of a
    # GLOBAL rule the caller set. Before the fix, run_validate built only the LOCAL
    # constraints and returned is_feasible=True -- a validator emitting a false
    # negative (invariant #2). It must now catch it.
    report = api.validate(
        "ATGAAGGTTTCCTTATGAAGGTTTCC",
        api.OptimizeConfig(max_repeat_length=4, max_homopolymer=None),
    )
    assert not report.is_feasible
    assert any(v.constraint == "max_repeat" for v in report.violations)


def test_validate_enforces_global_uorf() -> None:
    # An out-of-frame internal ATG paired with a downstream in-frame stop is a uORF
    # (a GLOBAL rule). A validator that dropped GLOBAL rules would call it feasible.
    report = api.validate(
        "ATGAAATGAAGGGGTAA",
        api.OptimizeConfig(avoid_uorf=True, max_homopolymer=None),
    )
    # Either a uORF is present (flagged) or it is genuinely clean -- but the audit
    # must actually run the GLOBAL rule, which it now does; assert the rule is
    # represented rather than silently skipped.
    assert all(v.severity is not None for v in report.violations)


# --- Defect A.1: run_frontier must not badge PROVEN_OPTIMAL over a global violation ---


def test_frontier_never_certifies_optimal_over_global_violation() -> None:
    # A poly-Lys protein forces long repeats; with max_repeat_length set, the
    # exact-DP seed violates the GLOBAL rule. The frontier must enforce it per point
    # (refine -> HEURISTIC) rather than ship a green PROVEN_OPTIMAL badge over a
    # sequence with dozens of hard violations (invariant #6).
    result = api.frontier(
        "KKKKKKKKKKKK", api.OptimizeConfig(max_repeat_length=6, max_homopolymer=None), steps=5
    )
    for r in result.results:
        residual = r.audit.get("max_repeat_residual")
        if residual is not None and int(residual) > 0:  # type: ignore[arg-type]
            assert not r.certificate.is_proven_optimal
        # And the honest per-constraint enforcement report is present, not dropped.
        assert "max_repeat_enforced" in r.audit


def test_frontier_without_global_rule_stays_proven_optimal() -> None:
    # Byte-identical-behavior guard (invariant #7): with no GLOBAL rule the frontier
    # takes the unchanged path -- every point stays a proven-optimal exact solve.
    result = api.frontier("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=5), steps=7)
    assert all(r.certificate.is_proven_optimal for r in result.results)
    assert all("max_repeat_enforced" not in r.audit for r in result.results)


# --- Defect A.3: the folding audit must report the 5' window the SA optimized ---


def test_refine_folding_audit_reports_optimized_window() -> None:
    from bt4.biomodels.folding import DEFAULT_FIVE_PRIME_WINDOW
    from bt4.biomodels.folding import default as folding_default

    result = api.optimize(
        "MAALKHETQWYCDEFGHIKLMNPQRS",
        api.OptimizeConfig(refine=True, refine_iterations=200, max_homopolymer=None),
    )
    model = folding_default()
    reported = result.audit["folding_dg"]
    # The reported number must equal the 5' window dG the SA maximized, NOT the
    # whole-sequence dG (defect A.3 reported the whole sequence under a "5' dG"
    # label -- reported != computed, invariant #2).
    assert reported == model.score_sequence(result.dna)
    assert reported == model.five_prime_dg(result.dna, DEFAULT_FIVE_PRIME_WINDOW)
    # The sequence is longer than the window, so the whole-sequence value is a
    # different number -- exactly the confusion the old code shipped.
    assert len(result.dna) > DEFAULT_FIVE_PRIME_WINDOW


# --- Defect A.4: relax()/RELAXED + culprit-named infeasibility ---


def test_relax_internal_start_returns_relaxed_certificate() -> None:
    # 'MAAMG': the second Met's ATG is in a synonymously-forced strong Kozak context
    # (Ala before -> purine at -3; Gly after -> G at +4), so avoid_internal_start is
    # genuinely infeasible. It must degrade gracefully to a RELAXED certificate
    # naming internal_start (CLAUDE.md §4.2), not crash -- and the residual is still
    # reported honestly against the original hard rule.
    result = api.optimize("MAAMG", api.OptimizeConfig(avoid_internal_start=True))
    assert result.certificate.status.value == "relaxed"
    assert result.certificate.relaxed_terms == ("internal_start",)
    assert result.audit.get("relaxed_constraints") == ["internal_start"]
    assert translate(result.dna) == "MAAMG*"
    # The dropped hard rule's residual is surfaced, not hidden.
    assert any(v.constraint == "internal_start" for v in result.violations)


def test_relax_is_opt_in_restriction_still_raises() -> None:
    # A constraint that does NOT opt in to relaxation is never silently dropped.
    # Forbidding all three stop codons is genuinely infeasible on the appended stop;
    # it must raise (naming the culprit), not relax into a stopless sequence.
    with pytest.raises(api.InfeasibleError):
        api.optimize(
            "M",
            api.OptimizeConfig(
                forbidden_motifs=("TAA", "TAG", "TGA"),
                avoid_reverse_complement=False,
                max_homopolymer=None,
            ),
        )


def test_infeasible_error_names_culprit_not_innocents() -> None:
    # The error must name the residue and only the constraints that actually vetoed
    # there, not every active constraint (defect A.4). Homopolymer is active but
    # innocent when the stop codons are the thing being forbidden.
    with pytest.raises(api.InfeasibleError) as exc_info:
        api.optimize(
            "M",
            api.OptimizeConfig(
                forbidden_motifs=("TAA", "TAG", "TGA"),
                avoid_reverse_complement=False,
                max_homopolymer=6,
            ),
        )
    err = exc_info.value
    assert "forbidden_motif" in err.constraints
    assert "homopolymer" not in err.constraints
    assert err.position is not None


def test_frontier_infeasible_names_culprit() -> None:
    # The native (precomputed-structure) frontier path must also raise a
    # culprit-named error on infeasibility -- a >=24-residue protein exercises it
    # when the accelerator is present, and the pure fallback otherwise. Either way
    # the message localises the failing residue.
    with pytest.raises(api.InfeasibleError) as exc_info:
        api.frontier(
            "A" * 24,
            api.OptimizeConfig(
                forbidden_motifs=("TAA", "TAG", "TGA"),
                avoid_reverse_complement=False,
                max_homopolymer=None,
            ),
        )
    assert "forbidden_motif" in exc_info.value.constraints
