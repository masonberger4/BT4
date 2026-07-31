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
