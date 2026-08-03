"""Tests for the optimality certificate."""

from __future__ import annotations

import pytest

from bt4.domain.certificate import OptimalityCertificate, OptimalityStatus


def test_proven_helper_is_optimal_with_zero_gap() -> None:
    cert = OptimalityCertificate.proven("exact_dp", detail="full state space")
    assert cert.is_proven_optimal
    assert cert.gap == 0.0
    assert cert.solver == "exact_dp"


def test_proven_cannot_carry_a_gap_or_relaxed_terms() -> None:
    with pytest.raises(ValueError):
        OptimalityCertificate(OptimalityStatus.PROVEN_OPTIMAL, "exact_dp", gap=0.1)
    with pytest.raises(ValueError):
        OptimalityCertificate(
            OptimalityStatus.PROVEN_OPTIMAL, "exact_dp", relaxed_terms=("gc",)
        )


def test_gap_must_be_in_unit_interval() -> None:
    with pytest.raises(ValueError):
        OptimalityCertificate(OptimalityStatus.GAP_BOUNDED, "cpsat_ilp", gap=1.5)


def test_relaxed_certificate_records_terms() -> None:
    cert = OptimalityCertificate(
        OptimalityStatus.RELAXED, "lagrangian", gap=0.05, relaxed_terms=("cpg", "gc")
    )
    assert not cert.is_proven_optimal
    assert cert.relaxed_terms == ("cpg", "gc")


def test_sampled_makes_no_optimality_claim() -> None:
    # Library / degenerate-design mode: a stochastic draw, never an optimum.
    cert = OptimalityCertificate(
        OptimalityStatus.SAMPLED, "library_sampler", detail="drawn from codon distribution"
    )
    assert cert.status is OptimalityStatus.SAMPLED
    assert cert.status.value == "sampled"
    assert not cert.is_proven_optimal
    assert cert.gap is None
    assert cert.relaxed_terms == ()
