"""Tests for the io serialization layer: FASTA text and versioned JSON."""

from __future__ import annotations

import json

from bt4.domain import (
    Metrics,
    ObjectiveVector,
    OptimalityCertificate,
    Result,
    Severity,
    Violation,
)
from bt4.io import result_to_dict, result_to_json, to_fasta


def _make_result() -> Result:
    """Build a small, fully populated Result by hand."""
    objective = ObjectiveVector({"cai_logw": -1.5, "gc_proximity": -0.1})
    metrics = Metrics(
        objective=objective,
        gc=0.55,
        length_nt=12,
        hard_violations=1,
        soft_violations=1,
    )
    certificate = OptimalityCertificate.proven("exact_dp", detail="full state space")
    violations = (
        Violation(
            constraint="homopolymer",
            severity=Severity.HARD,
            start=0,
            end=7,
            detail="run of 7 A",
        ),
        Violation(
            constraint="gc_window",
            severity=Severity.SOFT,
            start=3,
            end=9,
            detail="window GC low",
        ),
    )
    return Result(
        protein="MAA",
        dna="ATGGCTGCTTAA",
        metrics=metrics,
        certificate=certificate,
        violations=violations,
        audit={"cai": 0.8, "seed": 0},
    )


def test_to_fasta_wraps_at_60_columns() -> None:
    dna = "ACGT" * 37 + "AC"  # length 150
    assert len(dna) == 150
    text = to_fasta(dna, header="job1")
    lines = text.splitlines()
    assert lines[0] == ">job1"
    assert lines[0].startswith(">")
    body = lines[1:]
    assert all(len(line) <= 60 for line in body)
    assert "".join(body) == dna
    assert text.endswith("\n")


def test_to_fasta_empty_dna_is_header_only() -> None:
    text = to_fasta("", header="empty")
    assert text == ">empty\n"


def test_to_fasta_default_header() -> None:
    text = to_fasta("ATGTAA")
    assert text == ">bt4\nATGTAA\n"


def test_result_to_json_is_deterministic() -> None:
    result = _make_result()
    first = result_to_json(result)
    second = result_to_json(result)
    assert first == second


def test_result_to_json_round_trips_to_dict() -> None:
    result = _make_result()
    parsed = json.loads(result_to_json(result))
    assert isinstance(parsed, dict)
    assert parsed["dna"] == result.dna
    assert parsed["schema_version"] == "1"
    assert parsed["protein"] == result.protein
    assert parsed["length_nt"] == result.metrics.length_nt
    assert len(parsed["violations"]) == len(result.violations)


def test_certificate_status_serializes_to_value_string() -> None:
    result = _make_result()
    parsed = json.loads(result_to_json(result))
    assert parsed["certificate"]["status"] == "proven_optimal"
    assert parsed["certificate"]["solver"] == "exact_dp"
    assert parsed["certificate"]["gap"] == 0.0
    assert parsed["certificate"]["relaxed_terms"] == []


def test_severity_serializes_to_value_string() -> None:
    result = _make_result()
    data = result_to_dict(result)
    violations = data["violations"]
    assert isinstance(violations, list)
    severities = {v["severity"] for v in violations}
    assert severities == {"hard", "soft"}


def test_objective_terms_present_and_sorted() -> None:
    result = _make_result()
    data = result_to_dict(result)
    objective = data["objective"]
    assert isinstance(objective, dict)
    assert list(objective) == sorted(objective)
    assert objective["cai_logw"] == -1.5
    assert objective["gc_proximity"] == -0.1


def test_audit_passes_through() -> None:
    result = _make_result()
    data = result_to_dict(result)
    assert data["audit"] == {"cai": 0.8, "seed": 0}


def test_metrics_block_mirrors_recomputed_fields() -> None:
    result = _make_result()
    data = result_to_dict(result)
    metrics = data["metrics"]
    assert isinstance(metrics, dict)
    assert metrics["gc"] == 0.55
    assert metrics["length_nt"] == 12
    assert metrics["hard_violations"] == 1
    assert metrics["soft_violations"] == 1
