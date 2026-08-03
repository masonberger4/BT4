"""Tests for wiring the codon-pair-bias (CpbTerm) objective into the pipeline.

The honesty rule (CLAUDE.md §8): there is no bundled default codon-pair table, so
a run that asks for codon-pair bias (``cpb_weight != 0``) must supply a reference
CDS set or be refused -- never fall back to a fabricated table. When a reference
is given, the term is built from it, is exact in the trellis (PAIRWISE), enters
the objective vector, and its content hash enters the manifest.
"""

from __future__ import annotations

import pytest

from bt4 import api
from bt4.domain import translate
from bt4.pipeline.optimize import OptimizeConfig

_PROTEIN = "MKAILVDEQTRSFYNWHGP"
# A tiny reference CDS set (real reference sets would be far larger; these are
# valid ACGT, codon-aligned, and exercise the Coleman CPS build deterministically).
_REF_CDS = (
    "ATGGCCGCTAAAGGACAGACCTGA",
    "ATGAAACGTGACGAACAGACCTAA",
    "ATGGACGACCTGATGAAAGCCTAA",
)


def _clean_ref() -> tuple[str, ...]:
    return _REF_CDS


def test_cpb_weight_without_reference_raises() -> None:
    with pytest.raises(ValueError, match="cpb_weight is set but cpb_reference_cds"):
        api.optimize(_PROTEIN, OptimizeConfig(cpb_weight=1.0))


def test_cpb_with_reference_optimizes_and_scores() -> None:
    ref = _clean_ref()
    result = api.optimize(_PROTEIN, OptimizeConfig(cpb_weight=1.0, cpb_reference_cds=ref))
    assert translate(result.dna) == _PROTEIN + "*"
    # The codon-pair term is present in the recomputed objective vector.
    assert "codon_pair" in result.metrics.objective.values


def test_cpb_reference_enters_the_manifest() -> None:
    ref = _clean_ref()
    a = api.optimize(_PROTEIN, OptimizeConfig(cpb_weight=1.0, cpb_reference_cds=ref))
    manifest_a = a.audit["manifest"]
    assert isinstance(manifest_a, dict)
    inputs = manifest_a["inputs"]
    assert isinstance(inputs, dict)
    assert "codon_pair_cds_sha256" in inputs
    # A different reference CDS set yields a different manifest (invariant #9).
    ref2 = (*ref, "ATGCGTGCCGATTTAGGCTAA")
    b = api.optimize(_PROTEIN, OptimizeConfig(cpb_weight=1.0, cpb_reference_cds=ref2))
    manifest_b = b.audit["manifest"]
    assert isinstance(manifest_b, dict)
    assert manifest_b["inputs"]["codon_pair_cds_sha256"] != inputs["codon_pair_cds_sha256"]


def test_cpb_off_by_default_no_manifest_entry() -> None:
    result = api.optimize(_PROTEIN, OptimizeConfig())
    manifest = result.audit["manifest"]
    assert isinstance(manifest, dict)
    assert "codon_pair_cds_sha256" not in manifest["inputs"]
    assert "codon_pair" not in result.metrics.objective.values


def test_cpb_is_a_frontier_axis_when_weighted() -> None:
    ref = _clean_ref()
    frontier = api.frontier(
        _PROTEIN, OptimizeConfig(cpb_weight=1.0, cpb_reference_cds=ref), steps=7
    )
    # Every frontier point's objective vector carries the codon-pair axis.
    for r in frontier.results:
        assert "codon_pair" in r.metrics.objective.values
