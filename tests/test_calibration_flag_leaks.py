"""Regression tests for honesty-flag leaks found before promotion could arm them.

Every wrapped splice CNN and the RiboNN head currently ship ``calibrated=False``,
which masks a class of defect: code that *forwards* or *fails to record* a
calibration flag is harmless while the flag is always ``False``, and becomes wrong
the moment a maintainer promotes a backend through the attestation seam. An audit
of the flag's propagation found four such places; these tests pin the fixes so the
leaks cannot reopen.

The four, and why each is a defect rather than a preference:

1. **A flanked score is a different input regime.** The integration-fidelity gate
   is captured on the bare-CDS path, where the adapter pads its ~10 kb window with
   literal ``N``. Replacing that padding with real vector/UTR flanks is a regime the
   attestation never exercised -- ``score_in_context``'s own docstring says so --
   yet it forwarded the backend's ``calibrated`` onto the flanked result.
2. **The audit path does not read that result's flag.** ``audit_splice`` reads
   ``predictor.calibrated`` for every ``SpliceFlag``, every
   ``BackendCandidateAudit``, and the report-level ``all_calibrated`` (which drives
   Studio's banner), so ``_FlankedPredictor`` needs its own guard.
3. **An attestation claims the full weight map, so it may only promote the
   configuration that loads all of it.** A Pangolin gate run at one tissue touches
   3 of 12 weight files and 1 of 4 output channels, yet records all 12.
4. **A calibrated expression head steers delivery without entering the stamp.**
   Two different heads reranking one frontier returned different delivered DNA
   under byte-identical manifests -- invariant #9 broken, verified by execution.

None of these needs torch, TensorFlow or licensed weights: a stub predictor
reporting ``calibrated=True`` reproduces every one.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from bt4 import api
from bt4.biomodels.expression.base import ExpressionResult
from bt4.biomodels.splice import (
    AttestationError,
    PangolinSplicePredictor,
    SpliceResult,
    attest_backend,
    score_in_context,
    verified_predictor,
)
from bt4.biomodels.splice.pangolin import (
    DEFAULT_TISSUES,
    PINNED_WEIGHT_SHA256,
    FidelityReport,
)
from bt4.pipeline.optimize import OptimizeConfig
from bt4.pipeline.splice_audit import _FlankedPredictor

_PROTEIN = "MVSGDKLWYAAC"


@dataclass(frozen=True, slots=True)
class _CalibratedStub:
    """A splice backend that claims calibration -- stands in for a promoted CNN."""

    tag: str = "stub"

    @property
    def name(self) -> str:
        return f"stub-{self.tag}"

    @property
    def calibrated(self) -> bool:
        return True

    def score_sequence(self, dna: str) -> SpliceResult:
        n = len(dna)
        return SpliceResult(
            donor=tuple(0.5 for _ in range(n)),
            acceptor=tuple(0.0 for _ in range(n)),
            model_name=self.name,
            calibrated=True,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        return 0.0


@dataclass(frozen=True)
class _CalibratedHead:
    """An expression head that claims calibration and ranks by a tunable key."""

    tag: str
    flip: bool

    @property
    def name(self) -> str:
        return f"stub-head-{self.tag}"

    @property
    def calibrated(self) -> bool:
        return True

    @property
    def units(self) -> str:
        return "log-TE"

    def score_sequence(self, dna: str) -> ExpressionResult:
        score = float(sum(dna.encode())) % 97
        return ExpressionResult(
            score=-score if self.flip else score,
            model_name=self.name,
            calibrated=True,
            units=self.units,
        )


# --------------------------------------------------------------------------
# 1. score_in_context must not carry calibration into the flanked regime


def test_score_in_context_clears_calibrated_when_flanks_are_applied() -> None:
    """Real flanks are a regime the N-padded attestation never covered."""
    cds = "ATGGTAAGTGGCGATCGATCGATCGTAA"
    flanked = score_in_context(_CalibratedStub(), cds, "GGGCCCAAA", "TTTGGG")
    assert flanked.calibrated is False
    assert len(flanked.donor) == len(cds)


def test_score_in_context_preserves_calibrated_without_flanks() -> None:
    """With no context the call is a pass-through, so the flag must survive."""
    cds = "ATGGTAAGTGGCGATCGATCGATCGTAA"
    assert score_in_context(_CalibratedStub(), cds, "", "").calibrated is True


@pytest.mark.parametrize(
    ("upstream", "downstream"),
    [("GGGCCC", ""), ("", "TTTGGG"), ("GGGCCC", "TTTGGG")],
)
def test_either_flank_alone_is_enough_to_clear_it(upstream: str, downstream: str) -> None:
    """One-sided context still changes the regime, so it still clears the flag."""
    result = score_in_context(_CalibratedStub(), "ATGGTAAGTGGCTAA", upstream, downstream)
    assert result.calibrated is False


# --------------------------------------------------------------------------
# 2. _FlankedPredictor is the seam audit_splice actually reads


def test_flanked_predictor_reports_uncalibrated_when_wrapping_flanks() -> None:
    """``audit_splice`` reads ``predictor.calibrated``, so this property is load-bearing."""
    inner = _CalibratedStub()
    assert _FlankedPredictor(inner, "ACGT", "").calibrated is False
    assert _FlankedPredictor(inner, "", "ACGT").calibrated is False
    assert _FlankedPredictor(inner, "ACGT", "TTTT").calibrated is False


def test_flanked_predictor_passes_through_without_flanks() -> None:
    """Constructed with no context it is a transparent wrapper."""
    assert _FlankedPredictor(_CalibratedStub(), "", "").calibrated is True


def test_flanked_predictor_keeps_the_backend_name() -> None:
    """It is the same model -- only the regime changed, so the name must not."""
    inner = _CalibratedStub()
    assert _FlankedPredictor(inner, "ACGT", "").name == inner.name


# --------------------------------------------------------------------------
# 3. An attestation may only promote the configuration that loads all its weights


def _passing_attestation():
    """Build a passing Pangolin attestation over the full pinned weight map."""
    report = FidelityReport(passed=True, max_abs_deviation=1e-6, n_cases=5, tolerance=1e-3)
    return attest_backend("pangolin", report, dict(PINNED_WEIGHT_SHA256), bt4_version="0.4.0")


def test_default_tissue_set_promotes() -> None:
    """The configuration that loads every pinned file is promotable."""
    promoted = verified_predictor(PangolinSplicePredictor(), _passing_attestation())
    assert promoted.calibrated is True
    assert promoted.tissues == DEFAULT_TISSUES


@pytest.mark.parametrize("tissues", [("brain",), ("heart", "liver"), ("testis",)])
def test_partial_tissue_set_is_refused(tissues: tuple[str, ...]) -> None:
    """A subset configuration loads only part of the attested weights, so refuse it.

    Pangolin's tissue set picks both the weight files loaded and the output channel
    read, so a gate run at one tissue attests nothing about the others -- yet the
    attestation necessarily records the full 12-file map.
    """
    with pytest.raises(AttestationError, match="tissue"):
        verified_predictor(PangolinSplicePredictor(tissues=tissues), _passing_attestation())


def test_refusal_is_not_a_silent_downgrade() -> None:
    """A mismatch raises; it never returns a quietly-uncalibrated predictor."""
    predictor = PangolinSplicePredictor(tissues=("brain",))
    with pytest.raises(AttestationError):
        verified_predictor(predictor, _passing_attestation())
    assert predictor.calibrated is False


# --------------------------------------------------------------------------
# 4. A calibrated head that steers delivery must enter the provenance stamp


def test_two_calibrated_heads_cannot_share_a_stamp() -> None:
    """Invariant #9: one stamp must not map to two different delivered sequences."""
    frontier = api.frontier(_PROTEIN, OptimizeConfig(), 9)
    a = api.rerank_by_expression(frontier, _CalibratedHead("A", flip=False))
    b = api.rerank_by_expression(frontier, _CalibratedHead("B", flip=True))

    assert a.delivered() is not None
    assert b.delivered() is not None
    assert a.delivered().dna != b.delivered().dna, "heads must actually steer differently"
    assert a.manifest.stamp != b.manifest.stamp


def test_a_steering_head_changes_the_stamp_versus_no_rerank() -> None:
    """A reranked delivery must not carry the un-reranked run's stamp."""
    frontier = api.frontier(_PROTEIN, OptimizeConfig(), 9)
    reranked = api.rerank_by_expression(frontier, _CalibratedHead("A", flip=False))
    assert reranked.manifest.stamp != frontier.manifest.stamp
    assert reranked.manifest.extra["expression_model"] == "stub-head-A"


def test_uncalibrated_rerank_is_still_a_pure_no_op() -> None:
    """The placeholder steers nothing, so it must not perturb the stamp either.

    This is the other half of the honesty rule: an uncalibrated score neither moves
    delivery nor claims a place in provenance.
    """
    frontier = api.frontier(_PROTEIN, OptimizeConfig(), 9)
    annotated = api.rerank_by_expression(frontier)
    assert annotated.frontier.chosen == frontier.frontier.chosen
    assert annotated.manifest.stamp == frontier.manifest.stamp
    assert "expression_model" not in annotated.manifest.extra


def test_the_same_head_twice_is_reproducible() -> None:
    """Determinism (#7): the stamp is a function of the head, not of the call."""
    frontier = api.frontier(_PROTEIN, OptimizeConfig(), 9)
    one = api.rerank_by_expression(frontier, _CalibratedHead("A", flip=False))
    two = api.rerank_by_expression(frontier, _CalibratedHead("A", flip=False))
    assert one.manifest.stamp == two.manifest.stamp
    assert one.delivered().dna == two.delivered().dna
