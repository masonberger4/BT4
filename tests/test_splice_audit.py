"""Tests for the localize-and-flag splice audit (design step 4).

These use fake backends that return controlled per-position tracks, so the
localization (peak/NMS), the honest labeling (per-flag calibrated, all_calibrated,
combined-track kind), the two opposing sign conventions (per-flag added-risk vs
panel-level delta), the approximate cross-backend co-occurrence, and the
no-editing guarantee are all pinned without torch / a CNN checkout.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from bt4.biomodels.splice import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.audit import SpliceAuditReport, audit_splice
from bt4.biomodels.splice.base import SpliceResult, pooled_risk


@dataclass(frozen=True)
class _FakeSplice:
    """A splice backend returning caller-controlled donor/acceptor tracks."""

    name: str
    by_seq: dict[str, tuple[tuple[float, ...], tuple[float, ...]]]
    calibrated: bool = False
    top_k: int = 3

    def score_sequence(self, dna: str) -> SpliceResult:
        donor, acceptor = self.by_seq[dna]
        return SpliceResult(
            donor=donor, acceptor=acceptor, model_name=self.name, calibrated=self.calibrated
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        return pooled_risk(self.score_sequence(reference_dna), self.top_k) - pooled_risk(
            self.score_sequence(designed_dna), self.top_k
        )


def test_localization_is_peak_nms() -> None:
    # A broad above-threshold run yields ONE flag at its peak, not a cluster.
    donor = (0.0, 0.9, 0.95, 0.9, 0.0, 0.6, 0.0)  # run [1..3] peak@2; single @5
    acceptor = (0.0,) * 7
    ref = ("REF", ((0.0,) * 7, (0.0,) * 7))
    cand = ("CAND", (donor, acceptor))
    backend = _FakeSplice("fake", dict([ref, cand]))
    rep = audit_splice([backend], ["CAND"], "REF", threshold=0.5)
    flags = rep.candidates[0].by_backend[0].flags
    positions = [(f.position, f.kind) for f in flags]
    assert positions == [(2, "splice"), (5, "splice")]  # acceptor all-zero => combined => "splice"


def test_separated_donor_acceptor_kinds() -> None:
    donor = (0.0, 0.8, 0.0, 0.0)
    acceptor = (0.0, 0.0, 0.9, 0.0)
    seqs = {"R": ((0.0,) * 4, (0.0,) * 4), "C": (donor, acceptor)}
    backend = _FakeSplice("two-track", seqs)
    rep = audit_splice([backend], ["C"], "R", threshold=0.5)
    flags = rep.candidates[0].by_backend[0].flags
    assert [(f.position, f.kind) for f in flags] == [(1, "donor"), (2, "acceptor")]


def test_added_risk_sign_is_positive_worse_and_intra_backend() -> None:
    # Candidate scores higher than reference at the flagged position => added risk > 0.
    seqs = {"R": ((0.0, 0.30, 0.0), (0.0,) * 3), "C": ((0.0, 0.90, 0.0), (0.0,) * 3)}
    backend = _FakeSplice("b", seqs)
    rep = audit_splice([backend], ["C"], "R", threshold=0.5)
    flag = rep.candidates[0].by_backend[0].flags[0]
    assert flag.position == 1
    assert flag.added_risk_vs_reference == pytest.approx(0.90 - 0.30)  # positive = worse


def test_delta_splicing_is_larger_is_better() -> None:
    # Candidate has LESS pooled risk than the reference => delta_splicing > 0 (better).
    seqs = {"R": ((0.0, 0.99, 0.0), (0.0,) * 3), "C": ((0.0, 0.10, 0.0), (0.0,) * 3)}
    backend = _FakeSplice("b", seqs)
    rep = audit_splice([backend], ["C"], "R", threshold=0.5)
    ba = rep.candidates[0].by_backend[0]
    assert ba.delta_splicing > 0.0  # opposite sign convention to per-flag added_risk


def test_all_calibrated_false_and_flags_carry_backend_flag() -> None:
    seqs = {"R": ((0.0, 0.0), (0.0, 0.0)), "C": ((0.0, 0.9), (0.0, 0.0))}
    uncal = _FakeSplice("uncal", seqs, calibrated=False)
    rep = audit_splice([uncal], ["C"], "R", threshold=0.5)
    assert rep.all_calibrated is False
    assert all(f.calibrated is False for f in rep.candidates[0].by_backend[0].flags)


def test_cross_backend_co_occurrence_windowed() -> None:
    # Backend A flags at pos 10; backend B at pos 12; window 3 => they co-flag.
    def _peak(at: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
        return tuple(0.9 if i == at else 0.0 for i in range(20)), (0.0,) * 20

    a_seqs = {"R": ((0.0,) * 20, (0.0,) * 20), "C": _peak(10)}
    b_seqs = {"R": ((0.0,) * 20, (0.0,) * 20), "C": _peak(12)}
    a = _FakeSplice("A", a_seqs)
    b = _FakeSplice("B", b_seqs)
    rep = audit_splice([a, b], ["C"], "R", threshold=0.5, match_window=3)
    fa = rep.candidates[0].by_backend[0].flags[0]
    fb = rep.candidates[0].by_backend[1].flags[0]
    assert fa.also_flagged_by == ("B",)
    assert fb.also_flagged_by == ("A",)
    # Outside the window there is no co-occurrence.
    rep_narrow = audit_splice([a, b], ["C"], "R", threshold=0.5, match_window=1)
    assert rep_narrow.candidates[0].by_backend[0].flags[0].also_flagged_by == ()


def test_audit_never_edits_and_attaches_agreement() -> None:
    seqs = {
        "R": ((0.0, 0.0), (0.0, 0.0)),
        "C1": ((0.0, 0.9), (0.0, 0.0)),
        "C2": ((0.0, 0.1), (0.0, 0.0)),
    }
    backend = _FakeSplice("b", seqs)
    rep = audit_splice([backend], ["C1", "C2"], "R", threshold=0.5)
    assert [c.dna for c in rep.candidates] == ["C1", "C2"]  # unchanged, in order
    assert rep.agreement.backends == ("b",)
    assert rep.agreement.n_candidates == 2


def test_validation() -> None:
    seqs = {"R": ((0.0,), (0.0,)), "C": ((0.9,), (0.0,))}
    b = _FakeSplice("b", seqs)
    with pytest.raises(ValueError, match="at least one splice backend"):
        audit_splice([], ["C"], "R")
    with pytest.raises(ValueError, match="at least one candidate"):
        audit_splice([b], [], "R")
    with pytest.raises(ValueError, match="distinct"):
        audit_splice([b, _FakeSplice("b", seqs)], ["C"], "R")
    with pytest.raises(ValueError, match="match_window"):
        audit_splice([b], ["C"], "R", match_window=-1)
    with pytest.raises(ValueError, match="top_k"):
        audit_splice([b], ["C"], "R", top_k=0)


def test_determinism() -> None:
    seqs = {"R": ((0.0, 0.0, 0.0), (0.0,) * 3), "C": ((0.0, 0.9, 0.0), (0.0,) * 3)}
    b = _FakeSplice("b", seqs)
    a1 = audit_splice([b], ["C"], "R", threshold=0.5)
    a2 = audit_splice([b], ["C"], "R", threshold=0.5)
    assert a1 == a2


def test_baseline_backend_runs_and_is_advisory() -> None:
    # The real (uncalibrated) baseline works end-to-end and stays advisory.
    ref = "ATG" + "GCT" * 20 + "TAA"
    cand = "ATG" + "GCC" * 20 + "TAA"
    rep = audit_splice([ConsensusPwmSplicePredictor()], [cand], ref, threshold=0.5)
    assert rep.all_calibrated is False
    assert isinstance(rep, SpliceAuditReport)
    assert rep.candidates[0].by_backend[0].calibrated is False


# --------------------------------------------------------------------------- #
# Pipeline adapter over a step-3 CandidateSet.
# --------------------------------------------------------------------------- #


def test_audit_candidate_set_defaults_reference_to_delivered() -> None:
    from bt4.api import candidates
    from bt4.pipeline.optimize import OptimizeConfig
    from bt4.pipeline.splice_audit import audit_candidate_set, available_splice_backends

    cs = candidates("MKTAYIAKQRQISFVKSHFSRQLE", OptimizeConfig(), n=5)
    rep = audit_candidate_set(cs)
    # Reference defaults to the delivered candidate, whose added risk vs itself is 0.
    delivered_dna = cs.delivered().result.dna
    delivered_audit = next(c for c in rep.candidates if c.dna == delivered_dna)
    for ba in delivered_audit.by_backend:
        assert ba.delta_splicing == pytest.approx(0.0)
    # Only the honest baseline is available in CI (no CNN install).
    assert [b.name for b in available_splice_backends()] == ["consensus-pwm-baseline"]
    assert rep.all_calibrated is False


def test_audit_candidate_set_empty_raises() -> None:
    from dataclasses import replace

    from bt4.api import candidates
    from bt4.pipeline.optimize import OptimizeConfig
    from bt4.pipeline.splice_audit import audit_candidate_set

    cs = candidates("MKT", OptimizeConfig(), n=3)
    empty = replace(cs, candidates=(), chosen=0)
    with pytest.raises(ValueError, match="empty"):
        audit_candidate_set(empty)
