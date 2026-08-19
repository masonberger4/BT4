"""Tests for the cross-backend splice agreement report.

Covers the pure correlation primitives (``pearson`` / ``spearman``, including
ties and zero-variance) and :func:`backend_agreement`'s ranking / sign-agreement
logic, driven by controllable fake backends plus the real PWM baseline for one
realistic end-to-end check.
"""

from __future__ import annotations

import pytest

# The correlation primitives now live in the shared bt4.biomodels._stats; the
# splice agreement layer re-exports pearson/spearman for its public surface.
from bt4.biomodels._stats import _ranks
from bt4.biomodels.splice import (
    ConsensusPwmSplicePredictor,
    SpliceResult,
    backend_agreement,
    spearman,
)
from bt4.biomodels.splice.agreement import agreement_from_deltas, pearson
from bt4.biomodels.splice.base import pooled_risk


def test_pearson_basic() -> None:
    assert pearson([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == pytest.approx(1.0)
    assert pearson([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]) == pytest.approx(-1.0)
    # Zero variance is undefined -> honest 0.0.
    assert pearson([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0


def test_pearson_validation() -> None:
    with pytest.raises(ValueError):
        pearson([1.0, 2.0], [1.0])
    with pytest.raises(ValueError):
        pearson([1.0], [1.0])


def test_spearman_monotonic_nonlinear() -> None:
    # A perfectly monotone but nonlinear relation has Spearman 1.0...
    assert spearman([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0]) == pytest.approx(1.0)
    # ...while Pearson is < 1 for the same data.
    assert pearson([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0]) < 1.0
    assert spearman([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]) == pytest.approx(-1.0)


def test_spearman_with_ties() -> None:
    # Tie-averaged ranks keep Spearman well defined.
    rho = spearman([1.0, 1.0, 2.0, 3.0], [5.0, 5.0, 6.0, 7.0])
    assert rho == pytest.approx(1.0)


def test_ranks_are_tie_averaged() -> None:
    # Ties share the average of the ranks they span (1-based).
    assert _ranks([1.0, 1.0, 2.0]) == [1.5, 1.5, 3.0]
    assert _ranks([3.0, 1.0, 1.0]) == [3.0, 1.5, 1.5]
    assert _ranks([5.0, 5.0, 5.0, 5.0]) == [2.5, 2.5, 2.5, 2.5]


def test_spearman_asymmetric_ties_guards_averaging() -> None:
    # Asymmetric ties: x ties where y does not. A broken _ranks that assigned
    # sequential (non-averaged) ranks within a tie block would give 1.0 here;
    # correct tie-averaging gives sqrt(3)/2.
    assert spearman([1.0, 1.0, 2.0], [1.0, 2.0, 3.0]) == pytest.approx(0.8660254037844387)


class _FakeBackend:
    """A :class:`SplicePredictor` whose per-sequence risk is fully controllable.

    ``backend_agreement`` now computes Delta-splicing directly from
    ``score_sequence`` (the contract's ``pooled_risk(ref) - pooled_risk(cand)``,
    scoring the reference once), so a fake must supply real per-position scores,
    not just a delta. Each sequence maps to a single donor-site probability, so
    ``pooled_risk(seq) == max(0, logit(p))`` and the delta of a candidate vs the
    reference is exactly controllable through the two probabilities. The ``top_k``
    attribute is what the ``top_k=None`` path reads for this backend's depth.
    """

    def __init__(self, name: str, risk_by_seq: dict[str, float], top_k: int = 3) -> None:
        self._name = name
        self._risk_by_seq = risk_by_seq
        self.top_k = top_k

    @property
    def name(self) -> str:
        return self._name

    @property
    def calibrated(self) -> bool:
        return False

    def score_sequence(self, dna: str) -> SpliceResult:
        p = self._risk_by_seq[dna.upper()]
        return SpliceResult(
            donor=(p,),
            acceptor=(0.0,),
            model_name=self._name,
            calibrated=False,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        return pooled_risk(self.score_sequence(reference_dna), self.top_k) - pooled_risk(
            self.score_sequence(designed_dna), self.top_k
        )


# Reference and three candidates keyed by name; probabilities set per test.
_CANDS = ["X", "Y", "Z"]
_REF = "REF"


def test_agreement_perfect_rank() -> None:
    # Both backends rank candidates by risk identically (ref riskier than all).
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    b = _FakeBackend("b", {"REF": 0.95, "X": 0.55, "Y": 0.65, "Z": 0.75})
    report = backend_agreement([a, b], _CANDS, _REF)
    assert report.backends == ("a", "b")
    # All candidates are less risky than the reference -> all positive deltas.
    assert all(d > 0 for d in report.delta_by_backend["a"])
    # Risk increases X<Y<Z, so delta (ref - cand) decreases X>Y>Z for both.
    assert report.delta_by_backend["a"][0] > report.delta_by_backend["a"][2]
    assert report.rank_correlations[("a", "b")] == pytest.approx(1.0)
    assert report.sign_agreement == pytest.approx(1.0)
    assert report.n_candidates == 3


def test_agreement_opposite_rank() -> None:
    # a: ref riskiest, candidate risk increases X<Y<Z -> deltas positive, X>Y>Z.
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    # b: ref at background (risk 0), candidate risk decreases X>Y>Z -> deltas
    # negative, and ordered X<Y<Z (opposite ranking to a).
    b = _FakeBackend("b", {"REF": 0.5, "X": 0.8, "Y": 0.7, "Z": 0.6})
    report = backend_agreement([a, b], _CANDS, _REF)
    assert report.rank_correlations[("a", "b")] == pytest.approx(-1.0)
    # a says all positive, b says all negative -> never agree on sign.
    assert all(d > 0 for d in report.delta_by_backend["a"])
    assert all(d < 0 for d in report.delta_by_backend["b"])
    assert report.sign_agreement == pytest.approx(0.0)


def test_agreement_partial_sign() -> None:
    # a signs [+, -, +]; b signs [+, +, -] -> agree on candidate 0 only -> 1/3.
    a = _FakeBackend("a", {"REF": 0.7, "X": 0.6, "Y": 0.8, "Z": 0.55})
    b = _FakeBackend("b", {"REF": 0.7, "X": 0.6, "Y": 0.55, "Z": 0.8})
    report = backend_agreement([a, b], _CANDS, _REF)
    assert report.sign_agreement == pytest.approx(1.0 / 3.0)


def test_single_backend_has_full_sign_agreement_and_no_correlations() -> None:
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    report = backend_agreement([a], _CANDS, _REF)
    assert report.sign_agreement == 1.0
    assert report.rank_correlations == {}


def test_fewer_than_two_candidates_has_no_correlations() -> None:
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6})
    b = _FakeBackend("b", {"REF": 0.9, "X": 0.6})
    report = backend_agreement([a, b], ["X"], _REF)
    assert report.rank_correlations == {}  # rank correlation is undefined for n=1
    assert report.n_candidates == 1


def test_agreement_validation() -> None:
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    with pytest.raises(ValueError):
        backend_agreement([], _CANDS, _REF)
    with pytest.raises(ValueError):
        backend_agreement([a], [], _REF)
    # Distinct names required.
    dup = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    with pytest.raises(ValueError):
        backend_agreement([a, dup], _CANDS, _REF)


def test_top_k_override_pools_uniformly() -> None:
    # With an explicit top_k, the reference is scored once and every backend is
    # pooled at that depth. Ref riskier than all -> positive, decreasing deltas.
    a = _FakeBackend("a", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8}, top_k=1)
    report = backend_agreement([a], _CANDS, _REF, top_k=2)
    deltas = report.delta_by_backend["a"]
    assert deltas[0] > deltas[1] > deltas[2]
    assert all(d > 0 for d in deltas)


class _TwoTrackBackend:
    """A SpliceAI-shaped backend: BOTH donor and acceptor tracks populated.

    Each sequence maps to (acceptor_prob, donor_prob) placed as one site each, so
    pooled_risk pools two distinct sites -- exercising backend_agreement across a
    backend whose SpliceResult shape differs from the donor-only fake.
    """

    def __init__(
        self,
        name: str,
        risk_by_seq: dict[str, tuple[float, float]],
        top_k: int = 3,
    ) -> None:
        self._name = name
        self._risk_by_seq = risk_by_seq
        self.top_k = top_k

    @property
    def name(self) -> str:
        return self._name

    @property
    def calibrated(self) -> bool:
        return False

    def score_sequence(self, dna: str) -> SpliceResult:
        acc, don = self._risk_by_seq[dna.upper()]
        return SpliceResult(donor=(don,), acceptor=(acc,), model_name=self._name, calibrated=False)

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        return pooled_risk(self.score_sequence(reference_dna), self.top_k) - pooled_risk(
            self.score_sequence(designed_dna), self.top_k
        )


def test_agreement_is_shape_agnostic_across_backends() -> None:
    # A donor-only (Pangolin-shaped) backend and a both-tracks (SpliceAI-shaped)
    # backend compare fine: backend_agreement works at the pooled-risk level, so
    # the SpliceResult layout does not matter.
    donor_only = _FakeBackend("pangolin-shaped", {"REF": 0.9, "X": 0.6, "Y": 0.7, "Z": 0.8})
    two_track = _TwoTrackBackend(
        "spliceai-shaped",
        {"REF": (0.85, 0.9), "X": (0.55, 0.6), "Y": (0.65, 0.7), "Z": (0.75, 0.8)},
    )
    report = backend_agreement([donor_only, two_track], _CANDS, _REF)
    assert set(report.backends) == {"pangolin-shaped", "spliceai-shaped"}
    # Both rank candidates by risk the same way (X<Y<Z), so they agree.
    assert report.rank_correlations[("pangolin-shaped", "spliceai-shaped")] == pytest.approx(1.0)
    # The two-track backend pools two sites, so its deltas differ in magnitude
    # from the donor-only one -- agreement is about ranking, not identical values.
    assert report.delta_by_backend["spliceai-shaped"] != report.delta_by_backend["pangolin-shaped"]


def test_real_baseline_backend_runs() -> None:
    # A realistic end-to-end pass with the actual PWM baseline as one backend.
    baseline = ConsensusPwmSplicePredictor()
    carrier = "A" * 20
    donor = "AAGGTAAGA"
    reference = carrier + carrier
    clean = carrier + carrier
    risky = carrier + donor + "A" * 11
    fake = _FakeBackend("fake", {clean.upper(): 0.6, risky.upper(): 0.8, reference.upper(): 0.6})
    report = backend_agreement([baseline, fake], [clean, risky], reference)
    assert set(report.backends) == {"consensus-pwm-baseline", "fake"}
    assert len(report.delta_by_backend["consensus-pwm-baseline"]) == 2
    # The identical-to-reference candidate has delta 0 for the baseline.
    assert report.delta_by_backend["consensus-pwm-baseline"][0] == pytest.approx(0.0)
    # Adding a canonical donor raises risk -> negative baseline delta.
    assert report.delta_by_backend["consensus-pwm-baseline"][1] < 0.0


# --------------------------------------------------------------------------
# A backend whose scores never clear the pooling background


class _SubBackgroundPredictor:
    """Scores every position below 0.5, responding to sequence content — like Pangolin.

    Measured on the designed-CDS panel with the hash-verified weights, Pangolin's peak
    per-position score was 0.323 to 0.445 and varied more than twofold between a native
    CDS and its synonymous redesigns. Every one of those pooled to a risk of exactly
    zero. This stub reproduces that regime without needing the licensed weights: the
    scores differ per sequence, and none of them reaches the background.
    """

    name = "sub-background-stub"
    calibrated = False
    top_k = 3

    def score_sequence(self, dna: str) -> SpliceResult:
        # Content-dependent and strictly inside (0, 0.5): GC fraction scaled down.
        gc = sum(1 for base in dna.upper() if base in "GC") / max(len(dna), 1)
        peak = 0.05 + 0.4 * gc
        return SpliceResult(
            donor=(peak, peak * 0.9, peak * 0.8),
            acceptor=(0.0, 0.0, 0.0),
            model_name=self.name,
            calibrated=False,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        return pooled_risk(self.score_sequence(reference_dna), self.top_k) - pooled_risk(
            self.score_sequence(designed_dna), self.top_k
        )


def test_a_sub_background_backend_reports_zero_deltas_and_says_why() -> None:
    """The defect's shape: real per-sequence differences, reported as none."""
    predictor = _SubBackgroundPredictor()
    # Distinct GC fractions, so this stub's scores genuinely differ per candidate.
    candidates = ["ATGGCCGCCGCC", "ATGAAATTTAAA", "ATGGCCAAATTT"]
    report = backend_agreement([predictor], candidates, "ATGAAAGGGTTT")

    # The shipped risk deltas are all exactly zero...
    assert report.delta_by_backend[predictor.name] == (0.0, 0.0, 0.0)
    # ...and the report now says that is the pooling, not the model.
    assert report.n_sub_background[predictor.name] == len(candidates) + 1
    assert report.degenerate(predictor.name) is True
    # The peak that was discarded is carried, so "0.001" and "0.44" stay distinct.
    assert 0.0 < report.max_score_by_backend[predictor.name] < 0.5
    # And the background-free response does separate the candidates.
    responses = report.response_by_backend[predictor.name]
    assert len(set(responses)) == len(candidates)


def test_a_backend_above_background_is_not_reported_as_degenerate() -> None:
    """The diagnostic must not fire on the case it exists to distinguish from.

    The shipped PWM baseline scores a real donor consensus well above 0.5, so its
    pooled risk is a measurement and its zeros — where it has any — are real.
    """
    predictor = ConsensusPwmSplicePredictor()
    with_site = "ATGGTAAGTACCGGCGTAAGTGCC"
    report = backend_agreement([predictor], [with_site, "ATGAAACCCTTT"], "ATGAAACCCTTT")
    assert report.n_sub_background[predictor.name] < report.n_candidates + 1
    assert report.degenerate(predictor.name) is False


def test_degenerate_is_false_when_the_diagnostic_was_not_collected() -> None:
    """Deltas supplied precomputed carry no evidence — and absence is not soundness."""
    report = agreement_from_deltas({"a": (0.0, 0.0), "b": (0.0, 0.0)})
    assert report.n_sub_background == {}
    assert report.degenerate("a") is False
