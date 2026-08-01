"""Tests for the splice predictors behind the ``SplicePredictor`` contract.

Coverage:

* the runtime-checkable contract surface (``name`` / ``calibrated`` /
  ``score_sequence`` / ``delta_splicing``);
* the baseline's honesty (``calibrated is False``, name screams baseline) and
  determinism (dependency-free, seedless, repeatable);
* :func:`default` returning a working, labeled-uncalibrated model that never
  crashes (no calibrated model ships yet);
* ``delta_splicing(seq, seq) == 0.0`` and the fixed larger-is-better orientation
  (adding a canonical donor lowers the objective);
* an introduced strong canonical ``GT..`` donor raising pooled splice risk; and
* top-k / log-odds pooling **not saturating** -- two sites give a strictly
  larger pooled risk than one, unlike the saturating noisy-OR aggregation.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from bt4.biomodels.splice import (
    DEFAULT_TOP_K,
    ConsensusPwmSplicePredictor,
    SplicePredictor,
    SpliceResult,
    default,
    logit,
    pool_log_odds,
    pooled_risk,
)

# A strong canonical U2 donor motif MAG|GTRAG (invariant GT at +1,+2), carrying a
# single GT so each copy is exactly one donor site. An A-rich carrier forms no
# site of its own (no GT for a donor, no pyrimidine-tract + AG for an acceptor),
# so a plain carrier of matched length is a genuine "no donor" comparator.
_DONOR_MOTIF = "AAGGTAAGA"
_CARRIER = "A" * 20


def _noisy_or(probs: tuple[float, ...]) -> float:
    """Saturating noisy-OR aggregate ``1 - prod(1 - p)`` (the BT3 anti-pattern)."""
    product = 1.0
    for p in probs:
        product *= 1.0 - p
    return 1.0 - product


def test_baseline_contract_surface() -> None:
    model = ConsensusPwmSplicePredictor()
    assert isinstance(model, SplicePredictor)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)
    result = model.score_sequence("ATGGCCGGCTAA")
    assert isinstance(result, SpliceResult)
    assert isinstance(model.delta_splicing("ATGGCCGGCTAA", "ATGGCTGGGTAA"), float)


def test_baseline_is_labeled_uncalibrated() -> None:
    model = ConsensusPwmSplicePredictor()
    assert model.calibrated is False
    assert model.name == "consensus-pwm-baseline"
    assert "baseline" in model.name


def test_non_model_is_not_a_splice_predictor() -> None:
    assert not isinstance(object(), SplicePredictor)


def test_score_arrays_align_to_sequence_and_are_probabilities() -> None:
    model = ConsensusPwmSplicePredictor()
    seq = _CARRIER + _DONOR_MOTIF + _CARRIER
    result = model.score_sequence(seq)
    assert len(result.donor) == len(seq)
    assert len(result.acceptor) == len(seq)
    assert all(0.0 <= p <= 1.0 for p in result.donor)
    assert all(0.0 <= p <= 1.0 for p in result.acceptor)
    assert result.model_name == "consensus-pwm-baseline"
    assert result.calibrated is False


def test_score_result_is_frozen() -> None:
    result = ConsensusPwmSplicePredictor().score_sequence(_CARRIER + _DONOR_MOTIF)
    assert isinstance(result, SpliceResult)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.calibrated = True  # type: ignore[misc]


def test_baseline_determinism() -> None:
    seq = _CARRIER + _DONOR_MOTIF + _CARRIER
    a = ConsensusPwmSplicePredictor()
    b = ConsensusPwmSplicePredictor()
    first = a.score_sequence(seq)
    # Repeated calls and a fresh instance agree exactly (no global RNG, no time).
    assert a.score_sequence(seq) == first
    assert b.score_sequence(seq) == first
    ref = "A" * len(seq)
    assert a.delta_splicing(seq, ref) == b.delta_splicing(seq, ref)


def test_case_insensitive_and_validation() -> None:
    model = ConsensusPwmSplicePredictor()
    seq = _CARRIER + _DONOR_MOTIF
    assert model.score_sequence(seq.lower()) == model.score_sequence(seq)
    with pytest.raises(ValueError):
        model.score_sequence("ATGX")
    with pytest.raises(ValueError):
        model.score_sequence("")
    with pytest.raises(ValueError):
        ConsensusPwmSplicePredictor(top_k=0)


def test_delta_splicing_of_sequence_vs_itself_is_zero() -> None:
    model = ConsensusPwmSplicePredictor()
    seq = _CARRIER + _DONOR_MOTIF + _CARRIER
    assert model.delta_splicing(seq, seq) == 0.0
    # And for a plain sequence with no strong sites.
    assert model.delta_splicing(_CARRIER, _CARRIER) == 0.0


def test_introduced_canonical_donor_raises_splice_risk() -> None:
    model = ConsensusPwmSplicePredictor()
    with_donor = _CARRIER + _DONOR_MOTIF + _CARRIER
    without_donor = "A" * len(with_donor)
    risk_with = pooled_risk(model.score_sequence(with_donor), model.top_k)
    risk_without = pooled_risk(model.score_sequence(without_donor), model.top_k)
    # A canonical GT.. donor is a real site; a site-free carrier is not.
    assert risk_with > risk_without
    assert risk_without == 0.0

    # The donor position itself scores as a strong site (near-1 pseudo-prob).
    donor_scores = model.score_sequence(with_donor).donor
    assert max(donor_scores) > 0.9

    # Orientation: designing IN a donor (vs a reference without) is worse, so the
    # larger-is-better objective is negative.
    assert model.delta_splicing(with_donor, without_donor) < 0.0
    # Removing a donor (designed cleaner than reference) is better -> positive.
    assert model.delta_splicing(without_donor, with_donor) > 0.0


def test_top_k_log_odds_pooling_does_not_saturate() -> None:
    model = ConsensusPwmSplicePredictor()
    one_site = _CARRIER + _DONOR_MOTIF + _CARRIER + _CARRIER
    two_sites = _CARRIER + _DONOR_MOTIF + _CARRIER + _DONOR_MOTIF + _CARRIER

    risk_one = pooled_risk(model.score_sequence(one_site))
    risk_two = pooled_risk(model.score_sequence(two_sites))

    # Log-odds pooling is additive: two strong sites carry clearly more risk than
    # one (roughly double), and it never pegs at a ceiling.
    assert risk_two > risk_one
    assert risk_two > 1.5 * risk_one

    # Contrast with the saturating noisy-OR the BT3 model used (CLAUDE.md 10.14):
    # both sequences peg near 1.0, so it cannot tell one strong site from two.
    r_one = model.score_sequence(one_site)
    r_two = model.score_sequence(two_sites)
    nor_one = _noisy_or((*r_one.donor, *r_one.acceptor))
    nor_two = _noisy_or((*r_two.donor, *r_two.acceptor))
    assert nor_one > 0.99
    assert nor_two > 0.99
    assert abs(nor_two - nor_one) < 1e-3


def test_pool_log_odds_primitives() -> None:
    # Empty and all-background inputs pool to zero risk.
    assert pool_log_odds([]) == 0.0
    assert pool_log_odds([0.001, 0.01, 0.2]) == 0.0
    # Two above-background sites sum additively; a third is truncated by top_k=2.
    two = pool_log_odds([0.9, 0.9, 0.9], top_k=2)
    one = pool_log_odds([0.9], top_k=2)
    assert two == pytest.approx(2.0 * one)
    assert two == pytest.approx(2.0 * logit(0.9))
    # top_k must be positive.
    with pytest.raises(ValueError):
        pool_log_odds([0.9], top_k=0)


def test_default_returns_working_labeled_model() -> None:
    model = default()
    assert isinstance(model, SplicePredictor)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)
    assert isinstance(model.score_sequence("ATGGCCGGC"), SpliceResult)
    assert isinstance(model.delta_splicing("ATGGCCGGC", "ATGGCGGGA"), float)
    # No calibrated model ships yet -> default() is the honest, uncalibrated
    # baseline and must never crash getting there.
    assert isinstance(model, ConsensusPwmSplicePredictor)
    assert model.calibrated is False


def test_default_top_k_constant() -> None:
    assert DEFAULT_TOP_K == 3
    assert ConsensusPwmSplicePredictor().top_k == DEFAULT_TOP_K


def test_logit_is_finite_at_extremes() -> None:
    assert math.isfinite(logit(0.0))
    assert math.isfinite(logit(1.0))
    assert logit(0.5) == pytest.approx(0.0)
