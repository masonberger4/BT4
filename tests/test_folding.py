"""Tests for the folding models behind the ``FoldingModel`` contract.

Coverage:

* the runtime-checkable contract surface (``name``/``calibrated``/
  ``five_prime_dg``/``score_sequence``);
* the baseline's honesty (``calibrated is False``, name screams baseline) and
  determinism (dependency-free, seedless, repeatable);
* :func:`default` returning a working model without ViennaRNA installed (CI has
  no bindings) and never crashing;
* the fixed ``score_sequence`` orientation -- a more-structured 5' end scores
  worse (smaller) than a weakly-structured one;
* the ViennaRNA path, gated behind ``pytest.importorskip`` so CI without the
  bindings still passes.
"""

from __future__ import annotations

import dataclasses

import pytest

from bt4.biomodels.folding import (
    DEFAULT_FIVE_PRIME_WINDOW,
    BaselinePairingProxyModel,
    FoldingModel,
    FoldingResult,
    ViennaFoldingModel,
    default,
    five_prime_window,
)

# A strongly self-complementary 5' end (GC stem-loop) vs a weakly-structured one
# (no complementary partners at all).
STRUCTURED = "GGGGGGGGGGTTTTTCCCCCCCCCC"
WEAK = "AAAAAAAAAAAAAAAAAAAAAAAAA"


def test_baseline_contract_surface() -> None:
    model = BaselinePairingProxyModel()
    assert isinstance(model, FoldingModel)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)
    assert isinstance(model.five_prime_dg("ATGGCCGGC"), float)
    assert isinstance(model.score_sequence("ATGGCCGGC"), float)


def test_baseline_is_labeled_uncalibrated() -> None:
    model = BaselinePairingProxyModel()
    assert model.calibrated is False
    assert model.name == "baseline-pairing-proxy"
    assert "baseline" in model.name


def test_non_model_is_not_a_folding_model() -> None:
    assert not isinstance(object(), FoldingModel)


def test_baseline_determinism() -> None:
    seq = "ATGGCGCGCGCGCTTTTGCGCGCGCGCTAA"
    a = BaselinePairingProxyModel()
    b = BaselinePairingProxyModel()
    first = a.five_prime_dg(seq)
    # Repeated calls and a fresh instance agree exactly (no global RNG, no time).
    assert a.five_prime_dg(seq) == first
    assert b.five_prime_dg(seq) == first
    assert a.score_sequence(seq) == b.score_sequence(seq)


def test_baseline_orientation_more_structure_is_worse() -> None:
    model = BaselinePairingProxyModel()
    # A stable stem-loop yields a more-negative proxy deltaG...
    assert model.five_prime_dg(STRUCTURED) < model.five_prime_dg(WEAK)
    # ...and therefore a smaller (worse) score than the open sequence.
    assert model.score_sequence(STRUCTURED) < model.score_sequence(WEAK)
    # The unstructured window has no pairs at all.
    assert model.five_prime_dg(WEAK) == 0.0
    assert model.five_prime_dg(STRUCTURED) < 0.0


def test_baseline_window_slicing_and_validation() -> None:
    model = BaselinePairingProxyModel()
    # A 5-nt 5' window of the stem-loop is all G -- no pairs -- so it is higher
    # (less structured) than folding the whole self-complementary sequence.
    assert model.five_prime_dg(STRUCTURED, window=5) > model.five_prime_dg(STRUCTURED)
    assert model.five_prime_dg(STRUCTURED, window=5) == 0.0
    with pytest.raises(ValueError):
        model.five_prime_dg(STRUCTURED, window=0)
    with pytest.raises(ValueError):
        model.five_prime_dg(STRUCTURED, window=-4)
    with pytest.raises(ValueError):
        BaselinePairingProxyModel(five_prime_window=0)


def test_five_prime_window_helper() -> None:
    assert five_prime_window("atgGCC", None) == "ATGGCC"
    assert five_prime_window("ATGGCCTAA", 3) == "ATG"
    # A window longer than the sequence yields the whole sequence.
    assert five_prime_window("ATG", 999) == "ATG"
    with pytest.raises(ValueError):
        five_prime_window("ATG", 0)


def test_folding_result_is_frozen() -> None:
    model = BaselinePairingProxyModel()
    result = model.fold(STRUCTURED)
    assert isinstance(result, FoldingResult)
    assert result.dg == model.five_prime_dg(STRUCTURED)
    assert result.model_name == "baseline-pairing-proxy"
    assert result.calibrated is False
    assert result.structure is None
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.dg = 1.0  # type: ignore[misc]


def test_default_window_constant() -> None:
    assert DEFAULT_FIVE_PRIME_WINDOW == 48
    assert BaselinePairingProxyModel().five_prime_window == DEFAULT_FIVE_PRIME_WINDOW


def test_available_never_raises() -> None:
    # Must return a bool regardless of whether the bindings are present.
    assert isinstance(ViennaFoldingModel.available(), bool)


def test_default_returns_working_model() -> None:
    model = default()
    assert isinstance(model, FoldingModel)
    assert isinstance(model.name, str)
    assert isinstance(model.calibrated, bool)
    assert isinstance(model.five_prime_dg("ATGGCCGGC"), float)
    assert isinstance(model.score_sequence("ATGGCCGGC"), float)
    # Without the ViennaRNA bindings, default() must fall back to the honest,
    # uncalibrated baseline (and it must never crash getting there).
    if not ViennaFoldingModel.available():
        assert isinstance(model, BaselinePairingProxyModel)
        assert model.calibrated is False


def test_vienna_path() -> None:
    pytest.importorskip("RNA")
    model = ViennaFoldingModel()
    assert isinstance(model, FoldingModel)
    assert model.calibrated is True
    assert model.name == "viennarna-mfe"
    dg = model.five_prime_dg(STRUCTURED)
    assert isinstance(dg, float)
    # Real thermodynamics must agree with the fixed orientation: the stem-loop
    # is at least as structured (deltaG no larger) as the open sequence.
    assert model.score_sequence(STRUCTURED) <= model.score_sequence(WEAK)
    result = model.fold(STRUCTURED)
    assert isinstance(result, FoldingResult)
    assert result.structure is not None
    assert result.calibrated is True
