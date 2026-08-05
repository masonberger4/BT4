"""Tests for the wrapped RiboNN expression adapter.

These run **without** torch / pandas / the RiboNN checkout: the adapter is light at
import, the pinned-weight manifest is a bundled JSON, and everything the live
forward pass needs is gated behind :meth:`RiboNNExpressionModel.available`. They
pin the honesty guarantees:

* the adapter ships ``calibrated is False`` and only flips via a passing gate;
* the bundled weight manifest is 180 public SHA-256 hashes (90 human + 90 mouse);
* hash verification refuses missing / tampered weights before any load;
* it degrades honestly (``available() is False``) when the checkout is absent.
"""

from __future__ import annotations

import dataclasses
import re

import pytest

from bt4.biomodels.expression import (
    NullExpressionModel,
    RiboNNExpressionModel,
    default,
    load_pinned_sha256,
)
from bt4.biomodels.expression.ribonn import PINNED_WEIGHT_SHA256, _sha256_file

_HEX64 = re.compile(r"[0-9a-f]{64}")


def test_manifest_is_180_public_hashes() -> None:
    manifest = load_pinned_sha256()
    assert manifest == PINNED_WEIGHT_SHA256
    assert len(manifest) == 180
    assert sum(k.startswith("human/") for k in manifest) == 90
    assert sum(k.startswith("mouse/") for k in manifest) == 90
    for key, digest in manifest.items():
        assert key.endswith("/state_dict.pth")
        assert _HEX64.fullmatch(digest), key


def test_defaults_and_identity() -> None:
    model = RiboNNExpressionModel()
    assert model.species == "human"
    assert model.name == "ribonn[human]"
    assert model.calibrated is False  # never calibrated as shipped
    assert RiboNNExpressionModel(species="mouse").name == "ribonn[mouse]"


def test_argument_validation() -> None:
    with pytest.raises(ValueError, match="species"):
        RiboNNExpressionModel(species="rat")
    with pytest.raises(ValueError, match="top_k"):
        RiboNNExpressionModel(top_k=0)


def test_calibrated_mirrors_fidelity_verified() -> None:
    model = RiboNNExpressionModel()
    promoted = dataclasses.replace(model, fidelity_verified=True)
    assert promoted.calibrated is True  # the only way to flip it
    assert model.calibrated is False  # original untouched (frozen)


def test_unavailable_without_checkout(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    monkeypatch.delenv("BT4_RIBONN_WEIGHTS", raising=False)
    model = RiboNNExpressionModel()
    assert model.available() is False


def test_scoring_without_checkout_raises_clearly(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    # Non-empty UTRs so we get past the UTR guard and reach the missing-checkout error.
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    with pytest.raises(RuntimeError, match=r"BT4_RIBONN_DIR|clone not found"):
        model.score_sequence("ATGGCCTAA")


def test_scoring_rejects_empty_utr(monkeypatch: pytest.MonkeyPatch) -> None:
    # RiboNN's loader can't preprocess an all-empty UTR column (read back as NaN); the
    # adapter refuses up front with a clear message instead of a deep pandas crash.
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel()  # utr5/utr3 empty by default
    with pytest.raises(ValueError, match=r"UTR"):
        model.score_sequence("ATGGCCTAA")
    with pytest.raises(ValueError, match=r"UTR"):
        RiboNNExpressionModel(utr5="GCC").score_sequence("ATGGCCTAA")  # only 3' empty


def test_rejects_cds_not_ending_in_stop(monkeypatch: pytest.MonkeyPatch) -> None:
    # RiboNN hard-asserts CDS is 3N ending in a stop codon; the adapter refuses up
    # front (a pure input check, before touching the checkout).
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel()
    with pytest.raises(ValueError, match=r"3N ending in a stop codon"):
        model.score_sequence("ATGGCCGCC")  # 3N but no stop
    with pytest.raises(ValueError, match=r"3N ending in a stop codon"):
        model.score_sequence("ATGGCCTA")  # not length-3N


def test_verify_weights_rejects_missing(tmp_path) -> None:  # type: ignore[no-untyped-def]
    # An empty weights dir: the first pinned file is missing => refuse before load.
    (tmp_path / "human").mkdir()
    model = RiboNNExpressionModel(species="human")
    with pytest.raises(FileNotFoundError, match=r"state_dict\.pth"):
        model._verify_weights(tmp_path)


def test_verify_weights_rejects_tampered(tmp_path) -> None:  # type: ignore[no-untyped-def]
    # Populate every pinned human file with wrong bytes; the sha mismatch must fire
    # (not a false pass). We only need one to reach the ValueError branch, so give
    # every file identical junk and assert the digest-mismatch message.
    for rel in PINNED_WEIGHT_SHA256:
        if not rel.startswith("human/"):
            continue
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"not the real weights")
    model = RiboNNExpressionModel(species="human")
    with pytest.raises(ValueError, match=r"sha256|refusing to load"):
        model._verify_weights(tmp_path)


def test_sha256_file_roundtrip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    import hashlib

    blob = b"hello ribonn"
    path = tmp_path / "x.bin"
    path.write_bytes(blob)
    assert _sha256_file(path) == hashlib.sha256(blob).hexdigest()


def test_reduce_te_averages_ensemble_rows() -> None:
    # RiboNN returns one row per ensemble model, so a single input yields several rows
    # sharing a tx_id. The reducer must average them (mean over cell types, then over
    # the ensemble) into one scalar per input -- the bug was float() on a Series.
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out_df = pd.DataFrame(
        {
            "tx_id": ["bt4_0", "bt4_0", "bt4_1"],
            "predicted_TE_cellA": [1.0, 3.0, 10.0],
            "predicted_TE_cellB": [1.0, 3.0, 20.0],
        }
    )
    # bt4_0: row means [1.0, 3.0] -> ensemble mean 2.0; bt4_1: single row mean 15.0.
    assert _reduce_te_by_tx_id(out_df, ["bt4_0", "bt4_1"]) == [2.0, 15.0]


def test_reduce_te_missing_tx_id_raises() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out_df = pd.DataFrame({"tx_id": ["bt4_0"], "predicted_TE_x": [1.0]})
    with pytest.raises(ValueError, match=r"length cap|no prediction"):
        _reduce_te_by_tx_id(out_df, ["bt4_0", "bt4_1"])


def test_reduce_te_no_te_columns_raises() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out_df = pd.DataFrame({"tx_id": ["bt4_0"], "other_col": [1.0]})
    with pytest.raises(RuntimeError, match=r"predicted_TE"):
        _reduce_te_by_tx_id(out_df, ["bt4_0"])


def test_default_stays_null_placeholder() -> None:
    # No calibrated head ships, so default() is still the neutral placeholder.
    predictor = default()
    assert isinstance(predictor, NullExpressionModel)
    assert predictor.calibrated is False
