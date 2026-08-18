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
    ExpressionResult,
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


# --- Batched scoring (score_many / delta_logte_many) --------------------------
#
# These pin the batch surface without a live model: the real _predict_te needs
# torch + the checkout, so the wrapping/ordering tests stub it on the class (the
# adapter is a frozen, slotted dataclass, so patch the method on the class, which
# monkeypatch restores), and the guard-clause tests drive the real _predict_te far
# enough to hit its validation, exactly as the single-sequence tests above do.


def test_score_many_empty_returns_empty() -> None:
    # An empty candidate set is a no-op: no validation, no RiboNN invocation.
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    assert model.score_many([]) == []


def test_score_many_preserves_order_and_wraps(monkeypatch: pytest.MonkeyPatch) -> None:
    # score_many wraps each _predict_te value in an ExpressionResult, in input
    # order, carrying the honest calibrated=False and the CLR-residual units.
    scores = {"ATGTTTTGA": 0.25, "ATGGCCTAA": 1.5, "ATGAAATAA": -2.0}

    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        return [scores[d] for d in dnas]  # keyed by content => order-sensitive

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    dnas = ["ATGGCCTAA", "ATGAAATAA", "ATGTTTTGA"]  # not the dict's order
    results = model.score_many(dnas)

    assert [r.score for r in results] == [1.5, -2.0, 0.25]  # input order, not dict
    assert all(isinstance(r, ExpressionResult) for r in results)
    assert all(r.calibrated is False for r in results)
    assert all(r.model_name == "ribonn[human]" for r in results)
    assert all("CLR-residual" in r.units for r in results)


def test_score_sequence_delegates_to_batch(monkeypatch: pytest.MonkeyPatch) -> None:
    # score_sequence is now a one-element score_many; it must still return a single
    # ExpressionResult with the batched value.
    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        assert dnas == ["ATGGCCTAA"]
        return [3.14]

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)
    result = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT").score_sequence("ATGGCCTAA")
    assert result.score == 3.14
    assert result.calibrated is False


def test_score_many_validates_each_dna(monkeypatch: pytest.MonkeyPatch) -> None:
    # Every input is validate_dna'd before any RiboNN work; a bad base rejects the
    # whole batch (validation fires before _predict_te is reached).
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    with pytest.raises(ValueError, match="non-ACGT"):
        model.score_many(["ATGGCCTAA", "ATGXYZTAA"])


def test_score_many_rejects_bad_cds(monkeypatch: pytest.MonkeyPatch) -> None:
    # The per-input length-3N + stop-codon guard still fires inside the batch path.
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    with pytest.raises(ValueError, match=r"3N ending in a stop codon"):
        model.score_many(["ATGGCCTAA", "ATGGCCGCC"])  # second lacks a stop codon


def test_score_many_rejects_empty_utr(monkeypatch: pytest.MonkeyPatch) -> None:
    # The non-empty-UTR guard still fires for the batch path (empty by default).
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel()
    with pytest.raises(ValueError, match=r"UTR"):
        model.score_many(["ATGGCCTAA"])


def test_delta_logte_many_empty_returns_empty() -> None:
    # No designs => no deltas, and the reference is not even scored.
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    assert model.delta_logte_many([], "ATGGCCTAA") == []


def test_delta_logte_many_scores_reference_once(monkeypatch: pytest.MonkeyPatch) -> None:
    # The reference is a single shared baseline: scored ONCE (appended to the batch),
    # not once per design, in a single _predict_te invocation.
    calls: list[list[str]] = []

    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        calls.append(list(dnas))
        return [float(len(d)) for d in dnas]  # deterministic, per-input

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    designed = ["ATGAAATAA", "ATGTTTTTTTAA"]  # lengths 9 and 12
    reference = "ATGGGGTGA"  # length 9, distinct from every design

    deltas = model.delta_logte_many(designed, reference)

    assert len(calls) == 1  # a single RiboNN invocation for the whole set
    assert calls[0] == [*designed, reference]  # reference appended last...
    assert calls[0].count(reference) == 1  # ...exactly once
    assert deltas == [9.0 - 9.0, 12.0 - 9.0]  # TE(design) - TE(reference), in order


def test_delta_logte_delegates_to_batch(monkeypatch: pytest.MonkeyPatch) -> None:
    # delta_logte is now a one-design delta_logte_many; the scalar must match.
    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        assert dnas == ["ATGAAATAA", "ATGGGGTGA"]  # [designed, reference]
        return [5.0, 2.0]

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    assert model.delta_logte("ATGAAATAA", "ATGGGGTGA") == 3.0


def test_delta_logte_validates_design_before_reference(monkeypatch: pytest.MonkeyPatch) -> None:
    # The score_sequence/delta_logte delegations must stay behavior-preserving: when
    # both inputs are invalid DNA, the design's error surfaces first (as it did before
    # the batch refactor), so delta_logte_many validates designs before the reference.
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    with pytest.raises(ValueError, match=r"\['Z'\]"):  # the design's bad char, not the reference's
        model.delta_logte("ATGZZZTAA", "GGGQQQTAA")


def test_delta_logte_many_rejects_empty_utr(monkeypatch: pytest.MonkeyPatch) -> None:
    # Guard clauses still fire for the batched delta path too.
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel()
    with pytest.raises(ValueError, match=r"UTR"):
        model.delta_logte_many(["ATGGCCTAA"], "ATGGGGTGA")


def test_reduce_te_preserves_input_order_over_scrambled_df() -> None:
    # The reducer keys on tx_id and returns results in *input* (tx_id) order even
    # when RiboNN's output rows arrive scrambled and multi-row-per-input (ensemble),
    # so a batched score_many keeps its inputs aligned to their outputs.
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out_df = pd.DataFrame(
        {
            "tx_id": ["bt4_2", "bt4_0", "bt4_1", "bt4_2", "bt4_0", "bt4_1"],
            "predicted_TE_cA": [30.0, 10.0, 20.0, 32.0, 12.0, 22.0],
            "predicted_TE_cB": [30.0, 10.0, 20.0, 28.0, 8.0, 18.0],
        }
    )
    # Per-row cell-type means then the ensemble mean per tx_id:
    #   bt4_0: mean(10,10)=10, mean(12,8)=10  -> 10.0
    #   bt4_1: mean(20,20)=20, mean(22,18)=20 -> 20.0
    #   bt4_2: mean(30,30)=30, mean(32,28)=30 -> 30.0
    assert _reduce_te_by_tx_id(out_df, ["bt4_0", "bt4_1", "bt4_2"]) == [10.0, 20.0, 30.0]


# --- batch_size / num_workers passthrough -------------------------------------
#
# RiboNN's predict_using_nested_cross_validation_models exposes batch_size (its
# default 1024) and num_workers (its default 4). BT4 forwards both and defaults
# them down, because both upstream defaults are hostile in this adapter:
# num_workers>0 spawns workers that re-import without the mutated sys.path or the
# temporary cwd this adapter scores from (so they hang or fail wherever the start
# method is spawn -- Windows, macOS), and batch_size=1024 allocates 1024
# fixed-width (channels, 13318) float32 tensors at once. Neither can change a
# score: RiboNN pads to a fixed width, not to a batch's longest member, and its
# predict dataloader is built with shuffle=False.


def test_batch_knob_defaults_are_the_safe_ones() -> None:
    model = RiboNNExpressionModel()
    assert model.batch_size == 64  # below RiboNN's OOM-prone 1024
    assert model.num_workers == 0  # required wherever the start method is spawn


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"batch_size": 0}, "batch_size must be >= 1"),
        ({"batch_size": -1}, "batch_size must be >= 1"),
        ({"num_workers": -1}, "num_workers must be >= 0"),
    ],
)
def test_batch_knobs_are_validated(kwargs: dict[str, int], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        RiboNNExpressionModel(**kwargs)


def test_predict_layout_forwards_batch_knobs(tmp_path: object) -> None:
    # The load-bearing test for the Windows fix: whatever the adapter is configured
    # with must actually reach RiboNN's predict call, not be silently dropped.
    from pathlib import Path

    from bt4.biomodels.expression.ribonn import _run_predict_with_models_layout

    captured: dict[str, object] = {}

    def fake_predict(**kwargs: object) -> str:
        captured.update(kwargs)
        return "out"

    # Named "models" so the no-symlink branch runs (symlinks need Developer Mode
    # or an elevated prompt on Windows, which is exactly what we must not require).
    weights = Path(str(tmp_path)) / "models"
    (weights / "human").mkdir(parents=True)

    result = _run_predict_with_models_layout(
        fake_predict,
        weights,
        "in.tsv",
        "human",
        run_df="RUNS",
        top_k=3,
        batch_size=8,
        num_workers=0,
    )

    assert result == "out"
    assert captured["batch_size"] == 8
    assert captured["num_workers"] == 0
    assert captured["top_k_models_to_use"] == 3
    assert captured["species"] == "human"


def test_predict_layout_restores_cwd_even_on_failure(tmp_path: object) -> None:
    # The helper chdirs into the weights parent; a raising predict must not leave
    # the process in a different directory (it would break every later relative path).
    import os
    from pathlib import Path

    from bt4.biomodels.expression.ribonn import _run_predict_with_models_layout

    weights = Path(str(tmp_path)) / "models"
    (weights / "human").mkdir(parents=True)
    before = os.getcwd()

    def boom(**_: object) -> None:
        raise RuntimeError("upstream failed")

    with pytest.raises(RuntimeError, match="upstream failed"):
        _run_predict_with_models_layout(
            boom, weights, "in.tsv", "human", run_df=None, top_k=1,
            batch_size=1, num_workers=0,
        )
    assert os.getcwd() == before


def test_resolve_backend_threads_batch_knobs() -> None:
    from bt4.biomodels.expression import resolve_backend

    model = resolve_backend(
        "ribonn", utr5="GCCACC", utr3="GCTAAT", batch_size=16, num_workers=2
    )
    assert isinstance(model, RiboNNExpressionModel)
    assert (model.batch_size, model.num_workers) == (16, 2)
    assert model.calibrated is False  # a knob is not a calibration claim


# --- cell-type selection and fold semantics -----------------------------------
#
# RiboNN emits one row per input per outer fold (10), each already that fold's
# top_k-model mean, across 78 human / 68 mouse per-cell-type columns. Averaging
# every fold is right for a novel designed sequence and wrong for a natural one
# (nine folds trained on its label); averaging all 78 cell types is right for a
# generic design and wrong when the ground truth came from one cell line.


def _fake_out_df(pd: object) -> object:
    # Two inputs x two folds x two cell types. Per-row cell-type means are chosen
    # so every reduction below has an exactly-representable expected value.
    return pd.DataFrame(  # type: ignore[attr-defined]
        {
            "tx_id": ["bt4_0", "bt4_1", "bt4_0", "bt4_1"],
            "fold": [0, 0, 1, 1],
            "predicted_TE_HEK293T": [1.0, 3.0, 5.0, 7.0],
            "predicted_TE_HeLa": [3.0, 9.0, 11.0, 17.0],
        }
    )


def test_fold_reduction_keeps_fold_identity() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id_and_fold

    records = _reduce_te_by_tx_id_and_fold(_fake_out_df(pd), ["bt4_0", "bt4_1"])
    # sorted by (index, fold); te = mean over both cell-type columns
    assert [(r.index, r.fold, r.te) for r in records] == [
        (0, 0, 2.0),  # mean(1, 3)
        (0, 1, 8.0),  # mean(5, 11)
        (1, 0, 6.0),  # mean(3, 9)
        (1, 1, 12.0),  # mean(7, 17)
    ]


def test_fold_averaged_summary_is_the_mean_of_the_folds() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    # input 0: mean(2, 8) = 5; input 1: mean(6, 12) = 9
    assert _reduce_te_by_tx_id(_fake_out_df(pd), ["bt4_0", "bt4_1"]) == [5.0, 9.0]


def test_cell_type_selection_changes_the_number() -> None:
    # The point of the knob: selecting one cell line is a different number from the
    # mean of all of them, so it must not silently fall back to "all".
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out = _fake_out_df(pd)
    both = _reduce_te_by_tx_id(out, ["bt4_0", "bt4_1"])
    hek = _reduce_te_by_tx_id(out, ["bt4_0", "bt4_1"], ("HEK293T",))
    assert hek == [3.0, 5.0]  # input 0: mean(1, 5); input 1: mean(3, 7)
    assert hek != both


def test_unknown_cell_type_raises_and_lists_the_options() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    with pytest.raises(ValueError, match=r"no output for cell type\(s\) \['K562'\]"):
        _reduce_te_by_tx_id(_fake_out_df(pd), ["bt4_0", "bt4_1"], ("K562",))
    # and the message names what IS available, so the fix is obvious
    with pytest.raises(ValueError, match="HEK293T"):
        _reduce_te_by_tx_id(_fake_out_df(pd), ["bt4_0", "bt4_1"], ("K562",))


def test_missing_input_row_still_raises_the_length_cap_hint() -> None:
    pd = pytest.importorskip("pandas")
    from bt4.biomodels.expression.ribonn import _reduce_te_by_tx_id

    out = pd.DataFrame({"tx_id": ["bt4_0"], "fold": [0], "predicted_TE_a": [1.0]})
    with pytest.raises(ValueError, match="no prediction for input 1"):
        _reduce_te_by_tx_id(out, ["bt4_0", "bt4_1"])


def test_cell_types_are_validated() -> None:
    with pytest.raises(ValueError, match="must not contain blank names"):
        RiboNNExpressionModel(cell_types=("HEK293T", " "))
    with pytest.raises(ValueError, match="contains duplicates"):
        RiboNNExpressionModel(cell_types=("HEK293T", "HEK293T"))


def test_units_name_the_cell_type_selection(monkeypatch: pytest.MonkeyPatch) -> None:
    # The units label travels into reports and the manifest, so "all human cell
    # types" and "HEK293T only" must never share one.
    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        return [0.0] * len(dnas)

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)
    base = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    scoped = dataclasses.replace(base, cell_types=("HEK293T",))

    assert "all human cell types" in base.score_many(["ATGTAA"])[0].units
    assert scoped.score_many(["ATGTAA"])[0].units.endswith("(mean over HEK293T)")


def test_predict_folds_empty_is_a_no_op() -> None:
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    assert model.predict_folds([]) == []


def test_predict_folds_validates_dna(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("BT4_RIBONN_DIR", raising=False)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    with pytest.raises(ValueError, match="non-ACGT"):
        model.predict_folds(["ATGTAA", "ATGXAA"])


def test_predict_folds_shares_one_invocation(monkeypatch: pytest.MonkeyPatch) -> None:
    # predict_folds and _predict_te must be two views of ONE scoring path, so a
    # fold-resolved read cannot drift from the fold-averaged one.
    pd = pytest.importorskip("pandas")
    calls: list[list[str]] = []

    def fake_run(self: RiboNNExpressionModel, dnas: list[str]) -> tuple[object, list[str]]:
        calls.append(list(dnas))
        return _fake_out_df(pd), [f"bt4_{i}" for i in range(len(dnas))]

    monkeypatch.setattr(RiboNNExpressionModel, "_run_predict", fake_run)
    model = RiboNNExpressionModel(utr5="GCCACC", utr3="GCTAAT")
    dnas = ["ATGTAA", "ATGTGA"]

    per_fold = model.predict_folds(dnas)
    averaged = [r.score for r in model.score_many(dnas)]

    assert len(calls) == 2  # one invocation each; neither loops per sequence
    # The averaged view is exactly the mean of the fold-resolved view.
    for i, score in enumerate(averaged):
        folds = [r.te for r in per_fold if r.index == i]
        assert score == pytest.approx(sum(folds) / len(folds))
