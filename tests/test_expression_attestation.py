"""Tests for the expression attestation layer -- the only seam that flips ``calibrated``.

Two properties matter more than the rest:

1. **Licence-cleanliness is structural.** RiboNN's weights are Sanofi non-commercial, so
   its raw per-sequence outputs must never enter MIT-licensed BT4. The dataclass shape
   *is* that contract, and a drift must fail rather than be caught at review time.
2. **Scope is part of the claim.** An attestation earned within-protein, on human
   weights, scoring one cell line, certifies exactly that. A pooled run, a different
   species, or a different cell-type selection is a different quantity, and the promotion
   must refuse rather than quietly transfer the claim.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib

import pytest

from bt4.biomodels.expression import (
    MAX_ATTESTATION_COVERAGE_TOLERANCE,
    MAX_ATTESTATION_WIDTH_OVER_IQR,
    MIN_ATTESTATION_SPEARMAN,
    ExpressionAttestation,
    ExpressionAttestationError,
    NullExpressionModel,
    PanelRow,
    RiboNNExpressionModel,
    attest_expression,
    load_expression_attestation,
    panel_from_rows,
    verified_predictor,
)
from bt4.biomodels.expression.attestation import _ALLOWED_FIELDS, _pinned_weights
from bt4.pipeline.expression_gate import GateSettings, run_panel_gate


def _panel(n_groups: int = 10, n_variants: int = 6) -> object:
    """A panel whose within-protein driver is the variant index, not codon content."""
    rows = []
    for g in range(n_groups):
        for v in range(n_variants):
            body = ("AAA" * (v + 1)) + ("AAG" * (n_variants - v))
            rows.append(
                PanelRow(
                    group=f"P{g:02d}",
                    variant_id=f"g{g}v{v}",
                    cds="ATG" + body + "TAA",
                    measured=100.0 * g + v,
                    utr5="GCCACC",
                    utr3="GCTAAT",
                )
            )
    return panel_from_rows(rows)


def _passing_comparison() -> object:
    """An oracle head, gated within-group with a tolerance loose enough to pass."""
    panel = _panel()
    oracle = [row.measured for row in panel.rows]  # type: ignore[attr-defined]
    return run_panel_gate(
        panel,  # type: ignore[arg-type]
        "null",
        settings=GateSettings(
            within_group=True,
            recalibrate=True,
            coverage_tolerance=MAX_ATTESTATION_COVERAGE_TOLERANCE,
            bootstrap_resamples=200,
        ),
        head_scores=oracle,
    )


def _attestation() -> ExpressionAttestation:
    return attest_expression(
        _passing_comparison(),
        species="human",
        cell_types=("HEK293T",),
        readout="mean_ribosome_load",
        bt4_version="0.0.0-test",
    )


# --- licence-cleanliness is structural ----------------------------------------


def test_the_dataclass_shape_is_the_licence_clean_contract() -> None:
    assert {f.name for f in dataclasses.fields(ExpressionAttestation)} == _ALLOWED_FIELDS


def test_no_field_can_carry_a_raw_model_score() -> None:
    # RiboNN's per-sequence outputs are licence-encumbered. Nothing plural or
    # array-shaped that could hold them belongs in this record.
    forbidden = ("scores", "predictions", "per_sequence", "raw", "expected", "cds")
    for name in _ALLOWED_FIELDS:
        assert not any(token in name for token in forbidden), name


def test_from_dict_refuses_a_smuggled_field() -> None:
    payload = _attestation().to_dict()
    payload["per_sequence_scores"] = [1.0, 2.0, 3.0]
    with pytest.raises(ExpressionAttestationError, match="unexpected attestation field"):
        ExpressionAttestation.from_dict(payload)


def test_from_dict_refuses_a_truncated_record() -> None:
    payload = _attestation().to_dict()
    del payload["panel_sha256"]
    with pytest.raises(ExpressionAttestationError, match="missing field"):
        ExpressionAttestation.from_dict(payload)


# --- round trip and hashing ---------------------------------------------------


def test_round_trips_through_json(tmp_path: pathlib.Path) -> None:
    attestation = _attestation()
    path = tmp_path / "attestation.json"
    path.write_text(json.dumps(attestation.to_dict(), indent=2), encoding="utf-8")
    assert load_expression_attestation(path) == attestation


def test_content_hash_is_stable_and_value_sensitive() -> None:
    first = _attestation()
    assert first.content_hash() == _attestation().content_hash()  # no wall-clock
    changed = dataclasses.replace(first, panel_sha256="0" * 64)
    assert changed.content_hash() != first.content_hash()


def test_content_hash_changes_with_the_scope() -> None:
    # The scope IS the claim, so two attestations differing only in cell type must not
    # share a provenance stamp.
    base = _attestation()
    other = dataclasses.replace(base, cell_types=("HeLa",))
    assert other.content_hash() != base.content_hash()


# --- what may be attested -----------------------------------------------------


def test_a_passing_within_group_run_can_be_attested() -> None:
    attestation = _attestation()
    assert attestation.passed is True
    assert attestation.within_group is True
    assert attestation.readout == "mean_ribosome_load"
    assert len(attestation.weight_sha256) == 90  # the species' pinned files


def test_a_pooled_run_cannot_be_attested_however_good_it_looks() -> None:
    # The refusal that protects the whole exercise: pooled scoring credits
    # between-protein skill, which is not the regime BT4 deploys in.
    panel = _panel()
    oracle = [row.measured for row in panel.rows]  # type: ignore[attr-defined]
    pooled = run_panel_gate(
        panel,  # type: ignore[arg-type]
        "null",
        settings=GateSettings(
            within_group=False, recalibrate=True, coverage_tolerance=0.10,
            bootstrap_resamples=100,
        ),
        head_scores=oracle,
    )
    with pytest.raises(ExpressionAttestationError, match="refusing to attest a POOLED"):
        attest_expression(
            pooled, readout="mean_ribosome_load", bt4_version="0.0.0-test"
        )


def test_a_head_that_only_ties_a_baseline_cannot_be_attested() -> None:
    # The head is a perfect oracle AND the panel's within-protein driver is exactly GC3,
    # so the gc3 baseline ties it. The gate passes; the claim is still refused, because
    # a head that merely reproduces a feature BT4 already computes has earned nothing.
    from bt4.pipeline.expression_gate import _gc3

    rows = []
    for g in range(10):
        for v in range(6):
            body = ("AAA" * (v + 1)) + ("AAG" * (6 - v))
            cds = "ATG" + body + "TAA"
            rows.append(
                PanelRow(
                    group=f"P{g:02d}",
                    variant_id=f"g{g}v{v}",
                    cds=cds,
                    measured=100.0 * g + _gc3(cds),  # within-protein driver IS GC3
                    utr5="GCCACC",
                    utr3="GCTAAT",
                )
            )
    panel = panel_from_rows(rows)
    oracle = [row.measured for row in panel.rows]

    comparison = run_panel_gate(
        panel,
        "null",
        settings=GateSettings(
            within_group=True, recalibrate=True, coverage_tolerance=0.10,
            bootstrap_resamples=100,
        ),
        head_scores=oracle,
    )
    assert comparison.head.passed is True  # the thresholds are met ...
    assert comparison.beats_every_baseline is False  # ... but nothing was added
    with pytest.raises(ExpressionAttestationError, match="does not beat every baseline"):
        attest_expression(comparison, readout="mrl", bt4_version="0.0.0-test")


def test_a_failing_gate_cannot_be_attested() -> None:
    panel = _panel()
    blind = [1.0] * len(panel.rows)  # type: ignore[attr-defined]
    comparison = run_panel_gate(
        panel,  # type: ignore[arg-type]
        "null",
        settings=GateSettings(within_group=True, bootstrap_resamples=50),
        head_scores=blind,
    )
    with pytest.raises(ExpressionAttestationError, match="failing gate"):
        attest_expression(comparison, readout="mrl", bt4_version="0.0.0-test")


def test_an_unknown_backend_is_refused() -> None:
    with pytest.raises(ExpressionAttestationError, match="unknown expression backend"):
        attest_expression(
            _passing_comparison(), backend="nope", readout="mrl", bt4_version="0.0.0"
        )


# --- promotion: the single seam -----------------------------------------------


def test_verified_predictor_flips_calibrated_on_an_exact_scope_match() -> None:
    attestation = _attestation()
    model = RiboNNExpressionModel(
        species="human", utr5="GCCACC", utr3="GCTAAT", cell_types=("HEK293T",)
    )
    assert model.calibrated is False

    promoted = verified_predictor(model, attestation)
    assert promoted.calibrated is True
    assert model.calibrated is False  # frozen dataclass: the original is untouched


def test_promotion_refuses_a_different_cell_type_selection() -> None:
    # An attestation earned on HEK293T does not certify a head averaging all 78 cell
    # types -- those are different quantities.
    attestation = _attestation()
    all_cell_types = RiboNNExpressionModel(species="human", utr5="A", utr3="C")
    with pytest.raises(ExpressionAttestationError, match="different quantity"):
        verified_predictor(all_cell_types, attestation)

    other_line = RiboNNExpressionModel(
        species="human", utr5="A", utr3="C", cell_types=("HeLa",)
    )
    with pytest.raises(ExpressionAttestationError, match="attestation covers cell types"):
        verified_predictor(other_line, attestation)


def test_promotion_refuses_a_different_species() -> None:
    attestation = _attestation()
    mouse = RiboNNExpressionModel(
        species="mouse", utr5="A", utr3="C", cell_types=("HEK293T",)
    )
    with pytest.raises(ExpressionAttestationError, match="covers species 'human'"):
        verified_predictor(mouse, attestation)


def test_promotion_refuses_mismatched_weight_hashes() -> None:
    attestation = dataclasses.replace(
        _attestation(), weight_sha256=(("human/fake/state_dict.pth", "0" * 64),)
    )
    model = RiboNNExpressionModel(
        species="human", utr5="A", utr3="C", cell_types=("HEK293T",)
    )
    with pytest.raises(ExpressionAttestationError, match="weight hashes do not match"):
        verified_predictor(model, attestation)


def test_promotion_refuses_a_hand_edited_pooled_or_slack_attestation() -> None:
    # The floors are re-checked at promotion, so editing the JSON after the fact does
    # not buy a calibration.
    model = RiboNNExpressionModel(
        species="human", utr5="A", utr3="C", cell_types=("HEK293T",)
    )
    base = _attestation()

    for edited, match in (
        (dataclasses.replace(base, within_group=False), "POOLED"),
        (dataclasses.replace(base, passed=False), "did not pass"),
        (
            dataclasses.replace(base, min_spearman=MIN_ATTESTATION_SPEARMAN - 0.01),
            "below the floor",
        ),
        (
            dataclasses.replace(
                base, coverage_tolerance=MAX_ATTESTATION_COVERAGE_TOLERANCE + 0.01
            ),
            "looser than the floor",
        ),
        (
            dataclasses.replace(
                base, width_over_iqr=MAX_ATTESTATION_WIDTH_OVER_IQR + 0.1
            ),
            "vacuous",
        ),
    ):
        with pytest.raises(ExpressionAttestationError, match=match):
            verified_predictor(model, edited)


def test_promotion_refuses_a_non_attestable_predictor() -> None:
    with pytest.raises(ExpressionAttestationError, match="not an attestable"):
        verified_predictor(NullExpressionModel(), _attestation())


def test_the_default_head_is_still_the_uncalibrated_placeholder() -> None:
    # Adding the seam must not change what BT4 hands out by default.
    from bt4.biomodels.expression import default

    assert default().calibrated is False
    assert RiboNNExpressionModel().calibrated is False


def test_pinned_weights_are_filtered_by_species() -> None:
    human = _pinned_weights("human")
    mouse = _pinned_weights("mouse")
    assert len(human) == len(mouse) == 90
    assert not set(human) & set(mouse)
