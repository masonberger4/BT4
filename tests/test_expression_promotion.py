"""Tests for the seam that lets an earned expression calibration reach a user.

Until this seam existed, ``verified_predictor`` had **no caller in** ``src/``: a
maintainer could run the acceptance gate, earn a claim and commit a record, and BT4
would behave identically -- the head stayed ``calibrated=False`` and the candidate set
stayed in discovery order. These tests pin the two halves of fixing that:

1. **Opting in changes something, and only opting in does.** With the switch off (the
   default) nothing about BT4's behaviour moves, and no attestation is bundled today, so
   the shipped default is unchanged either way.
2. **Promotion refuses rather than downgrades.** Every part of the scope that changes the
   number -- species, cell-type selection, ``top_k``, UTR context, the pinned weights --
   is bound. A mismatch raises; it never quietly hands back an uncalibrated head to a
   caller who asked for a calibrated one.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib

import pytest

from bt4.biomodels.expression import (
    ATTESTATION_PATH_ENV_VAR,
    MAX_ATTESTATION_COVERAGE_TOLERANCE,
    USE_ATTESTED_ENV_VAR,
    ExpressionAttestation,
    ExpressionAttestationError,
    ExpressionResult,
    NullExpressionModel,
    PanelRow,
    RiboNNExpressionModel,
    attest_expression,
    attested_expression_backends,
    attested_promotion_enabled,
    bundled_expression_attestation,
    bundled_expression_attestation_path,
    panel_from_rows,
    promote_if_attested,
    resolve_backend,
    resolve_expression_attestation,
    utr_context_sha256,
)
from bt4.pipeline.candidates import assemble_and_rank_candidates
from bt4.pipeline.expression_gate import GateSettings, run_panel_gate

PANEL_UTR5 = "GCCACC"
PANEL_UTR3 = "GCTAAT"
CELL_TYPES = ("HEK293T",)


def _panel(
    n_groups: int = 10,
    n_variants: int = 6,
    *,
    cell_type: str = "",
    species: str = "",
    readout: str = "",
) -> object:
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
                    utr5=PANEL_UTR5,
                    utr3=PANEL_UTR3,
                    cell_type=cell_type,
                    species=species,
                    readout=readout,
                )
            )
    return panel_from_rows(rows)


def _comparison(panel: object | None = None, **kwargs: object) -> object:
    """A passing within-group gate over an oracle head, scored in ``CELL_TYPES``."""
    panel = panel if panel is not None else _panel()
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
        cell_types=CELL_TYPES,
        head_scores=oracle,
        **kwargs,  # type: ignore[arg-type]
    )


def _as_ribonn_run(comparison: object, species: str = "human") -> object:
    """Relabel a placeholder-scored comparison as a gate-scored RiboNN run.

    A **test double**, and the only way to exercise promotion without the Sanofi
    non-commercial weights. Production cannot reach this state, which is exactly why the
    two fields it overrides exist: `attest_expression` derives the backend from the head
    the gate actually constructed, and refuses a run whose scores were handed in rather
    than computed.
    """
    return dataclasses.replace(
        comparison,  # type: ignore[type-var]
        backend=f"ribonn[{species}]",
        scope=dataclasses.replace(
            comparison.scope, scoring_source="gate"  # type: ignore[attr-defined]
        ),
    )


def _attestation(**kwargs: object) -> ExpressionAttestation:
    return attest_expression(
        _as_ribonn_run(_comparison()),
        readout="mean_ribosome_load",
        bt4_version="0.0.0-test",
        **kwargs,  # type: ignore[arg-type]
    )


def _head(**kwargs: object) -> RiboNNExpressionModel:
    """A RiboNN head configured exactly as the attestation's scope."""
    fields: dict[str, object] = {
        "species": "human",
        "utr5": PANEL_UTR5,
        "utr3": PANEL_UTR3,
        "cell_types": CELL_TYPES,
    }
    fields.update(kwargs)
    return RiboNNExpressionModel(**fields)  # type: ignore[arg-type]


def _write(attestation: ExpressionAttestation, path: pathlib.Path) -> pathlib.Path:
    path.write_text(json.dumps(attestation.to_dict(), indent=2), encoding="utf-8")
    return path


# --- nothing ships, and nothing moves without the opt-in ----------------------


def test_no_attestation_is_bundled_today() -> None:
    # The honest shipped state: no expression head has passed its gate, so committing a
    # record -- even a plausible-looking one -- would be the fabricated artifact
    # CLAUDE.md 10.6 forbids. When one is legitimately earned this test changes with it.
    assert bundled_expression_attestation("ribonn") is None
    assert not bundled_expression_attestation_path("ribonn").exists()
    assert attested_expression_backends() == ()


def test_promotion_is_off_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(USE_ATTESTED_ENV_VAR, raising=False)
    assert attested_promotion_enabled() is False
    assert promote_if_attested(_head()).calibrated is False


def test_the_opt_in_alone_changes_nothing_without_an_attestation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Asking for promotion when nothing resolves is a no-op, not an error and not a
    # calibration: the switch cannot manufacture a claim out of an empty data directory.
    monkeypatch.delenv(ATTESTATION_PATH_ENV_VAR, raising=False)
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, "1")
    assert attested_promotion_enabled() is True
    assert promote_if_attested(_head()).calibrated is False
    assert resolve_backend("ribonn", utr5="A", utr3="C").calibrated is False


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
def test_env_var_truthy_spellings(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, value)
    assert attested_promotion_enabled() is True


@pytest.mark.parametrize("value", ["", "0", "false", "no", "off", "maybe"])
def test_env_var_falsy_spellings(monkeypatch: pytest.MonkeyPatch, value: str) -> None:
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, value)
    assert attested_promotion_enabled() is False


def test_without_the_opt_in_a_present_attestation_is_ignored(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    # The whole point of opt-in: a resolvable record on disk still changes nothing.
    monkeypatch.delenv(USE_ATTESTED_ENV_VAR, raising=False)
    monkeypatch.setenv(
        ATTESTATION_PATH_ENV_VAR, str(_write(_attestation(), tmp_path / "a.json"))
    )
    assert promote_if_attested(_head()).calibrated is False
    assert resolve_backend(
        "ribonn", utr5=PANEL_UTR5, utr3=PANEL_UTR3, cell_types=CELL_TYPES
    ).calibrated is False


# --- resolution ---------------------------------------------------------------


def test_a_maintainers_own_attestation_resolves_from_the_env_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    # An expression attestation is earned against a maintainer's own (often
    # unpublished) panel, so using one must not require committing it.
    attestation = _attestation()
    monkeypatch.setenv(
        ATTESTATION_PATH_ENV_VAR, str(_write(attestation, tmp_path / "a.json"))
    )
    assert resolve_expression_attestation("ribonn") == attestation
    assert promote_if_attested(_head(), enabled=True).calibrated is True


def test_a_mispointed_attestation_path_refuses_rather_than_falling_back(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    # A typo here would otherwise look exactly like "no attestation", leaving a user
    # believing a promotion happened for a reason nobody could see.
    monkeypatch.setenv(ATTESTATION_PATH_ENV_VAR, str(tmp_path / "nope.json"))
    with pytest.raises(ExpressionAttestationError, match="not a readable file"):
        resolve_expression_attestation("ribonn")


def test_an_explicit_attestation_beats_the_environment(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    monkeypatch.setenv(ATTESTATION_PATH_ENV_VAR, str(tmp_path / "nope.json"))
    explicit = _attestation()
    assert resolve_expression_attestation("ribonn", attestation=explicit) == explicit


# --- promotion binds the whole scope -----------------------------------------


def test_the_opt_in_promotes_a_matching_head() -> None:
    promoted = promote_if_attested(_head(), enabled=True, attestation=_attestation())
    assert promoted.calibrated is True
    assert isinstance(promoted, RiboNNExpressionModel)
    assert _head().calibrated is False  # frozen: the original is untouched


def test_resolve_backend_promotes_under_an_explicit_opt_in(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    monkeypatch.setenv(
        ATTESTATION_PATH_ENV_VAR, str(_write(_attestation(), tmp_path / "a.json"))
    )
    backend = resolve_backend(
        "ribonn",
        utr5=PANEL_UTR5,
        utr3=PANEL_UTR3,
        cell_types=CELL_TYPES,
        use_attested=True,
    )
    assert backend.calibrated is True


def test_the_placeholder_can_never_be_promoted() -> None:
    # Not an error, and not a calibration: the neutral head is simply not attestable, so
    # it passes through untouched however loudly a caller opts in.
    assert (
        promote_if_attested(
            NullExpressionModel(), enabled=True, attestation=_attestation()
        ).calibrated
        is False
    )
    assert resolve_backend("null", use_attested=True).calibrated is False


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"top_k": 7}, r"earned at top_k=5, predictor ensembles top_k=7"),
        ({"utr5": "AAAAAA"}, "not one this attestation covers"),
        ({"utr3": "TTTTTT"}, "not one this attestation covers"),
        ({"species": "mouse"}, "covers species 'human', predictor is 'mouse'"),
        ({"cell_types": ("HeLa",)}, r"covers cell types \['HEK293T'\]"),
        ({"cell_types": ()}, "different quantity"),
    ],
)
def test_promotion_refuses_every_scope_mismatch(
    overrides: dict[str, object], match: str
) -> None:
    # Each of these changes the number the gate measured, so the gated head is not this
    # head. The refusal is loud: a silent downgrade would hand a caller who asked for a
    # calibrated ranking an uncalibrated one with no way to notice.
    with pytest.raises(ExpressionAttestationError, match=match):
        promote_if_attested(
            _head(**overrides), enabled=True, attestation=_attestation()
        )


def test_promotion_refuses_an_attestation_for_another_backend() -> None:
    # Reachable only from a hand-edited record, since attest_expression derives the
    # backend from the run -- which is exactly why promotion re-checks it.
    edited = dataclasses.replace(_attestation(), backend="something_else")
    with pytest.raises(ExpressionAttestationError, match="attestation is for"):
        promote_if_attested(_head(), enabled=True, attestation=edited)


def test_a_record_whose_backend_and_species_disagree_is_refused() -> None:
    # "ribonn[mouse]" scored, "human" recorded: the record would name neither.
    comparison = _as_ribonn_run(_comparison(), species="mouse")
    with pytest.raises(ExpressionAttestationError, match="but recorded species"):
        attest_expression(comparison, readout="mrl", bt4_version="0.0.0-test")


def test_a_subclassed_head_is_promoted_rather_than_silently_skipped() -> None:
    # promote_if_attested's gate and verified_predictor's must agree: a class-name
    # comparison would skip a subclass here while the seam accepted it, which is a
    # silent downgrade in the one place that must never produce one.
    class _Subclassed(RiboNNExpressionModel):
        pass

    head = _Subclassed(
        species="human", utr5=PANEL_UTR5, utr3=PANEL_UTR3, cell_types=CELL_TYPES
    )
    promoted = promote_if_attested(head, enabled=True, attestation=_attestation())
    assert promoted.calibrated is True


def test_a_reordered_record_is_the_same_record() -> None:
    # content_hash is a *content* hash: reordering a JSON list must not move it, and
    # verified_predictor compares these tuples exactly.
    payload = _attestation().to_dict()
    payload["verified_against_panel"] = list(
        reversed(payload["verified_against_panel"])
    )
    payload["weight_sha256"] = list(reversed(payload["weight_sha256"]))
    assert (
        ExpressionAttestation.from_dict(payload).content_hash()
        == _attestation().content_hash()
    )


def test_promotion_refuses_mismatched_weight_hashes() -> None:
    edited = dataclasses.replace(
        _attestation(), weight_sha256=(("human/fake/state_dict.pth", "0" * 64),)
    )
    with pytest.raises(ExpressionAttestationError, match="weight hashes do not match"):
        promote_if_attested(_head(), enabled=True, attestation=edited)


def test_batch_and_worker_knobs_are_deliberately_not_bound() -> None:
    # RiboNN pads to a fixed width and does not shuffle when predicting, so neither knob
    # can change a score. Binding them would refuse a head that is provably the gated
    # one -- false precision, not extra rigour.
    promoted = promote_if_attested(
        _head(batch_size=8, num_workers=2), enabled=True, attestation=_attestation()
    )
    assert promoted.calibrated is True


# --- the scope recorded is the scope that ran --------------------------------


def test_the_record_carries_the_configuration_the_run_used() -> None:
    attestation = _attestation()
    assert attestation.cell_types == CELL_TYPES
    assert attestation.top_k == 5
    assert attestation.utr_context_sha256 == (
        utr_context_sha256(PANEL_UTR5, PANEL_UTR3),
    )
    assert attestation.scoring_source == "gate"
    assert attestation.backend == "ribonn"


def test_a_declared_cell_type_the_gate_did_not_score_is_refused() -> None:
    # The sharpest edge in the old procedure: run the gate across all 78 cell types,
    # then declare HEK293T, and every later check accepts the lie. It is now a refusal.
    panel = _panel()
    comparison = _as_ribonn_run(
        run_panel_gate(
            panel,  # type: ignore[arg-type]
            "null",
            settings=GateSettings(
                within_group=True,
                recalibrate=True,
                coverage_tolerance=MAX_ATTESTATION_COVERAGE_TOLERANCE,
                bootstrap_resamples=200,
            ),
            head_scores=[row.measured for row in panel.rows],  # type: ignore[attr-defined]
        )
    )
    # Built OUTSIDE the raises block: a fixture that raised the matching error would
    # make this pass while `attest_expression` did nothing at all. And it must succeed
    # when nothing is declared, so the refusal below is attributable to the declaration.
    assert attest_expression(
        comparison, readout="mrl", bt4_version="0.0.0-test"
    ).cell_types == ()

    with pytest.raises(
        ExpressionAttestationError, match=r"declared cell types \['HEK293T'\]"
    ):
        attest_expression(
            comparison,
            cell_types=("HEK293T",),
            readout="mrl",
            bt4_version="0.0.0-test",
        )


def test_a_declared_species_that_did_not_run_is_refused() -> None:
    with pytest.raises(ExpressionAttestationError, match="declared species"):
        _attestation(species="mouse")


def test_the_panel_verifies_the_scope_when_it_declares_one() -> None:
    panel = _panel(cell_type="HEK293T", species="human", readout="mean_ribosome_load")
    attestation = attest_expression(
        _as_ribonn_run(_comparison(panel)), bt4_version="0.0.0-test"
    )
    # All three were checked against the panel's own bytes, not taken on trust.
    assert attestation.verified_against_panel == ("cell_types", "readout", "species")
    assert attestation.readout == "mean_ribosome_load"


def test_a_scope_the_panel_cannot_confirm_is_recorded_as_merely_declared() -> None:
    # An honest gap, not a silent pass: a reader can tell which half of a scope had a
    # second check and which was the maintainer's word.
    assert _attestation().verified_against_panel == ()


def test_a_readout_the_panel_does_not_measure_is_refused() -> None:
    panel = _panel(readout="mean_ribosome_load")
    with pytest.raises(ExpressionAttestationError, match="not what the panel measures"):
        attest_expression(
            _as_ribonn_run(_comparison(panel)),
            readout="protein_yield",
            bt4_version="0",
        )


def test_a_run_with_no_readout_anywhere_is_refused() -> None:
    with pytest.raises(ExpressionAttestationError, match="must name the question"):
        attest_expression(_as_ribonn_run(_comparison()), bt4_version="0.0.0-test")


def test_the_gate_refuses_a_head_that_would_average_the_wrong_cell_types() -> None:
    # Checked before the scoring pass, because the gate is a run-once procedure: paying
    # for a full RiboNN pass and only then discovering the scope error is the failure.
    panel = _panel(cell_type="HEK293T")
    with pytest.raises(ValueError, match="scope error"):
        run_panel_gate(
            panel,  # type: ignore[arg-type]
            "null",
            settings=GateSettings(within_group=True, bootstrap_resamples=20),
            head_scores=[row.measured for row in panel.rows],  # type: ignore[attr-defined]
        )


def test_the_gate_refuses_a_species_the_panel_contradicts() -> None:
    panel = _panel(species="mouse")
    with pytest.raises(ValueError, match="panel declares species"):
        run_panel_gate(
            panel,  # type: ignore[arg-type]
            "null",
            settings=GateSettings(within_group=True, bootstrap_resamples=20),
            head_scores=[row.measured for row in panel.rows],  # type: ignore[attr-defined]
        )


def test_the_gate_never_scores_with_an_already_promoted_head(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    # The gate is what DECIDES promotion, so a standing opt-in must not colour the run
    # that judges the next one.
    monkeypatch.setenv(USE_ATTESTED_ENV_VAR, "1")
    monkeypatch.setenv(
        ATTESTATION_PATH_ENV_VAR, str(_write(_attestation(), tmp_path / "a.json"))
    )
    comparison = _comparison()
    assert comparison.backend_calibrated is False  # type: ignore[attr-defined]


def test_a_panel_that_mixes_assays_cannot_be_filed_under_one_readout() -> None:
    # Rows from the other assay were still scored, so one label would name neither.
    rows = list(_panel().rows)  # type: ignore[attr-defined]
    mixed = panel_from_rows(
        [
            dataclasses.replace(
                row, readout="mean_ribosome_load" if i % 2 else "protein_yield"
            )
            for i, row in enumerate(rows)
        ]
    )
    with pytest.raises(ExpressionAttestationError, match="mixes readouts"):
        attest_expression(
            _as_ribonn_run(_comparison(mixed)),
            readout="mean_ribosome_load",
            bt4_version="0.0.0-test",
        )


def test_an_explicit_opt_in_that_cannot_be_fulfilled_refuses(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # `enabled=True` is a request about this call. Answering it with a silently
    # uncalibrated head is the failure the layer exists to prevent.
    monkeypatch.delenv(ATTESTATION_PATH_ENV_VAR, raising=False)
    with pytest.raises(ExpressionAttestationError, match="no attestation resolves"):
        promote_if_attested(_head(), enabled=True)
    # The placeholder is simply not attestable, so it still passes through untouched.
    assert promote_if_attested(NullExpressionModel(), enabled=True).calibrated is False


def test_a_corrupt_attestation_refuses_rather_than_reading_as_absent(
    monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path
) -> None:
    # A bare JSONDecodeError escapes every handler that catches "unusable attestation",
    # which is how a hand-edited file took the desktop app down at startup.
    path = tmp_path / "a.json"
    path.write_text("{not json", encoding="utf-8")
    monkeypatch.setenv(ATTESTATION_PATH_ENV_VAR, str(path))
    with pytest.raises(ExpressionAttestationError, match="not valid JSON"):
        resolve_expression_attestation("ribonn")

    # And the "is this offerable" probe reports it as absent rather than propagating --
    # forced past the RiboNN-availability short-circuit, which in CI would otherwise
    # answer () before an attestation was ever consulted.
    monkeypatch.setattr(RiboNNExpressionModel, "available", lambda self: True)
    assert attested_expression_backends() == ()


def test_a_schema_1_record_is_refused_with_the_reason(tmp_path: pathlib.Path) -> None:
    payload = _attestation().to_dict()
    for field in ("top_k", "utr_context_sha256", "verified_against_panel", "scoring_source"):
        del payload[field]
    payload["schema_version"] = 1
    with pytest.raises(ExpressionAttestationError, match="cannot be filled in after"):
        ExpressionAttestation.from_dict(payload)


# --- the seam actually reaches the user --------------------------------------

_PROTEIN = "MKAYVQTL"


def _fake_scores(monkeypatch: pytest.MonkeyPatch) -> None:
    """Score without a checkout: rank by GC so the head's pick is not the solver's."""

    def fake_predict_te(self: RiboNNExpressionModel, dnas: list[str]) -> list[float]:
        return [sum(1 for ch in dna if ch in "GC") / len(dna) for dna in dnas]

    monkeypatch.setattr(RiboNNExpressionModel, "_predict_te", fake_predict_te)


def test_an_uncalibrated_head_still_only_annotates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _fake_scores(monkeypatch)
    cand_set = assemble_and_rank_candidates(_PROTEIN, steps=5, predictor=_head())
    assert cand_set.calibrated is False
    assert cand_set.order_basis == "discovery"
    assert all(not c.expression_calibrated for c in cand_set.candidates)


def test_a_promoted_head_reorders_and_repicks(monkeypatch: pytest.MonkeyPatch) -> None:
    # The payoff, and the thing that did not happen before this seam existed: a head
    # that earned its claim now steers the delivered candidate.
    _fake_scores(monkeypatch)
    promoted = promote_if_attested(_head(), enabled=True, attestation=_attestation())
    cand_set = assemble_and_rank_candidates(_PROTEIN, steps=5, predictor=promoted)

    assert cand_set.calibrated is True
    assert cand_set.order_basis == "expression_rank"
    assert cand_set.chosen == 0
    scores = [c.expression_score for c in cand_set.candidates]
    assert scores == sorted(scores, reverse=True)
    assert all(c.expression_calibrated for c in cand_set.candidates)

    # Invariant #9, pinned so it cannot pass for the wrong reason. Comparing against an
    # UNCALIBRATED run proves nothing here: `predictor_calibrated` alone already moves
    # the hash, so that assertion held with the attestation key deleted. Two runs that
    # are calibrated by DIFFERENT attestations isolate the claim.
    other = dataclasses.replace(_attestation(), panel_sha256="0" * 64)
    assert other.content_hash() != _attestation().content_hash()
    second = promote_if_attested(_head(), enabled=True, attestation=other)
    other_set = assemble_and_rank_candidates(_PROTEIN, steps=5, predictor=second)
    assert cand_set.manifest.config_hash != other_set.manifest.config_hash


def test_promotion_does_not_change_what_the_head_computes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # An attestation is a statement ABOUT a model, never an edit to it: the same DNA
    # must score identically before and after, or the claim would be about something
    # else. Only the honesty flag moves.
    _fake_scores(monkeypatch)
    plain = _head()
    promoted = promote_if_attested(plain, enabled=True, attestation=_attestation())
    dna = "ATGAAAGCGTATGTGCAAACCCTGTAA"
    before, after = plain.score_sequence(dna), promoted.score_sequence(dna)
    assert before.score == after.score
    assert (before.calibrated, after.calibrated) == (False, True)
    assert isinstance(after, ExpressionResult)
