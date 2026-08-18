"""Tests for running the splice gate over a panel against fixed baselines.

A threshold is only meaningful next to what a *dumb* predictor scores on the same
panel, and on this task the dumb predictors are unusually strong: ~99% of human introns
open ``GT`` and close ``AG``, and BT4 already ships a consensus PWM that reads that
motif for free. A wrapped CNN that cannot beat those has not earned a PyTorch
dependency, a hash-pinned weight set, or a non-commercial licence term.

These tests pin the properties that make the comparison hard to win dishonestly:

* **BT4's own baseline cannot certify itself** -- run as the head it ties the ``pwm``
  baseline exactly and the verdict is ``False``. The structural counterpart of the
  expression gate's "the null model provably cannot pass";
* a **misaligned** backend is diagnosed rather than scored as incompetent, because a
  one-base anchor disagreement and a hopeless model look identical in the metrics;
* a **combined-track** backend collapses to one stratum instead of being credited with
  an acceptor prediction it never makes;
* the comparison is **per stratum**, so beating the motif on donors cannot excuse
  losing to it on acceptors;
* a panel overlapping the models' **training chromosomes** can never be promotable;
* the defaults produce a report, not a pass.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from bt4.biomodels.splice import ConsensusPwmSplicePredictor, SpliceResult
from bt4.biomodels.splice.panel import SplicePanel, SpliceWindow, panel_from_windows
from bt4.cli.__main__ import main
from bt4.pipeline.splice_gate import (
    SPLICE_BASELINES,
    SpliceGateSettings,
    baseline_predictions,
    run_splice_panel_gate,
    score_splice_panel,
)

_NEG = "all other positions in the window"
_DONOR_SITE = "GTAAGT"
_ACCEPTOR_SITE = "TTTTTTTTTCAG"


def _hard_panel(groups: tuple[str, ...] = ("chr1", "chr3")) -> SplicePanel:
    """A panel the PWM baseline scores well but **not** perfectly.

    Real sites carry the consensus, and unannotated *decoys* carry it too -- so the
    motif-reading baselines generate false positives and a perfect oracle can strictly
    beat them. Without decoys every baseline ties at 1.0 and the comparison proves
    nothing.
    """
    rng = random.Random(11)
    windows = [
        _hard_window(f"w{index}", group, rng) for index, group in enumerate(groups)
    ]
    return panel_from_windows(windows, negative_construction=_NEG, annotation="synthetic")


def _hard_window(window_id: str, group: str, rng: random.Random) -> SpliceWindow:
    """Build one window: four real sites plus unannotated consensus decoys."""
    parts: list[str] = []
    donors: list[int] = []
    acceptors: list[int] = []
    position = 0

    def add(chunk: str) -> int:
        nonlocal position
        parts.append(chunk)
        start = position
        position += len(chunk)
        return start

    for repeat in range(4):
        add("".join(rng.choice("ACGT") for _ in range(50)))
        add("CAG")
        donors.append(add(_DONOR_SITE))
        add("".join(rng.choice("ACGT") for _ in range(40)))
        acceptors.append(add(_ACCEPTOR_SITE) + len(_ACCEPTOR_SITE) - 1)
        add("".join(rng.choice("ACGT") for _ in range(30)))
        # Decoys: the same consensus, deliberately NOT annotated.
        if repeat % 2 == 0:
            add("CAG")
            add(_DONOR_SITE)
            add("".join(rng.choice("ACGT") for _ in range(20)))
            add(_ACCEPTOR_SITE)
            add("".join(rng.choice("ACGT") for _ in range(20)))
    add("".join(rng.choice("ACGT") for _ in range(40)))
    return SpliceWindow(
        window_id, group, "".join(parts), donors=tuple(donors), acceptors=tuple(acceptors)
    )


def _oracle(panel: SplicePanel, *, combined: bool = False) -> list[SpliceResult]:
    """A perfect backend: 1.0 exactly at the annotated sites, 0.0 everywhere else."""
    results = []
    for window in panel.windows:
        donor = [0.0] * len(window.sequence)
        acceptor = [0.0] * len(window.sequence)
        for position in window.donors:
            donor[position] = 1.0
        for position in window.acceptors:
            (donor if combined else acceptor)[position] = 1.0
        results.append(
            SpliceResult(
                donor=tuple(donor),
                acceptor=tuple(acceptor),
                model_name="oracle",
                calibrated=False,
            )
        )
    return results


def _shifted(results: list[SpliceResult], offset: int) -> list[SpliceResult]:
    """Move every track ``offset`` bases downstream, simulating an anchor disagreement."""

    def move(track: tuple[float, ...]) -> tuple[float, ...]:
        n = len(track)
        return tuple(track[i - offset] if 0 <= i - offset < n else 0.0 for i in range(n))

    return [
        SpliceResult(
            donor=move(r.donor),
            acceptor=move(r.acceptor),
            model_name="misaligned",
            calibrated=False,
        )
        for r in results
    ]


# --------------------------------------------------------------------------
# The baseline BT4 already ships cannot certify itself


def test_the_pwm_backend_cannot_beat_the_pwm_baseline() -> None:
    """The structural refusal: BT4's own default is not evidence for itself.

    Run as the head it reproduces the ``pwm`` baseline exactly, so it can never clear
    the strictly-greater comparison -- no matter how good its absolute numbers look.
    """
    comparison = run_splice_panel_gate(_hard_panel(), "pwm")
    assert comparison.backend == "consensus-pwm-baseline"
    head = {s.name: s.pr_auc_skill for s in comparison.head.strata}
    pwm_report = next(r for name, r in comparison.baselines if name == "pwm")
    pwm = {s.name: s.pr_auc_skill for s in pwm_report.strata}
    assert head == pwm
    assert comparison.beats_every_baseline is False
    assert comparison.promotable is False


def test_a_perfect_backend_beats_every_baseline() -> None:
    """The comparison is winnable -- by a model that is genuinely better."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    assert comparison.beats_every_baseline is True
    assert comparison.held_out is True
    assert comparison.promotable is True
    for _, baseline, skill in comparison.best_baseline:
        assert baseline in SPLICE_BASELINES
        assert skill < 1.0


def test_the_decoys_really_do_cost_the_motif_baselines() -> None:
    """Guards the fixture: without imperfect baselines the comparison proves nothing."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    scores = {
        name: [s.pr_auc_skill for s in report.strata]
        for name, report in comparison.baselines
    }
    assert all(value < 1.0 for value in scores["pwm"])
    assert all(value < 1.0 for value in scores["gt_ag"])


def test_beating_one_stratum_is_not_enough() -> None:
    """Per-stratum, so donor strength cannot carry an acceptor failure."""
    panel = _hard_panel()
    results = _oracle(panel)
    blinded = [
        SpliceResult(
            donor=r.donor,
            acceptor=tuple(0.5 for _ in r.acceptor),
            model_name="donor-only",
            calibrated=False,
        )
        for r in results
    ]
    comparison = run_splice_panel_gate(panel, "pwm", results=blinded, combined_track=False)
    skill = {s.name: s.pr_auc_skill for s in comparison.head.strata}
    assert skill["donor"] == pytest.approx(1.0)
    assert skill["acceptor"] < 0.01
    assert comparison.beats_every_baseline is False


# --------------------------------------------------------------------------
# The two alignment traps


def test_a_misaligned_backend_is_diagnosed_not_scored() -> None:
    """A one-base anchor disagreement and a hopeless model look identical otherwise."""
    panel = _hard_panel()
    misaligned = _shifted(_oracle(panel), 2)
    comparison = run_splice_panel_gate(panel, "pwm", results=misaligned)
    assert all(s.pr_auc_skill < 0.01 for s in comparison.head.strata)
    assert comparison.alignment.modal_offset == 2
    assert "anchors DISAGREE" in comparison.alignment.note()
    assert any("DISAGREE" in note for note in comparison.notes)


def test_declaring_the_offset_recovers_the_perfect_score() -> None:
    """And with the anchor declared, the same backend scores exactly as it should."""
    panel = _hard_panel()
    misaligned = _shifted(_oracle(panel), 2)
    comparison = run_splice_panel_gate(panel, "pwm", results=misaligned, anchor_offset=2)
    assert all(s.pr_auc_skill == pytest.approx(1.0) for s in comparison.head.strata)
    assert comparison.alignment.modal_offset == 0
    assert "anchors agree" in comparison.alignment.note()


def test_an_aligned_backend_reports_a_clean_diagnostic() -> None:
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    assert comparison.alignment.modal_offset == 0
    assert comparison.alignment.fraction_at_zero == 1.0
    assert comparison.alignment.n_sites == panel.n_sites


def test_a_combined_track_backend_collapses_to_one_stratum() -> None:
    """Pangolin's binary head never separated donor from acceptor.

    Scoring its all-zero acceptor track as an acceptor prediction would report it as
    perfectly hopeless -- an artifact of the wrapper, not a finding about the model.
    """
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel, combined=True))
    assert comparison.strata == ("splice",)
    assert comparison.head.strata[0].pr_auc_skill == pytest.approx(1.0)
    assert any("combined P(splice) track" in note for note in comparison.notes)


def test_the_collapse_can_be_forced_off() -> None:
    """Detection is from the output, so a caller who knows better can override it."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(
        panel, "pwm", results=_oracle(panel, combined=True), combined_track=False
    )
    assert comparison.strata == ("acceptor", "donor")


def test_gt_ag_scores_the_combined_stratum_with_either_motif() -> None:
    """A baseline that could only see donors would be an unfair, meaningless control."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel, combined=True))
    gt_ag = next(report for name, report in comparison.baselines if name == "gt_ag")
    assert gt_ag.strata[0].pr_auc_skill > 0.0


# --------------------------------------------------------------------------
# Held-out-ness gates promotion


def test_a_training_chromosome_panel_can_never_be_promotable() -> None:
    """Both models trained on chr2; a metric computed there is optimistic."""
    panel = _hard_panel(groups=("chr2", "chr4"))
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    assert comparison.beats_every_baseline is True
    assert comparison.head.passed is True
    assert comparison.held_out is False
    assert comparison.promotable is False
    assert any("NOT HELD OUT" in note for note in comparison.notes)


# --------------------------------------------------------------------------
# Baselines


def test_constant_is_perfectly_calibrated_and_carries_no_information() -> None:
    """The permanent trap: excellent ECE, zero skill, visible in the same table."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    constant = next(report for name, report in comparison.baselines if name == "constant")
    for stratum in constant.strata:
        assert stratum.ece < 0.01
        assert stratum.pr_auc_skill == pytest.approx(0.0, abs=1e-9)
        assert stratum.brier_skill == pytest.approx(0.0, abs=1e-9)


def test_permutation_preserves_each_stratums_score_distribution() -> None:
    """The null must differ only in the pairing, never in the marginal."""
    panel = _hard_panel()
    results = _oracle(panel)
    comparison = run_splice_panel_gate(panel, "pwm", results=results)
    permutation = next(report for name, report in comparison.baselines if name == "permutation")
    for stratum in permutation.strata:
        assert stratum.pr_auc_skill < 0.1
        assert stratum.n_cases == {s.name: s.n_cases for s in comparison.head.strata}[stratum.name]


def test_permutation_is_deterministic_from_its_seed() -> None:
    """Invariant #7: the same panel and seed give the same null, every time."""
    panel = _hard_panel()
    results = _oracle(panel)
    first = run_splice_panel_gate(panel, "pwm", results=results)
    second = run_splice_panel_gate(panel, "pwm", results=results)
    assert [s.pr_auc for _, r in first.baselines for s in r.strata] == [
        s.pr_auc for _, r in second.baselines for s in r.strata
    ]


def test_a_different_seed_moves_the_permutation_baseline() -> None:
    """Scored with the graded PWM, not the oracle.

    A two-valued predictor's null is near-deterministic at this prevalence: shuffling
    8 positives into a 600-position track almost never lands one in the 1.0 tie block,
    so both seeds give the same average precision. That is correct behaviour and a
    useless test of the seed, so this uses a continuous score instead.
    """
    panel = _hard_panel()
    default = run_splice_panel_gate(panel, "pwm")
    reseeded = run_splice_panel_gate(panel, "pwm", settings=SpliceGateSettings(seed=99))

    def permutation(comparison: object) -> list[float]:
        return [
            s.pr_auc
            for name, report in comparison.baselines  # type: ignore[attr-defined]
            if name == "permutation"
            for s in report.strata
        ]

    assert permutation(default) != permutation(reseeded)


def test_an_unknown_baseline_is_refused() -> None:
    with pytest.raises(ValueError, match="unknown baseline"):
        run_splice_panel_gate(_hard_panel(), "pwm", baselines=("coin_flip",))
    panel = _hard_panel()
    cases, refs = [], []
    with pytest.raises(ValueError, match="unknown baseline"):
        baseline_predictions("coin_flip", panel, cases, refs, [], 0)


def test_dropping_a_baseline_is_possible_but_visible() -> None:
    panel = _hard_panel()
    comparison = run_splice_panel_gate(
        panel, "pwm", results=_oracle(panel), baselines=("constant",)
    )
    assert [name for name, _ in comparison.baselines] == ["constant"]


# --------------------------------------------------------------------------
# The gate cannot certify by default


def test_defaults_produce_a_report_not_a_pass() -> None:
    """A bare call yields numbers; the thresholds are the maintainer's to set."""
    settings = SpliceGateSettings()
    assert settings.min_pr_auc == 0.0
    assert settings.max_ece == 1.0
    comparison = run_splice_panel_gate(_hard_panel(), "pwm")
    assert comparison.head.min_pr_auc == 0.0


def test_thresholds_are_threaded_into_every_report() -> None:
    """The head and the baselines are judged by exactly the same bar."""
    panel = _hard_panel()
    settings = SpliceGateSettings(min_pr_auc=0.60)
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel), settings=settings)
    assert comparison.head.min_pr_auc == 0.60
    assert all(report.min_pr_auc == 0.60 for _, report in comparison.baselines)


def test_the_panel_provenance_reaches_every_report() -> None:
    """A number is never separated from the construction that produced it."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    assert comparison.panel_hash == panel.content_hash()
    assert comparison.head.negative_construction == _NEG
    assert comparison.head.panel_note == "synthetic"


def test_calibration_status_is_reported_beside_the_verdict_not_changed() -> None:
    """This gate flips nothing -- the fidelity flag is a different question entirely."""
    panel = _hard_panel()
    comparison = run_splice_panel_gate(panel, "pwm", results=_oracle(panel))
    assert comparison.backend_calibrated is False
    assert ConsensusPwmSplicePredictor().calibrated is False


def test_a_results_length_mismatch_is_refused() -> None:
    panel = _hard_panel()
    with pytest.raises(ValueError, match="entries for"):
        run_splice_panel_gate(panel, "pwm", results=_oracle(panel)[:1])


def test_score_splice_panel_returns_one_result_per_window() -> None:
    panel = _hard_panel()
    results = score_splice_panel(panel, "pwm")
    assert len(results) == len(panel.windows)
    assert all(len(r.donor) == len(w.sequence) for r, w in zip(results, panel.windows, strict=True))


def test_an_unknown_backend_is_refused() -> None:
    with pytest.raises(ValueError, match="unknown splice backend"):
        run_splice_panel_gate(_hard_panel(), "not-a-model")


# --------------------------------------------------------------------------
# The CLI surface


def _write_panel(tmp_path: Path) -> str:
    panel = _hard_panel()
    lines = ["window_id\tgroup\tsequence\tdonors\tacceptors"]
    for window in panel.windows:
        lines.append(
            "\t".join(
                (
                    window.window_id,
                    window.group,
                    window.sequence,
                    ",".join(str(p) for p in window.donors),
                    ",".join(str(p) for p in window.acceptors),
                )
            )
        )
    path = tmp_path / "panel.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def test_cli_reports_the_verdict_and_every_baseline(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """``bt4 splice-gate`` shows the head beside its controls, not alone."""
    assert main(
        [
            "splice-gate",
            _write_panel(tmp_path),
            "--negative-construction",
            _NEG,
            "--annotation",
            "synthetic",
        ]
    ) == 0
    out = capsys.readouterr().out
    for baseline in SPLICE_BASELINES:
        assert baseline in out
    assert "PROMOTABLE on this panel : False" in out
    assert "anchors agree" in out


def test_cli_requires_the_negative_construction(tmp_path: Path) -> None:
    """A PR-AUC threshold without a pinned denominator means nothing."""
    with pytest.raises(SystemExit):
        main(["splice-gate", _write_panel(tmp_path)])


def test_cli_refuses_the_network_backend(tmp_path: Path) -> None:
    """ASSP is excluded from the reproducible path, so it cannot support a gate."""
    with pytest.raises(SystemExit):
        main(
            [
                "splice-gate",
                _write_panel(tmp_path),
                "--negative-construction",
                _NEG,
                "--backend",
                "assp",
            ]
        )
