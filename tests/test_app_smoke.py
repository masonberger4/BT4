"""Headless smoke tests for the BT4 Studio desktop app.

These run under Qt's ``offscreen`` platform (no display needed). They drive the
engine synchronously through the worker's ``compute`` and call the window's
result/failure slots directly -- no real ``QThread`` is started and the event
loop is never entered, so the suite stays fast and hermetic.

The ASSP cross-check tests are driven from the committed **offline** fixtures
(``$BT4_ASSP_FIXTURE_DIR``), so the suite never makes a network call.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

pytest.importorskip("PySide6")

from PySide6 import QtGui, QtWidgets

import bt4
from bt4 import api
from bt4.app import studio
from bt4.app.studio import SequenceViewer, StudioWindow
from bt4.app.worker import (
    CandidatesResult,
    CandidatesWorker,
    CrossCheckWorker,
    LibraryWorker,
    OptimizeWorker,
)

ASSP_FIXTURES = Path(__file__).parent / "fixtures" / "assp"
# The sequence the committed ASSP fixtures were captured for (fixtures are keyed
# by sequence hash, so only this one resolves offline).
ASSP_SEQ = "ATGGCCGGCGATCGATCGATCGTAA"


@pytest.fixture(autouse=True, scope="module")
def _qapp() -> QtWidgets.QApplication:
    """A single offscreen QApplication shared by every test in the module."""
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    assert isinstance(app, QtWidgets.QApplication)
    return app


def test_worker_compute_runs() -> None:
    """The worker computes a frontier synchronously and it round-trips."""
    worker = OptimizeWorker("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), steps=5)
    result = worker.compute()

    assert len(result.results) >= 1
    delivered = result.delivered()
    assert delivered is not None
    assert delivered.protein == "MAALKHETQW"


def test_worker_reports_progress() -> None:
    """The worker emits per-point progress that ends at 100%."""
    seen: list[tuple[int, str]] = []
    worker = OptimizeWorker("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), steps=5)
    worker.progress.connect(lambda pct, label: seen.append((pct, label)))
    worker.compute()
    assert seen, "expected at least one progress update"
    assert seen[-1][0] == 100


def test_worker_cancel_before_any_point_raises() -> None:
    """Cancelling before the first point yields a ValueError, not a hang."""
    worker = OptimizeWorker("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), steps=9)
    worker.cancel()
    with pytest.raises(ValueError):
        worker.compute()


def test_window_populates() -> None:
    """Feeding a real frontier to the finished-slot populates every results widget."""
    app = QtWidgets.QApplication.instance()
    assert isinstance(app, QtWidgets.QApplication)

    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)

    window._on_finished(frontier)
    app.processEvents()

    delivered = frontier.delivered()
    assert delivered is not None
    assert delivered.dna in window.sequence_view.toPlainText()
    assert window.badge.text().strip() != ""
    assert window.metrics_table.rowCount() > 0
    assert window.export_fasta_btn.isEnabled()
    assert window.optimize_btn.isEnabled()


def test_failure_clears_stale_results() -> None:
    """A failure after a success clears the panel so nothing stale is exportable."""
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    assert window.export_fasta_btn.isEnabled()

    window._on_failed(ValueError("boom"))

    assert not window.export_fasta_btn.isEnabled()
    assert not window.export_json_btn.isEnabled()
    assert window.sequence_view.toPlainText() == ""
    assert window._delivered() is None
    assert window.optimize_btn.isEnabled()


def test_infeasible_is_handled_with_friendly_message() -> None:
    """The failure slot never raises and translates InfeasibleError to plain language."""
    window = StudioWindow()
    window._set_running(True)
    assert not window.optimize_btn.isEnabled()

    window._on_failed(api.InfeasibleError(["homopolymer", "restriction_site"]))

    assert window.optimize_btn.isEnabled()
    assert "satisfy these settings" in window.statusBar().currentMessage()


def test_protein_input_is_validated() -> None:
    """The protein box is cleaned (FASTA/case) and bad input is rejected, not crashed."""
    window = StudioWindow()

    window.protein_edit.setPlainText("maal khet qw")
    assert window._prepare_protein() == "MAALKHETQW"

    window.protein_edit.setPlainText(">seq1 description\nMAAL\nKHET")
    assert window._prepare_protein() == "MAALKHET"

    window.protein_edit.setPlainText("MAALKHET*")  # trailing stop
    assert window._prepare_protein() is None

    window.protein_edit.setPlainText("MAALBZ")  # B, Z are not amino acids
    assert window._prepare_protein() is None

    window.protein_edit.setPlainText("   ")  # empty
    assert window._prepare_protein() is None


def test_enzyme_names_are_case_insensitive() -> None:
    """Enzyme entries are canonicalized to catalog casing; unknown ones are rejected."""
    window = StudioWindow()

    window.enzymes_edit.setText("ecori, bamhi")
    assert window._prepare_enzymes() == ("EcoRI", "BamHI")

    window.enzymes_edit.setText("EcoRI, NotAnEnzyme")
    assert window._prepare_enzymes() is None

    window.enzymes_edit.setText("")
    assert window._prepare_enzymes() == ()


def test_enzyme_field_completes_the_last_entry_only() -> None:
    """The catalog has hundreds of enzymes, so the field has to be searchable.

    A stock completer matches the whole line, which breaks once the user has
    already listed one enzyme. This one completes only the token after the last
    comma and substitutes it back, leaving earlier entries intact.
    """
    window = StudioWindow()
    completer = window.enzyme_completer
    assert completer.model().rowCount() == len(api.available_enzymes())

    # Only the trailing token is completed against.
    assert completer.splitPath("EcoRI, bsa") == ["bsa"]
    assert completer.splitPath("ecori") == ["ecori"]

    # Substituting a completion preserves the earlier entries.
    window.enzymes_edit.setText("EcoRI, bsa")
    completer.setCompletionPrefix("BsaI")
    source_index = completer.completionModel().mapToSource(
        completer.completionModel().index(0, 0)
    )
    assert completer.pathFromIndex(source_index) == "EcoRI, BsaI"


def test_unknown_enzyme_warning_does_not_dump_the_catalog() -> None:
    """A 500-name wall of text hides the answer; near misses give it."""
    window = StudioWindow()
    window.enzymes_edit.setText("EcoR1")
    assert window._prepare_enzymes() is None

    box = window._msgbox
    assert box is not None
    detail = box.informativeText()
    assert "EcoRI" in detail, "the obvious correction should be offered"
    assert len(detail) < 600, "the whole catalog must not be pasted into the error"
    # Each suggestion carries its own site and is labelled a SPELLING match, so a
    # user cannot read a near-miss name as an equivalent recognition sequence.
    assert "SPELLING" in detail
    assert f"EcoRI ({api.resolve_enzyme('EcoRI')})" in detail
    assert "Forbidden motifs" in detail, "the escape hatch should be offered"


def test_sequence_viewer_highlights_and_locates_violations() -> None:
    """A violation span is highlighted and locatable by nucleotide position."""
    viewer = SequenceViewer(dark=False)
    dna = "ACGTACGTACGT"
    hard = api.Violation("max_repeat", api.Severity.HARD, 0, 4, "dispersed repeat")
    soft = api.Violation("cpg", api.Severity.SOFT, 6, 9, "elevated CpG")
    viewer.set_sequence(dna, (hard, soft))

    # Text is unchanged; both spans become extra-selection highlights.
    assert viewer.toPlainText() == dna
    assert len(viewer.extraSelections()) == 2

    # Position lookup drives the hover tooltip.
    assert viewer._violation_at(2) is hard
    assert viewer._violation_at(7) is soft
    assert viewer._violation_at(4) is None  # end is exclusive
    assert viewer._violation_at(11) is None


def test_sequence_viewer_drops_out_of_range_spans() -> None:
    """Spans outside the sequence are dropped defensively, not crashed on."""
    viewer = SequenceViewer(dark=True)
    dna = "ACGTAC"
    bad = api.Violation("bogus", api.Severity.HARD, 4, 99, "past the end")
    viewer.set_sequence(dna, (bad,))
    assert viewer.extraSelections() == []
    assert viewer._violation_at(5) is None


def test_hard_violation_wins_overlap() -> None:
    """When HARD and SOFT spans overlap a base, the tooltip resolves to HARD."""
    viewer = SequenceViewer(dark=False)
    dna = "ACGTACGTAC"
    soft = api.Violation("cpg", api.Severity.SOFT, 0, 8, "wide soft band")
    hard = api.Violation("homopolymer", api.Severity.HARD, 2, 5, "run")
    viewer.set_sequence(dna, (soft, hard))
    assert viewer._violation_at(3) is hard


def test_clean_result_hides_violation_legend() -> None:
    """A delivered sequence with no violations keeps the legend hidden."""
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)

    delivered = frontier.delivered()
    assert delivered is not None
    assert delivered.violations == ()  # feasible run: nothing to annotate
    assert not window.violations_legend.isVisible()
    assert len(window.sequence_view.extraSelections()) == 0


def test_candidates_worker_computes() -> None:
    """The candidates worker assembles a set and splice-audits it synchronously."""
    worker = CandidatesWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        steps=5,
        n=8,
        repeat_variants=2,
        include_cnns=False,
    )
    result = worker.compute()
    assert isinstance(result, CandidatesResult)
    assert result.candidate_set.candidates
    delivered = result.candidate_set.delivered()
    assert delivered is not None
    # The default expression head is the neutral placeholder: uncalibrated, so
    # the set stays in discovery order (never presented as a ranking).
    assert result.candidate_set.calibrated is False
    assert result.candidate_set.order_basis == "discovery"
    # A splice audit ran over a non-empty set and is advisory (no CNN installed).
    assert result.audit is not None
    assert result.audit.all_calibrated is False


def test_window_renders_candidates() -> None:
    """Feeding a candidate result populates the table and honest banners."""
    app = QtWidgets.QApplication.instance()
    assert isinstance(app, QtWidgets.QApplication)

    window = StudioWindow()
    worker = CandidatesWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        steps=5,
        n=8,
        repeat_variants=2,
        include_cnns=False,
    )
    result = worker.compute()

    window._on_cand_finished(result)
    app.processEvents()

    assert window.candidates_table.rowCount() == len(result.candidate_set.candidates)
    # Uncalibrated head => the banner says "discovery order, not a ranking".
    assert "not a ranking" in window.cand_banner.text().lower()
    # The splice audit is advisory and says so.
    assert "advisory" in window.splice_banner.text().lower()
    # The splice-flags column is filled for the delivered candidate.
    chosen = result.candidate_set.chosen
    assert window.candidates_table.item(chosen, 7) is not None
    # The rank button is re-enabled after a finished run.
    assert window.rank_btn.isEnabled()


def test_distinct_site_count_merges_close_positions() -> None:
    """Co-located flags (within the match window) count as one site."""
    from bt4.app.studio import _distinct_site_count

    assert _distinct_site_count([], 3) == 0
    assert _distinct_site_count([10], 3) == 1
    # 10 and 12 are within the window (one site); 40 is a second site.
    assert _distinct_site_count([10, 12, 40], 3) == 2
    # A chain each within the window collapses to a single site.
    assert _distinct_site_count([10, 13, 16, 19], 3) == 1


# --------------------------------------------------------------------------- #
# Expression head (RiboNN) wiring
# --------------------------------------------------------------------------- #


def test_ribonn_toggle_is_disabled_without_a_checkout() -> None:
    """With no RiboNN install, the toggle is disabled -- never a dead control."""
    window = StudioWindow()
    # CI never has the Sanofi non-commercial checkout/weights, so the backend
    # list is placeholder-only and the control disables itself with a tooltip
    # explaining exactly what is missing.
    assert "ribonn" not in api.available_expression_backends()
    assert window._ribonn_available is False
    assert not window.ribonn_check.isEnabled()
    assert "BT4_RIBONN_DIR" in window.ribonn_check.toolTip()
    # Its sub-controls follow the toggle rather than sitting enabled-but-useless.
    assert not window.utr5_edit.isEnabled()
    assert not window.ribonn_species_combo.isEnabled()


def test_prepare_predictor_defaults_to_the_placeholder() -> None:
    """With RiboNN unselected the candidate flow passes no predictor (=default)."""
    window = StudioWindow()
    ok, predictor = window._prepare_predictor()
    assert ok is True
    assert predictor is None


def test_prepare_predictor_refuses_missing_utrs() -> None:
    """Selecting RiboNN without UTR context is refused before the run starts.

    RiboNN cannot score an empty UTR column, so catching it here (rather than
    letting the engine raise mid-run) is what keeps the failure legible.
    """
    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)

    ok, predictor = window._prepare_predictor()
    assert ok is False
    assert predictor is None
    assert "UTR" in window.statusBar().currentMessage()


def test_prepare_predictor_builds_an_uncalibrated_ribonn() -> None:
    """With UTRs supplied the flow builds RiboNN -- still uncalibrated."""
    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window.utr5_edit.setText("acgtacgt")
    window.utr3_edit.setText("TTTTGGGG")

    ok, predictor = window._prepare_predictor()
    assert ok is True
    assert predictor is not None
    assert predictor.name == "ribonn[human]"
    # Constructing the adapter must not confer calibration (CLAUDE.md §10.6):
    # only a passing CDS-variant acceptance gate can.
    assert predictor.calibrated is False


def test_prepare_predictor_rejects_non_dna_utrs() -> None:
    """A UTR box holding non-DNA characters is refused with a plain message."""
    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window.utr5_edit.setText("ACGXTT")
    window.utr3_edit.setText("ACGT")

    ok, predictor = window._prepare_predictor()
    assert ok is False
    assert predictor is None


def test_candidates_worker_forwards_the_predictor() -> None:
    """A predictor handed to the worker reaches the candidate set's annotations."""
    predictor = api.resolve_expression_backend("null")
    worker = CandidatesWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        steps=5,
        n=6,
        repeat_variants=0,
        include_cnns=False,
        predictor=predictor,
    )
    result = worker.compute()
    delivered = result.candidate_set.delivered()
    assert delivered is not None
    assert delivered.expression_model == predictor.name
    # An uncalibrated head never steers delivery.
    assert result.candidate_set.calibrated is False
    assert result.candidate_set.order_basis == "discovery"


# --------------------------------------------------------------------------- #
# Expression attestation (the seam that lets an earned calibration reach the UI)
# --------------------------------------------------------------------------- #


def _attestation() -> api.ExpressionAttestation:
    """A passing attestation for a HEK293T panel (a test double -- see below)."""
    rows = [
        api.PanelRow(
            group=f"P{g:02d}",
            variant_id=f"g{g}v{v}",
            cds="ATG" + ("AAA" * (v + 1)) + ("AAG" * (6 - v)) + "TAA",
            measured=100.0 * g + v,
            utr5="GCCACC",
            utr3="GCTAAT",
        )
        for g in range(10)
        for v in range(6)
    ]
    panel = api.panel_from_rows(rows)
    comparison = api.expression_gate(
        panel,
        "null",
        settings=api.GateSettings(
            within_group=True, recalibrate=True, coverage_tolerance=0.10,
            bootstrap_resamples=200,
        ),
        cell_types=("HEK293T",),
        head_scores=[row.measured for row in panel.rows],
    )
    # A test double: relabelled as a gate-scored RiboNN run, which is the only way to
    # reach a promotable record without the Sanofi non-commercial weights. Production
    # refuses both overrides (see tests/test_expression_promotion.py).
    import dataclasses

    as_ribonn = dataclasses.replace(
        comparison,
        backend="ribonn[human]",
        scope=dataclasses.replace(comparison.scope, scoring_source="gate"),
    )
    return api.attest_expression(
        as_ribonn, readout="mean_ribosome_load", bt4_version="0.0.0-test"
    )


def test_attestation_toggle_is_disabled_and_explains_what_is_missing() -> None:
    """With no attestation resolvable, the toggle says why rather than doing nothing."""
    window = StudioWindow()
    assert window._expr_attestation is None
    assert not window.expr_attested_check.isEnabled()
    tip = window.expr_attested_check.toolTip()
    # The dead-control rule: name the missing thing AND how to supply it.
    assert "no expression attestation is available" in tip.lower()
    assert api.EXPRESSION_ATTESTATION_ENV_VAR in tip


def test_attestation_toggle_lights_up_and_names_its_scope(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A resolvable attestation enables the toggle and shows what it covers."""
    import json

    path = tmp_path / "ribonn_attestation.json"
    path.write_text(json.dumps(_attestation().to_dict()), encoding="utf-8")
    monkeypatch.setenv(api.EXPRESSION_ATTESTATION_ENV_VAR, str(path))

    window = StudioWindow()
    window._ribonn_available = True
    # CI has no RiboNN checkout, so the offerability probe says no; the attestation half
    # is what this test is about.
    window._expr_attestable = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window._update_ribonn_enabled(True)

    assert window._expr_attestation is not None
    assert window.expr_attested_check.isEnabled()
    tip = window.expr_attested_check.toolTip()
    # The scope IS the claim, so it is on the control, not buried in a doc.
    assert "HEK293T" in tip
    assert "mean_ribosome_load" in tip
    assert "REFUSED" in tip

    # And on the page, not only in a tooltip -- plus the species control is pinned to
    # the attestation while it is honoured, so the form cannot show one scope while the
    # run uses another.
    window.expr_attested_check.setChecked(True)
    assert "HEK293T" in window.expr_scope_label.text()
    assert window.ribonn_species_combo.currentText() == "human"
    assert not window.ribonn_species_combo.isEnabled()

    window.expr_attested_check.setChecked(False)
    assert window.expr_scope_label.text() == ""
    assert window.ribonn_species_combo.isEnabled()

    # And unticking the HEAD clears the claim too: a scope for a head that is no longer
    # in the run would be a calibration statement about nothing.
    window.expr_attested_check.setChecked(True)
    assert "HEK293T" in window.expr_scope_label.text()
    window.ribonn_check.setChecked(False)
    assert window.expr_scope_label.text() == ""


def test_attestation_toggle_promotes_the_head_it_covers(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Opting in promotes RiboNN -- but only when configured inside the scope."""
    import json

    path = tmp_path / "a.json"
    path.write_text(json.dumps(_attestation().to_dict()), encoding="utf-8")
    monkeypatch.setenv(api.EXPRESSION_ATTESTATION_ENV_VAR, str(path))

    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window.expr_attested_check.setEnabled(True)
    window.expr_attested_check.setChecked(True)
    window.utr5_edit.setText("GCCACC")
    window.utr3_edit.setText("GCTAAT")

    ok, predictor = window._prepare_predictor()
    assert ok is True
    assert predictor is not None
    assert predictor.calibrated is True
    # The head is built to the attestation's scope, not to whatever the form defaulted
    # to: averaging all 78 cell types is a different quantity from the one measured.
    assert predictor.cell_types == ("HEK293T",)
    assert predictor.species == "human"
    assert predictor.top_k == 5


def test_a_head_outside_the_attested_scope_is_refused_not_downgraded(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A UTR context the attestation never covered stops the run and says so.

    Silently handing back an uncalibrated head to someone who ticked "calibrated
    ranking" is the exact failure this layer exists to prevent, so the run is refused
    with the mismatch named.
    """
    import json

    path = tmp_path / "a.json"
    path.write_text(json.dumps(_attestation().to_dict()), encoding="utf-8")
    monkeypatch.setenv(api.EXPRESSION_ATTESTATION_ENV_VAR, str(path))

    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window.expr_attested_check.setEnabled(True)
    window.expr_attested_check.setChecked(True)
    window.utr5_edit.setText("AAAAAA")  # not the panel's context
    window.utr3_edit.setText("TTTTTT")

    ok, predictor = window._prepare_predictor()
    assert ok is False
    assert predictor is None
    assert "UTR context" in window.statusBar().currentMessage()


def test_losing_the_attestation_before_a_run_refuses_rather_than_downgrading(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The window opened with a record; by Run time it is gone. Refuse, do not downgrade.

    Building an uncalibrated head here would answer "give me a calibrated ranking" with
    an uncalibrated one and say nothing -- the failure this whole layer exists to
    prevent, and one the user has no way to notice from the run itself.
    """
    import json

    path = tmp_path / "a.json"
    path.write_text(json.dumps(_attestation().to_dict()), encoding="utf-8")
    monkeypatch.setenv(api.EXPRESSION_ATTESTATION_ENV_VAR, str(path))

    window = StudioWindow()
    window._ribonn_available = True
    window.ribonn_check.setEnabled(True)
    window.ribonn_check.setChecked(True)
    window.expr_attested_check.setEnabled(True)
    window.expr_attested_check.setChecked(True)
    window.utr5_edit.setText("GCCACC")
    window.utr3_edit.setText("GCTAAT")
    assert window._prepare_predictor()[0] is True  # it works while the record is there

    path.unlink()
    ok, predictor = window._prepare_predictor()
    assert ok is False
    assert predictor is None
    assert "not a readable file" in window.statusBar().currentMessage()

    # A corrupt record is refused the same way, and named.
    path.write_text("{not json", encoding="utf-8")
    ok, predictor = window._prepare_predictor()
    assert ok is False
    assert predictor is None
    assert "not valid JSON" in window.statusBar().currentMessage()


def test_the_banner_names_the_scope_when_the_ranking_is_calibrated() -> None:
    """A calibrated set flips the banner to a ranking -- and says what it covers."""
    window = StudioWindow()
    attestation = _attestation()
    window._expr_attestation = attestation
    # The banner names a scope only for the record that actually promoted this run's
    # head, not merely for one that happens to be loaded.
    window._run_attestation_sha = attestation.content_hash()

    worker = CandidatesWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        steps=5,
        n=6,
        repeat_variants=0,
        include_cnns=False,
    )
    uncalibrated = worker.compute().candidate_set
    # Re-label the same set as if a promoted head had produced it, so this test pins the
    # RENDERING rule without needing the licensed weights a real promotion would score
    # with (the promotion path itself is covered in test_expression_promotion.py).
    import dataclasses

    calibrated = dataclasses.replace(
        uncalibrated, calibrated=True, order_basis="expression_rank"
    )
    window._render_candidates(calibrated)
    text = window.cand_banner.text()
    assert "Ranked by predicted expression" in text
    assert "HEK293T" in text
    assert "mean_ribosome_load" in text
    assert "not a ranking" not in text

    window._render_candidates(uncalibrated)
    assert "not a ranking" in window.cand_banner.text().lower()


def test_the_banner_refuses_to_name_a_scope_that_did_not_promote_this_set() -> None:
    """A loaded attestation is not evidence about the head that produced a set.

    If the record on disk changed between opening the window and running, naming its
    scope would attribute one calibration's claim to another's numbers. The banner
    compares the promoting head's own recorded hash and declines rather than guessing.
    """
    import dataclasses

    window = StudioWindow()
    window._expr_attestation = _attestation()
    window._run_attestation_sha = "0" * 64  # a different record promoted this run

    worker = CandidatesWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        steps=5,
        n=4,
        repeat_variants=0,
        include_cnns=False,
    )
    calibrated = dataclasses.replace(
        worker.compute().candidate_set, calibrated=True, order_basis="expression_rank"
    )
    window._render_candidates(calibrated)
    text = window.cand_banner.text()
    assert "scope unavailable" in text
    assert "HEK293T" not in text


# --------------------------------------------------------------------------- #
# ASSP cross-check (offline fixtures; never a live call)
# --------------------------------------------------------------------------- #


def test_crosscheck_worker_renders_an_available_report(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A fixture-backed ASSP report renders with its honesty tags and its sites."""
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    report = CrossCheckWorker(ASSP_SEQ).compute()
    assert report.available is True
    assert report.network_derived is True
    assert report.calibrated is False

    window = StudioWindow()
    window._render_crosscheck(report)

    banner = window.assp_banner.text()
    assert "network-derived" in banner
    assert "UNCALIBRATED" in banner
    assert "not</b> part of the run manifest" in banner
    assert window.assp_table.rowCount() == len(report.sites)
    assert window.assp_btn.text() == "Validate with ASSP"


def test_crosscheck_unavailable_leaves_the_run_intact(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An ASSP outage degrades to a labeled banner and never fails the design."""
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)

    # No fixture exists for this sequence, so the transport reports unavailable.
    delivered = frontier.delivered()
    assert delivered is not None
    report = CrossCheckWorker(delivered.dna).compute()
    assert report.available is False
    assert report.reason

    window._on_crosscheck_finished(report)
    assert "Unavailable" in window.assp_banner.text()
    assert not window.assp_table.isVisible()
    # The delivered result is untouched and still exportable.
    assert window._delivered() is not None
    assert window.export_fasta_btn.isEnabled()


def test_crosscheck_is_cleared_when_the_delivered_sequence_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A report is about ONE sequence, so a new result must clear the panel.

    Otherwise the previous sequence's splice sites would be shown beside a
    different design, silently attributing them to it.
    """
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    window = StudioWindow()
    window._render_crosscheck(CrossCheckWorker(ASSP_SEQ).compute())
    assert window.assp_table.rowCount() > 0

    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    assert window.assp_table.rowCount() == 0
    assert "Not run" in window.assp_banner.text()


def test_crosscheck_for_another_sequence_is_discarded(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A report that finishes after the design changed is discarded, not shown.

    A cross-check describes exactly one sequence. Rendering a late report beside
    a different delivered design would attribute one sequence's predicted splice
    sites to another -- the misattribution the panel-clearing rule exists to
    prevent, arriving by the other ordering.
    """
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    delivered = frontier.delivered()
    assert delivered is not None

    stale = CrossCheckWorker(ASSP_SEQ).compute()
    assert stale.available is True
    assert stale.dna.upper() != delivered.dna.upper()

    window._on_crosscheck_finished(stale)
    assert window.assp_table.rowCount() == 0
    assert "discarded" in window.assp_banner.text().lower()
    assert "discarded" in window.statusBar().currentMessage().lower()
    # The design itself is untouched and still exportable.
    assert window._delivered() is not None
    assert window.export_fasta_btn.isEnabled()


def test_menu_actions_obey_the_same_run_gate() -> None:
    """Shortcuts must not be a back door around the "one flow at a time" gate.

    The Run-menu actions carry keyboard shortcuts, so gating only the buttons
    would let Ctrl+R start an optimization during an in-flight cross-check.
    """
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    assert window.act_optimize.isEnabled()
    assert window.act_assp.isEnabled()

    for setter in (
        window._set_running,
        window._set_candidates_running,
        window._set_library_running,
        window._set_crosscheck_running,
    ):
        setter(True)
        assert not window.act_optimize.isEnabled()
        assert not window.act_rank.isEnabled()
        assert not window.act_library.isEnabled()
        assert not window.act_assp.isEnabled()
        setter(False)
        assert window.act_optimize.isEnabled()

    # Exports follow their buttons through the menu too.
    assert window.act_export_fasta.isEnabled()
    window._on_failed(ValueError("boom"))
    assert not window.act_export_fasta.isEnabled()
    assert not window.act_export_json.isEnabled()


def test_crosscheck_escapes_service_supplied_text() -> None:
    """A remote service's text cannot inject markup into the honesty banner.

    The banner is RichText and ``reason`` carries whatever the service said, so
    unescaped markup could rewrite -- or hide -- the very labels that mark these
    numbers network-derived, uncalibrated, and advisory.
    """
    window = StudioWindow()
    hostile = "<b>calibrated</b><span style='display:none'>"
    report = api.SpliceCrossCheck(
        dna="ATGTAA",
        backend="assp",
        available=False,
        reason=hostile,
        calibrated=False,
        network_derived=True,
        threshold=0.5,
        top_k=3,
        pooled_risk=0.0,
        sites=(),
    )
    window._render_crosscheck(report)

    banner = window.assp_banner.text()
    assert hostile not in banner            # the raw markup did not survive
    assert "&lt;b&gt;calibrated&lt;/b&gt;" in banner  # it is shown as text
    assert "UNCALIBRATED" in banner         # ...and our own labels are intact
    assert "network-derived" in banner


def test_crosscheck_button_needs_a_delivered_sequence() -> None:
    """The ASSP control is enabled only once there is something to cross-check."""
    window = StudioWindow()
    assert not window.assp_btn.isEnabled()

    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    assert window.assp_btn.isEnabled()

    window._on_failed(ValueError("boom"))
    assert not window.assp_btn.isEnabled()


def test_export_json_never_carries_crosscheck_numbers(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A written export is byte-identical whether or not a cross-check was run.

    ASSP numbers are network-derived and excluded from the
    reproducible-from-manifest guarantee (CLAUDE.md §6, §10.15), so they must
    never reach an exported artifact. This drives the window's real export
    action -- not just the serializer -- on both sides of a rendered report.
    """
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)

    written: list[Path] = []

    def export_to(path: Path) -> str:
        monkeypatch.setattr(
            QtWidgets.QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (str(path), "")),
        )
        window._export_json()
        written.append(path)
        return path.read_text(encoding="utf-8")

    before = export_to(tmp_path / "before.json")

    # Render a real (fixture-backed) ASSP report into the window, then re-export.
    report = CrossCheckWorker(ASSP_SEQ).compute()
    window._render_crosscheck(report)
    assert window.assp_table.rowCount() > 0  # the report really is on screen
    assert report.pooled_risk != 0.0  # ...and carries numbers that could leak
    after = export_to(tmp_path / "after.json")

    assert len(written) == 2
    assert before == after
    assert "assp" not in after.lower()
    assert f"{report.pooled_risk:.3f}" not in after


def test_starting_a_second_flow_is_refused_not_just_greyed_out(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The one-flow-at-a-time rule is enforced in the code path, not only the UI.

    Greying out a button is a hint; the guard is what makes it an invariant. Each
    ``_start_*`` must refuse outright while any other flow is running.
    """
    monkeypatch.setenv("BT4_ASSP_FIXTURE_DIR", str(ASSP_FIXTURES))
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    window.protein_edit.setPlainText("MAALKHETQW")
    # Always consent, so a refusal cannot be mistaken for a declined dialog.
    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QtWidgets.QMessageBox.StandardButton.Yes),
    )

    # Pretend a cross-check is in flight, then try to start every other flow.
    window._set_crosscheck_running(True)
    for start in (
        window._start_optimize,
        window._start_candidates,
        window._start_library,
        window._start_crosscheck,
    ):
        start()
    assert window._thread is None
    assert window._cand_thread is None
    assert window._lib_thread is None
    assert window._cc_thread is None
    assert not window._optimize_running
    assert not window._cand_running
    assert not window._lib_running
    window._set_crosscheck_running(False)


# --------------------------------------------------------------------------- #
# Library (sampled) mode
# --------------------------------------------------------------------------- #


def test_library_worker_samples_and_renders() -> None:
    """A sampled library renders with its SAMPLED framing, never as a ranking."""
    app = QtWidgets.QApplication.instance()
    assert isinstance(app, QtWidgets.QApplication)

    worker = LibraryWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        n=6,
        temperature=1.0,
        seed=7,
    )
    result = worker.compute()
    assert len(result.results) == 6
    assert all(r.certificate.status.value == "sampled" for r in result.results)

    window = StudioWindow()
    window._on_library_finished(result)
    app.processEvents()

    assert window.library_table.rowCount() == 6
    banner = window.lib_banner.text().lower()
    assert "sampled, not optimized" in banner
    assert "not a ranking" in banner
    assert window.export_library_btn.isEnabled()
    # Selecting the first row shows that member's own sequence.
    assert window.library_view.toPlainText() == result.results[0].dna
    assert window.library_btn.isEnabled()


def test_second_library_draw_replaces_the_shown_member() -> None:
    """A second draw must not strand the first draw's sequence in the viewer.

    Repopulating the table in place leaves the old selection intact, so
    re-selecting row 0 emits no selection change -- the viewer has to be
    repainted explicitly or it would show one library's member beside another
    library's table.
    """
    window = StudioWindow()

    def draw(seed: int) -> api.LibraryResult:
        return LibraryWorker(
            "MAALKHETQW",
            api.OptimizeConfig(max_homopolymer=5),
            n=4,
            temperature=1.6,
            seed=seed,
        ).compute()

    first = draw(1)
    window._on_library_finished(first)
    assert window.library_view.toPlainText() == first.results[0].dna

    second = draw(999)
    assert second.results[0].dna != first.results[0].dna  # the draws differ
    window._on_library_finished(second)
    assert window.library_view.toPlainText() == second.results[0].dna


def test_library_is_deterministic_from_its_seed() -> None:
    """The same seed reproduces the same draw (invariant #7 reaches the UI)."""
    def draw() -> list[str]:
        result = LibraryWorker(
            "MAALKHETQW",
            api.OptimizeConfig(max_homopolymer=5),
            n=4,
            temperature=1.2,
            seed=11,
        ).compute()
        return [r.dna for r in result.results]

    assert draw() == draw()


def test_library_export_writes_every_member(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Exporting the library writes one labeled FASTA record per sampled member."""
    window = StudioWindow()
    result = LibraryWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        n=3,
        temperature=1.0,
        seed=3,
    ).compute()
    window._on_library_finished(result)
    window.jobname_edit.setText("job")

    path = tmp_path / "lib.fasta"
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: (str(path), "")),
    )
    window._export_library()

    text = path.read_text(encoding="utf-8")
    assert text.count(">") == 3
    # Each record names itself a sample, so a downstream reader cannot mistake a
    # sampled member for an optimized delivery.
    assert text.count("sampled") == 3
    flat = text.replace("\n", "")
    for member in result.results:
        assert member.dna in flat


def test_library_export_is_a_no_op_without_a_library(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """With no library drawn, the export writes nothing (no empty stub file)."""
    window = StudioWindow()
    path = tmp_path / "never.fasta"
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: (str(path), "")),
    )
    window._export_library()
    assert not path.exists()


def test_library_failure_resets_the_panel() -> None:
    """A failed draw clears the table so nothing stale stays exportable."""
    window = StudioWindow()
    result = LibraryWorker(
        "MAALKHETQW",
        api.OptimizeConfig(max_homopolymer=5),
        n=2,
        temperature=1.0,
        seed=1,
    ).compute()
    window._on_library_finished(result)
    assert window.export_library_btn.isEnabled()

    window._on_library_failed(ValueError("boom"))
    assert window.library_table.rowCount() == 0
    assert not window.export_library_btn.isEnabled()
    assert window._library is None
    assert window.library_btn.isEnabled()


# --------------------------------------------------------------------------- #
# Run gating and theming
# --------------------------------------------------------------------------- #


def test_only_one_engine_flow_runs_at_a_time() -> None:
    """Every start control is gated on "nothing is running", from one flag set."""
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    assert window.optimize_btn.isEnabled()
    assert window.rank_btn.isEnabled()
    assert window.library_btn.isEnabled()
    assert window.assp_btn.isEnabled()

    for setter in (
        window._set_running,
        window._set_candidates_running,
        window._set_library_running,
        window._set_crosscheck_running,
    ):
        setter(True)
        assert not window.optimize_btn.isEnabled()
        assert not window.rank_btn.isEnabled()
        assert not window.library_btn.isEnabled()
        assert not window.assp_btn.isEnabled()
        setter(False)
        assert window.optimize_btn.isEnabled()
        assert window.rank_btn.isEnabled()
        assert window.library_btn.isEnabled()
        assert window.assp_btn.isEnabled()

    # Cancel is live only during an optimization.
    assert not window.cancel_btn.isEnabled()
    window._set_running(True)
    assert window.cancel_btn.isEnabled()
    window._set_running(False)


def _action(window: StudioWindow, text: str) -> QtGui.QAction:
    """Find a menu action by its (mnemonic-stripped) label.

    Actions are parented to the window (not the menu bar), so search there.
    """
    for action in window.findChildren(QtGui.QAction):
        if action.text().replace("&", "") == text:
            return action
    raise AssertionError(f"no menu action named {text!r}")


def test_menu_actions_are_wired_to_the_real_flows() -> None:
    """Every action is reachable from the keyboard and triggers the actual slot."""
    window = StudioWindow()
    titles = [m.title().replace("&", "") for m in window.menuBar().findChildren(QtWidgets.QMenu)]
    assert titles == ["File", "Run", "View", "Help"]

    # Shortcuts are set, so the whole app is drivable without a mouse.
    assert _action(window, "Optimize").shortcut().toString() == "Ctrl+R"
    assert _action(window, "Sample library").shortcut().toString() == "Ctrl+L"

    # Triggering Optimize with an empty protein box reaches _start_optimize and
    # stops at its input validation -- proof the action is really connected.
    _action(window, "Optimize").trigger()
    assert "protein" in window.statusBar().currentMessage().lower()
    assert window.optimize_btn.isEnabled()  # nothing was started

    # The theme actions are an exclusive group with "System" checked by default.
    assert _action(window, "System").isChecked()
    _action(window, "Dark").trigger()
    assert window._dark is True
    assert _action(window, "Dark").isChecked()
    assert not _action(window, "System").isChecked()


def test_closing_mid_run_shuts_workers_down() -> None:
    """Closing the window while a flow runs must not destroy a running QThread.

    Qt aborts the process when a running thread is destroyed, and the threads are
    parented to the window -- so the close has to stop them first.
    """
    window = StudioWindow()
    window.protein_edit.setPlainText("MAALKHETQWCDEFGHIKLMNPQRSTVWY" * 2)
    window.lib_n_spin.setValue(24)
    window._start_library()
    assert window._lib_thread is not None

    thread = window._lib_thread
    window.close()
    assert not thread.isRunning()


def test_theme_switch_repaints_without_losing_results() -> None:
    """Switching theme restyles the surfaces and keeps the delivered result."""
    window = StudioWindow()
    frontier = api.frontier("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5), 5)
    window._on_finished(frontier)
    delivered = frontier.delivered()
    assert delivered is not None

    window._apply_theme("dark")
    assert window._dark is True
    assert window.sequence_view._dark is True
    assert window.sequence_view.toPlainText() == delivered.dna

    window._apply_theme("light")
    assert window._dark is False
    assert window.sequence_view._dark is False
    assert window.sequence_view.toPlainText() == delivered.dna
    assert window.export_fasta_btn.isEnabled()


def test_sequence_viewer_set_dark_repaints_highlights() -> None:
    """A theme change re-derives the violation bands from the held annotations."""
    viewer = SequenceViewer(dark=False)
    dna = "ACGTACGTACGT"
    hard = api.Violation("max_repeat", api.Severity.HARD, 0, 4, "dispersed repeat")
    viewer.set_sequence(dna, (hard,))
    light = viewer._colour(api.Severity.HARD).name()
    assert len(viewer.extraSelections()) == 1

    viewer.set_dark(True)
    # Same text, same annotations -- repainted in the dark palette.
    assert viewer.toPlainText() == dna
    assert len(viewer.extraSelections()) == 1
    assert viewer._violation_at(2) is not None
    assert viewer._colour(api.Severity.HARD).name() != light

    # Idempotent: setting the theme it already has changes nothing.
    viewer.set_dark(True)
    assert len(viewer.extraSelections()) == 1


def test_optimize_button_survives_optimize_then_rank() -> None:
    """An optimize then a candidate run must not leave Optimize stuck-disabled.

    Regression: the candidate flow gates Optimize on ``self._thread is None``; if
    the optimize thread ref is never cleared, Optimize stays disabled forever
    after an optimize-then-rank sequence.
    """
    window = StudioWindow()
    # Simulate an optimize having run and its thread having finished.
    window._thread = object()  # type: ignore[assignment]
    window._clear_optimize_thread()
    assert window._thread is None

    # A candidate run starts and finishes; Optimize must come back enabled.
    window._set_candidates_running(True)
    assert not window.optimize_btn.isEnabled()
    window._set_candidates_running(False)
    assert window.optimize_btn.isEnabled()


def test_candidates_failure_resets_panel() -> None:
    """A candidate-flow failure clears the table and re-enables the rank button."""
    window = StudioWindow()
    window._set_candidates_running(True)
    assert not window.rank_btn.isEnabled()

    window._on_cand_failed(api.InfeasibleError(["homopolymer"]))

    assert window.candidates_table.rowCount() == 0
    assert window.rank_btn.isEnabled()
    assert "satisfy these settings" in window.statusBar().currentMessage()


def test_tai_axis_is_offered_for_every_selectable_organism() -> None:
    """Every organism the user can pick now ships a tRNA table, so tAI is offered.

    E. coli was the last bundled organism without one -- tAI was unavailable exactly
    where translational selection is strongest. This asserts the gap is closed from
    the surface the user actually touches.
    """
    window = StudioWindow()
    for organism in api.available_organisms():
        window.organism_combo.setCurrentText(organism)
        window._update_tai_availability()
        assert window.tai_check.isEnabled(), organism


def test_tai_axis_disables_itself_when_a_table_is_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The guard must survive: a future organism without tRNA data must disable tAI.

    No bundled organism lacks a table any more, so the negative path is exercised
    against a stubbed availability list rather than being left uncovered.
    """
    window = StudioWindow()
    window.organism_combo.setCurrentText("homo_sapiens")
    window._update_tai_availability()
    window.tai_check.setChecked(True)

    monkeypatch.setattr(studio.api, "available_tai_organisms", lambda: ())
    window._update_tai_availability()
    assert not window.tai_check.isEnabled()
    assert not window.tai_check.isChecked()  # auto-unchecked when unavailable


def test_app_starts_on_the_engines_default_organism() -> None:
    """A freshly opened Studio must not default to a codon-commonness index.

    The combo is filled from ``available_organisms()``, which is sorted -- so
    without an explicit default it lands on *A. thaliana*, the one organism with
    no highly-expressed table. That silently opts the app out of the reference set
    every other surface defaults to.
    """
    window = StudioWindow()
    assert window.organism_combo.currentText() == api.OptimizeConfig().organism
    assert window.reference_combo.currentText() == api.default_reference_set(
        api.OptimizeConfig().organism
    )


def test_reference_set_survives_a_round_trip_through_a_one_set_organism() -> None:
    """Visiting an organism with only one reference set must not stick.

    Selecting *A. thaliana* forces the combo to ``genome_wide``. If that forced
    value were mistaken for a user preference, every later run on every other
    organism would quietly use codon commonness -- the same failure as the
    alphabetical default, two clicks away.
    """
    window = StudioWindow()
    window.organism_combo.setCurrentText("arabidopsis_thaliana")
    assert window.reference_combo.currentText() == "genome_wide"
    assert not window.reference_combo.isEnabled()

    window.organism_combo.setCurrentText("homo_sapiens")
    assert window.reference_combo.currentText() == "highly_expressed"
    assert window.reference_combo.isEnabled()


def test_an_explicit_reference_set_choice_is_remembered() -> None:
    """A real user pick, unlike a forced one, must carry across organisms."""
    window = StudioWindow()
    window.reference_combo.setCurrentText("genome_wide")
    window.organism_combo.setCurrentText("escherichia_coli")
    assert window.reference_combo.currentText() == "genome_wide"
    # ...and still survives a detour through the one-set organism.
    window.organism_combo.setCurrentText("arabidopsis_thaliana")
    window.organism_combo.setCurrentText("homo_sapiens")
    assert window.reference_combo.currentText() == "genome_wide"


def test_the_reference_set_tooltip_is_restored_not_stamped_once() -> None:
    """The one-set message must not outlive the organism it describes."""
    window = StudioWindow()
    explanatory = window.reference_combo.toolTip()
    window.organism_combo.setCurrentText("arabidopsis_thaliana")
    assert "Only the genome_wide" in window.reference_combo.toolTip()
    window.organism_combo.setCurrentText("homo_sapiens")
    assert window.reference_combo.toolTip() == explanatory


def test_delivered_tables_are_frozen_at_run_start() -> None:
    """Renders must label a result with the tables it was actually built from.

    The run is asynchronous; reading the combos when it finishes would attribute
    it to whatever the user happened to select in the meantime.
    """
    window = StudioWindow()
    window.protein_edit.setPlainText("MAALKHETQW")
    window.organism_combo.setCurrentText("escherichia_coli")
    window.reference_combo.setCurrentText("genome_wide")
    window._start_optimize()
    try:
        frozen = window._running_tables
        assert frozen == ("escherichia_coli", "genome_wide")

        # The user fiddles while the (already-launched) run is in flight...
        window.organism_combo.setCurrentText("homo_sapiens")
        assert window._running_tables == frozen
    finally:
        # closeEvent cancels the worker and joins the thread; letting a live
        # QThread be garbage-collected aborts the process.
        window.close()


# --------------------------------------------------------------------------- #
# Phase 1: presets, objective weights, epsilon-budgets, and the knobs that used
# to be hardcoded. The rule these pin down is "no dead control": every visible
# control must reach a real OptimizeConfig field.
# --------------------------------------------------------------------------- #


def test_no_preset_is_selected_at_launch() -> None:
    """BT4 is regime-agnostic, so the app must not pick a construct type for you."""
    window = StudioWindow()
    config, _steps = window._build_config(())
    assert window.preset_combo.currentIndex() == 0
    assert config.application_preset == ""
    assert config.gc_window_nt is None
    assert config.max_repeat_length is None


def test_new_design_controls_reach_the_config() -> None:
    """Every control added in this pass maps to its OptimizeConfig field."""
    window = StudioWindow()
    window.gc_window_spin.setValue(50)
    window.gc_window_min_spin.setValue(0.25)
    window.gc_window_max_spin.setValue(0.65)
    window.enzyme_sites_edit.setText("GANTC, CCWGG")
    window.weight_spins["gc_weight"].setValue(2.5)
    window.weight_spins["ramp_weight"].setValue(0.75)
    window.tandem_copies_spin.setValue(4)
    window.inverted_loop_spin.setValue(3)
    window.rc_check.setChecked(False)
    window.gc_min_spin.setValue(100)
    window.gc_max_spin.setValue(200)
    window.dinuc_combo.setCurrentIndex(1)  # CpG
    window.dinuc_min_spin.setValue(2)
    window.dinuc_max_spin.setValue(9)

    config, _steps = window._build_config(())

    assert config.gc_window_nt == 50
    assert (config.gc_window_min, config.gc_window_max) == (0.25, 0.65)
    assert config.restriction_extra_sites == ("GANTC", "CCWGG")
    assert config.gc_weight == 2.5
    assert config.ramp_weight == 0.75
    assert config.tandem_copies == 4
    assert config.inverted_loop == 3
    assert config.avoid_reverse_complement is False
    assert (config.gc_min, config.gc_max) == (100, 200)
    assert config.dinuc_budget == "CG"
    assert (config.dinuc_min, config.dinuc_max) == (2, 9)


def test_budget_controls_are_off_at_their_special_value() -> None:
    """The 'off' sentinel must mean None, not a literal -1 budget."""
    window = StudioWindow()
    config, _steps = window._build_config(())
    assert config.gc_min is None
    assert config.gc_max is None
    assert config.dinuc_budget is None
    assert config.dinuc_min is None
    assert config.dinuc_max is None


def test_codon_pair_weight_can_be_negative() -> None:
    """Codon-pair DE-optimization (attenuated-vaccine design) must be reachable."""
    window = StudioWindow()
    spin = window.weight_spins["cpb_weight"]
    assert spin.minimum() < 0, "a negative codon-pair weight must be selectable"


def test_choosing_a_preset_fills_the_visible_controls() -> None:
    """A preset writes into the controls so the user can see and edit what it did."""
    window = StudioWindow()
    keys = [window.preset_combo.itemData(i) for i in range(window.preset_combo.count())]
    index = keys.index("aav_transgene")
    window.preset_combo.setCurrentIndex(index)
    window._on_preset_chosen(index)

    assert window.gc_window_spin.value() == 50
    assert window.repeat_spin.value() == 20
    assert window.cpg_combo.currentText() == "deplete"
    assert window.splice_check.isChecked()
    assert window.uorf_check.isChecked()

    config, _steps = window._build_config(())
    assert config.application_preset == "aav_transgene"
    assert config.cpg_mode == "deplete"
    assert config.cpg_weight > 0


def test_a_user_edit_after_a_preset_wins() -> None:
    """A preset is a starting point, never a cage."""
    window = StudioWindow()
    keys = [window.preset_combo.itemData(i) for i in range(window.preset_combo.count())]
    index = keys.index("synthesis")
    window.preset_combo.setCurrentIndex(index)
    window._on_preset_chosen(index)
    window.homo_spin.setValue(4)  # the preset had set 6

    config, _steps = window._build_config(())
    assert config.max_homopolymer == 4


def test_choosing_none_preset_changes_nothing() -> None:
    window = StudioWindow()
    window.homo_spin.setValue(9)
    window._on_preset_chosen(0)  # the "(none)" entry
    assert window.homo_spin.value() == 9


def test_open_protein_fasta_loads_the_first_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Paste-only input was a named gap; a multi-record file says what it loaded."""
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">myprot description\nMAALKHETQWY\n>second\nMKV\n")
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *a, **k: (str(fasta), "")),
    )
    window = StudioWindow()
    window._open_protein()
    assert window.protein_edit.toPlainText() == "MAALKHETQWY"
    assert window.jobname_edit.text() == "myprot"
    assert "first of 2 records" in window.statusBar().currentMessage()


def test_validate_panel_reports_violations(monkeypatch: pytest.MonkeyPatch) -> None:
    """api.validate had no UI at all; this is the surface for an existing CDS."""
    monkeypatch.setattr(
        QtWidgets.QInputDialog,
        "getMultiLineText",
        staticmethod(lambda *a, **k: ("GGGGGGGGGGAAAA", True)),
    )
    seen: dict[str, str] = {}

    class _Box(QtWidgets.QMessageBox):
        def exec(self) -> int:
            seen["text"] = self.text()
            seen["detail"] = self.detailedText()
            return 0

    monkeypatch.setattr(QtWidgets, "QMessageBox", _Box)
    window = StudioWindow()
    window.homo_spin.setValue(4)
    window._validate_sequence()
    assert "NOT feasible" in seen["text"]
    assert "homopolymer" in seen["detail"]


def test_tracks_plot_includes_the_splice_track_and_labels_calibration() -> None:
    """The splice risk track was missing entirely; it must arrive labelled."""
    window = StudioWindow()
    result = api.optimize("MAALKHETQWYCDEFGHIKLM", api.OptimizeConfig(max_homopolymer=None))
    window._render_tracks(result)
    title = window.tracks_plot.plotItem.titleLabel.text
    assert "splice" in title
    # The shipped baseline is not calibrated, so the plot must say so.
    assert "UNCALIBRATED" in title


# --------------------------------------------------------------------------- #
# Phase 2/3: construct context in the app.
# --------------------------------------------------------------------------- #


def test_no_construct_context_by_default() -> None:
    """An empty context box must mean None, so a plain run is unchanged."""
    window = StudioWindow()
    config, _steps = window._build_config(())
    assert config.context is None
    assert config.context_provenance == "omit"


def test_construct_context_reaches_the_engine_and_changes_the_design() -> None:
    """The headline capability, driven from the app exactly as a user would."""
    window = StudioWindow()
    window.ctx_upstream_edit.setPlainText("ggcacca")
    window.ctx_downstream_edit.setPlainText("AAATTT")
    window.motifs_edit.setText("CCAGTG")
    window.rc_check.setChecked(False)
    window.homo_spin.setValue(0)

    config, _steps = window._build_config(())
    assert config.context is not None
    assert config.context.upstream == "GGCACCA"
    assert config.context.downstream == "AAATTT"

    result = api.optimize("VK", config)
    construct = config.context.assemble(result.dna)
    # The forbidden motif exists only across the junction, so avoiding it is proof
    # the design saw the flank.
    assert "CCAGTG" not in construct


def test_context_provenance_policy_is_selectable() -> None:
    window = StudioWindow()
    window.ctx_upstream_edit.setPlainText("GGCACCA")
    window.context_prov_combo.setCurrentIndex(1)
    config, _steps = window._build_config(())
    assert config.context_provenance == "hash"


def test_the_construct_context_box_is_not_the_ribonn_utr_box() -> None:
    """Two different 'UTR' inputs coexist; they must not be the same widget.

    The RiboNN boxes are model context (they annotate a score); the Design-panel
    boxes make the DESIGN construct-aware. Wiring one to the other would silently
    make a user's backbone change an expression annotation and nothing else.
    """
    window = StudioWindow()
    assert window.ctx_upstream_edit is not window.utr5_edit
    assert window.ctx_downstream_edit is not window.utr3_edit
    # The RiboNN tooltip points at the other box so the distinction is discoverable.
    assert "5' context" in window.utr5_edit.toolTip()


# --------------------------------------------------------------------------- #
# The metrics table is driven by the audit, not a fixed row list.
# --------------------------------------------------------------------------- #


def test_audit_rows_surface_partial_enforcement() -> None:
    """A rule the engine could only partly enforce must be visible in the GUI.

    The CLI has always printed this; the GUI showed nothing, so a user who set
    "Max repeat length" and got unremovable residuals saw only an unexplained
    hard-violation count.
    """
    rows = dict(
        studio._audit_rows(
            {"max_repeat_enforced": "partial", "max_repeat_residual": 19}
        )
    )
    assert "could not be removed" in rows["Max repeat"]
    assert "19" in rows["Max repeat"]


def test_audit_rows_report_clean_enforcement_too() -> None:
    rows = dict(
        studio._audit_rows({"uorf_enforced": "clean", "uorf_residual": 0})
    )
    assert rows["uORF"] == "clean (fully enforced)"


def test_audit_rows_surface_an_unknown_future_rule() -> None:
    """A rule added to the engine later must appear without editing the GUI.

    This is the property that matters: the poly(A) and windowed-GC rules would
    otherwise have been silently missing from the table while the CLI printed them.
    """
    rows = dict(
        studio._audit_rows(
            {"brand_new_rule_enforced": "partial", "brand_new_rule_residual": 3}
        )
    )
    assert "Brand new rule" in rows
    assert "3 could not be removed" in rows["Brand new rule"]


def test_audit_rows_flag_an_uncalibrated_folding_number() -> None:
    rows = dict(
        studio._audit_rows({"folding_dg": -41.0, "folding_calibrated": False})
    )
    assert "UNCALIBRATED" in rows["5' folding"]
    calibrated = dict(
        studio._audit_rows({"folding_dg": -12.5, "folding_calibrated": True})
    )
    assert "kcal/mol" in calibrated["5' folding"]
    assert "UNCALIBRATED" not in calibrated["5' folding"]


def test_audit_rows_name_relaxed_rules() -> None:
    rows = dict(studio._audit_rows({"relaxed_constraints": ["internal_start"]}))
    assert "internal_start" in rows["Relaxed rules"]


def test_audit_rows_are_empty_for_a_plain_run() -> None:
    assert studio._audit_rows({"cai": 1.0, "gc_percent": 55.0}) == []


def test_metrics_table_grows_with_the_audit() -> None:
    """End-to-end: the rendered table gains rows for what the run reported."""
    window = StudioWindow()
    plain = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    window._render_metrics(plain)
    baseline = window.metrics_table.rowCount()

    guarded = api.optimize(
        "MWWWWWWWWKLDE",
        api.OptimizeConfig(max_repeat_length=6, max_homopolymer=None),
    )
    window._render_metrics(guarded)
    assert window.metrics_table.rowCount() > baseline
    shown = {
        window.metrics_table.item(r, 0).text(): window.metrics_table.item(r, 1).text()
        for r in range(window.metrics_table.rowCount())
    }
    assert "could not be removed" in shown["Max repeat"]


def test_self_test_engine_runs_a_real_design() -> None:
    """``--self-test`` designs real sequences, not just a window.

    This is the check that runs inside the *frozen* bundle, where the v0.5.0 defect
    (a packaged data file absent from the freeze) lived, invisible to a suite that
    reads the same files straight off the source tree.
    """
    studio._self_test_engine()  # raises on failure; silent on success


def test_self_test_case_is_not_vacuous() -> None:
    """The rule under test must actually change the delivered sequence.

    The first version of this check asked only whether the delivered DNA lacked an
    EcoRI site -- and the unguarded optimum for its protein never contained one, so it
    passed identically whether or not the restriction rule was applied at all. That is
    the failure mode worth pinning: assert here that the pair is live, so the day a
    codon-table change makes it vacuous, this says so rather than the self-test quietly
    proving nothing.
    """
    plain = api.OptimizeConfig(organism="homo_sapiens", max_homopolymer=5, seed=0)
    guarded = api.OptimizeConfig(
        organism="homo_sapiens",
        max_homopolymer=5,
        seed=0,
        restriction_enzymes=(studio._SELF_TEST_ENZYME,),
    )
    unguarded_dna = api.optimize(studio._SELF_TEST_PROTEIN, plain).dna
    guarded_dna = api.optimize(studio._SELF_TEST_PROTEIN, guarded).dna
    assert studio._SELF_TEST_SITE in unguarded_dna
    assert studio._SELF_TEST_SITE not in guarded_dna
    assert unguarded_dna != guarded_dna


def test_self_test_engine_rejects_a_malformed_design(monkeypatch: pytest.MonkeyPatch) -> None:
    """A bad design fails the self-test.

    Pins the guard itself. A smoke check that passes whatever the engine returns would
    have let the broken bundle through just as quietly as no check at all.
    """

    class _Stub:
        dna = "ATG"  # far too short for the self-test protein

    monkeypatch.setattr(studio.api, "optimize", lambda *a, **k: _Stub())
    with pytest.raises(RuntimeError, match="expected 228"):
        studio._self_test_engine()


def test_self_test_engine_rejects_an_unremoved_site(monkeypatch: pytest.MonkeyPatch) -> None:
    """A guarded design still carrying the site fails the self-test."""
    real = api.optimize

    def _ignore_the_rule(protein: str, config: object = None) -> object:
        # Both solves come back unguarded -- i.e. the restriction rule did nothing.
        return real(protein, studio.api.OptimizeConfig(
            organism="homo_sapiens", max_homopolymer=5, seed=0
        ))

    monkeypatch.setattr(studio.api, "optimize", _ignore_the_rule)
    with pytest.raises(RuntimeError, match="survived a run"):
        studio._self_test_engine()


def test_self_test_engine_reports_a_vacuous_case(monkeypatch: pytest.MonkeyPatch) -> None:
    """If the unguarded optimum stops carrying the site, the self-test says so.

    The check must fail *loudly* when it can no longer prove anything, rather than
    passing on a comparison that has become trivially true.
    """
    real = api.optimize

    def _always_guarded(protein: str, config: object = None) -> object:
        return real(protein, studio.api.OptimizeConfig(
            organism="homo_sapiens",
            max_homopolymer=5,
            seed=0,
            restriction_enzymes=(studio._SELF_TEST_ENZYME,),
        ))

    monkeypatch.setattr(studio.api, "optimize", _always_guarded)
    with pytest.raises(RuntimeError, match="gone vacuous"):
        studio._self_test_engine()


def test_preset_reports_a_field_it_cannot_apply() -> None:
    """A preset field with no control is named as *not applied*, never silently dropped.

    `refine` used to be registered with a no-op setter, which kept it out of the
    "unmapped" list and so out of the message -- so choosing the IVT mRNA preset, which
    asks for refinement, produced a design without it and said nothing. The app's own
    rule is that a preset never sets something invisibly; a no-op setter is the one way
    to defeat it.
    """
    window = StudioWindow()
    ivt = next(
        (i for i in range(window.preset_combo.count())
         if "mRNA" in window.preset_combo.itemText(i)),
        None,
    )
    assert ivt is not None, "no IVT mRNA preset in the combo"
    window.preset_combo.setCurrentIndex(ivt)
    window._on_preset_chosen(ivt)
    message = window.statusBar().currentMessage()
    assert "Not applied" in message
    assert "refine" in message


def test_about_box_names_the_version() -> None:
    """A released app must be able to say which release it is.

    Without it, a user reporting "the app does X" gives no way to tell which build
    they are on -- and 0.5.0 changed what the engine delivers by default (the CAI
    reference set), so the answer matters.
    """
    window = StudioWindow()
    window._show_about()
    assert bt4.__version__ in window._msgbox.text()
