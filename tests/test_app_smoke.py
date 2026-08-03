"""Headless smoke tests for the BT4 Studio desktop app.

These run under Qt's ``offscreen`` platform (no display needed). They drive the
engine synchronously through the worker's ``compute`` and call the window's
result/failure slots directly -- no real ``QThread`` is started and the event
loop is never entered, so the suite stays fast and hermetic.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

pytest.importorskip("PySide6")

from PySide6 import QtWidgets

from bt4 import api
from bt4.app.studio import SequenceViewer, StudioWindow
from bt4.app.worker import OptimizeWorker


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


def test_tai_axis_tracks_organism_availability() -> None:
    """The tAI checkbox is enabled only for organisms with a bundled tRNA table."""
    window = StudioWindow()

    window.organism_combo.setCurrentText("homo_sapiens")
    window._update_tai_availability()
    assert window.tai_check.isEnabled()
    window.tai_check.setChecked(True)

    window.organism_combo.setCurrentText("escherichia_coli")
    window._update_tai_availability()
    assert not window.tai_check.isEnabled()
    assert not window.tai_check.isChecked()  # auto-unchecked when unavailable
