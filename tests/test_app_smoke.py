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
from bt4.app.studio import StudioWindow
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


def test_infeasible_is_handled() -> None:
    """The failure slot never raises and leaves the window usable."""
    window = StudioWindow()
    window._set_running(True)
    assert not window.optimize_btn.isEnabled()

    window._on_failed("no feasible codon")

    assert window.optimize_btn.isEnabled()
    assert "no feasible codon" in window.statusBar().currentMessage()
