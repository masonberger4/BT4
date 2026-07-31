"""BT4 Studio main window and application entry point.

``StudioWindow`` is the whole desktop UI: a left control panel (protein, target
organism, and per-run design knobs) and a right results panel (an honest
optimality-certificate badge, a recomputed-metrics table, an interactive
CAI/GC frontier scatter, and a monospaced sequence viewer, plus FASTA/JSON
export). Optimizations run on a background :class:`~bt4.app.worker.OptimizeWorker`
moved onto a ``QThread`` so the window never blocks, and any engine error is
surfaced as a non-fatal message rather than a crash.

This module talks to the engine exclusively through :mod:`bt4.api`; it never
imports the optimizer, pipeline, or biomodels directly. Everything is computed
locally and nothing leaves the machine.
"""

from __future__ import annotations

import sys
from itertools import pairwise

import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from bt4 import api
from bt4.app import theme
from bt4.app.worker import OptimizeWorker

__all__ = ["StudioWindow", "main"]

_METRIC_ROWS = (
    "CAI",
    "GC %",
    "Length (nt)",
    "Scored codons",
    "Hard violations",
    "Soft violations",
    "Optimality",
    "Solver",
)


def _is_dark() -> bool:
    """Guess whether the running app uses a dark theme from its window colour."""
    app = QtWidgets.QApplication.instance()
    if isinstance(app, QtWidgets.QApplication):
        window = app.palette().color(QtGui.QPalette.ColorRole.Window)
        return window.lightness() < 128
    return False


class StudioWindow(QtWidgets.QMainWindow):
    """The BT4 Studio main window: controls on the left, results on the right."""

    def __init__(self) -> None:
        """Build the full UI (no optimization is started here)."""
        super().__init__()
        self.setWindowTitle("BT4 Studio")
        self.resize(1180, 760)

        self._dark = _is_dark()
        self._last: api.FrontierResult | None = None
        self._thread: QtCore.QThread | None = None
        self._worker: OptimizeWorker | None = None
        self._msgbox: QtWidgets.QMessageBox | None = None

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_controls())
        splitter.addWidget(self._build_results())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([380, 800])
        self.setCentralWidget(splitter)

        self.setStatusBar(QtWidgets.QStatusBar())
        self.statusBar().showMessage("Ready. Enter a protein and click Optimize.")

        self._set_tab_order()
        self._reset_results()

    # ---- construction -----------------------------------------------------

    def _build_controls(self) -> QtWidgets.QWidget:
        """Build the left input/control panel."""
        box = QtWidgets.QGroupBox("Design")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)

        self.protein_edit = QtWidgets.QPlainTextEdit()
        self.protein_edit.setPlaceholderText("Paste a protein, e.g. MAALKHETQW")
        self.protein_edit.setAccessibleName("Protein sequence")
        self.protein_edit.setMinimumHeight(120)
        self._add_row(form, "Protein", self.protein_edit)

        self.organism_combo = QtWidgets.QComboBox()
        self.organism_combo.addItems(list(api.available_organisms()))
        self.organism_combo.setAccessibleName("Target organism")
        self._add_row(form, "Organism", self.organism_combo)

        self.jobname_edit = QtWidgets.QLineEdit()
        self.jobname_edit.setPlaceholderText("optional job name")
        self.jobname_edit.setAccessibleName("Job name")
        self._add_row(form, "Job name", self.jobname_edit)

        self.gc_spin = QtWidgets.QDoubleSpinBox()
        self.gc_spin.setRange(0.0, 1.0)
        self.gc_spin.setSingleStep(0.05)
        self.gc_spin.setValue(0.55)
        self.gc_spin.setAccessibleName("GC target fraction")
        self._add_row(form, "GC target", self.gc_spin)

        self.homo_spin = QtWidgets.QSpinBox()
        self.homo_spin.setRange(0, 20)
        self.homo_spin.setValue(6)
        self.homo_spin.setSpecialValueText("off")
        self.homo_spin.setAccessibleName("Maximum homopolymer run (0 = off)")
        self._add_row(form, "Max homopolymer", self.homo_spin)

        self.motifs_edit = QtWidgets.QLineEdit()
        self.motifs_edit.setPlaceholderText("comma-separated, e.g. GAATTC,GGATCC")
        self.motifs_edit.setAccessibleName("Forbidden motifs")
        self._add_row(form, "Forbidden motifs", self.motifs_edit)

        self.enzymes_edit = QtWidgets.QLineEdit()
        self.enzymes_edit.setPlaceholderText("comma-separated, e.g. EcoRI,BamHI")
        self.enzymes_edit.setAccessibleName("Restriction enzymes to avoid")
        self._add_row(form, "Restriction sites", self.enzymes_edit)

        self.cpg_combo = QtWidgets.QComboBox()
        self.cpg_combo.addItems(["off", "deplete", "elevate"])
        self.cpg_combo.setAccessibleName("CpG dinucleotide mode")
        self._add_row(form, "CpG", self.cpg_combo)

        self.steps_spin = QtWidgets.QSpinBox()
        self.steps_spin.setRange(1, 25)
        self.steps_spin.setValue(9)
        self.steps_spin.setAccessibleName("Frontier steps")
        self._add_row(form, "Frontier steps", self.steps_spin)

        self.beam_spin = QtWidgets.QSpinBox()
        self.beam_spin.setRange(0, 256)
        self.beam_spin.setValue(0)
        self.beam_spin.setSpecialValueText("exact")
        self.beam_spin.setAccessibleName("Beam width (0 = exact)")
        self._add_row(form, "Beam width", self.beam_spin)

        self.optimize_btn = QtWidgets.QPushButton("Optimize")
        self.optimize_btn.setAccessibleName("Run optimization")
        self.optimize_btn.clicked.connect(self._start_optimize)
        form.addRow(self.optimize_btn)

        return box

    def _add_row(
        self, form: QtWidgets.QFormLayout, text: str, widget: QtWidgets.QWidget
    ) -> None:
        """Add a labelled row whose label is the accessibility buddy of the widget."""
        label = QtWidgets.QLabel(text)
        label.setBuddy(widget)
        form.addRow(label, widget)

    def _build_results(self) -> QtWidgets.QWidget:
        """Build the right results panel."""
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)

        self.badge = QtWidgets.QLabel()
        self.badge.setObjectName("certBadge")
        self.badge.setWordWrap(True)
        self.badge.setAccessibleName("Optimality certificate")
        self.badge.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.badge)

        upper = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)

        self.metrics_table = QtWidgets.QTableWidget(0, 2)
        self.metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.metrics_table.verticalHeader().setVisible(False)
        self.metrics_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.metrics_table.horizontalHeader().setStretchLastSection(True)
        self.metrics_table.setAccessibleName("Delivered metrics")
        self.metrics_table.setMinimumWidth(280)
        upper.addWidget(self.metrics_table)

        self.plot = pg.PlotWidget()
        self._style_plot()
        self.plot.setAccessibleName("Pareto frontier: CAI versus GC")
        upper.addWidget(self.plot)
        upper.setStretchFactor(0, 0)
        upper.setStretchFactor(1, 1)
        layout.addWidget(upper, stretch=1)

        seq_label = QtWidgets.QLabel("Delivered coding DNA")
        layout.addWidget(seq_label)
        self.sequence_view = QtWidgets.QPlainTextEdit()
        self.sequence_view.setReadOnly(True)
        self.sequence_view.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        self.sequence_view.setFont(mono)
        self.sequence_view.setAccessibleName("Delivered sequence")
        self.sequence_view.setMaximumHeight(140)
        seq_label.setBuddy(self.sequence_view)
        layout.addWidget(self.sequence_view)

        exports = QtWidgets.QHBoxLayout()
        self.export_fasta_btn = QtWidgets.QPushButton("Export FASTA")
        self.export_fasta_btn.setAccessibleName("Export FASTA")
        self.export_fasta_btn.clicked.connect(self._export_fasta)
        self.export_json_btn = QtWidgets.QPushButton("Export JSON")
        self.export_json_btn.setAccessibleName("Export JSON")
        self.export_json_btn.clicked.connect(self._export_json)
        exports.addWidget(self.export_fasta_btn)
        exports.addWidget(self.export_json_btn)
        exports.addStretch(1)
        layout.addLayout(exports)

        return panel

    def _style_plot(self) -> None:
        """Apply theme-aware colours and axis labels to the frontier plot."""
        bg = "#242832" if self._dark else "#ffffff"
        fg = "#e6e9ed" if self._dark else "#1b1f24"
        self.plot.setBackground(bg)
        self.plot.setLabel("bottom", "CAI")
        self.plot.setLabel("left", "GC fraction")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        for edge in ("bottom", "left"):
            axis = self.plot.getAxis(edge)
            axis.setPen(fg)
            axis.setTextPen(fg)
            # Show raw CAI / GC-fraction values, not a rescaled "x0.001" SI prefix.
            axis.enableAutoSIPrefix(False)

    def _set_tab_order(self) -> None:
        """Wire a sensible keyboard tab order through the controls."""
        order = (
            self.protein_edit,
            self.organism_combo,
            self.jobname_edit,
            self.gc_spin,
            self.homo_spin,
            self.motifs_edit,
            self.enzymes_edit,
            self.cpg_combo,
            self.steps_spin,
            self.beam_spin,
            self.optimize_btn,
            self.export_fasta_btn,
            self.export_json_btn,
        )
        for first, second in pairwise(order):
            self.setTabOrder(first, second)

    # ---- run --------------------------------------------------------------

    def _read_config(self) -> tuple[str, api.OptimizeConfig, int]:
        """Read the controls into a protein, an OptimizeConfig, and a step count."""
        protein = "".join(self.protein_edit.toPlainText().split()).upper()
        homo = self.homo_spin.value()
        beam = self.beam_spin.value()
        motifs = tuple(
            m.strip().upper() for m in self.motifs_edit.text().split(",") if m.strip()
        )
        enzymes = tuple(
            e.strip() for e in self.enzymes_edit.text().split(",") if e.strip()
        )
        cpg = self.cpg_combo.currentText()
        cpg_weight = 0.0 if cpg == "off" else 1.0
        cpg_mode = "deplete" if cpg == "off" else cpg
        config = api.OptimizeConfig(
            organism=self.organism_combo.currentText(),
            gc_target=self.gc_spin.value(),
            max_homopolymer=homo if homo > 0 else None,
            forbidden_motifs=motifs,
            avoid_reverse_complement=True,
            restriction_enzymes=enzymes,
            cpg_weight=cpg_weight,
            cpg_mode=cpg_mode,
            beam=beam if beam > 0 else None,
            seed=0,
        )
        return protein, config, self.steps_spin.value()

    def _start_optimize(self) -> None:
        """Launch a frontier optimization on a background thread."""
        protein, config, steps = self._read_config()
        if not protein:
            self.statusBar().showMessage("Enter a protein sequence to optimize.")
            return

        self._set_running(True)
        self.statusBar().showMessage("Optimizing...")

        thread = QtCore.QThread(self)
        worker = OptimizeWorker(protein, config, steps)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(self._on_finished)
        worker.failed.connect(self._on_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        # Keep references so neither is garbage-collected mid-run.
        self._thread = thread
        self._worker = worker
        thread.start()

    def _set_running(self, running: bool) -> None:
        """Toggle controls while a run is in flight so the UI stays consistent."""
        self.optimize_btn.setEnabled(not running)
        self.optimize_btn.setText("Optimizing..." if running else "Optimize")

    # ---- slots ------------------------------------------------------------

    @QtCore.Slot(int, str)
    def _on_progress(self, value: int, label: str) -> None:
        """Show worker progress in the status bar."""
        self.statusBar().showMessage(f"{label} ({value}%)")

    @QtCore.Slot(object)
    def _on_finished(self, result: api.FrontierResult) -> None:
        """Populate the results panel from a finished frontier optimization."""
        self._last = result
        self._set_running(False)
        self._populate(result)

    @QtCore.Slot(str)
    def _on_failed(self, message: str) -> None:
        """Re-enable controls and show a non-blocking warning on any engine error."""
        self._set_running(False)
        self.statusBar().showMessage(f"Optimization failed: {message}")
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
        box.setWindowTitle("Optimization failed")
        box.setText("The optimization could not complete.")
        box.setInformativeText(message)
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        # show() (not exec()) so the app never blocks -- important headless too.
        self._msgbox = box
        box.show()

    # ---- rendering --------------------------------------------------------

    def _reset_results(self) -> None:
        """Clear the results panel to an honest empty state."""
        self.badge.setText("No optimization run yet.")
        self.badge.setStyleSheet(theme.badge_qss(""))
        self.metrics_table.setRowCount(0)
        self.sequence_view.setPlainText("")
        self.plot.clear()
        self.export_fasta_btn.setEnabled(False)
        self.export_json_btn.setEnabled(False)

    def _populate(self, result: api.FrontierResult) -> None:
        """Render badge, metrics, sequence, and frontier for a delivered result."""
        delivered = result.delivered()
        if delivered is None:
            self._reset_results()
            self.badge.setText("No feasible sequence was delivered.")
            self.badge.setStyleSheet(theme.badge_qss("relaxed"))
            self.statusBar().showMessage("No feasible result on the frontier.")
            return

        self._render_badge(delivered)
        self._render_metrics(delivered)
        self.sequence_view.setPlainText(delivered.dna)
        self._render_frontier(result)

        self.export_fasta_btn.setEnabled(True)
        self.export_json_btn.setEnabled(True)

        cai = float(delivered.audit["cai"])
        gc_pct = delivered.metrics.gc * 100.0
        self.statusBar().showMessage(
            f"delivered CAI {cai:.2f} at GC {gc_pct:.0f}% "
            f"- frontier of {len(result.results)} points"
        )

    def _render_badge(self, delivered: api.Result) -> None:
        """Show the delivered result's optimality certificate, honestly coloured."""
        cert = delivered.certificate
        title = cert.status.value.replace("_", " ").upper()
        text = f"{title} - {cert.solver}"
        if cert.detail:
            text += f"\n{cert.detail}"
        self.badge.setText(text)
        self.badge.setStyleSheet(theme.badge_qss(cert.status.value))

    def _render_metrics(self, delivered: api.Result) -> None:
        """Fill the metrics table from the delivered result's recomputed metrics."""
        metrics = delivered.metrics
        cert = delivered.certificate
        values = {
            "CAI": f"{float(delivered.audit['cai']):.4f}",
            "GC %": f"{metrics.gc * 100.0:.2f}",
            "Length (nt)": str(metrics.length_nt),
            "Scored codons": str(delivered.audit["n_scored_codons"]),
            "Hard violations": str(metrics.hard_violations),
            "Soft violations": str(metrics.soft_violations),
            "Optimality": cert.status.value,
            "Solver": cert.solver,
        }
        self.metrics_table.setRowCount(len(_METRIC_ROWS))
        for row, name in enumerate(_METRIC_ROWS):
            self.metrics_table.setItem(row, 0, QtWidgets.QTableWidgetItem(name))
            self.metrics_table.setItem(row, 1, QtWidgets.QTableWidgetItem(values[name]))
        self.metrics_table.resizeColumnsToContents()

    def _render_frontier(self, result: api.FrontierResult) -> None:
        """Draw the CAI/GC frontier scatter with the delivered point highlighted."""
        self.plot.clear()
        results = result.results
        xs = [float(r.audit["cai"]) for r in results]
        ys = [r.metrics.gc for r in results]
        if xs:
            self.plot.plot(
                xs,
                ys,
                pen=None,
                symbol="o",
                symbolSize=11,
                symbolBrush=pg.mkBrush("#4a83e0"),
                symbolPen=pg.mkPen("#ffffff"),
            )
        chosen = result.frontier.chosen
        if 0 <= chosen < len(xs):
            self.plot.plot(
                [xs[chosen]],
                [ys[chosen]],
                pen=None,
                symbol="star",
                symbolSize=24,
                symbolBrush=pg.mkBrush("#e0b64a"),
                symbolPen=pg.mkPen("#1b1f24"),
            )

    # ---- export -----------------------------------------------------------

    def _delivered(self) -> api.Result | None:
        """The currently delivered result, if any."""
        return self._last.delivered() if self._last is not None else None

    def _export_fasta(self) -> None:
        """Write the delivered sequence to a FASTA file chosen by the user."""
        delivered = self._delivered()
        if delivered is None:
            return
        header = self.jobname_edit.text().strip() or "bt4"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export FASTA", f"{header}.fasta", "FASTA (*.fasta *.fa);;All files (*)"
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(api.to_fasta(delivered.dna, header=header))
        self.statusBar().showMessage(f"Wrote FASTA to {path}")

    def _export_json(self) -> None:
        """Write the delivered result (with manifest) to a JSON file."""
        delivered = self._delivered()
        if delivered is None:
            return
        header = self.jobname_edit.text().strip() or "bt4"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export JSON", f"{header}.json", "JSON (*.json);;All files (*)"
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(api.result_to_json(delivered))
        self.statusBar().showMessage(f"Wrote JSON to {path}")


def main() -> int:
    """Create (or reuse) the QApplication, show BT4 Studio, and run the loop.

    Returns:
        The Qt event-loop exit code.
    """
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    if isinstance(app, QtWidgets.QApplication):
        app.setStyleSheet(theme.stylesheet(_is_dark()))
    window = StudioWindow()
    window.show()
    return app.exec()
