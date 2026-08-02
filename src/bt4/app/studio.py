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
        controls_scroll = QtWidgets.QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setWidget(self._build_controls())
        controls_scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        splitter.addWidget(controls_scroll)
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
        """Build the left input/control panel.

        Every control carries a hover tooltip explaining what its variable does,
        so a user can discover the design knobs without leaving the app.
        """
        box = QtWidgets.QGroupBox("Design")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)

        self.protein_edit = QtWidgets.QPlainTextEdit()
        self.protein_edit.setPlaceholderText("Paste a protein, e.g. MAALKHETQW")
        self.protein_edit.setAccessibleName("Protein sequence")
        self.protein_edit.setMinimumHeight(120)
        self._add_row(
            form, "Protein", self.protein_edit,
            "The amino-acid sequence to back-translate (one-letter codes; "
            "whitespace is ignored). A stop codon is appended automatically.",
        )

        self.organism_combo = QtWidgets.QComboBox()
        self.organism_combo.addItems(list(api.available_organisms()))
        self.organism_combo.setAccessibleName("Target organism")
        self._add_row(
            form, "Organism", self.organism_combo,
            "Codon-usage table used for CAI and codon choice; it decides which "
            "synonymous codons are preferred for the target organism.",
        )

        self.jobname_edit = QtWidgets.QLineEdit()
        self.jobname_edit.setPlaceholderText("optional job name")
        self.jobname_edit.setAccessibleName("Job name")
        self._add_row(
            form, "Job name", self.jobname_edit,
            "Optional label for this run, used only for your reference and in "
            "exported file headers.",
        )

        self.gc_spin = QtWidgets.QDoubleSpinBox()
        self.gc_spin.setRange(0.0, 1.0)
        self.gc_spin.setSingleStep(0.05)
        self.gc_spin.setValue(0.55)
        self.gc_spin.setAccessibleName("GC target fraction")
        self._add_row(
            form, "GC target", self.gc_spin,
            "Desired overall GC fraction (0-1); the GC-proximity objective pulls "
            "the sequence toward this value.",
        )

        self.homo_spin = QtWidgets.QSpinBox()
        self.homo_spin.setRange(0, 20)
        self.homo_spin.setValue(6)
        self.homo_spin.setSpecialValueText("off")
        self.homo_spin.setAccessibleName("Maximum homopolymer run (0 = off)")
        self._add_row(
            form, "Max homopolymer", self.homo_spin,
            "Longest allowed run of a single repeated base (e.g. AAAA). Long "
            "single-base runs cause synthesis errors and polymerase slippage. "
            "0 = off.",
        )

        self.gc_run_spin = QtWidgets.QSpinBox()
        self.gc_run_spin.setRange(0, 40)
        self.gc_run_spin.setValue(0)
        self.gc_run_spin.setSpecialValueText("off")
        self.gc_run_spin.setAccessibleName("Maximum GC-run length (0 = off)")
        self._add_row(
            form, "Max GC length", self.gc_run_spin,
            "Longest allowed run of consecutive strong (G or C) bases, mixed runs "
            "included (e.g. GCGCGC counts as 6). Long GC stretches form stable "
            "structure and are hard to synthesize. 0 = off.",
        )

        self.repeat_spin = QtWidgets.QSpinBox()
        self.repeat_spin.setRange(0, 40)
        self.repeat_spin.setValue(0)
        self.repeat_spin.setSpecialValueText("off")
        self.repeat_spin.setAccessibleName("Maximum repeat length (0 = off)")
        self._add_row(
            form, "Max repeat length", self.repeat_spin,
            "Longest allowed repeated substring anywhere in the sequence -- "
            "direct, inverted, or palindromic (reverse-complement aware). Longer "
            "internal repeats mis-prime PCR and seed recombination. This is a "
            "non-local rule enforced by a refinement pass, so the result is "
            "labeled heuristic (not proven-optimal), and any repeats it cannot "
            "remove are reported honestly. 0 = off.",
        )

        self.motifs_edit = QtWidgets.QLineEdit()
        self.motifs_edit.setPlaceholderText("comma-separated, e.g. GAATTC,GGATCC")
        self.motifs_edit.setAccessibleName("Forbidden motifs")
        self._add_row(
            form, "Forbidden motifs", self.motifs_edit,
            "Extra DNA substrings to ban, comma-separated; each motif's reverse "
            "complement is banned too.",
        )

        self._add_forbidden_presets(form)

        self.enzymes_edit = QtWidgets.QLineEdit()
        self.enzymes_edit.setPlaceholderText("comma-separated, e.g. EcoRI,BamHI")
        self.enzymes_edit.setAccessibleName("Restriction enzymes to avoid")
        self._add_row(
            form, "Restriction sites", self.enzymes_edit,
            "Restriction enzymes whose recognition sites (and reverse "
            "complements) must not appear, comma-separated (e.g. EcoRI,BamHI).",
        )

        self.cpg_combo = QtWidgets.QComboBox()
        self.cpg_combo.addItems(["off", "deplete", "elevate"])
        self.cpg_combo.setAccessibleName("CpG dinucleotide mode")
        self._add_row(
            form, "CpG", self.cpg_combo,
            "CpG dinucleotide handling: 'deplete' for fewer CpGs (stealth / "
            "reduced innate-immune signal), 'elevate' for more (immunostimulatory). "
            "'off' adds no CpG objective.",
        )

        self.minmax_combo = QtWidgets.QComboBox()
        self.minmax_combo.addItems(["off", "max", "min"])
        self.minmax_combo.setAccessibleName("%MinMax codon-commonness mode")
        self._add_row(
            form, "%MinMax", self.minmax_combo,
            "Codon-commonness bias: 'max' favours common codons, 'min' favours "
            "rare codons. 'off' adds no %MinMax objective.",
        )

        self.tandem_spin = QtWidgets.QSpinBox()
        self.tandem_spin.setRange(0, 12)
        self.tandem_spin.setValue(0)
        self.tandem_spin.setSpecialValueText("off")
        self.tandem_spin.setAccessibleName("Tandem-repeat unit length (0 = off)")
        self._add_row(
            form, "Tandem unit", self.tandem_spin,
            "Ban a tandem repeat of this unit length -- a short motif repeated "
            "back-to-back (e.g. unit 3 bans CAGCAGCAG). 0 = off.",
        )

        self.inverted_spin = QtWidgets.QSpinBox()
        self.inverted_spin.setRange(0, 20)
        self.inverted_spin.setValue(0)
        self.inverted_spin.setSpecialValueText("off")
        self.inverted_spin.setAccessibleName("Inverted-repeat stem length (0 = off)")
        self._add_row(
            form, "Hairpin stem", self.inverted_spin,
            "Ban a hairpin (inverted repeat) with arms of this length -- a stem "
            "that folds back on itself and can occlude ribosome loading. 0 = off.",
        )

        self.internal_start_check = QtWidgets.QCheckBox("avoid strong-Kozak internal ATG")
        self.internal_start_check.setAccessibleName("Avoid internal start codons")
        self._add_row(
            form, "Internal ATG", self.internal_start_check,
            "Forbid internal ATG codons sitting in a strong Kozak context "
            "(purine at -3, G at +4), which could drive spurious re-initiation.",
        )

        self.tai_check = QtWidgets.QCheckBox("add tAI axis (real tRNA)")
        self.tai_check.setAccessibleName("Add tRNA-adaptation-index objective")
        self._add_row(
            form, "tAI", self.tai_check,
            "Add a tRNA-adaptation-index objective axis, built from real tRNA "
            "gene copy numbers; available only for organisms with bundled tRNA "
            "data.",
        )

        self.steps_spin = QtWidgets.QSpinBox()
        self.steps_spin.setRange(1, 25)
        self.steps_spin.setValue(9)
        self.steps_spin.setAccessibleName("Frontier steps")
        self._add_row(
            form, "Frontier steps", self.steps_spin,
            "Resolution of the objective trade-off (Pareto) sweep: more steps "
            "trace a finer frontier but take longer.",
        )

        self.beam_spin = QtWidgets.QSpinBox()
        self.beam_spin.setRange(0, 256)
        self.beam_spin.setValue(0)
        self.beam_spin.setSpecialValueText("exact")
        self.beam_spin.setAccessibleName("Beam width (0 = exact)")
        self._add_row(
            form, "Beam width", self.beam_spin,
            "0 = exact dynamic programming (proven-optimal). A positive value "
            "caps the trellis beam for speed; the result is then labeled "
            "beam-truncated.",
        )

        self.optimize_btn = QtWidgets.QPushButton("Optimize")
        self.optimize_btn.setAccessibleName("Run optimization")
        self.optimize_btn.setToolTip(
            "Run the optimization and compute the Pareto frontier on a background "
            "thread; the UI stays responsive."
        )
        self.optimize_btn.clicked.connect(self._start_optimize)
        form.addRow(self.optimize_btn)

        return box

    def _add_forbidden_presets(self, form: QtWidgets.QFormLayout) -> None:
        """Add a row of checkboxes, one per named forbidden-sequence preset.

        Each checkbox forbids its preset's motifs and shows the preset's
        description as its hover tooltip (mirroring BT3's forbidden-sequence
        options, made legible).
        """
        container = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout(container)
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setSpacing(2)
        self.preset_checks: dict[str, QtWidgets.QCheckBox] = {}
        for preset in api.available_forbidden_presets():
            check = QtWidgets.QCheckBox(preset.label)
            check.setToolTip(preset.description)
            check.setAccessibleName(f"Forbid {preset.label}")
            check.setAccessibleDescription(preset.description)
            self.preset_checks[preset.key] = check
            vbox.addWidget(check)
        self._add_row(
            form, "Forbidden presets", container,
            "Tick named groups of sequences to forbid (their reverse complements "
            "too). Hover each option for what it bans.",
        )

    def _add_row(
        self,
        form: QtWidgets.QFormLayout,
        text: str,
        widget: QtWidgets.QWidget,
        tooltip: str = "",
    ) -> None:
        """Add a labelled row whose label is the accessibility buddy of the widget.

        When ``tooltip`` is given it is set on both the label and the widget, so
        hovering either shows the explanatory bubble.
        """
        label = QtWidgets.QLabel(text)
        label.setBuddy(widget)
        if tooltip:
            label.setToolTip(tooltip)
            widget.setToolTip(tooltip)
            widget.setAccessibleDescription(tooltip)
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

        tracks_label = QtWidgets.QLabel("Per-site composition tracks (GC / CpG density)")
        layout.addWidget(tracks_label)
        self.tracks_plot = pg.PlotWidget()
        self.tracks_plot.setMaximumHeight(150)
        self.tracks_plot.setAccessibleName("Per-site composition tracks: GC and CpG density")
        tracks_label.setBuddy(self.tracks_plot)
        self._style_tracks_plot()
        layout.addWidget(self.tracks_plot)

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

    def _style_tracks_plot(self) -> None:
        """Apply theme-aware colours and axis labels to the per-site tracks plot."""
        bg = "#242832" if self._dark else "#ffffff"
        fg = "#e6e9ed" if self._dark else "#1b1f24"
        self.tracks_plot.setBackground(bg)
        self.tracks_plot.setLabel("bottom", "window start (nt)")
        self.tracks_plot.setLabel("left", "fraction / density")
        self.tracks_plot.showGrid(x=True, y=True, alpha=0.2)
        self.tracks_plot.addLegend(offset=(-10, 10))
        self.tracks_plot.setYRange(0.0, 1.0)
        for edge in ("bottom", "left"):
            axis = self.tracks_plot.getAxis(edge)
            axis.setPen(fg)
            axis.setTextPen(fg)
            axis.enableAutoSIPrefix(False)

    def _set_tab_order(self) -> None:
        """Wire a sensible keyboard tab order through the controls."""
        order = (
            self.protein_edit,
            self.organism_combo,
            self.jobname_edit,
            self.gc_spin,
            self.homo_spin,
            self.gc_run_spin,
            self.repeat_spin,
            self.motifs_edit,
            *self.preset_checks.values(),
            self.enzymes_edit,
            self.cpg_combo,
            self.minmax_combo,
            self.tandem_spin,
            self.inverted_spin,
            self.internal_start_check,
            self.tai_check,
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
        gc_run = self.gc_run_spin.value()
        repeat = self.repeat_spin.value()
        beam = self.beam_spin.value()
        motifs = tuple(
            m.strip().upper() for m in self.motifs_edit.text().split(",") if m.strip()
        )
        presets = tuple(key for key, check in self.preset_checks.items() if check.isChecked())
        enzymes = tuple(
            e.strip() for e in self.enzymes_edit.text().split(",") if e.strip()
        )
        cpg = self.cpg_combo.currentText()
        cpg_weight = 0.0 if cpg == "off" else 1.0
        cpg_mode = "deplete" if cpg == "off" else cpg
        minmax = self.minmax_combo.currentText()
        minmax_weight = 0.0 if minmax == "off" else 1.0
        minmax_direction = "max" if minmax == "off" else minmax
        tandem = self.tandem_spin.value()
        inverted = self.inverted_spin.value()
        config = api.OptimizeConfig(
            organism=self.organism_combo.currentText(),
            gc_target=self.gc_spin.value(),
            max_homopolymer=homo if homo > 0 else None,
            max_gc_run=gc_run if gc_run > 0 else None,
            max_repeat_length=repeat if repeat > 0 else None,
            forbidden_motifs=motifs,
            forbidden_presets=presets,
            avoid_reverse_complement=True,
            restriction_enzymes=enzymes,
            cpg_weight=cpg_weight,
            cpg_mode=cpg_mode,
            minmax_weight=minmax_weight,
            minmax_direction=minmax_direction,
            tandem_unit=tandem if tandem > 0 else None,
            inverted_stem=inverted if inverted > 0 else None,
            avoid_internal_start=self.internal_start_check.isChecked(),
            tai_weight=1.0 if self.tai_check.isChecked() else 0.0,
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
        self.tracks_plot.clear()
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
        self._render_tracks(delivered)

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

    def _render_tracks(self, delivered: api.Result) -> None:
        """Plot the delivered sequence's per-site GC and CpG-density tracks.

        Honest reporting profiles (``api.tracks``): each point is a sliding-window
        statistic recomputed from the DNA, never a solver output. Both tracks are
        fraction-scaled (0-1), so they share the y-axis; %MinMax (a different
        scale) stays available via ``api.tracks`` / ``bt4 tracks``.
        """
        self.tracks_plot.clear()
        organism = self.organism_combo.currentText()
        try:
            tracks = api.tracks(delivered.dna, organism, nt_window=30)
        except ValueError:
            return
        gc = tracks.get("gc_fraction")
        cpg = tracks.get("cpg_density")
        if gc is not None and gc.values:
            self.tracks_plot.plot(
                list(range(len(gc.values))),
                list(gc.values),
                pen=pg.mkPen("#4a83e0", width=2),
                name="GC fraction",
            )
        if cpg is not None and cpg.values:
            self.tracks_plot.plot(
                list(range(len(cpg.values))),
                list(cpg.values),
                pen=pg.mkPen("#e0724a", width=2),
                name="CpG density",
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
