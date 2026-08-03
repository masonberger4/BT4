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

# Above this residue count, warn before starting: the exact frontier sweep can
# take a noticeable amount of time (the run is cancelable either way).
_LONG_PROTEIN_WARN = 500


def _is_dark() -> bool:
    """Guess whether the running app uses a dark theme from its window colour."""
    app = QtWidgets.QApplication.instance()
    if isinstance(app, QtWidgets.QApplication):
        window = app.palette().color(QtGui.QPalette.ColorRole.Window)
        return window.lightness() < 128
    return False


class SequenceViewer(QtWidgets.QPlainTextEdit):
    """Read-only monospaced DNA viewer that annotates constraint-violation spans.

    Each :class:`~bt4.api.Violation` on the delivered result is drawn as a
    coloured background over its ``[start, end)`` nucleotide span -- red for a
    HARD (feasibility) violation, amber for a SOFT (quality) one -- and a hover
    tooltip names the constraint, its severity, span, and detail. This is how
    residual GLOBAL violations (e.g. a dispersed max-repeat or uORF that
    refinement could not fully clear) become visible *where they occur*, not
    just as a count in the metrics table.

    The DNA text itself is never altered: highlights are Qt "extra selections"
    layered on top, so the exported sequence stays exactly the delivered one
    (invariant: what the app shows is reproducible from its stamp).
    """

    # Semi-transparent so the monospace bases stay legible over the band.
    _HARD_LIGHT = QtGui.QColor(192, 57, 43, 70)
    _HARD_DARK = QtGui.QColor(224, 90, 74, 95)
    _SOFT_LIGHT = QtGui.QColor(185, 119, 14, 70)
    _SOFT_DARK = QtGui.QColor(224, 150, 74, 95)

    def __init__(self, dark: bool) -> None:
        """Build an empty read-only viewer themed for ``dark`` or light mode."""
        super().__init__()
        self._dark = dark
        self._violations: tuple[api.Violation, ...] = ()
        self.setReadOnly(True)
        self.setMouseTracking(True)

    def set_sequence(
        self, dna: str, violations: tuple[api.Violation, ...] = ()
    ) -> None:
        """Show ``dna`` and highlight every in-range violation span.

        Spans outside ``[0, len(dna))`` are dropped defensively. HARD bands are
        painted last so they sit visually on top of any overlapping SOFT band.
        """
        self.setPlainText(dna)
        n = len(dna)
        kept = tuple(v for v in violations if 0 <= v.start < v.end <= n)
        self._violations = kept
        # SOFT first, HARD last => HARD backgrounds win where spans overlap.
        order = sorted(kept, key=lambda v: v.severity is api.Severity.HARD)
        selections: list[QtWidgets.QTextEdit.ExtraSelection] = []
        for violation in order:
            selection = QtWidgets.QTextEdit.ExtraSelection()
            fmt = QtGui.QTextCharFormat()
            fmt.setBackground(self._colour(violation.severity))
            selection.format = fmt
            cursor = self.textCursor()
            cursor.setPosition(violation.start)
            cursor.setPosition(violation.end, QtGui.QTextCursor.MoveMode.KeepAnchor)
            selection.cursor = cursor
            selections.append(selection)
        self.setExtraSelections(selections)

    def _colour(self, severity: api.Severity) -> QtGui.QColor:
        """The highlight colour for a severity, theme-aware."""
        if severity is api.Severity.HARD:
            return self._HARD_DARK if self._dark else self._HARD_LIGHT
        return self._SOFT_DARK if self._dark else self._SOFT_LIGHT

    def _violation_at(self, pos: int) -> api.Violation | None:
        """The narrowest violation covering nucleotide ``pos`` (HARD wins ties)."""
        hits = [v for v in self._violations if v.start <= pos < v.end]
        if not hits:
            return None
        return min(
            hits,
            key=lambda v: (v.severity is not api.Severity.HARD, v.end - v.start),
        )

    def event(self, event: QtCore.QEvent) -> bool:
        """Show a per-span tooltip naming the constraint under the cursor."""
        if event.type() == QtCore.QEvent.Type.ToolTip and isinstance(
            event, QtGui.QHelpEvent
        ):
            hit = self._violation_at(self.cursorForPosition(event.pos()).position())
            if hit is None:
                QtWidgets.QToolTip.hideText()
            else:
                severity = hit.severity.value.upper()
                text = f"{hit.constraint} [{severity}] - nt {hit.start}-{hit.end}"
                if hit.detail:
                    text += f"\n{hit.detail}"
                QtWidgets.QToolTip.showText(event.globalPos(), text, self)
            return True
        return super().event(event)


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
        self._cancel_requested = False

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

        self.cpb_check = QtWidgets.QCheckBox("codon-pair bias (needs reference CDS)")
        self.cpb_check.setAccessibleName("Codon-pair bias objective")
        self._add_row(
            form, "Codon-pair", self.cpb_check,
            "Add a codon-pair-bias (CPS) objective axis, preferring naturally "
            "over-represented adjacent codon pairs. There is no bundled default "
            "table (codon-pair scores are organism-specific), so you must point "
            "'Codon-pair CDS' at a reference coding-sequence FASTA.",
        )

        self.cpb_cds_edit = QtWidgets.QLineEdit()
        self.cpb_cds_edit.setPlaceholderText("path to a reference CDS FASTA")
        self.cpb_cds_edit.setAccessibleName("Codon-pair reference CDS FASTA path")
        self._add_row(
            form, "Codon-pair CDS", self.cpb_cds_edit,
            "Path to a FASTA of reference coding sequences the codon-pair table is "
            "built from (required when 'codon-pair bias' is ticked). The table's "
            "content hash enters the run manifest.",
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

        self.uorf_check = QtWidgets.QCheckBox("avoid out-of-frame uORFs")
        self.uorf_check.setAccessibleName("Avoid out-of-frame uORFs")
        self._add_row(
            form, "uORFs", self.uorf_check,
            "Suppress out-of-frame internal ATGs that pair with a downstream "
            "in-frame stop -- short uORFs that divert ribosomes from the main ORF "
            "and lower yield. Non-local, so enforced by a refinement pass (result "
            "labeled heuristic; any it cannot remove are reported). A structural "
            "flag, NOT a calibrated expression prediction.",
        )

        self.uorf_region_spin = QtWidgets.QSpinBox()
        self.uorf_region_spin.setRange(3, 100000)
        self.uorf_region_spin.setValue(100)
        self.uorf_region_spin.setAccessibleName("uORF 5' scan window (nt)")
        self._add_row(
            form, "uORF window", self.uorf_region_spin,
            "5'-proximal window (nucleotides) scanned for uORF-opening ATGs when "
            "'avoid out-of-frame uORFs' is on. uORFs matter most near the start "
            "(leaky scanning); a wider window is stricter but harder to satisfy. "
            "Default 100 (~first 33 codons).",
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
        self.cancel_btn = QtWidgets.QPushButton("Cancel")
        self.cancel_btn.setAccessibleName("Cancel the running optimization")
        self.cancel_btn.clicked.connect(self._cancel_optimize)
        self.cancel_btn.setEnabled(False)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.optimize_btn)
        buttons.addWidget(self.cancel_btn)
        form.addRow(buttons)

        # tAI is only meaningful for organisms with a bundled tRNA table; keep the
        # checkbox enabled/disabled in step with the chosen organism.
        self.organism_combo.currentTextChanged.connect(
            lambda *_: self._update_tai_availability()
        )
        self._update_tai_availability()

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
        self.sequence_view = SequenceViewer(self._dark)
        self.sequence_view.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        self.sequence_view.setFont(mono)
        self.sequence_view.setAccessibleName("Delivered sequence")
        self.sequence_view.setMaximumHeight(140)
        seq_label.setBuddy(self.sequence_view)
        layout.addWidget(self.sequence_view)

        # Legend for the inline violation highlights; hidden until a delivered
        # sequence actually carries a violation, so a clean run stays uncluttered.
        self.violations_legend = QtWidgets.QLabel()
        self.violations_legend.setWordWrap(True)
        self.violations_legend.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.violations_legend.setAccessibleName("Violation legend")
        self.violations_legend.hide()
        layout.addWidget(self.violations_legend)

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
            self.cpb_check,
            self.cpb_cds_edit,
            self.tandem_spin,
            self.inverted_spin,
            self.internal_start_check,
            self.uorf_check,
            self.uorf_region_spin,
            self.tai_check,
            self.steps_spin,
            self.beam_spin,
            self.optimize_btn,
            self.cancel_btn,
            self.export_fasta_btn,
            self.export_json_btn,
        )
        for first, second in pairwise(order):
            self.setTabOrder(first, second)

    # ---- run --------------------------------------------------------------

    def _update_tai_availability(self) -> None:
        """Enable the tAI axis only for organisms that ship a tRNA table."""
        organism = self.organism_combo.currentText()
        available = organism in api.available_tai_organisms()
        self.tai_check.setEnabled(available)
        if available:
            self.tai_check.setToolTip("Add a tRNA-adaptation-index (tAI) objective axis.")
        else:
            self.tai_check.setChecked(False)
            self.tai_check.setToolTip(
                f"No tRNA data is bundled for {organism}, so tAI is unavailable "
                "for this organism."
            )

    def _prepare_protein(self) -> str | None:
        """Read, clean, and validate the protein box; show a friendly note on error.

        Returns the upper-cased amino-acid string, or ``None`` (after explaining
        the problem) if the box is empty, holds a pasted FASTA header, ends in a
        stop, or contains non-amino-acid characters.
        """
        raw = self.protein_edit.toPlainText()
        if any(line.lstrip().startswith(">") for line in raw.splitlines()):
            # They pasted a FASTA record -- take the first sequence for them.
            records = api.parse_fasta(raw)
            if not records or not records[0][1]:
                self._warn(
                    "That looks like a FASTA file",
                    "Paste just the amino-acid sequence, or remove the '>' header line.",
                )
                self.protein_edit.setFocus()
                return None
            protein = records[0][1]
        else:
            protein = "".join(raw.split()).upper()

        if not protein:
            self._warn(
                "No protein yet",
                "Paste a protein sequence (for example MAALKHETQW) to optimize.",
            )
            self.protein_edit.setFocus()
            return None
        if protein.endswith("*") and set(protein[:-1]) <= api.AMINO_ACIDS:
            self._warn(
                "Remove the trailing stop",
                "This sequence ends in '*'. BT4 adds the stop codon itself -- paste "
                "the protein without the trailing stop.",
            )
            self.protein_edit.setFocus()
            return None
        bad = sorted({ch for ch in protein if ch not in api.AMINO_ACIDS})
        if bad:
            shown = ", ".join(repr(ch) for ch in bad)
            self._warn(
                "That isn't a valid protein",
                f"These characters aren't amino acids: {shown}. Use single-letter "
                "amino-acid codes (A, C, D, E, F, ...).",
            )
            self.protein_edit.setFocus()
            return None
        return protein

    def _prepare_enzymes(self) -> tuple[str, ...] | None:
        """Canonicalize the restriction-enzyme field case-insensitively.

        Returns the tuple of catalog-cased enzyme names, or ``None`` (after
        listing the valid names) if any entry is not in the catalog.
        """
        raw = [e.strip() for e in self.enzymes_edit.text().split(",") if e.strip()]
        catalog = {name.lower(): name for name in api.available_enzymes()}
        canonical: list[str] = []
        unknown: list[str] = []
        for entry in raw:
            name = catalog.get(entry.lower())
            (canonical if name else unknown).append(name or entry)
        if unknown:
            self._warn(
                "Unknown restriction enzyme",
                f"Not in the catalog: {', '.join(unknown)}.",
                f"Available enzymes: {', '.join(api.available_enzymes())}.",
            )
            self.enzymes_edit.setFocus()
            return None
        return tuple(canonical)

    def _confirm_long_run(self, protein: str) -> bool:
        """Warn before optimizing a long protein; return whether to proceed."""
        if len(protein) <= _LONG_PROTEIN_WARN:
            return True
        answer = QtWidgets.QMessageBox.question(
            self,
            "Long protein",
            f"This protein has {len(protein)} residues, so the frontier sweep may "
            "take a while. You can press Cancel while it runs.\n\nStart anyway?",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )
        return answer == QtWidgets.QMessageBox.StandardButton.Yes

    def _build_config(self, enzymes: tuple[str, ...]) -> tuple[api.OptimizeConfig, int]:
        """Read the design controls into an OptimizeConfig and a frontier step count."""
        homo = self.homo_spin.value()
        gc_run = self.gc_run_spin.value()
        repeat = self.repeat_spin.value()
        beam = self.beam_spin.value()
        motifs = tuple(
            m.strip().upper() for m in self.motifs_edit.text().split(",") if m.strip()
        )
        presets = tuple(key for key, check in self.preset_checks.items() if check.isChecked())
        cpg = self.cpg_combo.currentText()
        cpg_weight = 0.0 if cpg == "off" else 1.0
        cpg_mode = "deplete" if cpg == "off" else cpg
        minmax = self.minmax_combo.currentText()
        minmax_weight = 0.0 if minmax == "off" else 1.0
        minmax_direction = "max" if minmax == "off" else minmax
        cpb_weight, cpb_cds = self._read_cpb()
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
            cpb_weight=cpb_weight,
            cpb_reference_cds=cpb_cds,
            minmax_weight=minmax_weight,
            minmax_direction=minmax_direction,
            tandem_unit=tandem if tandem > 0 else None,
            inverted_stem=inverted if inverted > 0 else None,
            avoid_internal_start=self.internal_start_check.isChecked(),
            avoid_uorf=self.uorf_check.isChecked(),
            uorf_region_nt=self.uorf_region_spin.value(),
            tai_weight=1.0 if self.tai_check.isChecked() else 0.0,
            beam=beam if beam > 0 else None,
            seed=0,
        )
        return config, self.steps_spin.value()

    def _read_cpb(self) -> tuple[float, tuple[str, ...]]:
        """Read the codon-pair controls into a (weight, reference-CDS) pair.

        Returns ``(0.0, ())`` -- codon-pair bias off -- when the box is unticked,
        or when it is ticked but the reference CDS FASTA is missing or unreadable
        (a warning is recorded on ``self._cpb_warning`` for the status bar, since
        there is no bundled default table to fall back on).
        """
        self._cpb_warning: str | None = None
        if not self.cpb_check.isChecked():
            return 0.0, ()
        path = self.cpb_cds_edit.text().strip()
        if not path:
            self._cpb_warning = "Codon-pair bias needs a reference CDS FASTA; ignoring it."
            return 0.0, ()
        try:
            cds = tuple(seq for _header, seq in api.read_fasta(path))
        except (OSError, ValueError) as exc:
            self._cpb_warning = f"Could not read codon-pair CDS ({exc}); ignoring it."
            return 0.0, ()
        if not cds:
            self._cpb_warning = "Codon-pair CDS FASTA had no sequences; ignoring it."
            return 0.0, ()
        return 1.0, cds

    def _start_optimize(self) -> None:
        """Validate the inputs, then launch a frontier optimization off-thread."""
        protein = self._prepare_protein()
        if protein is None:
            return
        enzymes = self._prepare_enzymes()
        if enzymes is None:
            return
        if not self._confirm_long_run(protein):
            return
        config, steps = self._build_config(enzymes)

        self._cancel_requested = False
        self._set_running(True)
        # _build_config -> _read_cpb may set a non-fatal codon-pair warning; fold it
        # into the running status so it is actually seen (a bare status is replaced).
        message = "Optimizing..."
        if self._cpb_warning:
            message = f"Optimizing... ({self._cpb_warning})"
        self.statusBar().showMessage(message)

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

    def _cancel_optimize(self) -> None:
        """Ask the running worker to stop after its current frontier point."""
        if self._worker is not None:
            self._cancel_requested = True
            self._worker.cancel()
            self.cancel_btn.setEnabled(False)
            self.statusBar().showMessage("Cancelling...")

    def _warn(self, title: str, text: str, detail: str = "") -> None:
        """Show a non-blocking, plain-language warning and mirror it to the status bar."""
        self.statusBar().showMessage(text)
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
        box.setWindowTitle(title)
        box.setText(text)
        if detail:
            box.setInformativeText(detail)
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        # show() (not exec()) so the app never blocks -- important headless too.
        self._msgbox = box
        box.show()

    def _set_running(self, running: bool) -> None:
        """Toggle controls while a run is in flight so the UI stays consistent."""
        self.optimize_btn.setEnabled(not running)
        self.optimize_btn.setText("Optimizing..." if running else "Optimize")
        self.cancel_btn.setEnabled(running)

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
        if self._cancel_requested:
            self.statusBar().showMessage(
                "Cancelled -- showing the partial frontier computed so far."
            )

    @QtCore.Slot(object)
    def _on_failed(self, error: object) -> None:
        """Clear the results and show a plain-language message on any engine error."""
        self._set_running(False)
        # Never leave a stale, still-exportable result behind a failed run (P9):
        # a subsequent export would silently write the previous sequence.
        self._reset_results()
        if self._cancel_requested:
            self.statusBar().showMessage("Optimization cancelled.")
            return
        headline, detail = self._friendly_error(error)
        self.statusBar().showMessage(headline)
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
        box.setWindowTitle("Couldn't optimize")
        box.setText(headline)
        if detail:
            box.setInformativeText(detail)
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        # show() (not exec()) so the app never blocks -- important headless too.
        self._msgbox = box
        box.show()

    def _friendly_error(self, error: object) -> tuple[str, str]:
        """Translate an engine exception into (headline, detail) plain language."""
        if isinstance(error, api.InfeasibleError):
            names = (
                ", ".join(error.constraints) if error.constraints else "the active limits"
            )
            return (
                "No sequence can satisfy these settings.",
                f"Nothing works under {names}. Try increasing Max homopolymer, removing "
                "forbidden motifs or restriction sites, or relaxing the GC target and the "
                "tandem/hairpin limits.",
            )
        return ("The optimization couldn't complete.", str(error))

    # ---- rendering --------------------------------------------------------

    def _reset_results(self) -> None:
        """Clear the results panel to an honest empty state."""
        # Drop the delivered result too, so nothing stale stays exportable via
        # _delivered() after a cleared/failed run (P9).
        self._last = None
        self.badge.setText("No optimization run yet.")
        self.badge.setStyleSheet(theme.badge_qss(""))
        self.metrics_table.setRowCount(0)
        self.sequence_view.set_sequence("")
        self.violations_legend.hide()
        self.violations_legend.clear()
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
        self._render_sequence(delivered)
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

    def _render_sequence(self, delivered: api.Result) -> None:
        """Show the delivered DNA with inline violation highlights and a legend.

        The highlighted spans come straight from ``delivered.violations`` (the
        whole-sequence audit), so a residual GLOBAL violation left after
        refinement is marked exactly where it sits. The legend summarises the
        HARD/SOFT counts and their colours, and stays hidden on a clean run.
        """
        violations = delivered.violations
        self.sequence_view.set_sequence(delivered.dna, violations)
        if not violations:
            self.violations_legend.hide()
            self.violations_legend.clear()
            return

        hard = sum(1 for v in violations if v.severity is api.Severity.HARD)
        soft = len(violations) - hard
        parts: list[str] = []
        if hard:
            parts.append(
                f'<span style="color:#c0392b;">&#9632;</span> '
                f"{hard} hard (feasibility)"
            )
        if soft:
            parts.append(
                f'<span style="color:#b9770e;">&#9632;</span> '
                f"{soft} soft (quality)"
            )
        self.violations_legend.setText(
            "Highlighted spans — " + " &nbsp; ".join(parts) + " · hover for detail"
        )
        self.violations_legend.show()

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
    if "--self-test" in sys.argv:
        # Headless startup check used by the packaging CI. Constructing
        # StudioWindow already exercised the things that break a frozen bundle --
        # bundled-resource loading (available_organisms), the Qt platform plugin,
        # and pyqtgraph -- so flush pending events and exit 0 without entering the
        # blocking event loop. Lets CI assert the packaged app actually opens.
        app.processEvents()
        return 0
    return app.exec()
