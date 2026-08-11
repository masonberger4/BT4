"""BT4 Studio main window and application entry point.

``StudioWindow`` is the whole desktop UI: a left control panel (protein, target
organism, and per-run design knobs) and a right results area of three tabs --

* **Design** -- the delivered frontier result: an honest optimality-certificate
  badge, a recomputed-metrics table, an interactive CAI/GC frontier scatter, a
  monospaced sequence viewer with inline violation annotations, per-site
  composition tracks, FASTA/JSON export, and the opt-in **ASSP cross-check**;
* **Candidates & splice audit** -- design-flow steps 3-4, the expression-annotated
  finalist set with its advisory cryptic-splice audit;
* **Library (sampled)** -- Phase-5 degenerate-design mode, a stochastic sampler
  that is honestly labeled as such.

Every engine call runs on a background worker (:mod:`bt4.app.worker`) moved onto
a ``QThread`` so the window never blocks, and any engine error is surfaced as a
non-fatal message rather than a crash.

This module talks to the engine exclusively through :mod:`bt4.api`; it never
imports the optimizer, pipeline, or biomodels directly. Everything is computed
locally and **nothing leaves the machine** unless the user explicitly presses
"Validate with ASSP", whose network-derived numbers are labeled advisory and are
never folded into an export or a run manifest (CLAUDE.md §6, §6.6, §10.15).
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from itertools import pairwise

import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from bt4 import api
from bt4.app import theme
from bt4.app.worker import (
    CandidatesResult,
    CandidatesWorker,
    CrossCheckWorker,
    LibraryWorker,
    OptimizeWorker,
)

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

_LIBRARY_COLS = ("#", "CAI", "GC %", "Hard viol.", "Soft viol.", "Certificate")

_CROSSCHECK_COLS = ("Kind", "Position", "Score", "Class")

# Shown before the one control in BT4 Studio that leaves the machine. The app is
# offline-first, so sending a designed sequence to a third-party web service is an
# explicit, informed choice -- never a side effect of pressing Optimize.
_ASSP_CONSENT = (
    "This sends the delivered coding sequence to the public ASSP web service "
    "(Alternative Splice Site Predictor) over the network.\n\n"
    "Everything else in BT4 Studio is computed locally. ASSP's numbers are "
    "network-derived and uncalibrated: they are advisory only, are NOT part of the "
    "run manifest, and are never written into an export.\n\n"
    "Send the sequence?"
)


def _distinct_site_count(positions: list[int], match_window: int) -> int:
    """Count distinct splice sites, merging positions within ``match_window``.

    ``positions`` must be sorted ascending. Two flags whose anchors are no more
    than ``match_window`` nt apart are treated as the same localized site (the
    same co-occurrence rule the audit uses for cross-backend matching), so a site
    flagged by several backends is counted once rather than once per backend.
    """
    sites = 0
    last: int | None = None
    for pos in positions:
        if last is None or pos - last > match_window:
            sites += 1
        last = pos
    return sites


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

    def set_dark(self, dark: bool) -> None:
        """Switch the highlight palette and repaint the current annotations.

        The DNA and its violation spans are already held, so a theme change
        re-derives the bands rather than needing the caller to re-render.
        """
        if dark == self._dark:
            return
        self._dark = dark
        self.set_sequence(self.toPlainText(), self._violations)

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
        self._theme_choice = "system"
        self._last: api.FrontierResult | None = None
        self._msgbox: QtWidgets.QMessageBox | None = None
        self._cancel_requested = False
        # Live thread/worker references, kept only so neither is garbage-collected
        # mid-run. Button enablement is driven by the _*_running flags below, never
        # by these, so a missed clear can never strand a control (regression: the
        # optimize-then-rank stuck-Optimize bug).
        self._thread: QtCore.QThread | None = None
        self._worker: OptimizeWorker | None = None
        self._cand_thread: QtCore.QThread | None = None
        self._cand_worker: CandidatesWorker | None = None
        self._cc_thread: QtCore.QThread | None = None
        self._cc_worker: CrossCheckWorker | None = None
        self._lib_thread: QtCore.QThread | None = None
        self._lib_worker: LibraryWorker | None = None
        self._optimize_running = False
        self._cand_running = False
        self._cc_running = False
        self._lib_running = False
        self._library: api.LibraryResult | None = None

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        controls_scroll = QtWidgets.QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setWidget(self._build_controls())
        controls_scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        controls_scroll.setMinimumWidth(300)
        splitter.addWidget(controls_scroll)
        splitter.addWidget(self._build_results())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        # Let the results side absorb every resize: the control column is a fixed
        # form, so a narrow window shrinks the plots/tables rather than clipping
        # the knobs (which stay reachable through their scroll area).
        splitter.setCollapsible(0, False)
        splitter.setSizes([380, 800])
        self.setCentralWidget(splitter)

        self._build_menus()
        self.setStatusBar(QtWidgets.QStatusBar())
        self.statusBar().showMessage("Ready. Enter a protein and click Optimize.")

        self._set_tab_order()
        self._reset_results()
        self._reset_candidates()
        self._reset_library()

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

        self.splice_check = QtWidgets.QCheckBox("avoid strong splice-consensus motifs")
        self.splice_check.setAccessibleName("Avoid strong splice-site consensus motifs")
        self._add_row(
            form, "Splice sites", self.splice_check,
            "Forbid strong donor (GTRAGT) and acceptor (polypyrimidine + YAG) "
            "splice-consensus motifs on the sense strand, to suppress the most "
            "obvious cryptic splice sites. A structural heuristic (not a splice "
            "CNN); it never bans the bare GT/AG and makes no calibrated risk claim.",
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
        """Build the right results area: Design, Candidates, and Library tabs.

        The Design tab shows the delivered frontier result (badge, metrics,
        frontier scatter, sequence, tracks) plus the opt-in ASSP cross-check. The
        Candidates tab runs the expression/splice design flow
        (:func:`bt4.api.candidates` + :func:`bt4.api.splice_audit`) and shows the
        ranked, honestly-labeled candidate set with its advisory splice audit. The
        Library tab draws a sampled (not optimized) library via
        :func:`bt4.api.library`.
        """
        tabs = QtWidgets.QTabWidget()
        tabs.setAccessibleName("Results")
        tabs.addTab(self._build_design_tab(), "Design")
        tabs.addTab(self._build_candidates_tab(), "Candidates && splice audit")
        tabs.addTab(self._build_library_tab(), "Library (sampled)")
        self.tabs = tabs
        return tabs

    def _build_design_tab(self) -> QtWidgets.QWidget:
        """Build the delivered-frontier results tab."""
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

        layout.addWidget(self._build_crosscheck_group())

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

    def _build_crosscheck_group(self) -> QtWidgets.QWidget:
        """Build the opt-in ASSP cross-check control and its advisory readout.

        This is the single control in BT4 Studio that touches the network
        (CLAUDE.md §6.6), so it is explicit in every direction: the user presses
        it, consents to the sequence being sent, and gets a result banner that
        leads with **network-derived / uncalibrated / not in the manifest**. The
        report is never folded into the delivered :class:`~bt4.api.Result`, the
        metrics table, or an export -- exactly as the CLI prints it to stderr
        rather than into the stdout artifact (§10.15).
        """
        box = QtWidgets.QGroupBox("Cross-check (optional, network)")
        layout = QtWidgets.QVBoxLayout(box)

        row = QtWidgets.QHBoxLayout()
        self.assp_btn = QtWidgets.QPushButton("Validate with ASSP")
        self.assp_btn.setAccessibleName("Validate the delivered sequence with ASSP")
        self.assp_btn.setToolTip(
            "Send the delivered sequence to the online ASSP service (Alternative "
            "Splice Site Predictor) for an independent cryptic-splice opinion. "
            "Opt-in and out-of-loop: it never steers the optimizer. Its numbers are "
            "network-derived and uncalibrated -- advisory only, excluded from the "
            "run manifest, and never written into an export. Needs the bt4[assp] "
            "extra; an outage degrades gracefully and never fails a run."
        )
        self.assp_btn.clicked.connect(self._start_crosscheck)
        self.assp_btn.setEnabled(False)
        row.addWidget(self.assp_btn)
        row.addStretch(1)
        layout.addLayout(row)

        self.assp_banner = QtWidgets.QLabel()
        self.assp_banner.setWordWrap(True)
        self.assp_banner.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.assp_banner.setAccessibleName("ASSP cross-check summary")
        self.assp_banner.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        layout.addWidget(self.assp_banner)

        self.assp_table = QtWidgets.QTableWidget(0, len(_CROSSCHECK_COLS))
        self.assp_table.setHorizontalHeaderLabels(list(_CROSSCHECK_COLS))
        self.assp_table.verticalHeader().setVisible(False)
        self.assp_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.assp_table.horizontalHeader().setStretchLastSection(True)
        self.assp_table.setAccessibleName("ASSP predicted splice sites")
        self.assp_table.setMaximumHeight(150)
        self.assp_table.hide()
        layout.addWidget(self.assp_table)

        return box

    _CANDIDATE_COLS = (
        "#",
        "Source",
        "CAI",
        "GC %",
        "Expression",
        "Calibrated",
        "Hard viol.",
        "Splice flags",
    )

    def _build_candidates_tab(self) -> QtWidgets.QWidget:
        """Build the candidate-set + splice-audit tab (design-flow steps 3-4).

        A small control bar (candidate count, repeat-refined variants, an
        opt-in to run the installed splice CNNs) drives a background
        :class:`~bt4.app.worker.CandidatesWorker`; the result is an honestly
        labeled table plus advisory splice-audit banners. Every score here is
        uncalibrated today, and the UI says so rather than implying a ranking.
        """
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)

        intro = QtWidgets.QLabel(
            "Assemble the Pareto frontier (plus repeat-refined variants) into a "
            "finalist set, annotate each with the expression head, and run an "
            "advisory cryptic-splice audit. Uses the same design controls on the "
            "left; press <b>Rank &amp; audit</b> after (or instead of) Optimize."
        )
        intro.setWordWrap(True)
        intro.setTextFormat(QtCore.Qt.TextFormat.RichText)
        layout.addWidget(intro)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(QtWidgets.QLabel("Keep"))
        self.cand_n_spin = QtWidgets.QSpinBox()
        self.cand_n_spin.setRange(1, 100)
        self.cand_n_spin.setValue(24)
        self.cand_n_spin.setAccessibleName("Maximum candidates to keep")
        self.cand_n_spin.setToolTip(
            "Maximum candidates to keep after scoring. The delivered sequence is "
            "always retained; the cap is applied after scoring so it never drops it."
        )
        controls.addWidget(self.cand_n_spin)

        controls.addWidget(QtWidgets.QLabel("Repeat variants"))
        self.cand_repeat_spin = QtWidgets.QSpinBox()
        self.cand_repeat_spin.setRange(0, 20)
        self.cand_repeat_spin.setValue(4)
        self.cand_repeat_spin.setAccessibleName("Repeat-refined variants to attempt")
        self.cand_repeat_spin.setToolTip(
            "Repeat-refined variants to attempt when the delivered exact-DP seed "
            "violates a GLOBAL rule (max-repeat / uORF). 0 = frontier only."
        )
        controls.addWidget(self.cand_repeat_spin)

        self.splice_cnn_check = QtWidgets.QCheckBox("run installed splice CNNs")
        self.splice_cnn_check.setAccessibleName("Include installed splice CNN backends")
        self.splice_cnn_check.setToolTip(
            "Also run the wrapped SpliceAI / Pangolin CNNs in the splice audit "
            "when they are installed (out-of-loop, may be slow). Off = the honest "
            "PWM baseline only. Every backend is uncalibrated today, so the audit "
            "stays advisory either way."
        )
        controls.addWidget(self.splice_cnn_check)
        controls.addStretch(1)
        layout.addLayout(controls)

        layout.addWidget(self._build_expression_group())

        buttons = QtWidgets.QHBoxLayout()
        self.rank_btn = QtWidgets.QPushButton("Rank && audit")
        self.rank_btn.setAccessibleName("Rank and audit candidates")
        self.rank_btn.setToolTip(
            "Assemble, expression-rank, and splice-audit the candidate set on a "
            "background thread; the UI stays responsive."
        )
        self.rank_btn.clicked.connect(self._start_candidates)
        # No Cancel button here: unlike the point-wise-cancelable frontier sweep,
        # the assemble->audit flow exposes no cancellation hook, so a Cancel
        # control would be a dishonest no-op. The disabled "Ranking..." button is
        # the honest in-flight indicator instead.
        buttons.addWidget(self.rank_btn)
        buttons.addStretch(1)
        layout.addLayout(buttons)

        self.cand_banner = QtWidgets.QLabel()
        self.cand_banner.setObjectName("certBadge")
        self.cand_banner.setWordWrap(True)
        self.cand_banner.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.cand_banner.setAccessibleName("Candidate-set summary")
        self.cand_banner.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        layout.addWidget(self.cand_banner)

        self.candidates_table = QtWidgets.QTableWidget(0, len(self._CANDIDATE_COLS))
        self.candidates_table.setHorizontalHeaderLabels(list(self._CANDIDATE_COLS))
        self.candidates_table.verticalHeader().setVisible(False)
        self.candidates_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.candidates_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.candidates_table.horizontalHeader().setStretchLastSection(True)
        self.candidates_table.setAccessibleName("Ranked candidate set")
        layout.addWidget(self.candidates_table, stretch=1)

        self.splice_banner = QtWidgets.QLabel()
        self.splice_banner.setWordWrap(True)
        self.splice_banner.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.splice_banner.setAccessibleName("Splice-audit summary")
        self.splice_banner.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        layout.addWidget(self.splice_banner)

        return panel

    def _build_expression_group(self) -> QtWidgets.QWidget:
        """Build the opt-in expression-head controls for the candidate flow.

        BT4 ships one wrapped head -- RiboNN (Sanofi non-commercial), driven from
        the user's own checkout via ``$BT4_RIBONN_DIR``. The toggle is enabled
        only when :func:`bt4.api.available_expression_backends` reports it can
        actually run here, so the control is never a dead end, and the fixed
        5'/3' UTR context it needs is entered beside it: RiboNN refuses an empty
        UTR (its loader reads an all-empty column as NaN) and the UTRs carry most
        of its signal, so holding the real ones fixed while the CDS varies is the
        intended use.

        Selecting RiboNN changes **nothing** about delivery. It is
        ``calibrated=False`` -- reproducing it faithfully is not calibration for
        BT4's CDS-variant regime -- so :func:`bt4.api.candidates` keeps the set in
        discovery order and leaves the solver's pick delivered; the head only
        annotates (CLAUDE.md §6, §10.6).
        """
        box = QtWidgets.QGroupBox("Expression head (optional)")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)

        available = api.available_expression_backends()
        self._ribonn_available = "ribonn" in available

        self.ribonn_check = QtWidgets.QCheckBox("score candidates with RiboNN")
        self.ribonn_check.setAccessibleName("Use the wrapped RiboNN expression head")
        self.ribonn_check.setEnabled(self._ribonn_available)
        self.ribonn_check.toggled.connect(self._update_ribonn_enabled)
        if self._ribonn_available:
            ribonn_tip = (
                "Annotate every candidate with the wrapped RiboNN translation-"
                "efficiency model (your own checkout and weights, hash-verified "
                "before loading). Runs once over the whole set, off the GUI thread. "
                "RiboNN is UNCALIBRATED for BT4's CDS-variant regime, so it only "
                "annotates: the table stays in discovery order and the solver's "
                "pick stays delivered. Requires the 5'/3' UTR context below."
            )
        else:
            ribonn_tip = (
                "No usable RiboNN install found. Clone the RiboNN repository "
                "(Sanofi non-commercial), download its Zenodo weights into "
                "<repo>/models, install torch + pandas, and point $BT4_RIBONN_DIR "
                "at the checkout. Nothing is bundled: BT4 drives your own copy."
            )
        self._add_row(form, "RiboNN", self.ribonn_check, ribonn_tip)

        self.ribonn_species_combo = QtWidgets.QComboBox()
        self.ribonn_species_combo.addItems(["human", "mouse"])
        self.ribonn_species_combo.setAccessibleName("RiboNN species weight set")
        self._add_row(
            form, "Species", self.ribonn_species_combo,
            "Which RiboNN weight set to ensemble -- the human or mouse "
            "cross-validation runs. Match it to your target organism.",
        )

        self.utr5_edit = QtWidgets.QLineEdit()
        self.utr5_edit.setPlaceholderText("5' UTR (required for RiboNN)")
        self.utr5_edit.setAccessibleName("Fixed 5-prime UTR context")
        self._add_row(
            form, "5' UTR", self.utr5_edit,
            "The transcript's 5' UTR, held FIXED while the CDS varies. RiboNN "
            "scores a whole transcript and cannot accept an empty UTR; the 5' UTR "
            "and initiation region carry most of its signal, so an empty-UTR score "
            "would not be meaningful. Not part of the designed sequence -- it is "
            "context for the model only and is never exported.",
        )

        self.utr3_edit = QtWidgets.QLineEdit()
        self.utr3_edit.setPlaceholderText("3' UTR (required for RiboNN)")
        self.utr3_edit.setAccessibleName("Fixed 3-prime UTR context")
        self._add_row(
            form, "3' UTR", self.utr3_edit,
            "The transcript's 3' UTR, held FIXED while the CDS varies (as the 5' "
            "UTR). Context for the model only; never part of the exported CDS.",
        )

        self._update_ribonn_enabled(self.ribonn_check.isChecked())
        return box

    def _update_ribonn_enabled(self, checked: bool) -> None:
        """Enable the RiboNN sub-controls only while the head is selected."""
        for widget in (
            self.ribonn_species_combo,
            self.utr5_edit,
            self.utr3_edit,
        ):
            widget.setEnabled(checked and self._ribonn_available)

    def _prepare_predictor(self) -> tuple[bool, api.ExpressionPredictor | None]:
        """Build the expression head for a candidate run, or explain the problem.

        Returns:
            ``(ok, predictor)``. ``predictor`` is what to pass to
            :func:`bt4.api.candidates` -- ``None`` meaning the default neutral
            placeholder. ``ok`` is ``False`` when the user asked for RiboNN but
            has not supplied the UTR context it requires; a warning has already
            been shown and the run must not start. Refusing here, rather than
            letting the engine raise mid-run, keeps the failure legible.
        """
        if not (self.ribonn_check.isChecked() and self._ribonn_available):
            return True, None
        utr5 = "".join(self.utr5_edit.text().split()).upper()
        utr3 = "".join(self.utr3_edit.text().split()).upper()
        missing = [
            name for name, value in (("5'", utr5), ("3'", utr3)) if not value
        ]
        if missing:
            self._warn(
                "RiboNN needs UTR context",
                f"Enter the transcript's {' and '.join(missing)} UTR, or untick "
                "RiboNN to use the neutral placeholder.",
                "RiboNN scores a whole transcript: an empty UTR column breaks its "
                "preprocessing, and the UTRs carry most of its signal. They are "
                "held fixed while the CDS varies and are never exported.",
            )
            (self.utr5_edit if not utr5 else self.utr3_edit).setFocus()
            return False, None
        bad = sorted({ch for ch in utr5 + utr3 if ch not in "ACGT"})
        if bad:
            shown = ", ".join(repr(ch) for ch in bad)
            self._warn(
                "That isn't a UTR sequence",
                f"These characters aren't DNA bases: {shown}. Use A, C, G, and T.",
            )
            self.utr5_edit.setFocus()
            return False, None
        predictor = api.resolve_expression_backend(
            "ribonn",
            species=self.ribonn_species_combo.currentText(),
            utr5=utr5,
            utr3=utr3,
        )
        return True, predictor

    def _build_library_tab(self) -> QtWidgets.QWidget:
        """Build the Phase-5 library / degenerate-design tab.

        Library mode is a **sampler, not an optimizer** (CLAUDE.md §9, Phase 5):
        it draws each residue's synonymous codon from the organism's usage
        distribution (tempered), keeping only codons that pass every LOCAL
        constraint. Every member carries the ``SAMPLED`` certificate and makes no
        optimality or expression claim, GLOBAL rules are validated but not
        enforced during the draw, and the banner says all of that plainly rather
        than letting a table of sequences imply a ranking.
        """
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)

        intro = QtWidgets.QLabel(
            "Draw a diverse library instead of one optimum: each residue's "
            "synonymous codon is <b>sampled</b> from the organism's usage "
            "distribution, keeping only codons that satisfy the LOCAL constraints "
            "on the left. Useful for screening panels and degenerate designs."
        )
        intro.setWordWrap(True)
        intro.setTextFormat(QtCore.Qt.TextFormat.RichText)
        layout.addWidget(intro)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(QtWidgets.QLabel("Members"))
        self.lib_n_spin = QtWidgets.QSpinBox()
        self.lib_n_spin.setRange(1, 200)
        self.lib_n_spin.setValue(12)
        self.lib_n_spin.setAccessibleName("Library members to sample")
        self.lib_n_spin.setToolTip("How many sequences to draw.")
        controls.addWidget(self.lib_n_spin)

        controls.addWidget(QtWidgets.QLabel("Temperature"))
        self.lib_temp_spin = QtWidgets.QDoubleSpinBox()
        self.lib_temp_spin.setRange(0.05, 10.0)
        self.lib_temp_spin.setSingleStep(0.1)
        self.lib_temp_spin.setValue(1.0)
        self.lib_temp_spin.setAccessibleName("Sampling temperature")
        self.lib_temp_spin.setToolTip(
            "Sampling temperature. Toward 0 the draw approaches the most-favored "
            "codon at every position (low diversity); 1.0 is the organism's "
            "natural usage distribution; larger values approach uniform (high "
            "diversity, lower CAI)."
        )
        controls.addWidget(self.lib_temp_spin)

        controls.addWidget(QtWidgets.QLabel("Seed"))
        self.lib_seed_spin = QtWidgets.QSpinBox()
        self.lib_seed_spin.setRange(0, 2_000_000_000)
        self.lib_seed_spin.setValue(0)
        self.lib_seed_spin.setAccessibleName("Sampling seed")
        self.lib_seed_spin.setToolTip(
            "Master seed for the draw. The same seed reproduces the same library "
            "byte-for-byte, and the effective seed enters the run manifest."
        )
        controls.addWidget(self.lib_seed_spin)
        controls.addStretch(1)
        layout.addLayout(controls)

        buttons = QtWidgets.QHBoxLayout()
        self.library_btn = QtWidgets.QPushButton("Sample library")
        self.library_btn.setAccessibleName("Sample a sequence library")
        self.library_btn.setToolTip(
            "Draw the library on a background thread; the UI stays responsive."
        )
        self.library_btn.clicked.connect(self._start_library)
        buttons.addWidget(self.library_btn)
        self.export_library_btn = QtWidgets.QPushButton("Export library FASTA")
        self.export_library_btn.setAccessibleName("Export the library as FASTA")
        self.export_library_btn.setToolTip(
            "Write every sampled member to one multi-record FASTA file."
        )
        self.export_library_btn.clicked.connect(self._export_library)
        buttons.addWidget(self.export_library_btn)
        buttons.addStretch(1)
        layout.addLayout(buttons)

        self.lib_banner = QtWidgets.QLabel()
        self.lib_banner.setObjectName("certBadge")
        self.lib_banner.setWordWrap(True)
        self.lib_banner.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.lib_banner.setAccessibleName("Library summary")
        self.lib_banner.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse
        )
        layout.addWidget(self.lib_banner)

        self.library_table = QtWidgets.QTableWidget(0, len(_LIBRARY_COLS))
        self.library_table.setHorizontalHeaderLabels(list(_LIBRARY_COLS))
        self.library_table.verticalHeader().setVisible(False)
        self.library_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.library_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.library_table.horizontalHeader().setStretchLastSection(True)
        self.library_table.setAccessibleName("Sampled library members")
        self.library_table.itemSelectionChanged.connect(self._show_library_member)
        layout.addWidget(self.library_table, stretch=1)

        member_label = QtWidgets.QLabel("Selected member")
        layout.addWidget(member_label)
        self.library_view = SequenceViewer(self._dark)
        self.library_view.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth
        )
        self.library_view.setFont(
            QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        )
        self.library_view.setAccessibleName("Selected library member sequence")
        self.library_view.setMaximumHeight(110)
        member_label.setBuddy(self.library_view)
        layout.addWidget(self.library_view)

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

    def _build_menus(self) -> None:
        """Build the menu bar: File, Run, View, and Help.

        The menu is the keyboard-first route to everything the buttons do (each
        action carries a standard shortcut), and it is where the theme choice
        lives -- so the app is usable without a mouse and legible in either
        theme (CLAUDE.md §6.6: accessibility is a requirement, not a nice-to-have).
        """
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        self._add_action(
            file_menu, "Export &FASTA...", self._export_fasta, "Ctrl+E"
        )
        self._add_action(file_menu, "Export &JSON...", self._export_json, "Ctrl+J")
        self._add_action(
            file_menu, "Export &library FASTA...", self._export_library, "Ctrl+Shift+E"
        )
        file_menu.addSeparator()
        self._add_action(file_menu, "&Quit", self.close, "Ctrl+Q")

        run_menu = menubar.addMenu("&Run")
        self._add_action(run_menu, "&Optimize", self._start_optimize, "Ctrl+R")
        self._add_action(run_menu, "&Cancel", self._cancel_optimize, "Esc")
        run_menu.addSeparator()
        self._add_action(
            run_menu, "Rank && &audit candidates", self._start_candidates, "Ctrl+K"
        )
        self._add_action(
            run_menu, "&Sample library", self._start_library, "Ctrl+L"
        )
        run_menu.addSeparator()
        self._add_action(
            run_menu, "&Validate with ASSP (online)", self._start_crosscheck
        )

        view_menu = menubar.addMenu("&View")
        self.theme_group = QtGui.QActionGroup(self)
        self.theme_group.setExclusive(True)
        for label, key in (("&System", "system"), ("&Light", "light"), ("&Dark", "dark")):
            action = self._add_action(
                view_menu, label, lambda _=False, k=key: self._apply_theme(k)
            )
            action.setCheckable(True)
            action.setChecked(key == self._theme_choice)
            self.theme_group.addAction(action)

        help_menu = menubar.addMenu("&Help")
        self._add_action(help_menu, "&About BT4 Studio", self._show_about)

    def _add_action(
        self,
        menu: QtWidgets.QMenu,
        text: str,
        slot: Callable[..., object],
        shortcut: str = "",
    ) -> QtGui.QAction:
        """Add one menu action (optionally with a shortcut) and return it."""
        action = QtGui.QAction(text, self)
        if shortcut:
            action.setShortcut(QtGui.QKeySequence(shortcut))
        action.triggered.connect(slot)
        menu.addAction(action)
        return action

    def _show_about(self) -> None:
        """Show what BT4 Studio is -- and what it does and does not claim."""
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Information)
        box.setWindowTitle("About BT4 Studio")
        box.setText(
            "BT4 Studio designs a coding sequence by constrained, "
            "multi-objective optimization over a codon trellis."
        )
        box.setInformativeText(
            "Everything is computed locally and stays on this machine, except the "
            "opt-in ASSP cross-check, which you trigger explicitly. Results carry "
            "an optimality certificate and a content-addressed run manifest, so "
            "any exported design is reproducible from its stamp. The splice and "
            "expression models shipped today are UNCALIBRATED: they annotate and "
            "advise, and they never steer what is delivered."
        )
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        # show() (not exec()) so the app never blocks -- important headless too.
        self._msgbox = box
        box.show()

    def _apply_theme(self, choice: str) -> None:
        """Switch the app between the system, light, and dark palettes at runtime.

        Restyles everything that caches a colour: the application stylesheet, the
        two pyqtgraph plots, the certificate/summary badges, and the sequence
        viewers' violation highlights (which are re-applied from the still-live
        results, so nothing has to be recomputed).
        """
        self._theme_choice = choice
        self._dark = _is_dark() if choice == "system" else (choice == "dark")

        app = QtWidgets.QApplication.instance()
        if isinstance(app, QtWidgets.QApplication):
            app.setStyleSheet(theme.stylesheet(self._dark))

        self._style_plot()
        self._style_tracks_plot()
        for viewer in (self.sequence_view, self.library_view):
            viewer.set_dark(self._dark)

        # Badges carry their own inline QSS, so re-assert it against the new theme.
        delivered = self._delivered()
        self.badge.setStyleSheet(
            theme.badge_qss(delivered.certificate.status.value if delivered else "")
        )
        if self._last is not None and delivered is not None:
            self._render_frontier(self._last)
            self._render_tracks(delivered)

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
            self.splice_check,
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
            self.assp_btn,
            self.cand_n_spin,
            self.cand_repeat_spin,
            self.splice_cnn_check,
            self.ribonn_check,
            self.ribonn_species_combo,
            self.utr5_edit,
            self.utr3_edit,
            self.rank_btn,
            self.lib_n_spin,
            self.lib_temp_spin,
            self.lib_seed_spin,
            self.library_btn,
            self.export_library_btn,
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
            avoid_splice_sites=self.splice_check.isChecked(),
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

    def _wire_thread(
        self,
        worker: QtCore.QObject,
        on_finished: Callable[[object], None],
        on_failed: Callable[[object], None],
        on_cleared: Callable[[], None],
    ) -> QtCore.QThread:
        """Move ``worker`` onto a fresh thread and wire the standard signal set.

        Every background flow needs the same eight connections (start, progress,
        finished/failed to the window and to ``quit``, deletion, and a cleared
        callback), so they live in one place rather than being re-typed per flow.
        The thread is returned **unstarted**: the caller stores its references
        first, so neither object can be garbage-collected mid-run.
        """
        thread = QtCore.QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._on_progress)
        worker.finished.connect(on_finished)
        worker.failed.connect(on_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(on_cleared)
        return thread

    def _start_optimize(self) -> None:
        """Validate the inputs, then launch a frontier optimization off-thread."""
        if self._thread is not None:
            return  # an optimize is already in flight
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

        worker = OptimizeWorker(protein, config, steps)
        thread = self._wire_thread(
            worker, self._on_finished, self._on_failed, self._clear_optimize_thread
        )
        # Keep references so neither is garbage-collected mid-run.
        self._thread = thread
        self._worker = worker
        thread.start()

    def _clear_optimize_thread(self) -> None:
        """Drop the finished optimize thread/worker refs once the thread has quit."""
        self._thread = None
        self._worker = None

    def _cancel_optimize(self) -> None:
        """Ask the running worker to stop after its current frontier point."""
        if self._worker is not None:
            self._cancel_requested = True
            self._worker.cancel()
            self.cancel_btn.setEnabled(False)
            self.statusBar().showMessage("Cancelling...")

    # ---- candidates flow --------------------------------------------------

    def _start_candidates(self) -> None:
        """Validate inputs, then assemble + rank + splice-audit off-thread."""
        if self._cand_thread is not None:
            return  # a candidate run is already in flight
        protein = self._prepare_protein()
        if protein is None:
            return
        enzymes = self._prepare_enzymes()
        if enzymes is None:
            return
        ok, predictor = self._prepare_predictor()
        if not ok:
            return
        if not self._confirm_long_run(protein):
            return
        config, steps = self._build_config(enzymes)

        self._set_candidates_running(True)
        message = "Ranking candidates..."
        if self._cpb_warning:
            message = f"Ranking candidates... ({self._cpb_warning})"
        self.statusBar().showMessage(message)

        worker = CandidatesWorker(
            protein,
            config,
            steps=steps,
            n=self.cand_n_spin.value(),
            repeat_variants=self.cand_repeat_spin.value(),
            include_cnns=self.splice_cnn_check.isChecked(),
            predictor=predictor,
        )
        thread = self._wire_thread(
            worker,
            self._on_cand_finished,
            self._on_cand_failed,
            self._clear_candidates_thread,
        )
        self._cand_thread = thread
        self._cand_worker = worker
        thread.start()

    def _clear_candidates_thread(self) -> None:
        """Drop the finished thread/worker references once the thread has quit."""
        self._cand_thread = None
        self._cand_worker = None

    def _set_candidates_running(self, running: bool) -> None:
        """Toggle controls while a candidate run is in flight."""
        self._cand_running = running
        self.rank_btn.setText("Ranking..." if running else "Rank && audit")
        self._update_run_buttons()

    # ---- cross-check flow -------------------------------------------------

    def _start_crosscheck(self) -> None:
        """Ask for consent, then cross-check the delivered sequence off-thread.

        The only outbound network call BT4 Studio makes. It runs on the delivered
        sequence only (never mid-optimization), is confirmed by the user first,
        and can never fail a run: an outage returns ``available is False`` with a
        reason (CLAUDE.md §10.15).
        """
        if self._cc_thread is not None:
            return  # a cross-check is already in flight
        delivered = self._delivered()
        if delivered is None:
            return
        answer = QtWidgets.QMessageBox.question(
            self,
            "Send sequence to ASSP?",
            _ASSP_CONSENT,
            QtWidgets.QMessageBox.StandardButton.Yes
            | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )
        if answer != QtWidgets.QMessageBox.StandardButton.Yes:
            return

        self._set_crosscheck_running(True)
        self.statusBar().showMessage("Cross-checking with ASSP...")

        worker = CrossCheckWorker(delivered.dna, backend="assp")
        thread = self._wire_thread(
            worker,
            self._on_crosscheck_finished,
            self._on_crosscheck_failed,
            self._clear_crosscheck_thread,
        )
        self._cc_thread = thread
        self._cc_worker = worker
        thread.start()

    def _clear_crosscheck_thread(self) -> None:
        """Drop the finished cross-check thread/worker references."""
        self._cc_thread = None
        self._cc_worker = None

    def _set_crosscheck_running(self, running: bool) -> None:
        """Toggle controls while a cross-check is in flight."""
        self._cc_running = running
        self.assp_btn.setText("Cross-checking..." if running else "Validate with ASSP")
        self._update_run_buttons()

    @QtCore.Slot(object)
    def _on_crosscheck_finished(self, report: object) -> None:
        """Render a finished cross-check (available or gracefully unavailable)."""
        self._set_crosscheck_running(False)
        if not isinstance(report, api.SpliceCrossCheck):
            return
        self._render_crosscheck(report)

    @QtCore.Slot(object)
    def _on_crosscheck_failed(self, error: object) -> None:
        """Report a cross-check that could not even be attempted.

        ``run_splice_crosscheck`` degrades a backend failure into an unavailable
        *report*, so reaching here means something more basic went wrong. Either
        way this is advisory: the delivered result is left untouched and still
        exportable -- an optional cross-check never invalidates a run.
        """
        self._set_crosscheck_running(False)
        self.assp_table.hide()
        self.assp_table.setRowCount(0)
        self.assp_banner.setText(
            "<b>ASSP cross-check unavailable</b> &mdash; "
            f"{error}<br><i>Advisory only; the delivered result is unaffected.</i>"
        )
        self.statusBar().showMessage("ASSP cross-check unavailable.")

    # ---- library flow -----------------------------------------------------

    def _start_library(self) -> None:
        """Validate inputs, then sample a library off-thread."""
        if self._lib_thread is not None:
            return  # a library draw is already in flight
        protein = self._prepare_protein()
        if protein is None:
            return
        enzymes = self._prepare_enzymes()
        if enzymes is None:
            return
        config, _steps = self._build_config(enzymes)

        self._set_library_running(True)
        self.statusBar().showMessage("Sampling library...")

        worker = LibraryWorker(
            protein,
            config,
            n=self.lib_n_spin.value(),
            temperature=self.lib_temp_spin.value(),
            seed=self.lib_seed_spin.value(),
        )
        thread = self._wire_thread(
            worker,
            self._on_library_finished,
            self._on_library_failed,
            self._clear_library_thread,
        )
        self._lib_thread = thread
        self._lib_worker = worker
        thread.start()

    def _clear_library_thread(self) -> None:
        """Drop the finished library thread/worker references."""
        self._lib_thread = None
        self._lib_worker = None

    def _set_library_running(self, running: bool) -> None:
        """Toggle controls while a library draw is in flight."""
        self._lib_running = running
        self.library_btn.setText("Sampling..." if running else "Sample library")
        self._update_run_buttons()

    @QtCore.Slot(object)
    def _on_library_finished(self, result: object) -> None:
        """Render a finished library draw."""
        self._set_library_running(False)
        if not isinstance(result, api.LibraryResult):
            return
        self._library = result
        self._render_library(result)
        self.statusBar().showMessage(
            f"Sampled {len(result.results)} sequence(s); "
            f"{result.distinct} distinct."
        )

    @QtCore.Slot(object)
    def _on_library_failed(self, error: object) -> None:
        """Clear the library panel and explain the failure in plain language."""
        self._set_library_running(False)
        self._reset_library()
        headline, detail = self._friendly_error(error)
        self.statusBar().showMessage(headline)
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
        box.setWindowTitle("Couldn't sample a library")
        box.setText(headline)
        if detail:
            box.setInformativeText(detail)
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        self._msgbox = box
        box.show()

    @QtCore.Slot(object)
    def _on_cand_finished(self, result: object) -> None:
        """Render the delivered candidate set and its splice audit."""
        self._set_candidates_running(False)
        if not isinstance(result, CandidatesResult):
            return
        self._render_candidates(result.candidate_set)
        self._render_splice_audit(result.audit, result.candidate_set)
        delivered = result.candidate_set.delivered()
        if delivered is None:
            self.statusBar().showMessage("No candidates were assembled.")
        else:
            self.statusBar().showMessage(
                f"{len(result.candidate_set.candidates)} candidate(s) ranked "
                f"({result.candidate_set.order_basis})."
            )

    @QtCore.Slot(object)
    def _on_cand_failed(self, error: object) -> None:
        """Clear the candidate panel and show a plain-language message on error."""
        self._set_candidates_running(False)
        self._reset_candidates()
        headline, detail = self._friendly_error(error)
        self.statusBar().showMessage(headline)
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
        box.setWindowTitle("Couldn't rank candidates")
        box.setText(headline)
        if detail:
            box.setInformativeText(detail)
        box.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok)
        self._msgbox = box
        box.show()

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
        """Toggle controls while an optimization is in flight."""
        self._optimize_running = running
        self.optimize_btn.setText("Optimizing..." if running else "Optimize")
        self._update_run_buttons()

    def _update_run_buttons(self) -> None:
        """Enable each run control from the live flags (one source of truth).

        Only one engine flow may run at a time, so every start button is gated on
        "nothing is running". Enablement is derived from the ``_*_running`` flags
        rather than from thread references, so a missed reference clear can never
        strand a control (the optimize-then-rank regression). ASSP additionally
        needs a delivered sequence to cross-check.
        """
        busy = (
            self._optimize_running
            or self._cand_running
            or self._cc_running
            or self._lib_running
        )
        self.optimize_btn.setEnabled(not busy)
        self.cancel_btn.setEnabled(self._optimize_running)
        self.rank_btn.setEnabled(not busy)
        self.library_btn.setEnabled(not busy)
        self.assp_btn.setEnabled(not busy and self._delivered() is not None)

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
        self._reset_crosscheck()
        self._update_run_buttons()

    def _reset_crosscheck(self) -> None:
        """Clear the ASSP panel to its empty state.

        Called whenever the delivered sequence changes or is dropped: a
        cross-check report is *about one sequence*, so leaving the previous
        report on screen beside a new design would attribute someone else's
        splice sites to it.
        """
        self.assp_table.setRowCount(0)
        self.assp_table.hide()
        self.assp_banner.setText(
            "Not run. ASSP is an <b>optional online</b> second opinion on the "
            "delivered sequence &mdash; nothing is sent until you press the button."
        )

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
        # A cross-check describes one specific sequence; this is a new one.
        self._reset_crosscheck()

        self.export_fasta_btn.setEnabled(True)
        self.export_json_btn.setEnabled(True)
        self._update_run_buttons()

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

    def _render_crosscheck(self, report: api.SpliceCrossCheck) -> None:
        """Render an ASSP cross-check, leading with what it is and is not.

        Mirrors the CLI's stderr block (``_print_splice_crosscheck``): the tags
        come first -- network-derived, uncalibrated, advisory, not in the manifest
        -- then either the graceful "unavailable" reason or the pooled risk and
        the localized sites. Nothing here touches the delivered
        :class:`~bt4.api.Result`, so an export stays byte-identical whether or not
        a cross-check was ever run (CLAUDE.md §6, §10.15).
        """
        tags = ["network-derived" if report.network_derived else "local"]
        tags.append("calibrated" if report.calibrated else "UNCALIBRATED")
        parts = [
            f"<b>{report.backend}</b> &mdash; {', '.join(tags)}; advisory only, "
            "<b>not</b> part of the run manifest and never exported."
        ]
        if not report.available:
            parts.append(f"<b>Unavailable:</b> {report.reason}")
            parts.append(
                "<i>An opt-in cross-check outage never fails a run &mdash; the "
                "delivered result and its local audit are unaffected.</i>"
            )
            self.assp_table.setRowCount(0)
            self.assp_table.hide()
            self.assp_banner.setText("<br>".join(parts))
            self.statusBar().showMessage("ASSP unavailable (the run is unaffected).")
            return

        parts.append(
            f"Pooled risk <b>{report.pooled_risk:.3f}</b> "
            f"(top-{report.top_k} log-odds, uncalibrated) &middot; "
            f"{len(report.sites)} predicted site(s)."
        )
        self.assp_banner.setText("<br>".join(parts))

        self.assp_table.setRowCount(len(report.sites))
        for row, site in enumerate(report.sites):
            values = (
                site.kind,
                str(site.position),
                f"{site.score:.3f}",
                site.site_class or "-",
            )
            for col, text in enumerate(values):
                self.assp_table.setItem(row, col, QtWidgets.QTableWidgetItem(text))
        self.assp_table.resizeColumnsToContents()
        self.assp_table.setVisible(bool(report.sites))
        self.statusBar().showMessage(
            f"ASSP cross-check: {len(report.sites)} site(s) "
            "(network-derived, advisory)."
        )

    # ---- library rendering ------------------------------------------------

    def _reset_library(self) -> None:
        """Clear the library panel to an honest empty state."""
        self._library = None
        self.library_table.setRowCount(0)
        self.library_view.set_sequence("")
        self.lib_banner.setText(
            "No library yet. Press <b>Sample library</b> to draw one."
        )
        self.lib_banner.setStyleSheet(theme.badge_qss(""))
        self.export_library_btn.setEnabled(False)

    def _render_library(self, result: api.LibraryResult) -> None:
        """Fill the library table and banner from a finished draw.

        The banner leads with the honest framing: these sequences are **sampled,
        not optimized**, so their ``SAMPLED`` certificate asserts no optimality
        and the table order is a draw order, not a ranking. Diversity is reported
        as measured (distinct count and mean pairwise Hamming), and any residual
        GLOBAL violation -- not enforced during sampling -- is visible per member
        in the hard/soft columns and highlighted in the sequence viewer.
        """
        members = result.results
        self.library_table.setRowCount(len(members))
        for row, member in enumerate(members):
            values = (
                str(row),
                f"{float(member.audit['cai']):.4f}",
                f"{member.metrics.gc * 100.0:.2f}",
                str(member.metrics.hard_violations),
                str(member.metrics.soft_violations),
                member.certificate.status.value,
            )
            for col, text in enumerate(values):
                self.library_table.setItem(row, col, QtWidgets.QTableWidgetItem(text))
        self.library_table.resizeColumnsToContents()

        self.lib_banner.setText(
            "<b>Sampled, not optimized.</b> Every member carries the SAMPLED "
            "certificate: no optimality claim, no expression claim, and the row "
            "order is the draw order, not a ranking.<br>"
            f"{len(members)} member(s) &middot; {result.distinct} distinct "
            f"&middot; mean pairwise difference "
            f"{result.mean_pairwise_hamming:.1%}.<br>"
            "<i>Non-local rules (max-repeat, uORF) are not enforced during "
            "sampling; any residual violation is counted per member and "
            "highlighted in the sequence below.</i>"
        )
        # Colour the banner from the members' own certificate, so the badge can
        # never drift from the claim the engine actually made.
        status = members[0].certificate.status.value if members else ""
        self.lib_banner.setStyleSheet(theme.badge_qss(status))
        self.export_library_btn.setEnabled(bool(members))
        if members:
            self.library_table.selectRow(0)

    def _show_library_member(self) -> None:
        """Show the selected library member's DNA with its violation highlights."""
        if self._library is None:
            return
        row = self.library_table.currentRow()
        members = self._library.results
        if not 0 <= row < len(members):
            self.library_view.set_sequence("")
            return
        member = members[row]
        self.library_view.set_sequence(member.dna, member.violations)

    # ---- candidates rendering ---------------------------------------------

    def _reset_candidates(self) -> None:
        """Clear the candidate panel to an honest empty state."""
        self.candidates_table.setRowCount(0)
        self.cand_banner.setText(
            "No candidate set yet. Press <b>Rank &amp; audit</b> to assemble one."
        )
        self.cand_banner.setStyleSheet(theme.badge_qss(""))
        self.splice_banner.clear()

    def _render_candidates(self, cand_set: api.CandidateSet) -> None:
        """Fill the candidate table and summary banner from ``cand_set``.

        Every score is recomputed per candidate from its own DNA (invariant #2).
        The banner states the calibration/order basis honestly: an uncalibrated
        head only annotates (discovery order, solver-delivered pick), so the
        table order is NOT a ranking and the UI says so.
        """
        candidates = cand_set.candidates
        self.candidates_table.setRowCount(len(candidates))
        for row, cand in enumerate(candidates):
            result = cand.result
            hard = result.metrics.hard_violations
            expr = (
                f"{cand.expression_score:.4g} {cand.expression_units}".strip()
            )
            marker = " ★" if row == cand_set.chosen else ""
            values = (
                f"{row}{marker}",
                cand.source,
                f"{float(result.audit['cai']):.4f}",
                f"{result.metrics.gc * 100.0:.2f}",
                expr,
                "yes" if cand.expression_calibrated else "no",
                str(hard),
                "",  # filled by the splice audit if it runs
            )
            for col, text in enumerate(values):
                self.candidates_table.setItem(
                    row, col, QtWidgets.QTableWidgetItem(text)
                )
        self.candidates_table.resizeColumnsToContents()

        if not candidates:
            self.cand_banner.setText("No candidates were assembled.")
            self.cand_banner.setStyleSheet(theme.badge_qss("relaxed"))
            return

        delivered = cand_set.delivered()
        model = delivered.expression_model if delivered is not None else "?"
        if cand_set.calibrated:
            basis = (
                f"Ranked by predicted expression ({model}); the ★ delivered "
                "candidate is the top-scoring one."
            )
            style = "proven_optimal"
        else:
            basis = (
                f"Expression head <b>{model}</b> is <b>uncalibrated</b>, so this is "
                "<b>discovery order, not a ranking</b>; the ★ delivered "
                "candidate is the solver's pick. Scores annotate only."
            )
            style = "heuristic"
        counts = (
            f"{cand_set.n_frontier} frontier + {cand_set.n_repeat_refined} "
            f"repeat-refined assembled; {cand_set.n_dedup_dropped} duplicate(s) and "
            f"{cand_set.n_dropped_cap} over-cap dropped."
        )
        self.cand_banner.setText(
            f"{basis}<br>{counts}<br><i>{cand_set.repeat_note}</i>"
        )
        self.cand_banner.setStyleSheet(theme.badge_qss(style))

    def _render_splice_audit(
        self, audit: api.SpliceAuditReport | None, cand_set: api.CandidateSet
    ) -> None:
        """Fill the per-candidate splice-flag column and the audit banner.

        The audit is advisory: every shipped backend is uncalibrated, so the
        banner leads with that and the flag counts are labeled a heuristic
        localization, never a calibrated risk claim.
        """
        if audit is None:
            self.splice_banner.setText(
                "<i>Splice audit skipped (no candidates).</i>"
            )
            return

        # Per-candidate count = distinct localized sites, NOT the raw cross-backend
        # flag sum: a site flagged by two backends is one site, so positions within
        # the audit's match window are merged before counting (else the column
        # would grow just from running more backends).
        for entry in audit.candidates:
            if not 0 <= entry.index < self.candidates_table.rowCount():
                continue
            n_sites = _distinct_site_count(
                sorted(f.position for b in entry.by_backend for f in b.flags),
                audit.match_window,
            )
            self.candidates_table.setItem(
                entry.index, 7, QtWidgets.QTableWidgetItem(str(n_sites))
            )
        header = self.candidates_table.horizontalHeaderItem(7)
        if header is not None:
            header.setText("Splice flags")
        self.candidates_table.resizeColumnsToContents()

        backends = ", ".join(audit.backends) if audit.backends else "none"
        calib = "calibrated" if audit.all_calibrated else "UNCALIBRATED (advisory)"
        agree = audit.agreement
        parts = [
            f"Splice audit &mdash; backends: <b>{backends}</b> &middot; {calib}."
        ]
        if len(audit.backends) >= 2:
            pairs = "; ".join(
                f"{a}-{b}: rho={rho:.2f}"
                for (a, b), rho in agree.rank_correlations.items()
            )
            parts.append(
                f"Cross-backend agreement: sign {agree.sign_agreement:.0%}"
                + (f" &middot; {pairs}" if pairs else "")
            )
        parts.append(
            "Flags localize residual cryptic sites (heuristic threshold "
            f"{audit.threshold:g}); they do <b>not</b> assert a calibrated splice "
            "risk and nothing was edited."
        )
        self.splice_banner.setText("<br>".join(parts))

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

    def _export_library(self) -> None:
        """Write every sampled library member to one multi-record FASTA file.

        Each record is named ``<job>_<index>`` and carries the sampler's honest
        framing in its header, so a downstream reader cannot mistake a sampled
        member for an optimized delivery.
        """
        if self._library is None or not self._library.results:
            return
        job = self.jobname_edit.text().strip() or "bt4"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export library FASTA",
            f"{job}_library.fasta",
            "FASTA (*.fasta *.fa);;All files (*)",
        )
        if not path:
            return
        records = [
            api.to_fasta(member.dna, header=f"{job}_{index} sampled")
            for index, member in enumerate(self._library.results)
        ]
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("".join(records))
        self.statusBar().showMessage(
            f"Wrote {len(records)} sampled sequence(s) to {path}"
        )

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
