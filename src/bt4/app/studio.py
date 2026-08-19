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
from collections.abc import Callable, Mapping
from html import escape
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

#: The engine's default organism (``api.OptimizeConfig().organism``).
_DEFAULT_ORGANISM = api.OptimizeConfig().organism

#: Shown as the first entry of the application-preset combo. BT4 applies no preset
#: by default (it stays regime-agnostic), so this is the selected entry at launch.
_NO_PRESET = "(none)"

# One row per objective weight: (config field, label, min, max, default, tooltip).
# CAI defaults to 1.0 to match the engine; everything else is off (0.0) until the
# user asks for it. Codon-pair bias allows a negative weight -- that is the
# attenuated-vaccine (deoptimization) case the engine supports and the app could
# not previously reach.
_WEIGHT_CONTROLS: tuple[tuple[str, str, float, float, float, str], ...] = (
    (
        "cai_weight", "CAI", 0.0, 100.0, 1.0,
        "Pull toward the organism's preferred codons (the classic codon-adaptation "
        "index). A weak proxy for expression on its own -- which is why it is one "
        "axis among several here.",
    ),
    (
        "gc_weight", "GC proximity", 0.0, 100.0, 0.0,
        "Pull the overall GC fraction toward the 'GC target' above. At 0 the GC "
        "target does nothing in a single solve (the frontier still sweeps it).",
    ),
    (
        "tai_weight", "tAI", 0.0, 100.0, 1.0,
        "Pull toward codons matched to the organism's tRNA pool, from real tRNA "
        "gene copy numbers. Strength only -- the tAI checkbox above switches it on,\n"
        "and it needs bundled tRNA data for the chosen organism.",
    ),
    (
        "cpb_weight", "Codon-pair bias", -100.0, 100.0, 1.0,
        "Favour codon PAIRS that are over-represented in the reference CDS set. A "
        "NEGATIVE weight deoptimizes them, which is how attenuated-vaccine designs "
        "are made. Strength/sign only -- the codon-pair checkbox switches it on,\n"
        "and it requires a reference CDS FASTA.",
    ),
    (
        "ramp_weight", "5' ramp", 0.0, 100.0, 0.0,
        "Favour less-adapted codons in the first few dozen codons. A shaping PRIOR, "
        "not a mechanism: the 5' benefit tracks reduced RNA structure rather than "
        "codon rarity, so use folding refinement for that lever.",
    ),
    (
        "cpg_weight", "CpG", 0.0, 100.0, 1.0,
        "Strength of the CpG objective; the direction (deplete/elevate) is the CpG "
        "control above, which also switches it off.",
    ),
    (
        "minmax_weight", "%MinMax", 0.0, 100.0, 1.0,
        "Strength of the codon-commonness (%MinMax) prior; the direction is the "
        "%MinMax control above, which also switches it off.",
    ),
)

# Shown on the reference-set combo and its label. A module constant because
# _update_reference_sets has to restore it whenever the organism changes.
_REFERENCE_SET_TOOLTIP = (
    "Which genes the CAI weights were counted over. 'highly_expressed' uses the "
    "300 most abundant proteins (PaxDb) -- the reference set CAI was defined on, "
    "so w = 1 marks the codon translation prefers. 'genome_wide' counts every "
    "gene, marking the codon that is merely most common. Only the sets bundled "
    "for the chosen organism are listed; neither is a measured expression "
    "prediction."
)

_METRIC_ROWS = (
    "CAI",
    # CAI is meaningless without the gene set its weights were counted over, so
    # the reference set is reported directly beneath it rather than left to the
    # control panel (which the user may have changed since the run).
    "CAI reference set",
    "GC %",
    "Length (nt)",
    "Scored codons",
    "Hard violations",
    "Soft violations",
    "Optimality",
    "Solver",
)

# Human labels for the non-local rules that report an enforcement status. A rule
# absent from this map still shows up (its key is prettified), so a NEW global rule
# is surfaced automatically rather than silently missing from the GUI -- which is
# exactly how the poly(A) and windowed-GC rules would otherwise have been invisible
# here while the CLI printed them.
_RULE_LABELS: dict[str, str] = {
    "max_repeat": "Max repeat",
    "uorf": "uORF",
    "polya_signal": "Poly(A) signal",
    "gc_window": "GC window",
}


def _audit_rows(audit: Mapping[str, object]) -> list[tuple[str, str]]:
    """Build metrics rows from whatever the run's audit actually reported.

    Driven by the audit dict rather than a fixed list, so a rule added to the
    engine later shows up in the GUI without anyone remembering to update it here.
    Covers the non-local rules' enforcement status (with residual counts), the
    optional objective read-outs, and any relaxation the solver had to make.
    """
    rows: list[tuple[str, str]] = []

    if "tai" in audit:
        rows.append(("tAI", f"{float(audit['tai']):.4f}"))  # type: ignore[arg-type]
    for key, label in (("cg_count", "CpG count"), ("ta_count", "UpA count")):
        if key in audit:
            rows.append((label, str(audit[key])))

    # Non-local rules: "<rule>_enforced" plus its "<rule>_residual" partner.
    for key in sorted(audit):
        if not key.endswith("_enforced"):
            continue
        rule = key[: -len("_enforced")]
        label = _RULE_LABELS.get(rule, rule.replace("_", " ").capitalize())
        status = str(audit[key])
        residual = audit.get(f"{rule}_residual")
        if status == "clean":
            value = "clean (fully enforced)"
        else:
            value = f"{status} - {residual} could not be removed"
        rows.append((label, value))

    if audit.get("relaxed_constraints"):
        relaxed = audit["relaxed_constraints"]
        names = ", ".join(relaxed) if isinstance(relaxed, list) else str(relaxed)
        rows.append(("Relaxed rules", f"{names} - could not be satisfied"))

    if "folding_dg" in audit:
        calibrated = bool(audit.get("folding_calibrated", False))
        units = "kcal/mol" if calibrated else "arbitrary units, UNCALIBRATED"
        rows.append(("5' folding", f"{float(audit['folding_dg']):.3f} ({units})"))  # type: ignore[arg-type]

    return rows

# Above this residue count, warn before starting: the exact frontier sweep can
# take a noticeable amount of time (the run is cancelable either way).
_LONG_PROTEIN_WARN = 500

# How long a close waits for each still-running worker thread (see closeEvent).
_CLOSE_WAIT_MS = 5000

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


class EnzymeCompleter(QtWidgets.QCompleter):
    """Autocompletes the *last* comma-separated entry in the enzyme field.

    The REBASE-derived catalog holds hundreds of enzymes, so typing a name from
    memory is the weak point of an otherwise plain text field. A stock
    ``QCompleter`` matches against the entire line, which breaks as soon as the
    user has already listed one enzyme; this one splits on commas and completes
    only the token being typed, then substitutes it back in place.
    """

    def __init__(self, names: list[str], parent: QtCore.QObject | None = None) -> None:
        super().__init__(names, parent)
        self.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        self.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
        self.setCompletionMode(
            QtWidgets.QCompleter.CompletionMode.PopupCompletion
        )

    def splitPath(self, path: str) -> list[str]:
        """Complete against the text after the last comma only."""
        return [path.split(",")[-1].strip()]

    def pathFromIndex(self, index: QtCore.QModelIndex) -> str:
        """Rebuild the full field with the completed name in the last slot."""
        completed = str(super().pathFromIndex(index))
        widget = self.widget()
        if not isinstance(widget, QtWidgets.QLineEdit):
            return completed
        prefix = widget.text().rsplit(",", 1)[0] if "," in widget.text() else ""
        return f"{prefix}, {completed}" if prefix.strip() else completed


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
        # The organism/reference set the LAST DELIVERED run used. Renders read
        # these, never the live combos: a render can happen long after the run
        # (a theme switch re-renders), and by then the user may have changed the
        # controls -- reading them would relabel or recompute an old result
        # under a table it was not produced with.
        self._delivered_tables: tuple[str, str] = ("", "")
        # The reference set the USER picked, as distinct from one the engine
        # forced. Selecting an organism that has only genome_wide rewrites the
        # combo to it; without remembering the real choice separately, switching
        # back would keep genome_wide and silently opt every later run out of the
        # default -- the same "quietly hands out a codon-commonness index" failure
        # the organism default had, reachable in two clicks.
        self._reference_choice: str = ""
        # (organism, reference_set) of the run currently in flight; promoted to
        # _delivered_tables when it finishes.
        self._running_tables: tuple[str, str] = ("", "")
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
        controls = self._build_controls()
        controls_scroll.setWidget(controls)
        controls_scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        # Size the column from the form's OWN width requirement, plus room for the
        # vertical scrollbar. Horizontal scrolling is off (a form that scrolls
        # sideways is miserable to use), so a hard-coded minimum narrower than the
        # form silently clips its right edge -- checkbox labels and placeholders
        # were being cut mid-word with no way to reveal them. Clamped so the
        # controls can never eat the whole window on a small screen.
        scrollbar_allowance = controls_scroll.verticalScrollBar().sizeHint().width()
        natural_width = controls.sizeHint().width() + scrollbar_allowance + 4
        column_width = max(320, min(natural_width, 560))
        controls_scroll.setMinimumWidth(column_width)
        splitter.addWidget(controls_scroll)
        splitter.addWidget(self._build_results())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        # Let the results side absorb every resize: the control column keeps its
        # natural width so the knobs stay readable, and the plots/tables shrink.
        splitter.setCollapsible(0, False)
        splitter.setSizes([column_width, 900])
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
        # Without this the combo defaults to the alphabetically first organism,
        # which is A. thaliana -- the ONE organism with no highly-expressed table.
        # A freshly launched Studio would therefore hand out a codon-commonness
        # index, silently opting the app out of the default every other surface
        # uses. The engine's own default organism is the right starting point.
        if _DEFAULT_ORGANISM in api.available_organisms():
            self.organism_combo.setCurrentText(_DEFAULT_ORGANISM)
        self.organism_combo.setAccessibleName("Target organism")
        self._add_row(
            form, "Organism", self.organism_combo,
            "Codon-usage table used for CAI and codon choice; it decides which "
            "synonymous codons are preferred for the target organism.",
        )

        self.reference_combo = QtWidgets.QComboBox()
        self.reference_combo.setAccessibleName("CAI reference gene set")
        self._add_row(form, "Reference set", self.reference_combo, _REFERENCE_SET_TOOLTIP)

        self.jobname_edit = QtWidgets.QLineEdit()
        self.jobname_edit.setPlaceholderText("optional job name")
        self.jobname_edit.setAccessibleName("Job name")
        self._add_row(
            form, "Job name", self.jobname_edit,
            "Optional label for this run, used only for your reference and in "
            "exported file headers.",
        )

        self.ctx_upstream_edit = QtWidgets.QPlainTextEdit()
        self.ctx_upstream_edit.setPlaceholderText(
            "optional: 5'UTR / vector sequence before the CDS"
        )
        self.ctx_upstream_edit.setAccessibleName("Upstream construct context")
        self.ctx_upstream_edit.setMaximumHeight(60)
        self._add_row(
            form, "5' context", self.ctx_upstream_edit,
            "The sequence that will sit immediately BEFORE your CDS -- the 5'UTR, "
            "and as much of the vector backbone before it as you know. Give it and "
            "BT4 stops designing in a vacuum: every rule below is also checked "
            "across the junction (so a restriction site or repeat formed half in "
            "the leader and half in codon 1 is caught), and uORF pairing can see an "
            "ATG in the leader that reads into your CDS. Write N for unknown bases; "
            "the flank is cut at the N nearest the CDS. It is context only -- it is "
            "never part of the designed sequence and never exported.",
        )

        self.ctx_downstream_edit = QtWidgets.QPlainTextEdit()
        self.ctx_downstream_edit.setPlaceholderText("optional: sequence after the CDS")
        self.ctx_downstream_edit.setAccessibleName("Downstream construct context")
        self.ctx_downstream_edit.setMaximumHeight(60)
        self._add_row(
            form, "3' context", self.ctx_downstream_edit,
            "The sequence immediately AFTER your CDS (3'UTR / backbone). Used by the "
            "whole-construct audit so a defect formed at the 3' junction is visible "
            "too. Context only; never part of the designed sequence.",
        )

        self.context_prov_combo = QtWidgets.QComboBox()
        self.context_prov_combo.addItem("omit (record only its length)", "omit")
        self.context_prov_combo.addItem("hash (record a content hash)", "hash")
        self.context_prov_combo.setAccessibleName("Context provenance policy")
        self._add_row(
            form, "Context in stamp", self.context_prov_combo,
            "What the run's provenance stamp records about the sequence above. "
            "'omit' keeps only its length, so the stamp cannot fingerprint your "
            "backbone. 'hash' stores a content hash, which makes the run fully "
            "reproducible from its stamp but does identify the backbone to anyone "
            "who already has a copy to test against. Your call -- BT4 will not pick "
            "for you. Either way the sequence itself is never sent anywhere.",
        )

        self.preset_combo = QtWidgets.QComboBox()
        self.preset_combo.addItem(_NO_PRESET)
        for app_preset in api.available_presets():
            self.preset_combo.addItem(app_preset.label, app_preset.key)
        self.preset_combo.setAccessibleName("Application preset")
        self._add_row(
            form, "Preset", self.preset_combo,
            "Optional starting point for a kind of construct (AAV, lentiviral, "
            "plasmid, IVT mRNA, or plain synthesis limits). NO preset is applied by "
            "default. Choosing one fills the controls below in -- you can then change "
            "anything, and what you change wins. Hover an entry for why it sets what "
            "it sets; these are conventions, not calibrated expression predictions.",
        )
        for index, app_preset in enumerate(api.available_presets(), start=1):
            self.preset_combo.setItemData(
                index,
                f"{app_preset.description}\n\nWhy: {app_preset.rationale}",
                QtCore.Qt.ItemDataRole.ToolTipRole,
            )
        self.preset_combo.activated.connect(self._on_preset_chosen)

        self.gc_spin = QtWidgets.QDoubleSpinBox()
        self.gc_spin.setRange(0.0, 1.0)
        self.gc_spin.setSingleStep(0.05)
        self.gc_spin.setValue(0.55)
        self.gc_spin.setAccessibleName("GC target fraction")
        self._add_row(
            form, "GC target", self.gc_spin,
            "Desired overall GC fraction (0-1). This is a soft objective: Studio "
            "sweeps it as a frontier axis, so it always shapes the trade-off here. "
            "(In a single CLI/API solve it only applies when its weight is > 0.) "
            "It nudges toward this value; it is not a hard bound.",
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

        self.gc_window_spin = QtWidgets.QSpinBox()
        self.gc_window_spin.setRange(0, 200)
        self.gc_window_spin.setValue(0)
        self.gc_window_spin.setSpecialValueText("off")
        self.gc_window_spin.setAccessibleName("GC window length in nucleotides (0 = off)")
        self._add_row(
            form, "GC window (nt)", self.gc_window_spin,
            "Bound the GC fraction of EVERY window of this many nucleotides -- the "
            "rule synthesis vendors actually specify (e.g. 25-65% per 50 nt). A "
            "whole-sequence GC number can look fine while a local stretch is extreme, "
            "and the local extreme is what fails synthesis. Up to 12 nt this stays "
            "exact in the trellis (proven-optimal); wider windows are enforced by the "
            "refinement pass instead (labeled heuristic, residuals reported), because "
            "the exact search grows exponentially with the window. 0 = off.",
        )

        gc_bounds = QtWidgets.QWidget()
        gc_bounds_row = QtWidgets.QHBoxLayout(gc_bounds)
        gc_bounds_row.setContentsMargins(0, 0, 0, 0)
        self.gc_window_min_spin = QtWidgets.QDoubleSpinBox()
        self.gc_window_min_spin.setRange(0.0, 1.0)
        self.gc_window_min_spin.setSingleStep(0.05)
        self.gc_window_min_spin.setValue(0.25)
        self.gc_window_min_spin.setAccessibleName("Minimum GC fraction per window")
        self.gc_window_max_spin = QtWidgets.QDoubleSpinBox()
        self.gc_window_max_spin.setRange(0.0, 1.0)
        self.gc_window_max_spin.setSingleStep(0.05)
        self.gc_window_max_spin.setValue(0.65)
        self.gc_window_max_spin.setAccessibleName("Maximum GC fraction per window")
        gc_bounds_row.addWidget(QtWidgets.QLabel("min"))
        gc_bounds_row.addWidget(self.gc_window_min_spin)
        gc_bounds_row.addWidget(QtWidgets.QLabel("max"))
        gc_bounds_row.addWidget(self.gc_window_max_spin)
        self._add_row(
            form, "GC window bounds", gc_bounds,
            "Lowest and highest GC fraction any window may have (used only when "
            "'GC window (nt)' is on). The vendor-typical band is 0.25-0.65.",
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
            "Extra DNA substrings to ban, comma-separated (literal A/C/G/T only -- "
            "for a degenerate site use 'Extra sites (IUPAC)' below).",
        )

        self.rc_check = QtWidgets.QCheckBox("also ban each motif's reverse complement")
        self.rc_check.setChecked(True)
        self.rc_check.setAccessibleName("Ban reverse complements of forbidden motifs")
        self._add_row(
            form, "Both strands", self.rc_check,
            "On (the default), a forbidden motif is banned on both strands, which "
            "is what you want for anything double-stranded DNA is read on. Turn it "
            "off for a rule that is genuinely sense-strand-only, such as an mRNA "
            "regulatory motif that means nothing on the template strand.",
        )

        self._add_forbidden_presets(form)

        self.enzymes_edit = QtWidgets.QLineEdit()
        self.enzymes_edit.setPlaceholderText("comma-separated, e.g. EcoRI,BsaI")
        self.enzymes_edit.setAccessibleName("Restriction enzymes to avoid")
        enzyme_names = list(api.available_enzymes())
        self.enzyme_completer = EnzymeCompleter(enzyme_names, self)
        self.enzymes_edit.setCompleter(self.enzyme_completer)
        self._add_row(
            form, "Restriction sites", self.enzymes_edit,
            "Restriction enzymes whose recognition sites (and reverse "
            "complements) must not appear, comma-separated (e.g. EcoRI,BsaI). "
            f"Start typing to search the {len(enzyme_names)}-enzyme REBASE "
            "catalog. Only the recognition sequence is modelled -- not cut "
            "position, star activity, or methylation sensitivity.",
        )

        self.enzyme_sites_edit = QtWidgets.QLineEdit()
        self.enzyme_sites_edit.setPlaceholderText("comma-separated IUPAC, e.g. GANTC")
        self.enzyme_sites_edit.setAccessibleName("Extra recognition sites to avoid (IUPAC)")
        self._add_row(
            form, "Extra sites (IUPAC)", self.enzyme_sites_edit,
            "Ban recognition sites directly, written in IUPAC (e.g. GANTC, "
            "CCNNNNNNNGG); each site's reverse complement is banned too. Use this "
            "for an enzyme the catalog above does not carry: 'Forbidden motifs' "
            "accepts only literal A/C/G/T, so a degenerate site can only be given "
            "here.",
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

        self.tandem_copies_spin = QtWidgets.QSpinBox()
        self.tandem_copies_spin.setRange(2, 10)
        self.tandem_copies_spin.setValue(3)
        self.tandem_copies_spin.setAccessibleName("Tandem-repeat copy count")
        self._add_row(
            form, "Tandem copies", self.tandem_copies_spin,
            "How many back-to-back copies count as a banned tandem repeat (used "
            "only when 'Tandem unit' is on). 3 copies of a 3-nt unit bans "
            "CAGCAGCAG; 2 would also ban CAGCAG.",
        )

        self.inverted_spin = QtWidgets.QSpinBox()
        self.inverted_spin.setRange(0, 20)
        self.inverted_spin.setValue(0)
        self.inverted_spin.setSpecialValueText("off")
        self.inverted_spin.setAccessibleName("Inverted-repeat stem length (0 = off)")
        self._add_row(
            form, "Hairpin stem", self.inverted_spin,
            "Ban a hairpin (inverted repeat) with arms of this length -- a stem "
            "that folds back on itself and can occlude ribosome loading. Note this "
            "is an EXACT trellis rule, so its cost grows steeply with stem+loop; "
            "for a broad repeat limit prefer 'Max repeat length', which is "
            "reverse-complement aware and already covers hairpins. 0 = off.",
        )

        self.inverted_loop_spin = QtWidgets.QSpinBox()
        self.inverted_loop_spin.setRange(0, 12)
        self.inverted_loop_spin.setValue(0)
        self.inverted_loop_spin.setAccessibleName("Inverted-repeat maximum loop length")
        self._add_row(
            form, "Hairpin loop", self.inverted_loop_spin,
            "Longest loop allowed between the hairpin arms (used only when "
            "'Hairpin stem' is on). 0 means only arms that sit directly "
            "back-to-back are banned; a larger loop catches more real hairpins but "
            "widens the exact search sharply.",
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

        self.polya_check = QtWidgets.QCheckBox("avoid functional poly(A) signals")
        self.polya_check.setAccessibleName("Avoid functional polyadenylation signals")
        self._add_row(
            form, "Poly(A)", self.polya_check,
            "Forbid a poly(A) hexamer (AATAAA/ATTAAA) only when a U/GU-rich element "
            "follows it downstream -- the bipartite signal the cleavage machinery "
            "actually recognises (CPSF binds the hexamer, CstF the downstream "
            "element). A premature poly(A) site truncates the transcript. This is "
            "deliberately MORE permissive than the 'Poly-A signals' forbidden preset "
            "above, which bans every bare hexamer even where no site could form. "
            "Non-local, so enforced by a refinement pass (result labeled heuristic; "
            "any it cannot remove are reported). Structural, not a calibrated "
            "cleavage prediction.",
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

        self._add_objective_weights(form)
        self._add_budgets(form)

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

        # Both of these depend on the chosen organism: tAI needs a bundled tRNA
        # table, and the reference sets differ per organism (only A. thaliana
        # lacks a highly-expressed one). Repopulating from the engine rather than
        # from a hard-coded list means the app can never offer a table that is
        # not actually bundled.
        self.organism_combo.currentTextChanged.connect(
            lambda *_: self._update_organism_dependent_controls()
        )
        # Only a real user edit reaches this slot: _update_reference_sets blocks
        # signals while it repopulates.
        self.reference_combo.currentTextChanged.connect(self._on_reference_set_chosen)
        self._update_organism_dependent_controls()

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

    def _add_objective_weights(self, form: QtWidgets.QFormLayout) -> None:
        """Add a spin box per objective term so weights are continuous, not on/off.

        Without these every optional objective was pinned to 1.0, so a user could
        not trade (say) CpG depletion against CAI, and the negative ``cpb_weight``
        that expresses codon-pair *de*-optimization was unreachable from the app.
        """
        container = QtWidgets.QWidget()
        grid = QtWidgets.QFormLayout(container)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(2)
        self.weight_spins: dict[str, QtWidgets.QDoubleSpinBox] = {}
        for field, label, low, high, default, tip in _WEIGHT_CONTROLS:
            spin = QtWidgets.QDoubleSpinBox()
            spin.setRange(low, high)
            spin.setSingleStep(0.1)
            spin.setDecimals(2)
            spin.setValue(default)
            spin.setAccessibleName(f"{label} objective weight")
            spin.setToolTip(tip)
            spin.setAccessibleDescription(tip)
            self.weight_spins[field] = spin
            row_label = QtWidgets.QLabel(label)
            row_label.setBuddy(spin)
            row_label.setToolTip(tip)
            grid.addRow(row_label, spin)
        self._add_row(
            form, "Objective weights", container,
            "Relative pull of each objective in a single solve. The Pareto frontier "
            "sweeps its own weights across the active axes, so these mainly shape "
            "which axes are active and how a single solve trades them off. A weight "
            "of 0 switches an objective off; codon-pair bias accepts a NEGATIVE "
            "weight, which deoptimizes pairs (attenuated-vaccine design).",
        )

    def _add_budgets(self, form: QtWidgets.QFormLayout) -> None:
        """Add the hard whole-sequence count budgets (GC and a dinucleotide).

        These are *hard* windows enforced exactly by the amount-bucketed DP, as
        opposed to the soft GC-proximity objective -- the engine has always had
        them and the app could not reach them.
        """
        gc_box = QtWidgets.QWidget()
        gc_row = QtWidgets.QHBoxLayout(gc_box)
        gc_row.setContentsMargins(0, 0, 0, 0)
        self.gc_min_spin = QtWidgets.QSpinBox()
        self.gc_min_spin.setRange(-1, 100000)
        self.gc_min_spin.setValue(-1)
        self.gc_min_spin.setSpecialValueText("off")
        self.gc_min_spin.setAccessibleName("Minimum total GC base count")
        self.gc_max_spin = QtWidgets.QSpinBox()
        self.gc_max_spin.setRange(-1, 100000)
        self.gc_max_spin.setValue(-1)
        self.gc_max_spin.setSpecialValueText("off")
        self.gc_max_spin.setAccessibleName("Maximum total GC base count")
        gc_row.addWidget(QtWidgets.QLabel("min"))
        gc_row.addWidget(self.gc_min_spin)
        gc_row.addWidget(QtWidgets.QLabel("max"))
        gc_row.addWidget(self.gc_max_spin)
        self._add_row(
            form, "GC count budget", gc_box,
            "HARD bound on the total number of G/C bases in the whole sequence "
            "(a count, not a fraction), enforced exactly. This is the strict "
            "counterpart of the soft 'GC target' above. Not combinable with the "
            "refinement-enforced rules (max repeat length, uORFs, a wide GC window).",
        )

        dinuc_box = QtWidgets.QWidget()
        dinuc_row = QtWidgets.QHBoxLayout(dinuc_box)
        dinuc_row.setContentsMargins(0, 0, 0, 0)
        self.dinuc_combo = QtWidgets.QComboBox()
        self.dinuc_combo.addItem("off", "")
        self.dinuc_combo.addItem("CpG (CG)", "CG")
        self.dinuc_combo.addItem("UpA (TA)", "TA")
        self.dinuc_combo.setAccessibleName("Dinucleotide to budget")
        self.dinuc_min_spin = QtWidgets.QSpinBox()
        self.dinuc_min_spin.setRange(-1, 100000)
        self.dinuc_min_spin.setValue(-1)
        self.dinuc_min_spin.setSpecialValueText("off")
        self.dinuc_min_spin.setAccessibleName("Minimum dinucleotide count")
        self.dinuc_max_spin = QtWidgets.QSpinBox()
        self.dinuc_max_spin.setRange(-1, 100000)
        self.dinuc_max_spin.setValue(-1)
        self.dinuc_max_spin.setSpecialValueText("off")
        self.dinuc_max_spin.setAccessibleName("Maximum dinucleotide count")
        dinuc_row.addWidget(self.dinuc_combo)
        dinuc_row.addWidget(QtWidgets.QLabel("min"))
        dinuc_row.addWidget(self.dinuc_min_spin)
        dinuc_row.addWidget(QtWidgets.QLabel("max"))
        dinuc_row.addWidget(self.dinuc_max_spin)
        self._add_row(
            form, "Dinucleotide budget", dinuc_box,
            "HARD bound on how many times a dinucleotide occurs in the whole "
            "sequence -- a CpG cap for stealth, or a CpG floor for an "
            "immunostimulatory design; UpA is a comparable stability knob. Enforced "
            "exactly (a 2-mer straddling a codon boundary is still counted), and "
            "mutually exclusive with the GC count budget.",
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
        self.export_gb_btn = QtWidgets.QPushButton("Export GenBank")
        self.export_gb_btn.setAccessibleName("Export GenBank")
        self.export_gb_btn.setToolTip(
            "Write an annotated GenBank map. Every residual violation becomes a "
            "misc_feature span, so a defect the optimizer could not remove is "
            "visible where it occurs when you open the file in SnapGene, Benchling "
            "or ApE. Includes the 5'/3' construct context when you supplied it."
        )
        self.export_gb_btn.clicked.connect(self._export_genbank)
        exports.addWidget(self.export_fasta_btn)
        exports.addWidget(self.export_json_btn)
        exports.addWidget(self.export_gb_btn)
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

        self.splice_attested_check = QtWidgets.QCheckBox("honor fidelity attestation")
        self.splice_attested_check.setAccessibleName(
            "Honor a committed splice fidelity attestation"
        )
        attested = api.attested_backends_available()
        if attested:
            self.splice_attested_check.setToolTip(
                "Promote an attested backend to calibrated=True for this run "
                f"(attested and installed: {', '.join(attested)}).\n\n"
                "An attestation proves BT4's WRAPPER reproduces the published model "
                "bit-for-bit -- integration fidelity. It does NOT show the scores are "
                "calibrated probabilities for coding sequence; that is a separate, "
                "still-unmet gate, and these models are measured weakest on exonic "
                "variants (median prAUC 0.419 vs 0.773 intronic), which is BT4's "
                "entire regime."
            )
        else:
            self.splice_attested_check.setEnabled(False)
            self.splice_attested_check.setToolTip(
                "No installed backend carries a committed fidelity attestation, so "
                "this would do nothing. Pangolin ships one; it needs the GPL pangolin "
                "package and its weights installed (see "
                "docs/DESIGN_splice_cnn_calibration.md)."
            )
        controls.addWidget(self.splice_attested_check)
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
            "context for the RiboNN model only and is never exported. NOTE: this is "
            "separate from the '5' context' box on the Design panel, which is what "
            "makes the DESIGN construct-aware; this one only feeds the expression "
            "model's fixed UTR slots.",
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
        self._add_action(file_menu, "&Open protein FASTA...", self._open_protein, "Ctrl+O")
        file_menu.addSeparator()
        self.act_export_fasta = self._add_action(
            file_menu, "Export &FASTA...", self._export_fasta, "Ctrl+E"
        )
        self.act_export_json = self._add_action(
            file_menu, "Export &JSON...", self._export_json, "Ctrl+J"
        )
        self.act_export_genbank = self._add_action(
            file_menu, "Export &GenBank...", self._export_genbank, "Ctrl+G"
        )
        self.act_export_library = self._add_action(
            file_menu, "Export &library FASTA...", self._export_library, "Ctrl+Shift+E"
        )
        file_menu.addSeparator()
        self._add_action(file_menu, "&Quit", self.close, "Ctrl+Q")

        run_menu = menubar.addMenu("&Run")
        self.act_optimize = self._add_action(
            run_menu, "&Optimize", self._start_optimize, "Ctrl+R"
        )
        self.act_cancel = self._add_action(
            run_menu, "&Cancel", self._cancel_optimize, "Esc"
        )
        run_menu.addSeparator()
        self._add_action(
            run_menu, "&Validate a sequence...", self._validate_sequence, "Ctrl+Shift+V"
        )
        run_menu.addSeparator()
        self.act_rank = self._add_action(
            run_menu, "Rank && &audit candidates", self._start_candidates, "Ctrl+K"
        )
        self.act_library = self._add_action(
            run_menu, "&Sample library", self._start_library, "Ctrl+L"
        )
        run_menu.addSeparator()
        self.act_assp = self._add_action(
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
            self.reference_combo,
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
            self.export_gb_btn,
            self.assp_btn,
            self.cand_n_spin,
            self.cand_repeat_spin,
            self.splice_cnn_check,
            self.splice_attested_check,
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

    def _on_reference_set_chosen(self, name: str) -> None:
        """Record an explicit user pick of the reference set."""
        if name:
            self._reference_choice = name

    def _update_organism_dependent_controls(self) -> None:
        """Refresh every control whose valid values depend on the organism."""
        self._update_reference_sets()
        self._update_tai_availability()

    def _update_reference_sets(self) -> None:
        """Re-list the reference sets bundled for the chosen organism.

        The engine is the only authority on what is bundled, so the combo is
        rebuilt from :func:`api.available_reference_sets` on every organism
        change -- a stale entry would let the user ask for a table that does not
        exist and get an error at run time instead of an honest, absent option.
        The previous choice is kept when the new organism also has it.
        """
        try:
            sets = api.available_reference_sets(self.organism_combo.currentText())
        except ValueError:
            sets = ()
        self.reference_combo.blockSignals(True)
        self.reference_combo.clear()
        self.reference_combo.addItems(list(sets))
        # Restore the user's own choice when this organism has it; otherwise fall
        # to sets[0], which is the organism's default (highly-expressed wherever
        # it exists). Reading back the combo's current text instead would treat a
        # forced value as a preference.
        if self._reference_choice in sets:
            self.reference_combo.setCurrentText(self._reference_choice)
        self.reference_combo.blockSignals(False)
        # One option is not a choice; say why rather than showing a dead control.
        # BOTH branches set the tooltip: setting it only in the single-set branch
        # would leave the "only one is bundled" text stuck on afterwards, telling
        # the user the opposite of what the (now populated) combo offers.
        self.reference_combo.setEnabled(len(sets) > 1)
        self.reference_combo.setToolTip(
            f"Only the {sets[0]} reference set is bundled for this organism."
            if len(sets) == 1
            else _REFERENCE_SET_TOOLTIP
        )

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
        explaining which entries are unknown) if any is not in the catalog.

        The catalog holds hundreds of enzymes, so an error offers the closest
        matching names rather than listing every one -- a wall of 500+ names
        hides the answer instead of giving it.

        Each suggestion is shown **with its recognition site**, and labelled as
        matched by spelling. The suggestions come from a fuzzy *name* match with
        no notion of site similarity, so a near-miss usually cuts something
        entirely different (``NotI`` is ``GCGGCCGC``, ``NcoI`` is ``CCATGG``).
        Offering a bare list would invite the user to accept a substitute that
        does not ban the site they care about -- and the run would then report
        proven-optimal with zero violations while their real site remained.
        """
        raw = [e.strip() for e in self.enzymes_edit.text().split(",") if e.strip()]
        catalog = {name.lower(): name for name in api.available_enzymes()}
        canonical: list[str] = []
        unknown: list[str] = []
        for entry in raw:
            name = catalog.get(entry.lower())
            (canonical if name else unknown).append(name or entry)
        if unknown:
            hints = []
            for entry in unknown:
                close = api.enzyme_suggestions(entry)
                if close:
                    sites = ", ".join(f"{hit} ({api.resolve_enzyme(hit)})" for hit in close)
                    hints.append(f"{entry} -> {sites}")
            if hints:
                detail = (
                    "Closest catalog names by SPELLING, not by recognition site -- "
                    "check the sequence before substituting: " + "; ".join(hints) + ". "
                )
            else:
                detail = f"The catalog has {len(catalog)} enzymes. "
            detail += (
                "If BT4 does not carry your enzyme, ban its recognition sequence "
                "directly instead of substituting a similar name: put it in 'Extra "
                "sites (IUPAC)', which accepts degenerate codes such as GANTC "
                "('Forbidden motifs' takes only literal A/C/G/T)."
            )
            self._warn(
                "Unknown restriction enzyme",
                f"Not in the catalog: {', '.join(unknown)}.",
                detail,
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
        # Each optional objective has an on/off control (a combo or a checkbox) and
        # a strength spin. The control gates; the spin sets magnitude, so a weight
        # is continuous rather than pinned to 1.0 as it used to be.
        weight = {field: spin.value() for field, spin in self.weight_spins.items()}
        cpg = self.cpg_combo.currentText()
        cpg_weight = 0.0 if cpg == "off" else weight["cpg_weight"]
        cpg_mode = "deplete" if cpg == "off" else cpg
        minmax = self.minmax_combo.currentText()
        minmax_weight = 0.0 if minmax == "off" else weight["minmax_weight"]
        minmax_direction = "max" if minmax == "off" else minmax
        cpb_weight, cpb_cds = self._read_cpb()
        tandem = self.tandem_spin.value()
        inverted = self.inverted_spin.value()
        gc_window = self.gc_window_spin.value()
        extra_sites = tuple(
            s.strip().upper()
            for s in self.enzyme_sites_edit.text().split(",")
            if s.strip()
        )
        gc_min, gc_max = self.gc_min_spin.value(), self.gc_max_spin.value()
        dinuc = str(self.dinuc_combo.currentData() or "")
        dinuc_min, dinuc_max = self.dinuc_min_spin.value(), self.dinuc_max_spin.value()
        config = api.OptimizeConfig(
            organism=self.organism_combo.currentText(),
            reference_set=self.reference_combo.currentText() or None,
            gc_target=self.gc_spin.value(),
            cai_weight=weight["cai_weight"],
            gc_weight=weight["gc_weight"],
            ramp_weight=weight["ramp_weight"],
            max_homopolymer=homo if homo > 0 else None,
            max_gc_run=gc_run if gc_run > 0 else None,
            gc_window_nt=gc_window if gc_window > 0 else None,
            gc_window_min=self.gc_window_min_spin.value(),
            gc_window_max=self.gc_window_max_spin.value(),
            max_repeat_length=repeat if repeat > 0 else None,
            forbidden_motifs=motifs,
            forbidden_presets=presets,
            avoid_reverse_complement=self.rc_check.isChecked(),
            restriction_enzymes=enzymes,
            restriction_extra_sites=extra_sites,
            cpg_weight=cpg_weight,
            cpg_mode=cpg_mode,
            cpb_weight=cpb_weight,
            cpb_reference_cds=cpb_cds,
            minmax_weight=minmax_weight,
            minmax_direction=minmax_direction,
            tandem_unit=tandem if tandem > 0 else None,
            tandem_copies=self.tandem_copies_spin.value(),
            inverted_stem=inverted if inverted > 0 else None,
            inverted_loop=self.inverted_loop_spin.value(),
            avoid_splice_sites=self.splice_check.isChecked(),
            avoid_internal_start=self.internal_start_check.isChecked(),
            avoid_polya=self.polya_check.isChecked(),
            avoid_uorf=self.uorf_check.isChecked(),
            uorf_region_nt=self.uorf_region_spin.value(),
            tai_weight=weight["tai_weight"] if self.tai_check.isChecked() else 0.0,
            gc_min=gc_min if gc_min >= 0 else None,
            gc_max=gc_max if gc_max >= 0 else None,
            dinuc_budget=dinuc or None,
            dinuc_min=dinuc_min if dinuc and dinuc_min >= 0 else None,
            dinuc_max=dinuc_max if dinuc and dinuc_max >= 0 else None,
            beam=beam if beam > 0 else None,
            seed=0,
            application_preset=str(self.preset_combo.currentData() or ""),
            context=self._build_context(),
            context_provenance=str(self.context_prov_combo.currentData() or "omit"),
        )
        return config, self.steps_spin.value()

    def _build_context(self) -> api.ConstructContext | None:
        """Read the construct-context boxes, or ``None`` when both are empty.

        Returning ``None`` for empty input matters: it is what keeps a run without
        construct context byte-identical to one from before the feature existed.
        """
        upstream = "".join(self.ctx_upstream_edit.toPlainText().split())
        downstream = "".join(self.ctx_downstream_edit.toPlainText().split())
        if not upstream and not downstream:
            return None
        return api.ConstructContext(upstream=upstream, downstream=downstream)

    def _on_preset_chosen(self, index: int) -> None:
        """Fill the design controls in from the chosen application preset.

        The preset's values are written into the visible controls rather than
        applied invisibly at run time, so the user can see exactly what it changed
        and edit anything afterwards -- what they change wins. Choosing ``(none)``
        leaves every control alone (there is no default preset to restore).
        """
        key = self.preset_combo.itemData(index)
        if not key:
            return
        preset = api.resolve_preset(str(key))
        setters: dict[str, Callable[[object], None]] = {
            "gc_window_nt": lambda v: self.gc_window_spin.setValue(int(v)),  # type: ignore[arg-type]
            "gc_window_min": lambda v: self.gc_window_min_spin.setValue(float(v)),  # type: ignore[arg-type]
            "gc_window_max": lambda v: self.gc_window_max_spin.setValue(float(v)),  # type: ignore[arg-type]
            "max_homopolymer": lambda v: self.homo_spin.setValue(int(v)),  # type: ignore[arg-type]
            "max_gc_run": lambda v: self.gc_run_spin.setValue(int(v)),  # type: ignore[arg-type]
            "max_repeat_length": lambda v: self.repeat_spin.setValue(int(v)),  # type: ignore[arg-type]
            "cpg_weight": lambda v: self.weight_spins["cpg_weight"].setValue(float(v)),  # type: ignore[arg-type]
            "cpg_mode": lambda v: self.cpg_combo.setCurrentText(str(v)),
            "avoid_splice_sites": lambda v: self.splice_check.setChecked(bool(v)),
            "avoid_polya": lambda v: self.polya_check.setChecked(bool(v)),
            "avoid_uorf": lambda v: self.uorf_check.setChecked(bool(v)),
            "inverted_stem": lambda v: self.inverted_spin.setValue(int(v)),  # type: ignore[arg-type]
            "inverted_loop": lambda v: self.inverted_loop_spin.setValue(int(v)),  # type: ignore[arg-type]
            "refine": lambda v: None,  # no refinement control in the app yet
        }
        unmapped = sorted(set(preset.overrides) - set(setters))
        for field, value in preset.overrides.items():
            setter = setters.get(field)
            if setter is not None:
                setter(value)
        message = f"Applied preset '{preset.label}'. Change anything below - your edits win."
        if unmapped:
            # Never let a preset silently set something with no visible control.
            message += f" (Not shown in this panel: {', '.join(unmapped)}.)"
        self.statusBar().showMessage(message)

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
        # The weight spin carries the strength AND the sign: a negative weight is
        # codon-pair deoptimization (attenuated-vaccine design), which the engine
        # supports and the app previously could not express.
        return self.weight_spins["cpb_weight"].value(), cds

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        """Shut running workers down before the window (and its threads) die.

        Qt aborts the process if a ``QThread`` is destroyed while still running,
        and closing the window destroys the threads parented to it. So on close
        we cancel what is cancelable (the frontier sweep polls a flag between
        points) and give each live thread a bounded chance to finish.

        The wait is bounded on purpose: a solve that is deep inside a single
        long step cannot be interrupted, and hanging the close forever would be
        worse than the abort it avoids. If the budget runs out we accept the
        close anyway -- the same behaviour as before -- so this strictly reduces
        the window in which a close can abort, rather than pretending to
        eliminate it.
        """
        if self._worker is not None:
            self._worker.cancel()
        for thread in (
            self._thread,
            self._cand_thread,
            self._cc_thread,
            self._lib_thread,
        ):
            if thread is not None and thread.isRunning():
                thread.quit()
                thread.wait(_CLOSE_WAIT_MS)
        super().closeEvent(event)

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
        if self._busy():
            return  # another engine flow is already in flight
        protein = self._prepare_protein()
        if protein is None:
            return
        enzymes = self._prepare_enzymes()
        if enzymes is None:
            return
        if not self._confirm_long_run(protein):
            return
        config, steps = self._build_config(enzymes)
        # Frozen here, not in _on_finished: the run is asynchronous and the user
        # can change the combos while it is in flight, so reading them on
        # completion would label the result with tables it was not built from.
        self._running_tables = (
            config.organism,
            config.reference_set or api.default_reference_set(config.organism),
        )

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
        if self._busy():
            return  # another engine flow is already in flight
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
            use_attested=self.splice_attested_check.isChecked(),
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
        if self._busy():
            return  # another engine flow is already in flight
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
        """Render a finished cross-check, unless it is about a different sequence.

        A cross-check report describes exactly one sequence, and it carries that
        sequence (:attr:`SpliceCrossCheck.dna`). If the delivered design changed
        while the report was in flight, showing it would attribute one sequence's
        predicted splice sites to another -- so it is discarded and said so,
        rather than rendered. Comparing the report's own DNA to the live delivered
        DNA needs no extra bookkeeping and cannot drift out of sync.
        """
        self._set_crosscheck_running(False)
        if not isinstance(report, api.SpliceCrossCheck):
            return
        delivered = self._delivered()
        if delivered is None or report.dna.upper() != delivered.dna.upper():
            self._reset_crosscheck()
            self.assp_banner.setText(
                "Cross-check <b>discarded</b>: the delivered sequence changed "
                "while it was running, and a report describes one sequence only. "
                "Run it again on the current design."
            )
            self.statusBar().showMessage(
                "ASSP cross-check discarded (the delivered sequence changed)."
            )
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
            f"{escape(str(error))}<br>"
            "<i>Advisory only; the delivered result is unaffected.</i>"
        )
        self.statusBar().showMessage("ASSP cross-check unavailable.")

    # ---- library flow -----------------------------------------------------

    def _start_library(self) -> None:
        """Validate inputs, then sample a library off-thread."""
        if self._busy():
            return  # another engine flow is already in flight
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

        Only one engine flow may run at a time, so every start control is gated on
        "nothing is running". Enablement is derived from the ``_*_running`` flags
        rather than from thread references, so a missed reference clear can never
        strand a control (the optimize-then-rank regression). ASSP additionally
        needs a delivered sequence to cross-check.

        The **menu actions are gated here too**, not just the buttons: they carry
        keyboard shortcuts, so gating only the buttons would leave Ctrl+R a live
        back door into starting a second engine flow mid-run.
        """
        busy = self._busy()
        can_crosscheck = not busy and self._delivered() is not None
        for widget in (self.optimize_btn, self.act_optimize):
            widget.setEnabled(not busy)
        for widget in (self.cancel_btn, self.act_cancel):
            widget.setEnabled(self._optimize_running)
        for widget in (self.rank_btn, self.act_rank):
            widget.setEnabled(not busy)
        for widget in (self.library_btn, self.act_library):
            widget.setEnabled(not busy)
        for widget in (self.assp_btn, self.act_assp):
            widget.setEnabled(can_crosscheck)

    def _busy(self) -> bool:
        """Whether any engine flow is currently running.

        The single source of truth for "only one flow at a time". Both the UI
        gate (:meth:`_update_run_buttons`) and each flow's own entry guard read
        it, so the invariant is enforced in the code path itself and not only by
        greying out a control.
        """
        return (
            self._optimize_running
            or self._cand_running
            or self._cc_running
            or self._lib_running
        )

    def _set_export_enabled(self, enabled: bool) -> None:
        """Enable/disable the delivered-result exports (button and menu together)."""
        for widget in (
            self.export_fasta_btn,
            self.export_json_btn,
            self.export_gb_btn,
            self.act_export_fasta,
            self.act_export_json,
            self.act_export_genbank,
        ):
            widget.setEnabled(enabled)

    def _set_library_export_enabled(self, enabled: bool) -> None:
        """Enable/disable the library export (button and menu together)."""
        for widget in (self.export_library_btn, self.act_export_library):
            widget.setEnabled(enabled)

    # ---- slots ------------------------------------------------------------

    @QtCore.Slot(int, str)
    def _on_progress(self, value: int, label: str) -> None:
        """Show worker progress in the status bar."""
        self.statusBar().showMessage(f"{label} ({value}%)")

    @QtCore.Slot(object)
    def _on_finished(self, result: api.FrontierResult) -> None:
        """Populate the results panel from a finished frontier optimization."""
        self._last = result
        self._delivered_tables = self._running_tables
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
        self._set_export_enabled(False)
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

        self._set_export_enabled(True)
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
        """Fill the metrics table from the delivered result's recomputed metrics.

        The fixed core rows are followed by rows derived from the run's **audit**,
        so what the GUI shows tracks what the engine actually reported. That matters
        for honesty, not tidiness: a non-local rule that could only be *partly*
        enforced says so on the CLI, and before this the GUI showed nothing at all --
        a user who set "Max repeat length" and got 19 unremovable residuals saw only
        an unexplained hard-violation count.
        """
        metrics = delivered.metrics
        cert = delivered.certificate
        rows: list[tuple[str, str]] = [
            ("CAI", f"{float(delivered.audit['cai']):.4f}"),
            (
                "CAI reference set",
                str(delivered.audit.get("codon_reference_set", "unknown")),
            ),
            ("GC %", f"{metrics.gc * 100.0:.2f}"),
            ("Length (nt)", str(metrics.length_nt)),
            ("Scored codons", str(delivered.audit["n_scored_codons"])),
            ("Hard violations", str(metrics.hard_violations)),
            ("Soft violations", str(metrics.soft_violations)),
            ("Optimality", cert.status.value),
            ("Solver", cert.solver),
        ]
        rows.extend(_audit_rows(delivered.audit))

        self.metrics_table.setRowCount(len(rows))
        for row, (name, value) in enumerate(rows):
            self.metrics_table.setItem(row, 0, QtWidgets.QTableWidgetItem(name))
            item = QtWidgets.QTableWidgetItem(value)
            # A partly-enforced rule is the one row a user must not skim past.
            if "could not" in value or "UNCALIBRATED" in value:
                # Same red the sequence viewer uses for a HARD span, so "partly
                # enforced" reads as the warning it is rather than as ordinary text.
                item.setForeground(
                    SequenceViewer._HARD_DARK if _is_dark() else SequenceViewer._HARD_LIGHT
                )
            # The full text is also the tooltip: an enforcement note is the one
            # row a user must be able to read, and it is the longest.
            item.setToolTip(value)
            self.metrics_table.setItem(row, 1, item)
        self.metrics_table.resizeColumnsToContents()
        # Let the value column take the remaining width rather than being clipped
        # -- truncating "partial - N could not be removed" would hide exactly the
        # disclosure the row exists to make.
        header = self.metrics_table.horizontalHeader()
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)

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
        """Plot the delivered sequence's per-site GC, CpG-density and splice tracks.

        The first two are honest reporting profiles (``api.tracks``): each point is
        a sliding-window statistic recomputed from the DNA, never a solver output.
        The splice track is different in kind -- it comes from a *model*, so the
        plot title carries that model's calibration state, and with the shipped PWM
        baseline it says UNCALIBRATED. All three are fraction-scaled (0-1) so they
        share a y-axis; %MinMax (a different scale) stays available via
        ``api.tracks`` / ``bt4 tracks``.
        """
        self.tracks_plot.clear()
        organism, reference_set = self._delivered_tables
        try:
            tracks = api.tracks(
                delivered.dna,
                organism or self.organism_combo.currentText(),
                reference_set=reference_set or None,
                nt_window=30,
                splice=True,
            )
        except ValueError:
            return
        gc = tracks.get("gc_fraction")
        cpg = tracks.get("cpg_density")
        splice = tracks.get("splice_site")
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
        if splice is not None and splice.values:
            self.tracks_plot.plot(
                list(range(len(splice.values))),
                list(splice.values),
                pen=pg.mkPen("#8a5ad0", width=1),
                name="splice site",
            )
        if tracks.splice_model:
            state = "calibrated" if tracks.splice_calibrated else "UNCALIBRATED"
            self.tracks_plot.setTitle(
                f"per-site tracks - splice: {tracks.splice_model} [{state}]"
            )
        else:
            self.tracks_plot.setTitle("")

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
        # This banner is RichText, and `reason` can carry a remote service's own
        # error text. Escape everything that came from outside before it is
        # interpolated: untrusted markup must never be able to rewrite (or hide)
        # the very labels that mark these numbers advisory and uncalibrated.
        parts = [
            f"<b>{escape(report.backend)}</b> &mdash; {', '.join(tags)}; advisory "
            "only, <b>not</b> part of the run manifest and never exported."
        ]
        if not report.available:
            parts.append(f"<b>Unavailable:</b> {escape(report.reason or '')}")
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
        self._set_library_export_enabled(False)

    @staticmethod
    def _cai_header(results: list[api.Result] | tuple[api.Result, ...]) -> str:
        """Return the CAI column header, naming the reference set behind it.

        The Design tab spells the reference set out on its own metrics row for
        the reason given there -- a CAI without it is a number with no question
        attached. These tables have one number per row and no room for a second
        column, so the label goes in the header instead. Read from each result's
        own audit, never from the control panel, so it cannot drift from the
        numbers beneath it.
        """
        labels = {str(r.audit.get("codon_reference_set", "")) for r in results}
        return f"CAI ({labels.pop()})" if len(labels) == 1 and any(labels | {""}) else "CAI"

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
        self.library_table.setHorizontalHeaderLabels(
            tuple(
                self._cai_header(members) if col == "CAI" else col
                for col in _LIBRARY_COLS
            )
        )
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
        self._set_library_export_enabled(bool(members))
        self.library_table.clearSelection()
        if members:
            self.library_table.selectRow(0)
        # Repaint the viewer explicitly rather than relying on selectRow to emit
        # itemSelectionChanged: repopulating the table in place leaves the old
        # selection intact, so re-selecting row 0 after a second draw emits
        # nothing and would strand the PREVIOUS draw's sequence on screen --
        # attributing one library's member to another.
        self._show_library_member()

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
        self.candidates_table.setHorizontalHeaderLabels(
            tuple(
                self._cai_header([c.result for c in candidates])
                if col == "CAI"
                else col
                for col in self._CANDIDATE_COLS
            )
        )
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
        calib = (
            "fidelity-attested (NOT a calibration claim)"
            if audit.all_calibrated
            else "UNCALIBRATED (advisory)"
        )
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
        if audit.all_calibrated:
            # `calibrated=True` here means the wrapper reproduces the published model,
            # not that its numbers are calibrated probabilities. Saying so is the whole
            # point of surfacing the opt-in in a GUI rather than only an env var.
            parts.append(
                "An attestation proves BT4's <b>wrapper</b> is faithful to the "
                "published model &mdash; not that its scores are calibrated "
                "probabilities for coding sequence. These models are measured weakest "
                "on exonic variants (median prAUC 0.419 vs 0.773 intronic), which is "
                "BT4's entire regime."
            )
        self.splice_banner.setText("<br>".join(parts))

    # ---- export -----------------------------------------------------------

    def _delivered(self) -> api.Result | None:
        """The currently delivered result, if any."""
        return self._last.delivered() if self._last is not None else None

    def _open_protein(self) -> None:
        """Load a protein FASTA into the input box (the design target in §6.6).

        The app was paste-only, so a user with a sequence in a file had to route it
        through a text editor. Multi-record files load the first record and say so,
        rather than silently designing for a sequence the user did not pick.
        """
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open protein FASTA", "", "FASTA (*.fasta *.fa *.faa);;All files (*)"
        )
        if not path:
            return
        try:
            records = list(api.read_fasta(path))
        except (OSError, ValueError) as exc:
            self._warn("Could not read FASTA", str(exc), "")
            return
        if not records:
            self._warn("Empty FASTA", f"{path} contained no records.", "")
            return
        header, sequence = records[0]
        self.protein_edit.setPlainText(sequence)
        if not self.jobname_edit.text().strip() and header:
            self.jobname_edit.setText(header.split()[0])
        note = f"Loaded {header or 'record 1'} ({len(sequence)} residues) from {path}"
        if len(records) > 1:
            note += f" - first of {len(records)} records; the rest were not loaded"
        self.statusBar().showMessage(note)

    def _validate_sequence(self) -> None:
        """Audit an existing coding sequence against the current design controls.

        ``api.validate`` was fully built and had no UI at all, so a user holding a
        CDS -- their own, or one a collaborator sent -- had to drop to the CLI to
        check it. The audit uses the same controls as a design run, and reports
        every rule they set (LOCAL and non-local alike).
        """
        dna, ok = QtWidgets.QInputDialog.getMultiLineText(
            self, "Validate a coding sequence", "Paste an ACGT coding sequence:", ""
        )
        if not ok:
            return
        cleaned = "".join(dna.split()).upper()
        if not cleaned:
            return
        enzymes = self._prepare_enzymes()
        if enzymes is None:  # an unknown enzyme name was already reported
            return
        config, _steps = self._build_config(enzymes)
        try:
            report = api.validate(cleaned, config)
        except ValueError as exc:
            self._warn("Could not validate", str(exc), "")
            return
        hard = report.metrics.hard_violations
        soft = report.metrics.soft_violations
        verdict = "FEASIBLE" if report.is_feasible else "NOT feasible"
        lines = [
            f"{verdict} under the current design controls.",
            f"{len(report.dna)} nt, GC {report.metrics.gc * 100:.1f}%",
            f"{hard} hard / {soft} soft violation(s).",
        ]
        detail = "\n".join(
            f"[{v.severity.name}] {v.constraint} @ {v.start}-{v.end}: {v.detail}"
            for v in report.violations[:100]
        )
        if len(report.violations) > 100:
            detail += f"\n... and {len(report.violations) - 100} more"
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Validation report")
        box.setText("\n".join(lines))
        if detail:
            box.setDetailedText(detail)
        box.exec()

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

    def _export_genbank(self) -> None:
        """Write the delivered result as an annotated GenBank map.

        This is the export that puts BT4's honesty where the user actually looks:
        every residual violation becomes a ``misc_feature`` span, so a defect the
        optimizer could not remove is visible in SnapGene / Benchling rather than
        only in the JSON audit. Uses the run's construct context when one was
        supplied, so flanks and junction findings are placed correctly.
        """
        delivered = self._delivered()
        if delivered is None:
            return
        header = self.jobname_edit.text().strip() or "bt4"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export GenBank", f"{header}.gb", "GenBank (*.gb *.gbk);;All files (*)"
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(
                api.write_genbank(
                    delivered, context=self._build_context(), locus=header
                )
            )
        self.statusBar().showMessage(f"Wrote annotated GenBank to {path}")


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
