# BT4 — next-session build brief

A ready-to-use prompt/brief for the next building session. Paste it (or point the
session at this file) to resume work. **Read [`../CLAUDE.md`](../CLAUDE.md) first —
it is the constitution and it overrides anything here that has drifted.**

---

## Where BT4 is right now

Phases 0–2 complete; **Phase 3 groundwork landed**; **Phase 5 opened**. Shipped
and green on `main`:

- **Honest exact-DP core** — codon trellis with true per-constraint context and a
  real optimality certificate; beam as an explicit knob.
- **Objectives:** CAI, tAI (8 organisms, real GtRNAdb data), GC-proximity, 5′
  ramp, CpG, %MinMax, and **codon-pair bias** (built from a user reference CDS).
  Returned as a multi-objective **Pareto frontier**.
- **Constraints:** homopolymer, **max-GC-run**, **max-repeat-length** (dispersed,
  RC-aware, GLOBAL/refinement-enforced), tandem/inverted repeats, forbidden motifs
  + **named presets**, restriction sites (IUPAC, RC), strong-Kozak internal-ATG,
  and **out-of-frame uORF** (GLOBAL/refinement-enforced, structural).
- **Budget backends:** OR-Tools CP-SAT and an honest Lagrangian/exact-bucketed
  budget DP — now **context-aware**, so **CpG/UpA whole-sequence count budgets**
  (`dinuc_budget`/`dinuc_min`/`dinuc_max`) ship as exact, proven-optimal
  dinucleotide-count budgets (completed the last Phase 2 item).
- **Library / degenerate-design mode (Phase 5):** `api.library` / `bt4 library` —
  an honest deterministic codon-distribution sampler with a `SAMPLED` certificate
  (not an optimizer; local-constraint-respecting; no optimality/expression claim).
- **Phase 3 groundwork:** `FoldingModel` (ViennaRNA + labeled baseline),
  `SplicePredictor` (labeled PWM baseline), the SA refinement engine (with a
  global-constraint gate, invariant #5), per-site tracks plotted in BT4 Studio.
- **`ExpressionPredictor` contract scaffolded** (`biomodels/expression/`) with a
  neutral, honestly-uncalibrated placeholder and a frontier-rerank hook that never
  steers delivery unless the predictor is calibrated.
- **Native primitives:** `bt4_native` now also ships `max_gc_run` and
  `longest_repeat` (byte-identical Python fallbacks + equivalence tests);
  `max_gc_run` backs the GC-run `ok_suffix` veto.
- **Surfaces:** stable `bt4.api`, the `bt4` CLI, BT4 Studio (PySide6, tooltips on
  every control), an optional FastAPI service, content-hashed provenance manifests.

## What's left (priority order — see CLAUDE.md §9 for the authoritative list)

1. **Validated splice model** (Phase 3) — a SpliceAI/Pangolin-class per-nucleotide
   CNN trained on real GENCODE, Δsplicing objective, **held-out-chromosome gate**
   (PR-AUC/MCC/ECE), hash-pinned artifact out of git. Slots behind the existing
   `SplicePredictor` contract; until it passes its gate the PWM baseline stays
   `calibrated=False`. Needs GPU + data.
2. **Learned expression head** (Phase 4) — trained on real MPRA / ribosome-load
   data, hash-pinned, calibrated + uncertainty (conformal). Slots behind the
   scaffolded `ExpressionPredictor` (set `calibrated=True` only after the gate).
   **Blocker/watch-out:** the honest hard part is *matched-regime data* — most
   uORF/ribosome-load data is 5′UTR, but the tool controls the CDS (see the
   uORF-calibration analysis). Do **not** relabel a hand-weighted composite as
   "calibrated" (§10.5/§10.6).
3. **Full Rust trellis port** (Phase 1 perf) — the DP inner loop still runs in
   pure Python. The `bt4_native` primitive set has grown (`gc_count`,
   `max_homopolymer_run`, `reverse_complement`, `max_gc_run`, `longest_repeat`),
   but porting the **DP inner loop** to `bt4_native` (PyO3/maturin, numpy
   fallback, `abi3` wheels) remains. Runtime is a first-class concern (§7); keep
   the perf regression test green. *(Lesson from this wave: a whole-sequence
   O(n²) native call on a per-move hot path is a pessimization in the
   no-extension case — measure before wiring, and keep the pure-Python path
   fast.)*
4. **Packaged installers** (Phase 4) — PyInstaller/Briefcase for macOS/Windows/
   Linux; polish Studio theming/accessibility; optional external-validation report.
5. **Phase 5 (continued)** — library/degenerate-design mode has **landed**;
   remaining Phase 5 is **more organisms with authoritative provenance**,
   restriction-enzyme catalog growth, and tissue/condition-specific tables. A
   natural follow-up: add a library-mode control to BT4 Studio (kept out of this
   wave to avoid touching the app layer).

## Working agreements (do not violate)

- **Honesty is structural.** Never present an unenforced constraint, an
  unvalidated number, or a heuristic result as if it were real. New model →
  `calibrated=False` until it passes a held-out gate. Never fabricate a data table
  — refuse and say why (as tAI / codon-pair bias do).
- **Adding a constraint/objective/model = a new file + a registry/export entry +
  its honesty property test** (`ok_suffix⇔validate` / `delta==score` / calibration
  gate). Never an engine edit. Keep the strict layering (import-linter enforces it).
- **Keep CLAUDE.md current in the same change** (§10.11). Update README when a
  user-facing surface changes.
- **Single-trunk + CI.** Branch, open a PR, merge on green. The full local gate:
  ```
  python -m ruff check src tests scripts
  python -m mypy <changed files>            # whole-package mypy shows env-only
                                            # cpsat/_accel noise; CI's dep-free
                                            # quality job is the source of truth
  lint-imports
  QT_QPA_PLATFORM=offscreen python -m pytest tests/ -p no:cacheprovider
  ```
- **Sandbox limits (must be done by a human):** deleting remote branches, pushing
  git tags, and cutting releases are blocked here (HTTP 403). Leave those for the
  maintainer; don't work around them.

## Suggested first move

The self-contained, no-external-data option is the **full Rust trellis port**
(item 3) — it directly serves the runtime goal and needs no GPU or data. The
splice model (item 1) and expression head (item 2) are the highest-value items
but both need **real held-out data + GPU**; confirm scope and data access with
the maintainer before starting either, and never ship an uncalibrated model as
if it were validated (§10.6). Packaged installers (item 4) can be advanced in
the sandbox up to the point where signing/tag-pushing/release-cutting is needed
(those are human-only here, HTTP 403).
