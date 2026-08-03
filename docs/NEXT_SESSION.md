# BT4 — next-session build brief

A ready-to-use prompt/brief for the next building session. Paste it (or point the
session at this file) to resume work. **Read [`../CLAUDE.md`](../CLAUDE.md) first —
it is the constitution and it overrides anything here that has drifted.**

---

## Where BT4 is right now

Phases 0–1 complete; **Phase 2 essentially complete**; **Phase 3 groundwork
landed**. Shipped and green on `main`:

- **Honest exact-DP core** — codon trellis with true per-constraint context and a
  real optimality certificate; beam as an explicit knob.
- **Objectives:** CAI, tAI (8 organisms, real GtRNAdb data), GC-proximity, 5′
  ramp, CpG, %MinMax, and **codon-pair bias** (built from a user reference CDS).
  Returned as a multi-objective **Pareto frontier**.
- **Constraints:** homopolymer, **max-GC-run**, **max-repeat-length** (dispersed,
  RC-aware, GLOBAL/refinement-enforced), tandem/inverted repeats, forbidden motifs
  + **named presets**, restriction sites (IUPAC, RC), strong-Kozak internal-ATG,
  and **out-of-frame uORF** (GLOBAL/refinement-enforced, structural).
- **Budget backends:** OR-Tools CP-SAT and an honest Lagrangian/exact-budget DP.
- **Phase 3 groundwork:** `FoldingModel` (ViennaRNA + labeled baseline),
  `SplicePredictor` (labeled PWM baseline), the SA refinement engine (with a
  global-constraint gate, invariant #5), per-site tracks plotted in BT4 Studio.
- **`ExpressionPredictor` contract scaffolded** (`biomodels/expression/`) with a
  neutral, honestly-uncalibrated placeholder and a frontier-rerank hook that never
  steers delivery unless the predictor is calibrated.
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
3. **Rust trellis port** (Phase 1 perf) — port the DP inner loop to `bt4_native`
   (PyO3/maturin) with a numpy fallback; `abi3` wheels. Runtime is a first-class
   concern (§7); keep the perf regression test green.
4. **CpG / UpA whole-sequence count budgets** (Phase 2 remaining) — non-local
   count budgets via the budget/Lagrangian or refinement path.
5. **Packaged installers** (Phase 4) — PyInstaller/Briefcase for macOS/Windows/
   Linux; polish Studio theming/accessibility; optional external-validation report.
6. **Phase 5** — library/degenerate-design mode (sample the codon distribution,
   not one MFC target), more organisms with provenance, tissue-specific tables.

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

Pick item 1 or 3 (both are self-contained and high-value). If data/GPU aren't
available for the splice model, do the **Rust trellis port** (item 3) — it needs
no external data and directly serves the runtime goal. Confirm scope with the
maintainer before a large data-dependent effort.
