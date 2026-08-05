# BT4 — next-session build brief

A ready-to-use brief for the next building session. Paste it (or point the session
at this file) to resume work. **Read [`../CLAUDE.md`](../CLAUDE.md) first — it is
the constitution and it overrides anything here that has drifted.**

---

## Where BT4 is right now

Phases 0–2 complete; **Phase 3 groundwork landed**; **Phase 4 in progress**;
**Phase 5 opened**. All merged and green on `main`.

- **Honest exact-DP core** — codon trellis with true per-constraint context and a
  real optimality certificate; beam as an explicit knob. The **full DP inner loop
  is now ported to Rust** (`bt4_native.trellis_solve`, regime-gated with a
  byte-identical pure-Python twin + equivalence test; the win is the amortized
  Pareto frontier, ~2.7–5.5× faster, byte-identical output/certificates).
- **Objectives:** CAI, tAI (8 organisms, real GtRNAdb data), GC-proximity, 5′
  ramp, CpG, %MinMax, and **codon-pair bias** (built from a user reference CDS).
  Returned as a multi-objective **Pareto frontier**.
- **Constraints:** homopolymer, **max-GC-run**, **max-repeat-length** (dispersed,
  RC-aware, GLOBAL/refinement-enforced), tandem/inverted repeats, forbidden motifs
  + **named presets**, restriction sites (IUPAC, RC), strong-Kozak internal-ATG,
  and **out-of-frame uORF** (GLOBAL/refinement-enforced, structural).
- **Budget backends:** OR-Tools CP-SAT and an honest Lagrangian/exact-bucketed
  budget DP — context-aware, so **CpG/UpA whole-sequence count budgets**
  (`dinuc_budget`/`dinuc_min`/`dinuc_max`) ship as exact, proven-optimal
  dinucleotide-count budgets (completed the last Phase 2 item).
- **Refinement:** the incremental SA engine (`optimize/anneal_refine.py`) now has
  **block/segment moves + parallel tempering** (opt-in `block_size`/`block_prob`
  and `replicas`/`temps`/`swap_every`), so it can cross barriers a single-codon
  chain cannot — **without weakening invariant #5** on the delivered result (every
  replica gates against its own whole-sequence hard-violation count; all knobs
  default off and reproduce the prior single-chain trajectory byte-for-byte). A
  repeat pinned to synonymously-immovable bases is an honestly-disclosed
  feasibility floor (`max_repeat_residual`, enforcement `"partial"`), not a defect.
- **Splice (Phase 3):** `SplicePredictor` with a labeled PWM baseline **plus both
  wrapped CNN backends** — `PangolinSplicePredictor` (GPL-3.0) and
  `SpliceAiSplicePredictor` (PolyForm Strict code + CC BY-NC weights) — each lazily
  imported from the user's own install, hash-pinned, reproducing upstream
  bit-for-bit, `calibrated=False` until its fidelity gate, with a two-backend
  agreement harness and a license-clean fidelity-attestation layer.
- **Expression (Phase 4):** `ExpressionPredictor` contract + a neutral
  `NullExpressionModel` placeholder + a frontier-rerank hook that never steers
  delivery unless calibrated, **plus the wrapped `RiboNNExpressionModel`**
  (`biomodels/expression/ribonn.py`, Sanofi non-commercial, driven from the user's
  own checkout via `$BT4_RIBONN_DIR`, hash-pinned against a bundled 180-entry
  manifest, CLR-residual TE units, `delta_logte` the CDS-attributable signal). Ships
  `calibrated=False`; `default()` still returns the placeholder. The model-agnostic
  acceptance gate (`biomodels/expression/gate.py`, `verify_expression_gate`) is the
  promotion path. **`RiboNNExpressionModel` has now had its first real end-to-end
  runs against the licensed weights** (a maintainer's machine — non-commercial
  weights, never bundled or CI-run), which validated the adapter and fixed two
  live-only integration bugs (see the archive marker below).
- **Native primitives:** `bt4_native` ships `gc_count`, `max_homopolymer_run`,
  `reverse_complement`, `max_gc_run`, `longest_repeat`, and `trellis_solve` — each
  with a byte-identical Python fallback + equivalence test.
- **Surfaces:** stable `bt4.api`; the `bt4` CLI; BT4 Studio (PySide6, tooltips on
  every control, per-site risk tracks, and a sequence viewer with inline
  HARD/SOFT violation annotations); an optional FastAPI service; content-hashed
  provenance manifests; per-site reporting tracks (`api.tracks` / `bt4 tracks`).
- **Packaging:** now supports **Python 3.10+** (was 3.11+), so the RiboNN backend
  installs into the same environment as its pinned `torch==1.13.1` stack.

## What's left (see CLAUDE.md §9 for the authoritative list)

The two big self-contained engine items from the last brief — the **Rust trellis
port** and **block/tempering refinement** — have both **landed**. What remains is
mostly data-gated or human-only.

1. **Promote the splice CNNs to `calibrated=True`** (Phase 3 tail, human-only) —
   capture reference panels and **record the fidelity gates**
   (`verify_pangolin_fidelity` / `verify_spliceai_fidelity`, then a committed
   `FidelityAttestation` via the attestation layer), so `default()` can prefer a
   calibrated backend. Needs the licensed weights — a maintainer step, not
   fabricated. Also the opt-in **ASSP** cross-check with offline fixtures.
2. **Promote RiboNN to `calibrated=True`** (Phase 4 tail, data-gated, human-only) —
   assemble a **license-clean, regime-matched CDS-variant TE panel** and run
   `verify_expression_gate` (Spearman + split-conformal coverage on a group-disjoint
   split). Reproducing RiboNN faithfully is *not* calibration for BT4's CDS-variant
   regime (its own ablation puts only ~31% of per-nt signal in the CDS). Do **not**
   relabel a hand-weighted composite as "calibrated" (§10.5/§10.6).
   - **RiboNN perf/UX follow-up (self-contained, no data):** scoring is dominated by
     fixed per-call overhead (weight hashing + model load + Windows DataLoader
     worker-process spawn), and it is mostly *per RiboNN invocation*, not per
     sequence. `score_sequence` runs one CDS per call and `delta_logte` two, so
     scoring a whole frontier one-at-a-time pays that cost N times. Worth adding: a
     **public batch scoring method** (amortize a whole frontier in one call) and a
     `num_workers=0` path (usually faster on Windows for small inference jobs).
3. **Packaged installers** (Phase 4) — PyInstaller/Briefcase for macOS/Windows/
   Linux; polish Studio theming/accessibility; optional external-validation report.
   Advance up to the point where signing / tag-pushing / release-cutting is needed
   (human-only here — HTTP 403 in the sandbox).
4. **Phase 5 (continued)** — library/degenerate-design mode has landed; remaining is
   **more organisms with authoritative provenance**, restriction-enzyme catalog
   growth, and tissue/condition-specific tables. A natural follow-up: a library-mode
   control in BT4 Studio.

## Working agreements (do not violate)

- **Honesty is structural.** Never present an unenforced constraint, an unvalidated
  number, or a heuristic result as if it were real. New model → `calibrated=False`
  until it passes a held-out / fidelity gate. Never fabricate a data table — refuse
  and say why (as tAI / codon-pair bias do).
- **Adding a constraint/objective/model = a new file + a registry/export entry + its
  honesty property test** (`ok_suffix⇔validate` / `delta==score` / calibration
  gate). Never an engine edit. Keep the strict layering (import-linter enforces it).
- **Keep CLAUDE.md current in the same change** (§10.11). Update README when a
  user-facing surface changes.
- **Single-trunk + CI.** Branch, open a PR, merge on green. The full local gate:
  ```
  python -m ruff check src tests scripts
  python -m mypy                            # whole-package; CI's dep-free quality
                                            # job is the source of truth
  lint-imports
  QT_QPA_PLATFORM=offscreen python -m pytest tests/ -p no:cacheprovider
  ```
- **Sandbox limits (must be done by a human):** deleting remote branches, pushing
  git tags, cutting releases, and anything needing the licensed splice/expression
  weights are blocked here. Leave those for the maintainer; don't work around them.

## Suggested first move

The **expression/splice design flow is spec'd** in
[`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md) and its build
order is fixed. Start there:

- **Batched RiboNN scoring** (design step 1) — a `score_many` / `delta_logte_many`
  on the RiboNN adapter that runs the whole candidate set in one RiboNN invocation
  (amortizes the large fixed per-call overhead), plus an optional `num_workers=0`
  path. Pure engineering, no data, testable without the GUI. **The recommended first
  PR.**
- **Finish the calibration tails** — record the splice fidelity gates and run the
  expression acceptance gate. Both need licensed weights / matched-regime data and
  are human-only; don't fabricate a panel.
- **Phase 5 breadth** — more organisms with authoritative provenance is
  self-contained and always welcome (a good non-expression alternative).

The learned-expression calibration (item 2) is the one item gated on real
matched-regime data — never ship an uncalibrated composite as validated
(§10.5/§10.6). Packaged installers (item 3) can be advanced up to signing/release
(human-only here).

---

## Session archive marker (for continuity)

**Prior building session delivered (all merged to `main`, green):** both wrapped
splice adapters — `PangolinSplicePredictor` (PR #33) and `SpliceAiSplicePredictor`
(PR #34) — plus the two-backend agreement harness, and design plans for the Rust
trellis port and the block/tempering refinement moves. (Those two plans have since
been implemented and merged — see the constitution.)

**This session delivered (all merged to `main`, green):**
- **Python 3.10 support** (PR #42) — `requires-python >=3.10`, 3.10 classifier,
  ruff/mypy target 3.10, CI quality matrix adds 3.10. Pure compatibility widening
  (the core uses no 3.11-only features). It unblocks installing BT4 into the same
  environment as the RiboNN backend, whose pinned `torch==1.13.1` ships only CPython
  ≤3.10 wheels.
- **RiboNN adapter fixes** (PR #43) — the **first real end-to-end runs against the
  licensed RiboNN weights** (maintainer's Windows machine) surfaced two integration
  bugs that only appear once the live forward pass runs, both fixed:
  1. RiboNN returns its ensemble as **one row per cross-validation model**, so the
     per-input realignment now **groups by `tx_id` and averages** (mean over cell
     types *and* the ensemble) via the tested helper `_reduce_te_by_tx_id` — a plain
     `set_index` left duplicate labels and `float(Series)` raised `TypeError`.
  2. Scoring now **requires non-empty `utr5`/`utr3`** and refuses empty ones up front
     (RiboNN's loader reads an all-empty UTR column as `NaN` and its `.str`
     preprocessing crashes; the UTRs carry most of RiboNN's signal anyway).
  Both are property-tested against a synthetic RiboNN output table; `calibrated`
  stays `False`.
- **Docs sync** (PR #44 — this brief + CLAUDE.md §7/§9 + README) — recorded the
  Python-version change and the RiboNN validation.
- **First real RiboNN score through BT4** — `RiboNNExpressionModel(species="human")`
  scored a sequence end-to-end on the maintainer's machine (≈1.624 CLR-residual TE,
  `calibrated=False`), validating the whole path against the licensed weights.
- **Honest competitive positioning** (PR #45 — `docs/COMPARISON.md`, linked from the
  README) — a sourced review of BT4 vs IDT/Twist/GeneArt/GenScript/ATUM/DNA Chisel
  from a three-agent research sweep. Verdict: BT4 is near-unique on *rigor* (exact
  multi-objective optimization + certificates, byte-reproducible provenance,
  validated ML with honest calibration), but **not** yet differentiated on "expresses
  better" — CAI barely predicts expression (Kudla 2009, Welch 2009) and the CDS is
  only ~31% of the per-nt signal (RiboNN), and BT4's splice/expression models are
  still `calibrated=False`.
- **Expression/splice design of record** (PR #46 —
  `docs/DESIGN_expression_splice_flow.md`) — the agreed A→D pipeline for wiring
  cryptic-splice screening and RiboNN ranking into the design flow and BT4 Studio,
  spec'd before any code. **Locked decisions:** the RiboNN-ranked library is sourced
  from **frontier + repeat-fix variants**; the splice CNN **localize-and-flags now**
  (`calibrated=False`) and gains an **auto-edit** path only once a backend passes its
  fidelity gate. Heavy CNNs stay out-of-loop; the interacting constraints (repeats +
  uORF + splice-motifs + forbidden) are resolved in one **jointly-gated** refinement
  (invariant #5), not a linear chain.

**Environment notes for the RiboNN backend (learned on real hardware this session):**
RiboNN's own stack needs `numpy<2` (torch 1.13.1 ABI), `setuptools<81` (its older
`pytorch_lightning` calls `pkg_resources`), and the Zenodo `weights.zip` extracted to
a directory literally named `models/` under `$BT4_RIBONN_DIR` (so RiboNN's hard-coded
`models/<species>/<run_id>/state_dict.pth` path resolves without a Windows symlink).

**Deliberately NOT done, and why:** no splice/expression model was promoted to
`calibrated=True` — that needs the licensed weights + captured/matched-regime panels
(human-only), and the constitution forbids shipping an uncalibrated model as
validated. No engine code was written for the expression/splice flow — it was
**spec'd first** (`DESIGN_expression_splice_flow.md`) by explicit choice, so the
build has a fixed target.

**To resume — the next build is the expression/splice flow, and its order is
already fixed in [`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md)
§"Implementation phasing":**
1. **Batched RiboNN scoring** — a `score_many` / `delta_logte_many` on the RiboNN
   adapter that runs the whole candidate set in **one** RiboNN invocation (amortizes
   the large fixed per-call overhead), plus an optional `num_workers=0` path.
   Independently useful, testable without the GUI. **← start here.**
2. Strong splice-consensus donor/acceptor **motif constraint** (LOCAL; new file +
   registry entry + `ok_suffix⇔validate` test; honestly a heuristic, not a CNN).
3. **Candidate-set assembly + rerank** over frontier + repeat-fix library
   (`rerank_by_expression` applied across the set, calibrated-gated selection).
4. **Splice CNN localize-and-flag** audit (batched SpliceAI+Pangolin over the set).
5. **BT4 Studio UI** — UTR fields, the two toggles, the annotated frontier + ranked
   table with uncalibrated badges, on the background thread.
6. **(Gated, human-data)** splice **auto-edit** + RiboNN **auto-select**, each
   unlocked only when its backend passes its fidelity/acceptance gate.

Read [`../CLAUDE.md`](../CLAUDE.md) (§6, §9 Phases 3–4, §10.6, invariants #5/#7),
then [`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md) and
[`COMPARISON.md`](COMPARISON.md), then this brief. Phase 5 organism breadth remains
a good self-contained alternative if you want a non-expression PR.
