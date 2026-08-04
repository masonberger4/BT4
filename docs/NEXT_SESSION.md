# BT4 — next-session build brief

A ready-to-use prompt/brief for the next building session. Paste it (or point the
session at this file) to resume work. **Read [`../CLAUDE.md`](../CLAUDE.md) first —
it is the constitution and it overrides anything here that has drifted.**

---

## Where BT4 is right now

Phases 0–2 complete; **Phase 3 groundwork landed**; **Phase 5 opened**. All
merged and green on `main` (the two wrapped splice CNN backends from this session,
PRs #33 and #34, are now merged):

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
  `SplicePredictor` (labeled PWM baseline **plus both wrapped CNN backends, now
  merged to `main`** — Pangolin (GPL-3.0, PR #33) and SpliceAI (PolyForm Strict
  code + CC BY-NC weights, PR #34), lazily imported, hash-pinned,
  `calibrated=False` until their fidelity gates, with a two-backend agreement
  harness), the SA refinement engine (with a global-constraint gate, invariant
  #5), per-site tracks plotted in BT4 Studio.
- **`ExpressionPredictor` contract scaffolded** (`biomodels/expression/`) with a
  neutral, honestly-uncalibrated placeholder and a frontier-rerank hook that never
  steers delivery unless the predictor is calibrated.
- **Native primitives:** `bt4_native` now also ships `max_gc_run` and
  `longest_repeat` (byte-identical Python fallbacks + equivalence tests);
  `max_gc_run` backs the GC-run `ok_suffix` veto.
- **Surfaces:** stable `bt4.api`, the `bt4` CLI, BT4 Studio (PySide6, tooltips on
  every control, and a sequence viewer that renders **inline violation
  annotations** — each `Violation` highlighted over its `[start, end)` span,
  HARD red / SOFT amber, with a hover tooltip and legend, so residual GLOBAL
  violations show *where* they occur), an optional FastAPI service, content-hashed
  provenance manifests.

## What's left (see CLAUDE.md §9 for the authoritative list)

> **Live priority (this list is kept in original numbering for continuity):** item
> 1 (wrap SpliceAI + Pangolin) is now **essentially done** — both adapters and the
> agreement harness landed this session (PRs #33/#34); only a maintainer tail
> remains (recording the fidelity gates + the ASSP cross-check). The genuinely
> next self-contained work is **item 3 (Rust trellis port)** and **item 4
> (block/tempering refinement)** — see "Suggested first move" below.

1. **Wrap published SpliceAI + Pangolin as calibrated splice backends** (Phase 3)
   — **Decision: no self-training.** Wrap the already-validated **Pangolin** and
   **SpliceAI** as *inference-only* backends behind the **existing**
   `SplicePredictor` contract; the Δsplicing framing and top-k/log-odds pooling are
   already in `biomodels/splice/base.py`. **✅ Both adapters + the agreement harness
   have landed:** `PangolinSplicePredictor` (PR #33) and `SpliceAiSplicePredictor`
   (PR #34) — both merged to `main` — plus `backend_agreement` +
   `scripts/compare_splice_backends.py`, with the PWM baseline still the
   `calibrated=False` default.
   - **License corrections (both were wrong in the earlier brief).** Pangolin is
     **GPL-3.0** (not MIT). SpliceAI is stricter still: **code = PolyForm Strict
     1.0.0, weights = CC BY-NC 4.0** (noncommercial) — the `setup.py` "GPLv3"
     string is contradicted by the authoritative LICENSE files. Both follow BT4's
     GPL-ViennaRNA pattern: **lazily import the user's own installed package +
     weights, bundle neither code nor weights** — BT4 stays MIT. Install each
     yourself (github.com/tkzeng/Pangolin, github.com/Illumina/SpliceAI; each
     ships its weights). SpliceAI's CC BY-NC weights make that backend
     noncommercial-only.
   - **What landed (both):** lazy heavy-dep imports (so `import bt4` stays light),
     weights **SHA-256 hash-pinned** (published digests, re-verified by hand) and
     checked *before* load, out-of-loop scoring, and per-adapter fidelity gates
     (`verify_pangolin_fidelity` / `verify_spliceai_fidelity`). Pangolin emits one
     combined `P(splice)` → `SpliceResult.donor` (acceptor zero); SpliceAI's 3-way
     softmax maps cleanly to `donor` + `acceptor` (both populated). Each was
     verified to reproduce its upstream model **bit-for-bit** against the real
     weights, but ships `calibrated=False` (no reference panel bundled). The
     agreement harness needed no change — with both CNNs installed it compares two
     real, independently-trained models.
   - **What remains for item 1:** (a) capturing reference panels and **recording
     the fidelity gates** to promote either backend to `calibrated=True` (and
     having `default()` prefer a calibrated one) — a maintainer step needing the
     licensed weights, not fabricated here; (b) the opt-in **ASSP** cross-check
     with offline fixtures. Keep the honest scope note (predicts splice-*site
     presence*; lower Δ = lower *predicted* cryptic-splice risk, a strong prior,
     not validated expression gain). **Needs no GPU/training data.**
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
4. **Refinement reach: block/segment moves + parallel tempering** (Phase 3,
   self-contained, no data) — the SA engine (`optimize/anneal_refine.py`) only
   proposes **single-codon** moves and its global gate strictly forbids any
   increase in the hard-violation count, so it can leave a dispersed
   **max-repeat / uORF** in place when clearing it needs a *coordinated
   multi-codon* change or a temporary count increase (it can already traverse
   count-`==` plateaus, so pure lateral multi-position effects are reachable). Add
   **block/segment moves** (propose synonymous swaps at several positions at once)
   and a **parallel-tempering** schedule (hot replicas may accept uphill moves and
   swap into the cold chain) so refinement can cross those barriers **without
   weakening invariant #5 on the delivered result**. Keep it honest: a repeat
   pinned to synonymously-immovable bases (Met `ATG` / Trp `TGG`, or a
   base-locked degenerate position) is a genuine feasibility floor — still report
   it as a residual, don't pretend a block move can remove it. (An alternative
   worth weighing: encode max-repeat as a global constraint in an ILP/CP-SAT
   solve for an exact answer — currently it lives only in the refinement layer.)
   See CLAUDE.md §7 and §9 Phase 3.
5. **Packaged installers** (Phase 4) — PyInstaller/Briefcase for macOS/Windows/
   Linux; polish Studio theming/accessibility; optional external-validation report.
6. **Phase 5 (continued)** — library/degenerate-design mode has **landed**;
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

**Item 1's wrapped splice backends are both done** — Pangolin (PR #33) and
SpliceAI (PR #34), each hash-pinned and `calibrated=False` until its fidelity
gate, with the two-backend agreement harness. The best next moves, all
self-contained (no GPU, no external data):

- **Full Rust trellis port** (item 3) — the DP inner loop still runs in pure
  Python. A grounded plan exists (position-independent-regime `trellis_solve`
  primitive + byte-identical Python twin + inline-DP oracle).
- **Block/segment + parallel-tempering refinement moves** (item 4) — a grounded
  plan exists (opt-in replica/temper/block kwargs on `anneal_refine`, strict
  global gate on every replica so invariant #5 holds).
- **Finish item 1's tail:** record the fidelity gates to promote a backend to
  `calibrated=True` — a maintainer step needing the licensed weights + a captured
  panel; don't fabricate one. The opt-in **ASSP** cross-check also remains.

The **learned expression head** (item 2) is the one item that still needs real
matched-regime data — defer it and never ship an uncalibrated composite as
validated (§10.5/§10.6). Packaged installers (item 5) can be advanced up to the
point where signing/tag-pushing/release-cutting is needed (human-only here, HTTP
403).

---

## Session archive marker (for continuity)

**Last building session delivered (all merged to `main`, green):**
- Phase 2 **completed** — CpG/UpA whole-sequence count budgets (PR #26).
- Phase 5 **opened** — library / degenerate-design mode (PR #25).
- Phase 1 perf — native `max_gc_run` + `longest_repeat` primitives (PR #27),
  including a fix for an O(n²) `longest_repeat` fast-path that regressed the
  pure-Python `MaxRepeatConstraint.validate` hot path.
- Docs — status sync (PR #28), the single-codon-SA refinement limitation note
  (PR #30), and this **splice decision** (wrap SpliceAI/Pangolin, no self-train).

**This session delivered (both PRs merged to `main`, green): both wrapped splice
adapters.** `PangolinSplicePredictor` (**PR #33**, merged; wraps the user's
installed GPL-3.0 Pangolin) and `SpliceAiSplicePredictor` (**PR #34**, stacked on
#33 then retargeted to `main` and merged; wraps the user's installed SpliceAI —
code PolyForm Strict 1.0.0, weights CC BY-NC 4.0) — neither bundled, hash-pinned
weights verified before load, `calibrated=False` until per-adapter fidelity gates
— plus the `backend_agreement` two-backend harness +
`scripts/compare_splice_backends.py`, the `bt4[splice-pangolin]` /
`bt4[splice-spliceai]` extras, and the **license corrections** (Pangolin is
GPL-3.0 not MIT; SpliceAI is PolyForm+CC BY-NC not GPL). Both were CI-green before
merge; each adapter was verified to reproduce its upstream model **bit-for-bit**
against the real weights. (This archive doc-sync itself is a small follow-up PR.)

**Also produced this session (design plans, not code):** grounded,
execution-ready implementation plans for the **Rust trellis port** (item 3) and
the **block/segment + parallel-tempering refinement moves** (item 4), from a
parallel design fan-out — so the next session can pick either up quickly (see
"Suggested first move").

**Deliberately NOT done, and why:** no bespoke splice CNN and no expression head
were trained — both would need real held-out data (+ GPU for a from-scratch CNN),
and the constitution forbids shipping an uncalibrated model as validated. Both
CNN backends ship `calibrated=False` because **no reference panel is bundled**
(capturing one needs the licensed weights and reproduces licensed outputs) — the
promotion is a maintainer step, not fabricated here. The data-blocked
**expression head** (item 2) remains. **To resume: read `../CLAUDE.md` (§6 Splice,
§9 Phase 3) then this brief; pick up at the Rust trellis port (item 3), the
block/tempering refinement moves (item 4), or recording the splice fidelity
gates.**
