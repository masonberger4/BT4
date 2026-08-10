# Design: expression-aware, splice-screened design flow

> **Status: IMPLEMENTED through step 5 (only step 6 remains, calibration-gated).**
> This document is the design-of-record for wiring cryptic-splice screening and
> RiboNN expression ranking into BT4's design flow (and BT4 Studio). The
> "Implementation phasing" section below tracks which steps have shipped: steps
> 1–5 (batched RiboNN scoring, the splice-consensus motif constraint, candidate
> assembly + rerank, the localize-and-flag splice audit, and the BT4 Studio UI)
> are **landed**; step 6 (auto-edit / auto-select) is deferred until a backend is
> `calibrated=True`. For live status read
> [`NEXT_SESSION.md`](NEXT_SESSION.md); [`../CLAUDE.md`](../CLAUDE.md) is the
> authoritative constitution — keep this in sync (§10.11).

## Purpose

Produce a coding sequence that is (a) exactly optimized on BT4's cheap additive
objectives, (b) free of forbidden motifs and *strong* splice-consensus motifs,
(c) audited by real splice CNNs for residual cryptic sites, and (d) ranked for
predicted translation efficiency by RiboNN — **efficiently and honestly**, without
ever putting a heavy CNN in the optimizer's inner loop or presenting an
uncalibrated model's choice as validated.

## Governing constraints (why the flow is shaped this way)

1. **Cheap terms in the loop, heavy models at the end.** The exact DP / SA
   refinement makes O(context) incremental moves, thousands of them. RiboNN and
   SpliceAI/Pangolin are ~10 kb-context CNNs with large fixed per-call overhead —
   structurally out-of-loop (CLAUDE.md §6). Their outputs are also **non-local and
   non-additive**, so they cannot be a DP `delta()` or a cheap SA `delta_score`.
2. **Interacting constraints must be resolved jointly, not in a linear chain.**
   Fixing a repeat can re-introduce a splice or forbidden motif and vice versa. A
   strict "avoid → fix splice → fix repeats" ordering plays whack-a-mole. BT4's
   answer is one **jointly-gated refinement** that counts *all* hard violations
   and never lets any rise (invariant #5), not sequential passes.
3. **Motif-avoidance ≠ CNN-risk-reduction.** The bare donor `GT` / acceptor `AG`
   dinucleotides are ubiquitous and cannot be forbidden; only the longer
   *consensus* (e.g. a `GTRAGT`-type donor, a polypyrimidine + `AG` acceptor) can.
   Real splicing depends on context far beyond the core motif, which is exactly
   why the CNNs exist. So the in-loop motif constraint reduces *obvious* risk; the
   CNN audit does the real work.
4. **Calibration honesty (§10.6).** SpliceAI/Pangolin and RiboNN currently ship
   `calibrated=False`. Until each passes its gate, its output **advises**; the
   human confirms delivery. The in-loop splice-*motif* constraint is a structural
   heuristic (always-on, honestly labeled), not a calibrated-risk claim.
5. **Determinism (invariant #7).** Every candidate — including the repeat-fix
   library — is generated from a seed and is byte-reproducible.

## The pipeline

```
A. FRONTIER — exact DP, in-loop                                    [seconds]
   Additive objectives (CAI / tAI / codon-pair / 5' ramp) + LOCAL constraints:
   forbidden motifs, restriction sites, GC-run, homopolymer, strong-Kozak,
   + NEW strong splice-consensus donor/acceptor motif constraints.
   → a Pareto frontier whose points are already free of forbidden and
     strong-splice-consensus motifs.

B. JOINT REFINEMENT — out-of-loop, gated                          [fast]
   For each carried candidate, resolve the GLOBAL rules TOGETHER — dispersed
   repeats + out-of-frame uORF + splice-consensus motifs — in one
   `anneal_refine` pass whose global gate never lets ANY hard-violation count
   rise (invariant #5). Repeat resolution is non-unique, so a DETERMINISTIC
   enumeration of break choices / seeds yields a small LIBRARY of feasible
   variants per candidate.

C. SPLICE CNN AUDIT — out-of-loop, batched, once                  [minutes]
   Run SpliceAI + Pangolin over the candidate set ONCE (batched) to LOCALIZE
   residual cryptic sites the motif constraint missed. Backend agreement =
   confidence. Present each site with its span, score, and which backends agree.
   (See "Localize-and-flag vs auto-edit" below.)

D. RiboNN RANK — out-of-loop, batched, once                       [minutes]
   Batch-score the whole candidate set (frontier points + repeat-fix library)
   with RiboNN `delta_logte` vs the reference (UTRs held fixed). Present each
   candidate annotated with ΔTE and its splice flags. The user selects the
   delivered sequence.
```

### Candidate set (locked decision): **frontier + repeat-variants**

The library RiboNN ranks is sourced from **both** the Pareto frontier (CAI/GC/tAI
trade-off diversity — more TE-relevant and more discriminable) **and** the
repeat-fix variants from Stage B (structural diversity). Repeat-fix variants alone
are near-clones and would give RiboNN little to discriminate; the frontier widens
the spread so the ranking is meaningful. The combined set is deduplicated and
capped at N (each member is a heavy-CNN input).

### Splice CNN role (locked decision): **localize-and-flag now → auto-edit once calibrated**

- **Now (SpliceAI/Pangolin `calibrated=False`):** Stage C **localizes and flags**
  residual cryptic sites (span, score, backend agreement) but does **not** edit
  the sequence. The exported sequence is exactly what the user selected; splice
  findings are shown as annotations (mirroring BT4 Studio's existing inline
  violation highlighting). This respects §10.6 — an uncalibrated model advises,
  it does not silently rewrite the deliverable.
- **Later (a backend passes its fidelity gate → `calibrated=True`):** the *same*
  Stage C gains an **auto-edit** path — a targeted synonymous edit at each flagged
  locus that re-enters the **Stage-B jointly-gated refinement** (so removing a
  splice site can never raise the repeat/uORF/forbidden count), followed by one
  confirming CNN pass. The auto-edit is gated on `calibrated=True` exactly as
  `rerank_by_expression` gates auto-selection on the expression head being
  calibrated. No rework: the calibration flip is the only switch between
  advisory-flag and auto-edit.

### Delivery selection

With RiboNN `calibrated=False`, Stage D **highlights** RiboNN's suggested top and
flags splice-risky candidates, but the **user picks** the delivered sequence. When
RiboNN is calibrated, `rerank_by_expression` auto-selects the top (Stage D
collapses to one click). Again: the gate is the only switch.

## Why not a CNN in the inner loop (recorded for posterity)

Scoring RiboNN or SpliceAI/Pangolin per DP move or per SA move would require
thousands of forward passes of a 10 kb-context CNN — infeasible — and their
non-additive outputs cannot be a DP `delta()` anyway. The elegant substitute:
BT4's fast additive terms (CAI, tAI, codon-pair, ramp, and ViennaRNA ΔG in
refinement) **are the cheap correlates** of what the CNNs capture, so the loop
optimizes the proxies and the CNNs adjudicate a small finalist set. Even the
commercial incumbents (GenScript's "generate >10k then sort," Twist/IDT
generate-then-screen) use generate-and-rank, never inner-loop expression scoring.

## Computational efficiency summary

- **Prune before you spend:** score only the frontier + capped library, never
  every intermediate.
- **Batch the heavy models:** one RiboNN call and one SpliceAI+Pangolin pass over
  the whole candidate set + reference. Cost is dominated by fixed overhead, so
  scoring N candidates ≈ the wall-clock of scoring one (this is why the batched
  `score_many` is the first implementation step).
- **Cache by `(sequence, model-version)` hash:** reruns and overlapping designs
  are free and fit BT4's provenance model.
- **Background thread + `num_workers=0`:** BT4 Studio already runs engine work on
  a `QThread`; `num_workers=0` removes most of RiboNN's Windows worker-spawn
  overhead for small batches.

## Honest caveats (must survive into the UI and exports)

- The CDS is a **minority of the expression signal** (~31% per RiboNN's ablation;
  5′UTR/initiation folding dominates — see [`COMPARISON.md`](COMPARISON.md)). A
  CDS-only optimizer, however exact, optimizes a minority of the signal.
- The in-loop splice-motif constraint reduces *obvious* risk only; it is **not** a
  CNN-equivalent guarantee.
- SpliceAI/Pangolin/RiboNN are `calibrated=False` today; their contributions are
  advisory until gated, and every export stamps that status in the manifest.

## Implementation phasing (maps to future PRs)

1. **Batched RiboNN scoring** — a `score_many` / `delta_logte_many` on the adapter
   that runs the whole candidate set in one RiboNN invocation. Independently
   useful, testable without the GUI. *(First PR.)* ✅ **Landed** —
   `RiboNNExpressionModel.score_many` / `.delta_logte_many` reuse the batched
   `_predict_te` path; `delta_logte_many` scores the shared reference once;
   `score_sequence` / `delta_logte` delegate to them; `calibrated` stays `False`.
   The optional `num_workers=0` path was left out (RiboNN's predict entry point
   exposes no worker-count parameter; batching already amortizes the one-time
   worker spawn).
2. **Strong splice-consensus motif constraint** — a LOCAL constraint (new file +
   registry entry + `ok_suffix⇔validate` test), donor/acceptor consensus only,
   honestly labeled a heuristic. Wired through config/CLI/app. ✅ **Landed** —
   `bt4.constraints.SpliceSiteMotifConstraint` (`avoid_splice_sites`), IUPAC donor
   `GTRAGT` + acceptor `YYYYYYNYAGG`, **sense strand only** (no RC), never the bare
   `GT`/`AG`, `calibrated`-free structural heuristic; wired through
   `OptimizeConfig`, `--avoid-splice-sites`, the service schema, and BT4 Studio.
3. **Candidate-set assembly + `rerank` over frontier+library** — an API surface
   that assembles the frontier + deterministic repeat-fix library and batch-ranks
   it (`rerank_by_expression` applied across the set, calibrated-gated selection).
   ✅ **Landed** — `bt4.api.candidates` / `assemble_and_rank_candidates` →
   `CandidateSet`. Frontier + repeat-refined variants (gated on the seed actually
   violating a GLOBAL rule), de-duplicated, batch-scored via the new
   `BatchExpressionPredictor` contract when available, calibrated-gated delivery
   (discovery order + solver-delivered `chosen` when uncalibrated; reorder + top
   pick when calibrated). The solver-delivered sequence is pinned so it is
   invariant to `n`; de-dup/cap counts and the predictor identity are reported.
4. **Splice CNN audit (localize-and-flag)** — batched SpliceAI+Pangolin over the
   candidate set, returning per-site flags + backend agreement; no editing.
   ✅ **Landed** — `bt4.api.splice_audit` / `biomodels.splice.audit_splice` →
   `SpliceAuditReport`. Peak/NMS localization; per-flag `added_risk_vs_reference`
   (positive=worse, intra-backend) kept distinct from panel `delta_splicing`
   (larger=better); the pooled `backend_agreement` is the authoritative
   cross-backend signal, with `also_flagged_by` a raw ±window positional
   co-occurrence (not kind-agreement — Pangolin's combined track can't disagree on
   kind). Advisory only — `all_calibrated=False`, per-flag `calibrated`, **no
   editing**. Raw-sequence core in `biomodels/splice/audit.py`; `CandidateSet`
   adapter + `available_splice_backends()` in `pipeline/splice_audit.py`.
5. **BT4 Studio UI** — UTR fields, the two toggles, the annotated frontier +
   ranked table with uncalibrated badges, all on the background thread. ✅
   **Landed** — the **Candidates & splice audit** tab
   (`app/worker.py::CandidatesWorker`, `app/studio.py`) runs `api.candidates` →
   `api.splice_audit` on a background `QThread` and renders the ranked,
   honestly-labeled candidate set (delivered pick starred; per-member source /
   CAI / GC / expression+units / calibration / hard-violation / **distinct**
   splice-site counts) with two advisory banners. Same calibrated-gating honesty
   as the API: an uncalibrated head is shown as **discovery order, not a ranking**
   with the solver's pick starred and scores annotating only; the splice banner
   leads with **UNCALIBRATED (advisory)** when `all_calibrated` is `False`. Every
   table metric is recomputed per candidate from its own DNA (invariant #2); an
   opt-in toggle routes the installed SpliceAI/Pangolin CNNs into the audit.
6. **(Gated, future)** splice **auto-edit** and RiboNN **auto-select**, each
   unlocked only when its backend passes its fidelity/acceptance gate.

## Future / out of scope for v1

- **Expression-guided generational search** (mutate-the-best across batched RiboNN
  generations) — tractable *because* of batched scoring, but costly and steering
  by an (currently) uncalibrated model, so it would be a clearly-labeled
  **exploratory** mode, never the default delivery path, and only after RiboNN is
  calibrated.
- **Beyond-CDS design** (5′UTR / initiation-region) — where the literature says
  most of the expression signal lives; a larger scope change, noted here because
  it is the real ceiling on a CDS-only tool's expression claim.

## Cross-references

- [`../CLAUDE.md`](../CLAUDE.md) §6 (splice/expression contracts, out-of-loop
  scoring), §9 (Phase 3/4 status, acceptance/fidelity gates), §10.6 (no
  uncalibrated model presented as validated), invariant #5 (refinement never
  raises the hard-violation count), invariant #7 (determinism).
- [`COMPARISON.md`](COMPARISON.md) — why generate-and-rank (not inner-loop
  expression) is the field norm, and the CDS-minority-of-signal caveat.
- [`NEXT_SESSION.md`](NEXT_SESSION.md) — the batched-scoring perf follow-up that is
  step 1 here.
