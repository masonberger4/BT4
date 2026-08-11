# Changelog

All notable changes to BT4 are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it cuts
its first tagged release.

## [Unreleased]

### Added
- **BT4 Studio surfaces the engine-ready backends, gains library mode, and gets
  its Phase-4 polish** (`bt4.app`) — the two models that already existed behind
  `bt4.api` but had no UI are now wired in, plus the sampler and the accessibility
  work called for in CLAUDE.md §6.6. All of it is pure plumbing over the stable
  API (no engine change, no calibration claim):
  - **RiboNN in the Candidates tab.** An opt-in *Expression head* group (toggle,
    species, and the fixed 5'/3' UTR context the model requires) routes a
    `RiboNNExpressionModel` into `api.candidates`. The toggle is enabled **only**
    when `available_expression_backends()` reports the user's own checkout and
    weights actually resolve, so it is never a dead control, and it explains what
    is missing otherwise. Missing/non-DNA UTRs are refused *before* the run starts
    rather than raising mid-flight. RiboNN stays `calibrated=False`, so the banner
    still reads **discovery order, not a ranking** and the solver's pick stays
    delivered (§10.6).
  - **Validate with ASSP.** The one control that leaves the machine. It asks for
    consent first (naming the service and what is sent), runs
    `api.splice_crosscheck` on a background thread, and renders the report led by
    its tags — *network-derived, UNCALIBRATED, advisory, **not** part of the run
    manifest and never exported* — with the localized sites in a table. An outage
    degrades to a labeled "unavailable" banner and never fails a run (§10.15). The
    panel is cleared whenever the delivered sequence changes, so one sequence's
    splice sites can never be shown beside another's, and an export is
    byte-identical whether or not a cross-check ran (regression-tested).
  - **Library (sampled) tab.** `api.library` with members / temperature / seed
    controls, a per-member table, the selected member's sequence with its
    violation highlights, and a multi-record FASTA export whose every record is
    named `sampled`. The banner leads with **sampled, not optimized** — the
    `SAMPLED` certificate colours the badge directly, so it cannot drift from the
    claim the engine made — and reports measured diversity (distinct count, mean
    pairwise difference).
  - **Phase-4 polish.** A File/Run/View/Help menu bar with standard shortcuts
    makes every action keyboard-reachable; **View → System / Light / Dark**
    switches theme at runtime (restyling the stylesheet, both plots, the badges,
    and the sequence viewers' violation bands from the still-live results, via a
    new `SequenceViewer.set_dark`); tab order covers every new control and each
    carries an accessible name plus an explanatory tooltip.
  - **One source of truth for run gating.** All four flows (optimize, rank+audit,
    cross-check, library) share a `_wire_thread` helper and a `_update_run_buttons`
    gate driven by explicit running-flags rather than thread references — so a
    missed reference clear can no longer strand a control (the previous
    optimize-then-rank stuck-button class of bug is now structurally impossible).
  - New shared `_EngineWorker` base in `bt4.app.worker` (signal trio + the
    never-raise contract) with `CrossCheckWorker` and `LibraryWorker` alongside
    the existing two.
- **Public expression-backend registry** (`bt4.biomodels.expression.available_backends`
  / `resolve_backend`, re-exported as `bt4.api.available_expression_backends` /
  `resolve_expression_backend`) — the mirror of the splice resolver, so a frontend
  selects an expression head by name through the stable API instead of importing
  `biomodels` across a layer (§3, §10.9). `available_backends()` never raises and
  lists `"ribonn"` only when it can genuinely run; resolution is lazy (no torch
  import, no weight load) and confers **no** calibration.
- **Opt-in, out-of-loop ASSP splice cross-check** (`bt4.api.splice_crosscheck` /
  `bt4.pipeline.run_splice_crosscheck`, `bt4.biomodels.splice.AsspSplicePredictor`)
  — a **network** validator that runs the online ASSP service (Alternative Splice
  Site Predictor, Wang & Marín 2006) over an already-delivered sequence behind the
  existing `SplicePredictor` contract, closing the last non-human-gated gap in the
  splice subsystem (CLAUDE.md §6, §10.15). BT3's fatal splice bug was scraping this
  exact service **in the optimizer's inner loop as its only splice path**; BT4
  inverts every property of that mistake, structurally:
  - **Opt-in and out-of-the-inner-loop.** Requested explicitly by name and gated
    behind the `bt4[assp]` extra (httpx, lazily imported); it runs only as a final
    audit / validation pass on the delivered sequence, never per optimizer move, and
    is **never** returned by `splice.default()` or `available_splice_backends()`.
  - **Never blocking.** Rate-limited with exponential backoff and cached by
    sequence hash; if the service is unreachable or returns a garbled body the raw
    predictor raises an `AsspError`, but `run_splice_crosscheck` catches it and
    reports "unavailable" — a cross-check outage can never fail an optimization. The
    same graceful path covers a wrapped CNN's missing deps.
  - **Network-derived and non-reproducible.** `network_derived` is `True` and
    `calibrated` is `False`; ASSP numbers are excluded from the
    reproducible-from-manifest guarantee and reported as a separate advisory section
    (the CLI prints them to **stderr**, never into the stdout FASTA/JSON artifact or
    a `Result` manifest).
  - **Wired through the CLI** — `bt4 validate --splice-backend assp` and `bt4
    optimize --check-splice assp` (both flags also accept `pwm` / `pangolin` /
    `spliceai` for an offline or installed-CNN cross-check).
  - **CI never makes a live call.** The adapter is driven from committed **offline
    fixtures** (`tests/fixtures/assp/`, `FixtureAsspTransport`, selected via
    `$BT4_ASSP_FIXTURE_DIR`). Honest caveat: the live wire format is *unverified
    against the service* (unreachable during development), so the fixtures are
    *synthetic ASSP-format reports*, not real captures — the same "no bundled panel
    ships" posture as the wrapped CNNs.

### Fixed
- **RiboNN adapter: correct ensemble aggregation and honest empty-UTR guard.** The
  first end-to-end runs against real RiboNN weights surfaced two integration bugs.
  (1) RiboNN returns the ensemble as **one row per cross-validation model**, so a
  single input yields several rows sharing a `tx_id`; the adapter's `set_index`
  realignment then made `float(ordered[tx_id])` operate on a Series and raised
  `TypeError: cannot convert the series to <class 'float'>`. Realignment now groups
  by `tx_id` and averages (the ensemble mean, also averaging over cell types) via a
  new tested helper `_reduce_te_by_tx_id` — a no-op when rows are already unique.
  (2) Scoring with the default **empty** `utr5`/`utr3` crashed deep inside RiboNN's
  data loader (pandas reads an all-empty UTR column as `NaN` and its `.str`
  preprocessing fails); the adapter now refuses up front with a clear message, since
  the UTRs carry most of RiboNN's signal and an empty-UTR score is not meaningful.

### Changed
- **Python 3.10 is now supported** (was 3.11+). `requires-python` is lowered to
  `>=3.10`, the 3.10 classifier is added, ruff/mypy target 3.10, and CI's quality
  matrix now runs 3.10 alongside 3.11–3.13. The pure core uses no 3.11-only
  features, so this is a compatibility widening with no behavior change. It notably
  lets the wrapped **RiboNN** expression backend be installed into the same
  environment as its own dependency stack, whose pinned `torch==1.13.1` ships only
  CPython ≤3.10 wheels.

### Added
- **BT4 Studio "Candidates & splice audit" tab** — step 5 (final) of the
  expression/splice design flow, surfacing `api.candidates` → `api.splice_audit`
  in the desktop app. A background `CandidatesWorker` (mirroring the known-good
  `OptimizeWorker` `QThread` lifecycle) runs both on a worker thread and hands the
  window the candidate set + splice audit in one signal. The tab renders the
  ranked, honestly-labeled candidate table (delivered pick starred; per-member
  source / CAI / GC / expression+units / calibration / hard-violation / **distinct**
  splice-site counts) with two advisory banners: an *uncalibrated* expression head
  is shown as **discovery order, not a ranking** (solver's pick starred, scores
  annotating only; a calibrated head switches to ranked-by-expression), and the
  splice banner leads with **UNCALIBRATED (advisory)** whenever `all_calibrated` is
  `False`, reporting cross-backend agreement and stating the flags localize sites
  heuristically and edit nothing. Every metric is recomputed per candidate from its
  own DNA (invariant #2); an opt-in toggle routes the installed SpliceAI/Pangolin
  CNNs into the audit. The results area is now a `QTabWidget` (Design | Candidates &
  splice audit); the Design tab is unchanged. No Cancel control on this tab (the
  assemble→audit flow is not point-cancelable), and the cross-flow Optimize/Rank
  gating clears the worker-thread reference so neither button can deadlock.
- **Localize-and-flag splice audit** (`bt4.api.splice_audit` /
  `bt4.biomodels.splice.audit_splice`) — step 4 of the expression/splice design
  flow (`docs/DESIGN_expression_splice_flow.md` Stage C). An **out-of-loop,
  advisory** audit that runs the available `SplicePredictor` backends over a step-3
  candidate set to **localize** residual cryptic splice sites (one flag per
  contiguous above-threshold run, at its peak — non-maximal suppression) and attach
  the whole-panel **backend agreement** (pooled rank + sign) as the authoritative
  cross-backend confidence signal — built from the Delta-splicing values the audit
  already computed (a new shared `agreement_from_deltas` helper), so each backend
  scores every sequence **once**, never twice (§7). **It never edits** the sequences — a targeted synonymous auto-edit at flagged loci is a
  deliberately deferred, calibrated-gated future step. Honesty (CLAUDE.md §6/§10.6):
  every shipped backend is `calibrated=False` today, so `all_calibrated` is `False`
  and every `SpliceFlag` carries its **emitting backend's** `calibrated` flag; the
  site `threshold` is a **heuristic display knob** (not a validated cutoff) and the
  PWM baseline's per-position `score` is an uncalibrated **arbitrary-units**
  pseudo-score. Per-flag `added_risk_vs_reference` is **positive = worse** and
  strictly *intra-backend*, kept distinct from the panel-level `delta_splicing`
  (larger = better). Cross-backend `also_flagged_by` is a **raw positional
  co-occurrence** (±`match_window` nt, sized to the backends' anchor offsets),
  explicitly **not** a kind-level agreement (Pangolin reports one combined
  `P(splice)` and so can never disagree on kind — its flags are labelled `"splice"`,
  never donor-specific). New `biomodels/splice/audit.py` (raw-sequence core, imports
  only `domain` + the splice backends) + `pipeline/splice_audit.py` (the
  `CandidateSet` adapter + `available_splice_backends()`, which adds the wrapped
  SpliceAI/Pangolin CNNs when installed). Deterministic (#7). API-level surface (the
  BT4 Studio annotation UI is step 5).
- **Candidate-set assembly + expression rerank** (`bt4.api.candidates` /
  `assemble_and_rank_candidates`) — step 3 of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). Assembles the finalist set an
  expression head ranks: the **Pareto frontier** plus, when a GLOBAL rule is active
  *and* the delivered exact-DP seed actually violates it, a small **deterministic
  library of repeat-refined variants** (distinct refinement seeds over the delivered
  seed). The set is de-duplicated and scored by an `ExpressionPredictor` — in **one
  batched call** when the backend implements the new `BatchExpressionPredictor`
  contract (`score_many`, e.g. RiboNN), else per sequence — and delivered under the
  same **calibrated-gating** honesty rule as `rerank_by_expression`: an uncalibrated
  head (the default placeholder, and the shipped RiboNN adapter) only *annotates* —
  the set stays in **discovery order** (`order_basis="discovery"`) with the
  solver-delivered sequence `chosen` — while a calibrated head reorders by predicted
  expression (`order_basis="expression_rank"`, total order `(score desc, index asc)`)
  and re-picks the top (CLAUDE.md §10.5/§10.6). Hardened for correctness/honesty: the
  **delivered (`chosen`) sequence is invariant to `n`** (uncalibrated, the
  solver-delivered sequence is pinned first in discovery order; calibrated, the
  head's top pick is the top of the top-n keep — the cap is applied *after* scoring
  so a calibrated reranker never loses its best candidate);
  every member is a full `Result` (round-trips, metrics recomputed, certificate,
  residual GLOBAL violations disclosed); variants are labelled `repeat_refined` (the
  *process*, not a guaranteed fix); and de-dup/cap counts, the batch-path flag, and
  the predictor identity (folded into the manifest, invariant #9) are all reported.
  New `BatchExpressionPredictor` Protocol in `bt4.biomodels.expression`; `_refine`
  gains an optional `seed` (default unchanged). API-level surface (UI wiring is
  step 5). No calibration claim — ranking is a reporting no-op until a head is
  calibrated.
- **Strong splice-consensus motif constraint** (`bt4.constraints.SpliceSiteMotifConstraint`,
  `avoid_splice_sites`) — step 2 of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). A new **LOCAL, exact-in-the-trellis**
  hard constraint that forbids the *strong* splice-consensus **donor** (`GTRAGT`,
  the intronic +1..+6 core) and **acceptor** (`YYYYYYNYAGG`, a polypyrimidine tract
  + `NYAG|G`) motifs on the mRNA **sense strand only** (splicing is strand-specific,
  so — unlike restriction/repeat motifs — there is **no** reverse-complement
  banning). It is an honest **structural heuristic**, not a splice model: it reduces
  only the most *obvious* cryptic-splice risk and makes no calibrated claim; the
  wrapped SpliceAI/Pangolin CNNs do the real audit out of loop (CLAUDE.md §6,
  §10.6). It **never** bans the ubiquitous bare `GT`/`AG` (governing rule 3). The
  default patterns (Shapiro & Senapathy 1987; Zhang 1998; human/mammalian
  major-spliceosome only) are deliberately specific (~1/2048 donor, ~1/8192
  acceptor) so the hard veto rarely over-constrains a design, and are configurable
  via `donor_motifs`/`acceptor_motifs`. `ok_suffix⇔validate` and `context_len`
  sufficiency (5 donor / 10 acceptor) are property-tested (invariant #3). Wired
  through `OptimizeConfig`, the `bt4` CLI (`--avoid-splice-sites`), the `service`
  schema, and BT4 Studio (a checkbox with an explanatory tooltip); off by default.
- **Batched RiboNN scoring** (`RiboNNExpressionModel.score_many` /
  `.delta_logte_many`) — the first step of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). RiboNN's cost is dominated by fixed
  *per-invocation* overhead (weight hashing + model load + its DataLoader worker
  spawn), so scoring a whole candidate set one sequence at a time paid that cost
  N times. The new public batch methods route the entire set through the existing
  batched `_predict_te` path (one temporary TSV, one `predict` invocation — RiboNN's
  `top_k`-model ensemble runs inside that single call), so scoring a Pareto frontier
  costs roughly the wall-clock of scoring a single sequence.
  `delta_logte_many` additionally scores the shared **reference once** (appended to
  the batch), not once per design. Both preserve per-input validation (valid DNA,
  length-3N ending in a stop codon, non-empty `utr5`/`utr3`) and the `tx_id`
  realignment; results come back **in input order**. `score_sequence` /
  `delta_logte` now delegate to the batch methods (single source of truth). A
  `num_workers=0` DataLoader path was investigated and **deliberately left out**:
  RiboNN's `predict_using_nested_cross_validation_models` exposes no worker-count
  parameter, so requesting 0 workers would mean patching RiboNN internals (against
  the "wrap, never reimplement" contract), and batching already amortizes the
  one-time worker spawn across the set. `calibrated` stays **`False`** — no
  calibration claim. Tested without torch / pandas / the RiboNN checkout (batch
  ordering, ensemble averaging per `tx_id`, reference-scored-once, and the
  empty-UTR / bad-CDS guards still firing).
- **Wrapped RiboNN expression backend** (`bt4.biomodels.expression.RiboNNExpressionModel`)
  — the Phase-4 learned expression head behind the `ExpressionPredictor` contract
  (CLAUDE.md §6/§9). It runs the published **RiboNN** translation-efficiency CNN
  (Zheng, Persyn, Wang et al., *Nat Biotechnol* 2025; Sanofi / Cenik Lab)
  inference-only as an out-of-loop frontier reranker. **License:** RiboNN's code
  and weights are each **Sanofi non-commercial** (academic/non-commercial only) —
  compatible with BT4's open-source non-commercial scope and, like SpliceAI's
  CC BY-NC weights, **never bundled**: the adapter drives the user's own RiboNN
  clone (lazily importing the repo's `src`, pointed at via `$BT4_RIBONN_DIR`) and
  their Zenodo weights. Every weight it loads is verified against a bundled
  180-entry SHA-256 manifest (`data/ribonn_sha256.json`, 90 human + 90 mouse —
  public content hashes only) **before** `torch.load`. The score is in RiboNN's
  native **CLR-residual TE** units (never exponentiated); `delta_logte(designed,
  reference)` gives the UTR-fixed, CDS-attributable Δ (negative = a CDS change
  predicted to *reduce* expression), analogous to Pangolin's `delta_splicing`.
  Ships **`calibrated=False`** (`default()` still returns `NullExpressionModel`):
  faithful reproduction is not calibration for BT4's CDS-variant regime, so
  promotion requires a passing `verify_expression_gate` on a regime-matched panel
  (human-only, data-gated). New `bt4[expression-ribonn]` extra (torch + pandas),
  lazily imported so `import bt4` stays light.
- **Model-agnostic expression acceptance-gate harness**
  (`bt4.biomodels.expression.gate`) — the honest gate a learned expression head
  must pass to earn `calibrated=True` (CLAUDE.md §6/§8/§10.6, Phase 4). For a
  log-TE regression head it reports **Spearman** (primary), **Pearson**, **R²**,
  and **split-conformal coverage** at a target level (default 90%), evaluated on a
  **group-disjoint split** (homology cluster / chromosome) so no group leaks
  across calibration and test — the distribution-shift-aware check that a head
  validated only on natural-gene TE has *not* earned calibration for BT4's
  CDS-variant regime. `passed` requires both the Spearman threshold **and**
  conformal coverage near target (point accuracy *and* honest uncertainty). The
  gate never flips anything: thresholds are inputs set at gate time, and the
  neutral `NullExpressionModel` provably cannot pass (its zero-variance scores
  give Spearman 0). New `ExpressionEvalCase` / `ExpressionGateReport` and a
  `run_expression_gate(predictor, samples)` wrapper. Fully dependency-free and
  tested without torch or any real model, mirroring how the splice
  fidelity/attestation machinery shipped before a calibrated backend.
- **Shared dependency-free statistics** (`bt4.biomodels._stats`) — `pearson`,
  `spearman` (moved from `splice.agreement`, which now re-exports them), plus
  `r2_score`, `conformal_quantile` (finite-sample split-conformal), and
  `empirical_coverage`. Single well-tested home for the estimators the splice
  agreement report and the expression gate both use.
- **License-clean splice fidelity-attestation layer**
  (`bt4.biomodels.splice.attestation`) — the honest promotion path for the wrapped
  Pangolin / SpliceAI backends (CLAUDE.md §6, §10). A `FidelityAttestation` records
  **only** a passing integration-fidelity gate's derived scalars (`passed`,
  `max_abs_deviation`, `n_cases`, `tolerance`) plus the public pinned weight
  SHA-256s and the tool version — **never** a `FidelityCase` raw per-position score
  (those are the license-encumbered model outputs). The shape is enforced
  structurally (`_ALLOWED_FIELDS` + an honesty test asserting no raw-score field is
  serializable), and `from_dict` refuses any unexpected key. `attest_backend`
  refuses to record a failing or too-loose gate; `verified_predictor(predictor,
  attestation)` is the single seam that flips a backend to `calibrated=True`, and
  only when the attestation passed, clears the `MAX_ATTESTATION_TOLERANCE` floor,
  and its weight SHAs exactly match the adapter's `PINNED_WEIGHT_SHA256` (a
  refusal, never a silent downgrade). A deterministic, timestamp-free
  `content_hash` makes an attestation a provenance-manifest stamp. This layers the
  committed-record / private-execution / user-opt-in / baseline-fallback options;
  no attestation ships, so `default()` still returns the honest PWM baseline. Both
  Pangolin (GPL) and SpliceAI (CC BY-NC) are eligible to certify under BT4's
  open-source, non-commercial scope.

## [0.4.0] - 2026-08-04

First tagged release since 0.3.1, capturing the Phase 1 performance and Phase 3
refinement/splice wave: the full Rust trellis port, richer refinement moves, the
wrapped SpliceAI splice backend, and the last Phase 2 budget item.

### Added
- **Full Rust trellis port** (`bt4_native.trellis_solve`) — the exact-DP inner
  loop of `bt4.optimize.exact_dp.solve_exact` now runs in Rust (Phase 1, CLAUDE.md
  §7), following the existing native-primitive pattern: a PyO3 `#[pyfunction]` with
  a byte-identical pure-Python twin (`bt4._accel._py_trellis_solve`) and a
  Hypothesis equivalence test pinning the two. The DP is callback-driven, so Rust
  never calls back into Python: a **regime gate** restricts the native path to
  position-independent objectives (no `POSITIONAL` term — `CpbTerm` was made
  context-based so PAIRWISE terms stay position-independent), Python **precomputes**
  the reachable-context transition graph and the pre-summed per-transition deltas
  (fixing the float summation order, so the lexicographic tie-break is bit-for-bit
  identical), and the layer DP runs in Rust; it **falls back to the pure-Python
  DP** whenever the regime does not hold, the extension is absent, or a
  context-count cap is exceeded. A single solve is not accelerated (the Python
  precompute costs ~a whole pure DP), so `run_optimize` stays on the pure path; the
  win is the **Pareto frontier**, which builds the transition graph once and reuses
  it across every scalarization point (only the cheap deltas recomputed) with the
  DP in Rust — a measured ~2.7–5.5x `run_frontier` speedup with **byte-identical**
  DNA, objective scalars, and certificates.
- **Block/segment moves + parallel tempering in the SA refinement engine**
  (`bt4.optimize.anneal_refine`, Phase 3 — CLAUDE.md §7, §9). The engine gained
  four opt-in knobs: `block_size` / `block_prob` (coordinated multi-position
  synonymous swaps) and `replicas` / `temps` / `swap_every` (a parallel-tempering
  replica ladder with standard replica-exchange Metropolis swaps). These widen the
  refinement's *reach* so it can cross a barrier that only clears when several
  codons move **together** — a dispersed max-repeat or out-of-frame uORF the
  single-codon chain could leave in place — **without weakening invariant #5**:
  block candidates pass the same local (union-of-windows `ok_suffix`) and global
  (whole-sequence recount) feasibility gates, every replica gates against its own
  current hard-violation count, every visited configuration keeps a global count
  `<=` the seed's, and the delivered result is ranked lower-global-count-first then
  higher-score. All four default off, and with them off the engine reproduces the
  prior single-chain trajectory **byte-for-byte** (invariant #7). Block moves
  always full-`score` re-score (never `delta_score`), since summing per-position
  deltas is only valid for additive disjoint-context terms. The honest
  **feasibility floor** is preserved: a repeat pinned to synonymously-immovable
  bases (Met `ATG` / Trp `TGG`) is unremovable by any move and is still reported as
  a residual, never claimed clean. New Hypothesis tests pin the never-raise-global
  guarantee under block+tempering, determinism/round-trip with replicas and blocks,
  the default-knobs no-op, and the immovable-repeat feasibility floor.
- **Wrapped SpliceAI splice backend** (`bt4.biomodels.splice.SpliceAiSplicePredictor`)
  — the second *wrapped published* splice CNN behind the `SplicePredictor`
  contract, the cross-check to Pangolin (Phase 3, CLAUDE.md §6). It runs the
  published **SpliceAI** model (Jaganathan et al. 2019) inference-only, and its
  3-way per-position softmax (null/acceptor/donor) maps *cleanly* onto
  `SpliceResult.acceptor` and `.donor` (both populated, unlike Pangolin's single
  combined track). **License (verified): SpliceAI code is PolyForm Strict 1.0.0
  and its weights are CC BY-NC 4.0 (noncommercial) — even more restrictive than
  Pangolin's GPL, so no-bundle is mandatory**; the adapter lazily imports the
  user's own installed `spliceai` package + weights, SHA-256 hash-pinning them
  (verified before load). Ships **`calibrated=False`** (`verify_spliceai_fidelity`
  is the gate; no reference panel bundled), so `default()` still returns the PWM
  baseline. With both CNNs installed, the agreement harness now compares two real,
  independently-trained splice models (no harness code change needed — it already
  compares at the pooled-Δ level). New `bt4[splice-spliceai]` extra
  (TensorFlow), lazily imported so `import bt4` stays light.
- **Wrapped Pangolin splice backend** (`bt4.biomodels.splice.PangolinSplicePredictor`)
  — the first *wrapped published* splice model behind the existing
  `SplicePredictor` contract (Phase 3, CLAUDE.md §6). It runs the already-validated
  **Pangolin** CNN (Zeng & Li 2022) as an inference-only backend, feeding its
  per-nucleotide `P(splice)` into the shipped Δsplicing / top-k-log-odds framing.
  **License-clean:** Pangolin is **GPL-3.0** (the earlier roadmap's "MIT" was
  wrong), so — exactly as BT4 wraps GPL ViennaRNA — the adapter **lazily imports
  the user's own installed `pangolin` package and weights and bundles neither**;
  BT4 stays MIT. Weights are **SHA-256 hash-pinned** (the published v1.0.2 digests)
  and verified *before* they are unpickled, keeping runs
  reproducible-from-manifest. The adapter reproduces upstream Pangolin's scores
  **bit-for-bit** yet ships **`calibrated=False`** (no reference panel is bundled;
  `verify_pangolin_fidelity` is the promotion gate), so `default()` keeps returning
  the honest PWM baseline. Heavy deps behind the new `bt4[splice-pangolin]` extra,
  lazily imported so `import bt4` stays light.
- **Two-backend splice agreement harness** — `bt4.biomodels.splice.backend_agreement`
  reports each available backend's Δsplicing ranking, pairwise **Spearman rank
  agreement**, and sign agreement across candidates (the first-class uncertainty
  signal of CLAUDE.md §6/§8); it reports, it does not judge. Exposed as the
  standalone runner `scripts/compare_splice_backends.py` (`--fasta`, `--json`),
  which degrades to the baseline alone — and says so — when neither CNN backend
  (Pangolin nor SpliceAI) is installed.
- **CpG / UpA whole-sequence count budget** (`dinuc_budget` + `dinuc_min` /
  `dinuc_max`; CLI `--cpg-min/--cpg-max` and `--upa-min/--upa-max`) — the last
  Phase 2 item. A dinucleotide count does not decompose per-codon (a 2-mer
  straddles the codon boundary), so the amount-bucketed budget DP
  (`bt4.optimize.lagrangian`) now takes a **context-aware** per-codon amount
  (`bt4.objectives.dinucleotide.dinucleotide_amount`) attributing each occurrence
  to the codon holding its END base, with a new `budget_context` folded into the
  trellis state so a straddling count stays exact. Enforced by the same **exact
  bucketed DP** as the GC budget, with a `proven_optimal` certificate and every
  local constraint still honored. Mutually exclusive with the GC budget, and (like
  it) not combinable with `refine` / `max_repeat_length` / `avoid_uorf`. Wired
  through `OptimizeConfig`, the CLI, and the `service` request schema.
- **Library / degenerate-design mode (opens Phase 5).** `api.library(protein,
  config, n, *, seed, temperature)` and `bt4 library PROTEIN --n N` sample a
  *library* of coding sequences by drawing from each residue's synonymous-codon
  distribution (organism usage frequencies raised to `1/temperature`), keeping
  only codons that pass every LOCAL constraint. This is an honest **stochastic
  sampler, not an optimizer**: every member round-trips and carries metrics
  recomputed from its own DNA, the library is fully deterministic from its seed,
  and each result carries the new **`OptimalityStatus.SAMPLED`** certificate,
  which makes no optimality or expression claim. GLOBAL constraints
  (`max_repeat_length`, `avoid_uorf`) are not enforced during sampling but are
  validated and any residual violation reported honestly per member. New modules
  `bt4.optimize.sample` (deterministic constrained sampler, `domain`-only) and
  `bt4.pipeline.library` (`LibraryResult` + `run_library`).
- **Two more `bt4_native` hot-loop primitives** (`max_gc_run`, `longest_repeat`),
  each with a byte-for-byte pure-Python fallback in `bt4._accel` and a Hypothesis
  equivalence property test that pins the Rust and Python paths together (and, for
  `longest_repeat`, cross-checks `longest_repeat(seq) > m` iff
  `MaxRepeatConstraint(m).validate(seq)` flags a hard violation). This is honest
  incremental native acceleration — **not** a full trellis inner-loop port, which
  still remains (CLAUDE.md §7, §9 Phase 1).

### Changed
- **`GcRunConstraint.ok_suffix` now calls the (optionally Rust-accelerated)
  `bt4._accel.max_gc_run`** on its bounded trailing window, with no change to
  observable behavior (the pure-Python fallback is the same scan as before). The
  `longest_repeat` primitive is added and cross-checked against
  `MaxRepeatConstraint`, but is **deliberately not** placed on the per-SA-move
  `MaxRepeatConstraint.validate` hot path: the whole-sequence longest-repeat is
  O(n²), which is *slower* than the constraint's existing O(n·k) k-mer scan when
  the native extension is absent — so wiring it there would regress the common
  pure-Python path (CLAUDE.md §7, "everything incremental"). Every existing
  `ok_suffix ⇔ validate` and constraint test passes unchanged.

## [0.3.1] - 2026-08-01

BT4 Studio first-run polish: the desktop app now guides a non-technical user
through mistakes with plain-language messages instead of raw Python errors, and
never leaves a stale result behind a failed run.

### Added
- **Cancel button + live progress** for BT4 Studio. The frontier sweep now
  reports per-point progress (`solving frontier point 3 of 9`) and can be stopped
  mid-run; cancelling returns the partial frontier computed so far. `api.frontier`
  / `run_frontier` gained optional `on_progress` and `should_cancel` hooks.
- A one-time **warning before optimizing a very long protein** (it may take a
  while, and the run is cancelable).
- `bt4.api` now re-exports `InfeasibleError`, `validate_protein`, `AMINO_ACIDS`,
  and `available_tai_organisms` so frontends can validate input and translate
  failures without reaching past the API layer.

### Changed
- **Plain-language input handling in BT4 Studio.** Pasting a FASTA record strips
  its header automatically; an empty box, a trailing `*` stop, or non-amino-acid
  characters get a clear, specific message (not a Python `repr`); restriction-
  enzyme names are matched case-insensitively and unknown ones list the valid
  catalog; an infeasible constraint set explains which knobs to relax instead of
  saying "no feasible codon". The **tAI** checkbox is now labelled correctly and
  enabled only for organisms that ship a tRNA table.

### Fixed
- **A failed run no longer leaves a stale, exportable result on screen** — the
  results panel (and the delivered result behind Export) is cleared on failure,
  so Export can't silently write the previous sequence.
- `scripts/sensitivity.py` detected tAI availability via the pre-0.3.0 organism-
  list quirk and silently returned `None` for every organism after that quirk was
  fixed; it now uses `api.available_tai_organisms()`.

## [0.3.0] - 2026-08-01

First release with a **downloadable, double-clickable BT4 Studio app** for
Windows / macOS / Linux, plus a wave of Phase 2/3 objectives, constraints, and
solver backends.

### Added
- **5' translation-ramp objective** (`RampTerm`) -- a heuristic that prefers
  slower codons in the first N codons (`ramp_weight` / `ramp_codons`).
- **CpG / dinucleotide objective** (`DinucleotideTerm`) to deplete (stealth) or
  elevate (immunostimulatory) CpG content (`cpg_weight` / `cpg_mode`).
- **Codon-pair bias** (`CpbTerm` + `build_codon_pair_table`): a pairwise objective
  built from a reference CDS set, solved exactly in the trellis via a new
  `objective_context` on the DP (the state now carries the previous codon).
- **OR-Tools CP-SAT backend** (`bt4.optimize.cpsat.solve_cpsat`, `bt4[ilp]`
  extra): solves the additive objective under a global **GC budget** (`gc_min` /
  `gc_max`) with a proven-optimal / gap-bounded certificate. New `ilp` CI job.
- CLI flags for all of the above (`--ramp-weight`, `--cpg-weight`, `--cpg-mode`,
  `--gc-min`, `--gc-max`) and a CpG control in BT4 Studio.

### Changed
- **Idiot-proof, double-clickable app packaging.** The PyInstaller spec now emits
  a *single* file per desktop OS instead of a one-folder zip: a one-file
  `BT4-Studio-Windows.exe`, a one-file `BT4-Studio-Linux-x86_64`, and (on macOS) a
  `.app` that CI wraps in a drag-to-Applications `BT4-Studio-macOS.dmg`. Verified
  end-to-end on Linux: the one-file build launches BT4 Studio and runs its event
  loop. The README's install section is rewritten for non-technical users
  (download-one-file table + how to click past the unsigned-app OS warnings), with
  the from-source/CLI install moved to a "for developers" section.
- **Release pipeline is now re-drivable and self-healing.** `release.yml` accepts
  a `workflow_dispatch` `ref` input to rebuild an existing tag's source and
  idempotently (re)attach the per-OS app + wheel/sdist to its release — the
  honest, non-destructive way to repair a release that has no assets. The publish
  step now also fails loudly instead of publishing an empty, asset-less release.
  See [`packaging/README.md`](packaging/README.md#repairing-a-release).
- **CI now launches the packaged app.** A `bt4-studio --self-test` hook builds the
  main window (loading the bundled data + Qt/pyqtgraph) and exits without the
  event loop; the release workflow runs it against the freshly built bundle on
  each OS, so a bundle that builds but crashes on first launch fails CI instead of
  shipping. The macOS `.app` also now carries its real version in `Info.plist`,
  the codon/tRNA data dir is a regular package (reliable frozen-bundle resource
  loading), and the Windows asset rename/upload no longer depends on a fragile
  cross-shell absolute path. A full non-technical [`docs/INSTALL.md`](docs/INSTALL.md)
  guide was added.

### Fixed
- The only tagged release (`v0.2.0`) had **no downloadable app**: its publish step
  ran the pre-idempotency workflow and `gh release create` failed on "release
  already exists" (the tag/release were made in the UI first), so the built
  bundles never attached. The pipeline is now idempotent and re-drivable, and the
  docs no longer point users at an empty Releases page.
- **`available_organisms()` listed bogus organisms.** It matched every `*.tsv`,
  so the tAI tRNA tables leaked in as `homo_sapiens.trna`, `mus_musculus.trna`,
  and `saccharomyces_cerevisiae.trna` — visible in the app's organism dropdown and
  `bt4 organisms`, and unloadable as codon tables. The tRNA tables are now
  excluded (they remain available via `available_tai_organisms()`).

## [0.2.0] - 2026-07-31

Richer biology and surfaces on top of the exact-DP core.

### Added
- **Restriction-site constraint** (`bt4.constraints.RestrictionSiteConstraint`,
  `available_enzymes`): an IUPAC-aware matcher and a catalog of common enzymes
  (EcoRI, BamHI, NotI, ...), always avoiding each site's reverse complement.
  Wired into `OptimizeConfig.restriction_enzymes`, the CLI (`--enzyme`,
  `bt4 enzymes`), and BT4 Studio.
- **More organisms**: representative *E. coli* K-12 and *S. cerevisiae*
  codon-usage tables (auto-discovered; clearly labeled representative).
- **`bt4 build-table`** and `bt4.io` FASTA parsing: recompute an authentic codon
  table from a user-supplied CDS FASTA (Laplace-smoothed so the result always
  loads), with a content-hashed provenance sidecar.
- **`bt4.service`**: an optional FastAPI HTTP API (`/optimize`, `/frontier`,
  `/validate`, `/organisms`, `/health`) that calls only `bt4.api`.
- **Benchmark harness** (`scripts/benchmark.py`) and a golden/regression test
  suite pinning current optimizer output.

### Fixed
- BT4 Studio frontier plot now shows raw CAI / GC-fraction axis values instead of
  a rescaled "x0.001" SI-prefix label.

## [0.1.0] - 2026-07-31

First tagged release: an honest exact-DP codon optimizer with a CLI and the BT4
Studio desktop app.

### Added
- **Exact codon-trellis DP solver** (`bt4.optimize`) over the true per-constraint
  context, with an explicit `beam` speed knob and a machine-readable
  `OptimalityCertificate` (`proven_optimal` / `beam_truncated`).
- **Objective terms** (`bt4.objectives`): `CaiTerm` (log relative-adaptiveness)
  and `GcProximityTerm`, both additive with `delta == score` property tests.
- **Constraints** (`bt4.constraints`): `HomopolymerConstraint` and
  `ForbiddenMotifConstraint` (with automatic reverse complements), with
  `ok_suffix ⇔ validate` agreement property tests.
- **Pipeline + stable API** (`bt4.pipeline`, `bt4.api`): `optimize()`,
  `frontier()` (a CAI/GC Pareto frontier), and `validate()`, with metrics
  recomputed from the delivered DNA and a content-hashed provenance manifest.
- **`ObjectiveTerm` / `Constraint` protocols** and the `Scope` enum in the pure
  `domain` layer (the shared vocabulary the optimizer speaks).
- **`bt4` CLI**: `optimize`, `validate`, `organisms`, and `--version`.
- **BT4 Studio** (`bt4.app`): a native PySide6 desktop app calling `bt4.api` on a
  background thread — constraint controls, an honest optimality-certificate
  badge, a recomputed-metrics table, an interactive CAI/GC frontier plot, a
  sequence viewer, and FASTA/JSON export. Offline; nothing leaves the machine.
- **IO** (`bt4.io`): FASTA and versioned, deterministic JSON export.
- **Packaging & distribution**: a `packaging` extra, a PyInstaller spec
  (`packaging/bt4-studio.spec`) that builds a standalone BT4 Studio bundle, and a
  `Release` workflow that publishes the sdist + wheel and per-OS app bundles on a
  version tag.
- **Community health**: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`,
  and GitHub issue / pull-request templates; a landing-page `README.md` with a
  screenshot of BT4 Studio.

### Fixed
- **Wheel/sdist builds** were broken by a `pyproject.toml` `force-include` that
  double-added the codon data files (`pip install .` failed with "a second file
  is being added to the wheel archive at the same path"); replaced with
  `artifacts` so the data and `py.typed` marker ship exactly once.
- **Two import-linter layering violations** (`optimize → constraints`,
  `objectives → biomodels`) that surfaced once `bt4.app` existed — resolved by
  lifting the protocols into `domain` and decoupling `CaiTerm` from the codon
  table, keeping every pure layer importing only `domain`.

### Notes
- Richer objectives (tAI, codon-pair, 5′ ramp), ILP / relaxation backends, and
  the validated splice / folding / expression models are on the roadmap and are
  **not** yet shipped — see [`CLAUDE.md`](./CLAUDE.md) §9.
