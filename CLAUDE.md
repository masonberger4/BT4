# CLAUDE.md — BT4

Guidance for Claude Code (and humans) building **BT4**, the from-scratch
successor to BT3. This file is the constitution of the new repository: read it
before writing code, and keep it current as the architecture evolves.

> Status: **Phases 0-2 complete; Phase 3 groundwork landed; Phase 5 opened.** The
> pure
> `domain` layer, provenance manifest, packaging, layering contract, and CI are
> in place, and on top of them an **honest exact-DP core** now ships: a codon
> trellis with true per-constraint context and a real optimality certificate,
> the `CaiTerm`/`GcProximityTerm` objectives and `Homopolymer`/`ForbiddenMotif`
> constraints (with their `delta==score` and `ok_suffix⇔validate` property
> tests), a **CAI/GC Pareto frontier**, the stable `bt4.api`, a `bt4` CLI, and
> the first cut of **BT4 Studio** (the PySide6 desktop app). Phase 2 has since
> added codon-pair/ramp/CpG/%MinMax objectives, tandem/inverted-repeat plus
> **max-GC-run** and dispersed **max-repeat-length** constraints, named
> **forbidden-sequence presets**, and **two budget backends** — OR-Tools CP-SAT
> and an honest **Lagrangian/exact-bucketed budget DP** that keeps local
> constraints under a global budget; that budget DP is now **context-aware**, so
> the final Phase 2 item — **CpG/UpA whole-sequence count budgets** — ships as an
> exact, proven-optimal dinucleotide-count budget. **tAI has since landed on real
> GtRNAdb tRNA data for eight organisms**, and Phase 3 groundwork is in: the
> `FoldingModel` (ViennaRNA + labeled baseline) and `SplicePredictor` (labeled PWM
> baseline) contracts, the incremental SA refinement engine (with a
> global-constraint gate), and per-site risk tracks plotted in BT4 Studio.
> **Phase 5 has opened** with an honest **library / degenerate-design mode** (a
> deterministic codon-distribution sampler with a `SAMPLED` certificate, not an
> optimizer). Two more native `bt4_native` primitives (`max_gc_run`,
> `longest_repeat`) also landed. **Both wrapped published splice backends** have
> now landed too: `PangolinSplicePredictor` (Pangolin, GPL-3.0, PyTorch) and
> `SpliceAiSplicePredictor` (SpliceAI, PolyForm Strict code + CC BY-NC weights,
> TensorFlow), each driving the user's own installed package (*not* bundled;
> lazily imported like ViennaRNA), hash-pinned and reproducing upstream scores
> bit-for-bit, shipped `calibrated=False` pending their integration-fidelity
> gates, with a two-backend agreement harness that makes agreement between two
> real CNNs an uncertainty signal. The **full Rust trellis port** has also landed
> (`bt4_native.trellis_solve`, regime-gated with a byte-identical pure-Python
> fallback, amortized across the Pareto frontier). The **wrapped RiboNN expression
> head** has landed too (`RiboNNExpressionModel`, Sanofi non-commercial, driven
> from the user's own checkout, hash-pinned, `calibrated=False` until its
> CDS-variant gate) and has now had its **first real end-to-end runs against the
> licensed weights**, which validated the adapter and fixed two live-only
> integration bugs (ensemble row-per-model aggregation; a required-non-empty-UTR
> guard). BT4 now also supports **Python 3.10** (was 3.11+), so RiboNN installs into
> the same environment as its `torch==1.13.1` stack. Still ahead: recording the
> fidelity/acceptance gates to promote splice + expression to `calibrated=True`, and
> packaged installers — see §9. This document was written
> after a full review of the BT3 codebase and *every* BT3 branch (`master`,
> `almost-there`, `gemini`, `streamlit`, and the merged
> `claude/ultracode-app-redesign` line); the lessons are folded in below.
>
> **Keep this current.** CLAUDE.md is the constitution: when a phase lands, a
> contract changes, or the architecture evolves, update this file *in the same
> change* — a stale constitution is a BT3 anti-pattern (§10.11).

---

## 1. What BT4 is

BT4 back-translates a **protein** into a **coding DNA / mRNA** sequence that is
**optimized for real expression outcome** in a target organism (default *Homo
sapiens*) **subject to biological constraints** (GC content, homopolymers,
tandem/inverted repeats, forbidden & restriction motifs, internal ATG / Kozak /
uORF, CpG budget, cryptic splice sites, 5′ mRNA folding).

The load-bearing idea BT4 inherits from BT3 and keeps as bedrock:

> **This is constrained combinatorial optimization, not greedy per-codon
> substitution.**

Everything else BT4 rebuilds. The one-paragraph thesis that shapes the whole
design:

- The objective is a **vector**, not a scalar. CAI is *one weak, cheap prior*
  among many (tAI, codon-pair bias, 5′ ramp, folding ΔG, splice Δrisk, CpG/UpA
  budget, learned expression). BT4 returns a **Pareto frontier** and tells you
  where on it your answer landed.
- Objectives and constraints are **partitioned by locality**. Additive +
  bounded-context terms are solved *exactly* in a codon trellis whose state
  carries the *true* union of every constraint's declared context — **no global
  context cap that silently over-merges**. Genuinely non-local terms (folding,
  whole-sequence splice, learned expression) live in a **refinement /
  relaxation layer** with incremental delta-scoring — never pretending to fit a
  window.
- The solver is **honest about optimality**. Exact DP when the state space is
  small; ILP / CP-SAT or Lagrangian relaxation with a real **optimality-gap
  certificate** when it isn't; beam search only as an explicit speed knob. Every
  result states *how optimal it is and what was relaxed*.
- The ML is **real and validated, or it refuses to claim otherwise**. Shipped,
  hash-pinned models — **ViennaRNA folding** (calibrated today), the **published
  SpliceAI/Pangolin** wrappers for Δsplicing inference (shipped `calibrated=False`
  until their integration-fidelity gate is recorded), and a *planned* optional
  learned expression head — each gated on real held-out data or, for a wrapped
  published model, an integration-fidelity check against its known outputs. A
  model reports `calibrated=True` only once its gate passes; until then it is
  loudly labeled uncalibrated (§10.6).

### The governing principle

**Never present an unenforced constraint, an unvalidated number, or a
heuristic result as if it were real.** BT3's own origin story is the cautionary
tale: the original tool had constraint "checks" that were `return dna` no-ops
and a splice model that predicted from *position with no sequence content*. BT4
makes honesty **structural** (enforced by tests and certificates), not
aspirational (documented in prose).

---

## 2. What we learned from BT3 (all branches)

BT3's history is a real asset — study it, don't re-walk it.

| Branch | What it was | What BT4 takes from it |
|---|---|---|
| `almost-there` | Original design at its most mature: **backtracking + CAI heuristic pruning** (log-weight suffix bounds ≈ branch-and-bound), a real **splice-site optimizer**, **committed benchmark FASTA datasets** (L10–L1000 × N01–N10), a devcontainer, and a broad test suite. | Keep the benchmark corpus idea and the devcontainer. The B&B suffix bound is the seed of BT4's exact-DP admissible bound. |
| `gemini` | Performance-focused MVP fork ("600 bp", "performance seq gen"). | Runtime is a first-class concern from day one — but measured, not vibed. |
| `streamlit` | A **Streamlit UI** with FASTA upload, per-constraint toggles, CAI upload. | BT4 ships a real UI/service story; don't let it rot on a side branch. |
| `master` (rebuild) | From-scratch **beam-DP trellis**, layered pure core, local ML splice subsystem, honest metrics, provenance, multi-organism tables, GenBank/JSON IO, batch/stream/parallel API, FastAPI. | The layering and the two contracts (Constraint, Predictor) are excellent — generalize them. |

**Cross-branch anti-patterns BT4 must avoid:**

- The original splice optimizer **scraped the live ASSP web service** to detect
  splice sites — fragile, non-reproducible, network-bound. BT4 splice risk is
  **local, versioned, calibrated** — never a live scrape.
- Real **benchmark datasets existed on `almost-there` and were dropped** in the
  rebuild. BT4 keeps a committed, reproducible benchmark harness in-tree.
- Work **fragmented across parallel divergent branches** with weak merge
  discipline. BT4 runs **single-trunk** with CI gates (see §7).

---

## 3. Architecture (strict acyclic layering)

BT4 keeps BT3's best structural idea — a **pure core with heavy deps lazily
imported behind contracts** — and generalizes it. `import bt4` must stay
lightweight: **no torch, no ViennaRNA, no OR-Tools at base import.**

```
bt4/
  domain/        pure, stdlib-only. Depends on NOTHING.
                 genetic_code (translate, synonymous codons), sequence
                 validation, Result/Metrics/Violation/Severity,
                 ObjectiveVector, Frontier, OptimalityCertificate.
  biomodels/     data-driven, provenanced biological models:
      codon/       CodonUsageTable (CAI, w=f/f_max), tAITable, CodonPairTable —
                   each a provenanced dataset (source, build, CDS count, SHA-256).
      folding/     ViennaRNA adapter behind a FoldingModel contract.
      splice/      SplicePredictor contract + backends (consensus/PSSM/CNN).
      expression/  optional ExpressionPredictor head.
  objectives/    ObjectiveTerm contract: declares scope (LOCAL/PAIRWISE/
                 POSITIONAL/GLOBAL), context length, incremental delta, and a
                 whole-sequence score. CAI/tAI/CPB/GC/ramp/CpG/folding/…
  constraints/   Constraint contract (see §4). Shared IUPAC + reverse-complement
                 motif engine. Registry. One file per constraint.
  optimize/      Solver contract. Backends: exact_dp, beam_dp, cpsat_ilp,
                 lagrangian, anneal_refine. A planner picks/combines them and
                 emits an OptimalityCertificate.
  pipeline/      composes biomodels + objectives + constraints + solver into a
                 run. Owns the two-stage (exact core → non-local refinement)
                 orchestration and the Pareto sweep.
  io/            fasta, genbank, json_out (versioned self-describing schema),
                 run_manifest.
  provenance/    config_hash + full manifest (git SHA, table content hashes,
                 model SHAs, seed, tool version). Deterministic, no timestamp.
  api/           stable, print-free: optimize(), frontier(), validate(),
                 optimize_batch/stream/many. Results carry manifest + certificate.
  cli/           the ONLY layer allowed to print.
  app/           BT4 Studio — the first-class native desktop app (PySide6; see
                 §6.6). Calls `bt4.api` on a background QThread; offline-first.
                 In-tree, tested, NOT a stranded side branch (the lesson from
                 BT3's `streamlit` branch).
  service/       OPTIONAL headless HTTP API (async FastAPI): job queue, SSE
                 streaming, bounded resources, auth, OpenAPI tied to the result
                 schema version. For automation/remote use; the desktop app does
                 not depend on it.
```

**Layering rules — enforced by `import-linter` in CI, not by good intentions:**

- `domain` imports nothing.
- `biomodels` / `objectives` / `constraints` / `optimize` import only `domain`.
- `pipeline` composes them; `api` composes `pipeline`; `cli` / `service` / `ui`
  import only `api`.
- Every heavy dependency is lazily imported behind a contract **and** behind an
  optional extra: `bt4` (core, numpy), `bt4[ilp]`, `bt4[fold]`, `bt4[ml]`,
  `bt4[service]`, `bt4[dev]`.
- **Public surfaces are public.** BT3 leaked `_REGISTRY` / `_spec_to_dict`
  across cli/web/api. In BT4 no private symbol crosses a layer; registries and
  loaders are public and content-addressed from day one.

---

## 4. The contracts that hold the system together

Three contracts, generalizing BT3's two. Adding a term/constraint/model is a new
file + a registry entry — **never an engine edit.**

### 4.1 `ObjectiveTerm` (new — generalizes BT3's scattered soft penalties)

```
scope() -> LOCAL | PAIRWISE | POSITIONAL | GLOBAL
context_len() -> int                      # exact trailing context needed
delta(state, next_codon, pos) -> float    # incremental, used by DP and SA
score(dna) -> float                       # whole-sequence; MUST equal Σ deltas
```

### 4.2 `Constraint` (evolved from BT3's three-method contract)

```
scope() -> LOCAL | PAIRWISE | GLOBAL
context_len() -> int
ok_suffix(prefix, next_codon) -> bool     # hard veto; uses only context_len chars
penalty(prefix, next_codon) -> float      # soft cost (hard => 0)
validate(dna) -> Iterable[Violation]      # whole-sequence audit
relax() -> SoftConstraint                 # graceful degradation, not a dead-end
```

GLOBAL constraints additionally expose `delta_penalty(dna, edit)` for the
refinement layer. **No more `window=0, penalty=0` hacks** to encode "inert in the
DP" (that's how BT3 smuggled splice past the optimizer). Scope is declared.

### 4.3 `SplicePredictor` / `FoldingModel` / `ExpressionPredictor` (unified)

```
score_sequence(dna) -> per-position or scalar
default() -> a hash-verified registered model, else a safe baseline; NEVER crashes
```

Backends swap behind the contract; consumers never change. Models are
content-addressed and SHA-256-verified, and honestly labeled by calibration
status — `calibrated=True` only after a validation / integration-fidelity gate
passes; baselines and un-gated wrapped CNNs report `calibrated=False`. BT4's own
bundled models are pickle-free; a wrapped published backend's weights (e.g.
Pangolin's torch state dicts, SpliceAI's Keras `.h5`) are hash-verified *before*
any unpickling/load, so unverified bytes are never loaded.

### 4.4 `Solver`

```
solve(problem) -> (Result, OptimalityCertificate)
```

The certificate states one of: `proven_optimal | gap_bounded(gap) |
beam_truncated | context_capped | relaxed(terms)`. A run **cannot silently lose
optimality**.

---

## 5. Invariants — do not break (each is a test, not a sentence)

1. **Round-trip:** `translate(result.dna) == protein (+ stop)`. Property-tested
   over thousands of random proteins.
2. **Reported == computed:** every metric in `result.metrics` is recomputed from
   `result.dna` by the owning model — never trusted from an accumulator.
3. **`ok_suffix` ⇔ `validate` agreement — enforced, not prose.** For every
   registered constraint, a Hypothesis test builds feasible sequences respecting
   `ok_suffix` and asserts `validate` yields zero hard violations. A
   registry-level check asserts each hard constraint's `context_len` actually
   suffices for its own `ok_suffix`. *(BT3 documented this and then shipped a
   12-nt cap that broke it silently — this is our #1 regression to prevent.)*
4. **Delta == score:** an `ObjectiveTerm`'s accumulated `delta`s equal its
   whole-sequence `score`. *(BT3's GC penalty re-scored the same 50-nt window
   ~17× per codon — a non-additive term masquerading as additive.)*
5. **No new hard violations from refinement:** SA/refinement moves are
   synonymous-only and rejected if they raise the hard-violation count.
6. **Certificate honesty:** if it says `proven_optimal`, an independent exact
   solve on a bounded instance agrees. Any cap/prune/relaxation is reflected.
7. **Determinism:** seed everything; identical input + manifest ⇒ byte-identical
   output. End-to-end determinism job in CI.
8. **Stop-codon feasibility:** a caller-pinned/overridden stop is re-validated
   through `ok_suffix`. *(BT3 let an override bypass feasibility on the final
   codon.)*
9. **Provenance completeness:** the manifest hashes actual **table contents and
   model weights** (plus git SHA + seed), not config field names. Two runs with a
   swapped custom TSV must produce different manifests. *(BT3's `config_hash`
   over `to_dict()` gave a false "identical provenance" for a custom organism.)*

---

## 6. Scientific scope

### Objective terms (partitioned by locality)

| Term | Nature | Solved where |
|---|---|---|
| CAI (Sharp & Li) | additive per-site | exact DP |
| tAI (dos Reis) | additive per-site | exact DP |
| Codon-pair bias (CPS) | pairwise (prev codon) | DP, state extended by prev codon |
| 5′ ramp / %MinMax | positional, first ~30–50 codons | DP, position-aware weights |
| GC / GC-window | cumulative, bounded window | ILP linear constraint **or** incremental running-GC counter (never re-scored per window) |
| CpG / UpA budget | whole-sequence count | ILP budget or Lagrangian dual |
| 5′ folding ΔG | non-local, range-bounded near start | refinement layer (ViennaRNA) |
| Cryptic-splice Δrisk | non-local, whole-sequence | refinement layer (learned model) |
| Learned expression | non-local, learned | frontier reranking |

### Constraints (each a new file + registry entry)

GC, homopolymer, **max GC-run** (longest run of consecutive G/C, the "max GC
length" knob — LOCAL, exact), tandem **and reverse-complement/inverted** repeats,
**max repeat length** (any dispersed direct/inverted/palindromic repeat, anywhere,
RC-aware — genuinely GLOBAL, so refinement-enforced and honestly reported, *never*
merged into the trellis, §10.1), forbidden-motif (with named, documented
**forbidden-sequence presets**), **restriction sites (REBASE-style, auto
reverse-complement + IUPAC)**, internal-ATG → proper **Kozak context** + **uORF
pairing** (out-of-frame ATG..in-frame stop), **CpG content** (elevate for vaccines
/ deplete for stealth), **5′ folding**, splice. Hard constraints degrade
gracefully via `relax()` and report *which* constraints conflict rather than
aborting with "no feasible codon".

### ML models (real or explicitly labeled baseline)

- **Splice — wrap published SpliceAI + Pangolin (no self-training).** *Decision:*
  BT4 does **not** train its own splice CNN. A bespoke model would be both *worse*
  (SpliceAI/Pangolin were trained on the full GTEx/GENCODE corpus with serious
  compute) *and* **unvalidated** — which the honesty gate forbids shipping as
  calibrated anyway. Instead BT4 wraps the already-validated, published
  **SpliceAI** (Illumina; **PolyForm Strict 1.0.0** code + **CC BY-NC 4.0**
  weights, TensorFlow) and **Pangolin** (Zeng & Li 2022; **GPL-3.0**, PyTorch) as
  *inference-only* backends behind the existing
  `SplicePredictor` contract, feeding their dense per-nucleotide site scores into
  the already-shipped **Δsplicing** framing `P(site|designed) − P(site|reference)`
  with **top-k / log-odds pooling** (never saturating noisy-OR, §10.14). *(License
  correction — the earlier roadmap called Pangolin "MIT"; the upstream repo is in
  fact **GPL-3.0**, and **SpliceAI is more restrictive still — PolyForm Strict
  1.0.0 code + CC BY-NC 4.0 (noncommercial) weights** (its `setup.py` "GPLv3"
  string is contradicted by the authoritative LICENSE files). All are handled the
  license-clean way BT4 already handles GPL ViennaRNA: the adapter **lazily
  imports the user's own installed package and weights and never bundles or
  reimplements them**, so BT4 stays MIT. Pangolin reports **one combined
  per-position `P(splice)`** — the adapter puts it in `SpliceResult.donor` with
  `acceptor` all-zero so union-pooling counts each site once — whereas SpliceAI's
  **3-way softmax (null/acceptor/donor) maps cleanly to `donor` + `acceptor`,
  both populated**.)* Both run **out
  of the inner loop** — a ~10 kb-context CNN is far too slow to score per SA move —
  as a **final audit / frontier reranker**, exactly where the contract already
  places splice. Weights are **hash-pinned** (SHA-256, kept out of git,
  content-hash in the manifest), so — unlike the BT3 ASSP scrape — local inference
  stays **reproducible-from-manifest**. Running *both* and reporting their
  **agreement/disagreement** is a first-class uncertainty signal (§8), not
  redundancy. `calibrated=True` is set only after an **integration-fidelity
  check** — the adapter reproduces the published model's scores on a panel of
  known real sites / non-sites — *not* a from-scratch held-out-chromosome training
  gate. The heavy PyTorch / TF deps live behind optional extras
  (`bt4[splice-pangolin]` / `bt4[splice-spliceai]`, lazily imported like
  ViennaRNA). **Status: both wrapped CNN backends have landed** —
  `PangolinSplicePredictor` (`biomodels/splice/pangolin.py`) and
  `SpliceAiSplicePredictor` (`biomodels/splice/spliceai.py`) — with the
  two-backend agreement harness (`biomodels/splice/agreement.py`,
  `scripts/compare_splice_backends.py`) that makes agreement between two real,
  independently-trained CNNs reachable. Each reproduces its upstream model's
  scores **bit-for-bit** (verified against the published weights) yet ships
  **`calibrated=False`** — no reference panel is bundled (capturing one needs the
  licensed weights and reproduces licensed outputs), so `default()` keeps
  returning the PWM baseline and a maintainer promotes to calibrated only after
  recording the fidelity gate. **Honest scope:** these predict splice-*site
  presence*, and a lower Δ means lower *predicted cryptic-splice risk* — a strong
  prior, but not the same as validated expression gain (the same CAI-as-weak-proxy
  caution); SpliceAI's CC BY-NC weights additionally make that backend
  noncommercial-only.
- **Folding — ViennaRNA** MFE / partition-function ΔG (real thermodynamics, not
  a hand-weighted proxy); objective *and* constraint (avoid RBS/Kozak-occluding
  hairpins).
- **Expression — optional learned head** behind the predictor contract, trained
  on massively-parallel expression / ribosome-load data; used to rerank the
  frontier with reported calibration and uncertainty. Clearly labeled by
  training provenance.
- **ASSP — optional online cross-check backend** (Alternative Splice Site
  Predictor). This is a *kept, supported feature*, behind the same
  `SplicePredictor` contract as every other backend, but with guardrails that
  make it honest rather than the BT3 liability it was:
  - **Opt-in and out-of-the-inner-loop.** ASSP is a network service, so it is
    never used to score per-move inside the optimizer. It runs as a **final
    audit / validation pass** on a delivered sequence (`bt4 validate
    --splice-backend assp`, or `--check-splice assp` on `optimize`) and as a
    **ground-truth comparator** in the eval/benchmark harness.
  - **Never silent, never blocking.** Gated behind an explicit flag *and* the
    `bt4[assp]` extra; polite rate-limiting + backoff; responses cached by
    sequence hash; if the service is unreachable it degrades gracefully and says
    so — it can never fail an optimization. CI uses stored offline fixtures, no
    live calls.
  - **Labeled non-reproducible.** Any ASSP-derived number is stamped
    network-derived and excluded from the reproducible-from-manifest guarantee
    (the local calibrated model remains the default, reproducible splice path).

Keep BT3's excellent honesty machinery and **feed it real biology**: shared
`encode_window` encoder (no train/serve skew), chromosome-grouped splits, PR-AUC
/ MCC / ECE / Brier (never bare accuracy), acceptance gates vs strawmen. **Until
a model passes its gate on real held-out data, `--refine` loudly labels its
output consensus-baseline only and does not claim biological calibration.**

### 6.6 The app — BT4 Studio (a first-class native desktop app)

BT4 ships a **real, usable desktop application**, not just a CLI and a raw HTTP
API. BT3's UI lived and died on a `streamlit` side branch; BT4 Studio is an
in-tree, tested, supported surface. It is **offline-first and native**, which
fits BT4's reproducibility ethos (no server to stand up, results computed
locally, nothing leaves the machine unless the user explicitly clicks the ASSP
cross-check).

- **`app/` = the desktop app, built in PySide6 (Qt for Python).** Native,
  cross-platform (macOS/Windows/Linux), packaged into installers via
  PyInstaller/Briefcase. It calls the stable **`bt4.api`** directly and runs each
  optimization on a **background `QThread`** so the UI never blocks; it never
  imports `optimize`/`biomodels` directly (goes through `api/`). Interactive
  plots via **pyqtgraph** (frontier scatter, per-site risk tracks). Light/dark
  Qt theming.
- **`service/` = an *optional* headless HTTP API (async FastAPI)** for
  automation, remote/batch use, and future non-desktop frontends — job queue,
  SSE streaming, bounded resources, auth. The desktop app does **not** depend on
  it; it exists so the same engine is scriptable over HTTP.

What the app does (design targets):

- **Input:** paste a protein or upload FASTA; pick target organism; name the job.
- **Design controls:** toggle each constraint with real controls — GC target,
  max homopolymer, **max GC length** (GC-run) and **max repeat length**
  spinboxes, CpG mode (elevate/deplete), %MinMax, tandem/hairpin repeats,
  internal-ATG, tAI, forbidden motifs, a **checkbox per named forbidden-sequence
  preset**, and the **restriction-enzyme picker** (REBASE catalog); choose
  objective **weights** *or* **ε-budgets** per term. **Every control carries a
  hover tooltip** explaining what its variable does.
- **Run & watch:** live progress from the background worker; the
  **optimality-certificate badge** (proven-optimal vs gap-bounded vs
  beam-truncated) shown honestly, not hidden.
- **Visualize:** an **interactive Pareto-frontier plot** with the delivered
  point marked and the trade-off it made spelled out; **per-site risk tracks**
  along the sequence (splice / 5′-folding / GC / CpG heatmaps); a metrics table;
  a sequence viewer with constraint annotations.
- **Cross-check:** a one-click **"validate with ASSP"** button (§6, opt-in,
  clearly labeled network-derived) on the delivered sequence.
- **Export:** FASTA / GenBank / JSON **plus the run manifest**, so anything the
  app shows is reproducible from its stamp (except explicitly network-derived
  ASSP numbers).

Accessibility, light/dark theming, and responsive layout are in scope — "nice"
is a requirement, not a nice-to-have.

---

## 7. Engineering

- **Stack:** Python 3.10+ for orchestration/API; a compiled **Rust core (PyO3 /
  maturin)** for the trellis inner loop and incremental scorers, with a
  pure-numpy fallback when the extension isn't built. Rust over Cython for
  memory safety, GIL-free DP, and rayon parallelism (Pareto sweep + batch).
  ILP/CP-SAT via OR-Tools; folding via ViennaRNA bindings; ML via torch (lazy).
- **Performance:** everything **incremental** — each refinement move is
  O(context), not O(L). *(BT3's SA recomputed full CAI + full validation per
  move — quadratic.)* Length-scaled iteration counts, block/segment moves,
  adaptive / parallel-tempering schedules (block moves and tempering also widen
  the refinement's *reach* — escaping single-codon barriers such as coordinated
  multi-codon repeat removal, §9 Phase 3 — not just throughput). **Runtime and
  peak-memory scaling are
  regression-tested in CI** with an asserted curve and a wall-clock ceiling —
  BT3 had no performance test despite runtime being BT2's original weakness.
- **Testing:** Hypothesis property tests for every §5 invariant; **golden tests**
  (fixed protein panel → pinned DNA + metrics, regenerable by one command);
  **optimality tests** (beam/relaxed must match exact DP on small instances; the
  certificate must not lie); the committed **benchmark harness**
  (`scripts/benchmark.py`) comparing BT4 vs input vs GeneOptimizer/IDT/Twist on
  CAI/tAI/GC/CpG/repeats/splice/folding over a pinned panel (revive
  `almost-there`'s benchmark corpus).
- **CI from commit #1** (BT3 had none): GitHub Actions matrix (3.10/3.11/3.12/3.13)
  running ruff, `mypy --strict`, import-linter (layering), pytest + coverage
  gate, determinism job, and per-extra jobs (`[ml]`/`[fold]`/`[ilp]`). Merge
  blocked on failure. `abi3` wheels for the Rust core across platforms.
- **Provenance & packaging:** every result emits a **run manifest** (config hash,
  table provenance with SHA-256, model SHAs, solver certificate, seed, git
  commit, tool version) — reproducible from the stamp alone. Single-sourced
  version surfaced via `--version`. A real license (not "TBD"), a CHANGELOG,
  published wheel + sdist on tagged releases, model weights as hash-referenced
  release artifacts **out of git**.

---

## 8. Data honesty & validation

- Codon / tAI / CPB tables are **provenanced datasets** with sidecar manifests
  (source DB, accession, genome build, CDS count, extraction date, SHA-256); the
  table's *content hash* enters the provenance stamp.
- Reproducible **table-build pipeline**: `bt4 build-table --cds genes.fasta
  --organism X [--highly-expressed ref]` recomputes fractions from a declared CDS
  set. A test asserts a bundled table reproduces a *published* CAI for a
  benchmark gene set. *(BT3 asserted "representative" with no such check.)*
- **External ground truth, not just self-consistency.** Beyond round-trip and
  reported==computed, add a held-out set of real highly-expressed genes and test
  that BT4's output distributions (codon/GC/CpG) fall within the natural range
  rather than degenerating to a GC-skewed CAI=1 extreme. Committed, regenerable
  comparison tables vs published tools.
- **Uncertainty propagation:** report sensitivity of CAI/GC to table choice and
  solver budget; conformal intervals on model predictions. No point estimate
  presented as ground truth.

---

## 9. Phased roadmap

- **Phase 0 — Foundations & honesty scaffolding.** ✅ **Done.** Repo, strict
  layering + import-linter, CI (lint/type/test/determinism), pure `domain` with
  `ObjectiveVector` / `Frontier` / `OptimalityCertificate`, content-hashed
  provenance manifest, single-sourced version, license. Port BT3's genetic code +
  round-trip property test.
- **Phase 1 — Exact, honest single-objective core.** ✅ **Done.** Exact codon-
  trellis DP with **true per-constraint context (no global cap)**, exact DP +
  beam as an explicit knob, certificate emission, and the `ok_suffix⇔validate`
  and `delta==score` property tests all shipped. The trellis currently runs in
  pure Python with the Rust `bt4_native` primitives (`gc_count`,
  `max_homopolymer_run`, `reverse_complement`, `max_gc_run`, and `longest_repeat`)
  available and CI-built — each with a byte-identical pure-Python fallback and an
  equivalence property test. `max_gc_run` backs the GC-run `ok_suffix` veto on its
  bounded window; `longest_repeat` (a reverse-complement-aware longest-repeat
  query) is cross-checked against `MaxRepeatConstraint` but **deliberately kept
  off** the per-SA-move `validate` hot path, because whole-sequence O(n²)
  longest-repeat is slower there than the constraint's existing O(n·k) k-mer scan
  when the native extension is absent (§7). The **full trellis inner loop is now
  ported to Rust** as `bt4_native.trellis_solve` (with a byte-identical pure-Python
  twin `_py_trellis_solve` and an equivalence property test, exactly like the other
  primitives). Because the DP is callback-driven (`scalar_delta` + each
  constraint's `ok_suffix`) and Rust must never call back into Python, the port
  carries an honest **regime gate**: it runs only when the objective is
  position-independent (no `POSITIONAL` term — the codon-pair term was made
  context-based so PAIRWISE stays position-independent too), Python **precomputes**
  the reachable-context transition graph and the pre-summed per-transition deltas
  (so the float summation order — and thus the lexicographic tie-break — stay
  bit-for-bit identical), and the layer DP runs in Rust; the code **falls back to
  the pure-Python DP** whenever the regime does not hold, the extension is absent,
  or a context-count cap is exceeded. Because the Python-side precompute costs
  about as much as the whole pure DP, a *single* solve is not sped up (and
  `run_optimize` deliberately stays on the pure path); the win is the **Pareto
  frontier**, where the transition graph is built **once** and reused across every
  scalarization grid point (only the cheap deltas are recomputed) with the DP in
  Rust — a measured ~2.7–5.5x frontier speedup, byte-identical output and
  certificates. *This alone already beats BT3* (honest optimality, correct
  incremental GC).
- **Phase 2 — Multi-objective, richer biology & first app.** ✅ **Complete.**
  Delivered: the **multi-objective Pareto-frontier API** (a unit-simplex
  scalarization sweep over *every* active objective axis - CAI and GC always,
  plus ramp/CpG/%MinMax when weighted - so 3+ objectives trace a real trade-off
  surface, with the CAI/GC two-objective case unchanged; each point stays an
  exact proven-optimal solve, the frontier a bounded sample of the surface); the
  **BT4 Studio desktop MVP
  (PySide6)** (§6.6) calling `bt4.api` on a background thread, with the frontier
  plot and honest certificate badge; **restriction-site constraints** (REBASE-
  style catalog, IUPAC-aware, auto reverse-complement); **more organisms**
  (representative *E. coli* and *S. cerevisiae* tables, auto-discovered);
  **`bt4 build-table`** to recompute an authentic codon table from a user CDS
  FASTA; a committed **benchmark harness** (`scripts/benchmark.py`); the optional
  **`service/` FastAPI HTTP API**; a **5′ translation-ramp** term; **codon-pair
  bias** (`CpbTerm`, built from a reference CDS) exact in the trellis via the
  extended-state DP (objective context); a **CpG/dinucleotide** term
  (deplete/elevate); a first **OR-Tools CP-SAT backend** that solves the
  additive objective under a **global GC budget** with a gap certificate;
  **tandem- and inverted-repeat constraints** (`constraints/repeats.py`, LOCAL,
  `ok_suffix⇔validate`-tested); a **%MinMax** objective — an additive
  codon-commonness DP term kept honestly separate from the true non-additive
  sliding-window %MinMax reporting metric (`min_max_profile`); and a
  **Lagrangian-relaxation backend** (`optimize/lagrangian.py`) that dualizes a
  global GC budget into the *exact DP* — so unlike CP-SAT it keeps **local
  constraints and pairwise objective terms** honored under the budget — reporting
  an honest gap-bounded certificate from its subgradient dual bound. The pipeline
  now routes a GC budget to CP-SAT for the pure-additive case and to the
  Lagrangian backend whenever a local constraint or pairwise term is present. The
  **CpG/UpA whole-sequence *count* budgets** — the last remaining Phase 2 item —
  have now **landed** (`dinuc_budget`/`dinuc_min`/`dinuc_max`): a dinucleotide
  count does *not* decompose per-codon (a 2-mer straddles the codon boundary), so
  the amount-bucketed budget DP takes a **context-aware** per-codon `amount` that
  attributes each occurrence to the codon holding its END base (exactly as
  `DinucleotideTerm.delta` does), with `budget_context=1` folded into the state so
  a straddling count stays exact — enforced by the same **exact bucketed DP** as
  the GC budget with a `PROVEN_OPTIMAL` certificate, mutually exclusive with the
  GC budget and (like it) incompatible with refinement-enforced rules. Also
  landed: an **internal-ATG strong-Kozak constraint**
  (`constraints/kozak.py`, LOCAL, `ok_suffix⇔validate`-tested) that forbids
  internal `ATG` in a strong Kozak context (purine at -3, G at +4) — with the
  genuinely non-local **uORF pairing** part *honestly deferred* to Phase 3, not
  faked with a padded window; an **enriched benchmark harness** reporting CpG,
  tandem-repeat, and %MinMax alongside GC/homopolymer/CAI (still naive-vs-BT4, no
  fabricated competitor numbers); and BT4's first **performance/scaling
  regression test** (`tests/test_performance.py`, §7) asserting sub-quadratic
  exact-DP runtime under a wall-clock ceiling. **`CpbTerm` (codon-pair bias) is
  now wired** (`cpb_weight` + a user-supplied `cpb_reference_cds`; the pipeline
  builds the `CodonPairTable` at run time, PAIRWISE so exact in the trellis and a
  frontier axis, its reference-CDS content hash in the manifest) - honestly
  refusing when no reference CDS is given, since no default codon-pair table is
  bundled (§8). Remaining (Phase 2): none — the CpG/UpA whole-sequence *count*
  budgets, previously the last deferred item, have now shipped as an exact
  amount-bucketed DP budget (see the dinucleotide-budget note above). *(tAI and
  uORF pairing have also landed - see their bullets.)* The **published comparison vs
  GeneOptimizer/IDT/Twist** has landed as `scripts/compare_tools.py` over a real,
  cited, CC BY 4.0 panel (Ranaghan et al. 2021, KRas4B) - every metric recomputed
  by BT4's own functions, BT4 never claimed "better", each tool's output
  attributed and the ATUM/DNA2.0 truncation flagged. Also shipped: a
  `scripts/sensitivity.py` uncertainty report (CAI/GC/tAI spread across codon-
  table and solver-budget choices, §8) and a windowed CpG/dinucleotide reporting
  profile (`objectives/dinuc_profile.py`, honestly separate from the additive
  `DinucleotideTerm`, mirroring the %MinMax split). A companion
  `scripts/compare_reproducibility.py` adds the **run-to-run variability** view
  over Ranaghan Table 4 (three proteins x three *anonymized* algorithms x ten
  repeat runs, CC BY 4.0) - kept deliberately separate from the named-tool board
  (anonymized tools, repeat runs = a determinism axis, not a scoreboard), with
  BT4 shown as a zero-spread deterministic reference. Also landed (BT3
  `almost-there` parity): a **max GC-run** constraint (`constraints/gc_run.py`,
  the "max GC length" - LOCAL, exact in the trellis, `ok_suffix⇔validate`-tested);
  a **max repeat length** constraint (`constraints/max_repeat.py`) that bans any
  dispersed direct/inverted/palindromic repeat longer than the limit anywhere in
  the sequence (RC-aware, `k = max_length + 1` sufficiency proof) - genuinely
  GLOBAL, so it is **never merged into the exact DP** (that would silently
  over-merge, §10.1); instead the exact-DP seed is polished by the SA refinement
  engine (extended with a `global_constraints` gate that re-counts whole-sequence
  hard violations per move and never lets the count rise, invariant #5), with any
  residual repeats reported honestly and the certificate degraded to heuristic
  only when refinement actually runs (single-codon SA cannot always clear a
  repeat — see the Phase 3 refinement note on coordinated multi-codon moves and
  immovable-codon feasibility floors; such residuals are disclosed, never
  hidden); and named, documented **forbidden-sequence
  presets** (`constraints/forbidden.py`, e.g. poly-A signals, TATA box, telomere
  repeat, BT3 synthesis artifacts). All three are wired through `OptimizeConfig`,
  the `bt4` CLI (`--max-gc-run`, `--max-repeat-length`, `--forbid-preset`, `bt4
  presets`), the `service` schema, and **BT4 Studio** - which now also renders a
  **hover tooltip on every design control** explaining what each variable does,
  the two new spinboxes, and a checkbox per forbidden preset.
- **Phase 3 — Non-local models & refinement done right.** 🔶 **Groundwork
  landed.** The **`FoldingModel` contract** (`biomodels/folding/`) ships with a
  lazy **ViennaRNA** backend (`calibrated=True`, behind the `bt4[fold]` extra) and
  an honestly-labeled uncalibrated **baseline** (`calibrated=False`, a Nussinov
  base-pair proxy in *arbitrary units* — never presented as real ΔG), and a
  `default()` that never crashes. The **incremental SA refinement engine**
  (`optimize/anneal_refine.py`) does synonymous-only single-codon moves with an
  O(context) feasibility check that **provably never raises the hard-violation
  count** (invariant #5, brute-forced over 56k swaps), a deterministic seeded
  trajectory (#7), and an honest **HEURISTIC** certificate. A **`--refine`** path
  wires the folding ΔG into an SA pass over the exact-DP seed, loudly flagging
  `folding_calibrated=False` when the baseline is used (per §6/§10.6).
  **Reaching farther without weakening #5 — block moves + parallel tempering have
  landed.** `anneal_refine` now takes opt-in `block_size`/`block_prob` (coordinated
  multi-position synonymous swaps) and `replicas`/`temps`/`swap_every` (a
  parallel-tempering replica ladder with standard replica-exchange Metropolis
  swaps). Block moves make a barrier that only clears when two codons move
  *together* (a dispersed max-repeat / uORF) reachable in one step, and hot
  replicas can accept uphill moves and swap their already-feasible configuration
  into the cold chain — so the search crosses barriers the single-codon chain
  could not. Invariant #5 is preserved *structurally*: block candidates pass the
  same local (union-of-windows `ok_suffix`) and global (whole-sequence recount)
  gates, every replica gates against **its own** current count, every visited
  configuration keeps a global count `<=` the seed's, and the delivered result is
  ranked **lower-global-count-first then higher-score** (`_better`). All four knobs
  default off, and with them off the engine reproduces the prior single-chain
  trajectory **byte-for-byte** (#7). Block moves **always full-`score` re-score**
  (never `delta_score`) — summing per-position deltas is valid only for additive
  disjoint-context terms, which folding/splice/expression are not.
  **Feasibility floor (honestly disclosed, not hidden):** a repeat pinned to
  synonymously-**immovable** bases — e.g. Met `ATG` / Trp `TGG`, or a degenerate
  position all synonymous codons share — is unremovable by *any* synonymous scheme,
  block move or hot replica alike, and is a genuine feasibility floor, not an engine
  defect. Every such case is reported as a disclosed residual
  (`max_repeat_residual`, enforcement `"partial"`), never a silent "clean" claim.
  **Both wrapped
  published splice backends have landed**: `PangolinSplicePredictor`
  (`biomodels/splice/pangolin.py`, Pangolin — GPL-3.0, PyTorch, one combined
  `P(splice)` → `donor`) and `SpliceAiSplicePredictor`
  (`biomodels/splice/spliceai.py`, SpliceAI — PolyForm Strict code + CC BY-NC
  weights, TensorFlow, 3-way softmax → `donor` + `acceptor`). Each drives the
  user's own installed package (*not* bundled or reimplemented; lazily imported
  like GPL ViennaRNA), scored out-of-loop as an audit/reranker, with its
  **weights hash-pinned** (SHA-256 of the published files, verified *before* they
  are loaded) so runs stay reproducible-from-manifest. Each reproduces its
  upstream model's per-nucleotide scores **bit-for-bit** yet ships
  `calibrated=False` (no reference panel is bundled — the integration-fidelity
  gates `verify_pangolin_fidelity` / `verify_spliceai_fidelity` are the promotion
  path), so `default()` keeps returning the PWM baseline. The **two-backend
  agreement harness** (`backend_agreement` + `scripts/compare_splice_backends.py`)
  reports pairwise Spearman rank / sign agreement across whichever backends are
  available — with both CNNs installed, agreement between two real,
  independently-trained models — the first-class uncertainty signal of §6. The
  **license-clean fidelity-attestation layer** (`biomodels/splice/attestation.py`)
  now provides the honest promotion path: a `FidelityAttestation` records **only**
  a passing gate's derived scalars (`passed`, `max_abs_deviation`, `n_cases`,
  `tolerance`) plus the public pinned weight SHA-256s and the tool version —
  **never** a `FidelityCase` raw per-position score (those *are* the license-
  encumbered model outputs), a shape enforced structurally (`_ALLOWED_FIELDS` +
  a "no raw-score field" test). It layers four ways: a **committed** attestation
  (re-verifiable by anyone holding the same weights, with a deterministic
  `content_hash` for the manifest), **private execution** where the weights live,
  a **user opt-in** via `verified_predictor(predictor, attestation)` (which flips
  `calibrated=True` only when the attestation passed, clears the
  `MAX_ATTESTATION_TOLERANCE` floor, and its weight SHAs exactly match the
  adapter's `PINNED_WEIGHT_SHA256` — a refusal, never a silent downgrade), and the
  baseline **fallback**. Because BT4 is open-source and **non-commercial**, both
  Pangolin (GPL) and SpliceAI (CC BY-NC) are eligible to certify. Still ahead:
  running an actual gate to emit a committed attestation (human-only — needs the
  licensed weights + a captured panel) and the opt-in **ASSP** cross-check (§6)
  with offline fixtures (block/segment moves and parallel tempering have now
  landed — see the refinement note above).
  The **per-site risk tracks** now ship as honest reporting profiles through
  `api.tracks()` and `bt4 tracks` (sliding-window GC / CpG density / %MinMax,
  each recomputed from the sequence, never fed to the solver) **and are plotted
  in BT4 Studio** as per-site tracks. BT4 Studio's **sequence viewer now renders
  inline violation annotations**: each `Violation` from the delivered result's
  whole-sequence audit is highlighted over its `[start, end)` span (red for HARD
  / feasibility, amber for SOFT / quality) with a hover tooltip naming the
  constraint, severity, span, and detail, plus a HARD/SOFT legend — so a
  *residual GLOBAL violation* left after refinement (e.g. a dispersed max-repeat
  or uORF) is shown **where it occurs**, not just as a metrics-table count. The
  highlights are Qt extra-selections layered over the text, so the exported
  sequence stays exactly the delivered one. The
  **`SplicePredictor` contract**
  (`biomodels/splice/`) has now landed with an honestly-labeled uncalibrated
  **consensus/PWM baseline** (`calibrated=False`, top-k/log-odds Δsplicing pooling
  -- never noisy-OR, §10.14) and a `default()` that never crashes; the wrapped
  **Pangolin** and **SpliceAI** CNN backends (above) now slot behind exactly this
  contract, and until a backend passes its gate the baseline remains the
  never-calibrated default. The **uORF-pairing
  constraint** (`constraints/uorf.py`, `avoid_uorf`/`uorf_region_nt`) has now
  landed as the genuinely non-local half of internal-ATG handling: an
  out-of-frame internal ATG paired with a downstream in-its-frame stop is a short
  uORF. It is `Scope.GLOBAL` (ATG and stop arbitrarily far apart), so it is
  refinement-enforced through the same `anneal_refine` global-constraint gate as
  max-repeat (drives the count down, never raises it, invariant #5), **never
  merged into the exact DP** (§10.1), with residual uORFs reported honestly - a
  purely *structural* rule complementing the LOCAL strong-Kozak
  `InternalStartConstraint`, making **no** calibrated-expression claim. Also
  landed: the **`ExpressionPredictor` contract scaffold** (`biomodels/expression/`)
  for the Phase 4 learned head. Expression is non-local/learned, so it is not an
  `ObjectiveTerm` and never runs in the optimizer loop; instead a validated head
  will **rerank the frontier** as a post-solve pass
  (`pipeline/rerank.py::rerank_by_expression`, exposed via `bt4.api`). `default()`
  returns a **neutral placeholder** (`NullExpressionModel`, `calibrated=False`,
  every score `0.0`) because expression has no structural anchor the way
  folding/splice do - a hand-weighted CAI+GC+ΔG composite dressed as "expression"
  would be the §10.5/§10.6 trap. The rerank hook **only re-picks the delivered
  point when the predictor is calibrated**; with the placeholder it is a pure
  reporting no-op (an uncalibrated score never steers delivery). The
  **model-agnostic acceptance gate** a head must pass to *earn* `calibrated=True`
  has now landed (`biomodels/expression/gate.py`): for a **log-TE regression**
  head it reports Spearman (primary) / Pearson / R² plus **split-conformal
  coverage** on a **group-disjoint split** (homology/chromosome), so a head
  validated only on natural-gene TE cannot claim calibration for the CDS-variant
  regime BT4 optimizes; `passed` needs both the Spearman threshold and coverage
  near target, thresholds are inputs set at gate time, and the `NullExpressionModel`
  provably cannot pass. The shared estimators live in `biomodels/_stats.py`
  (`pearson`/`spearman`/`r2_score`/`conformal_quantile`/`empirical_coverage`),
  reused by the splice agreement report. **Still human-only:** obtaining a
  license-clean, regime-matched CDS-variant panel (e.g. wrapping the published
  RiboNN log-TE CNN, Sanofi non-commercial — eligible under BT4's non-commercial
  scope, handled non-vendored/hash-pinned like SpliceAI) and running the gate; the
  `calibrated` flip is earned on data, never assigned. **The wrapped RiboNN adapter
  has now landed** (`biomodels/expression/ribonn.py`): `RiboNNExpressionModel`
  drives the user's own RiboNN checkout (Sanofi **non-commercial** code + Zenodo
  weights — *not* vendored; lazily imports the repo's `src`, pointed at via
  `$BT4_RIBONN_DIR`), verifying every weight it loads against a bundled 180-entry
  SHA-256 manifest (`data/ribonn_sha256.json`, 90 human + 90 mouse — public content
  hashes only) *before* `torch.load`. It scores in RiboNN's native **CLR-residual
  TE** units (never exponentiated — the model applies no log/exp itself, per the
  paper's centered-log-ratio target) and exposes `delta_logte(designed, reference)`
  — the UTR-fixed, CDS-attributable Δ that encodes the *avoid-expression-limiting-
  sequences* framing (a negative Δ flags a CDS change predicted to reduce
  expression), analogous to Pangolin's `delta_splicing`. It ships **`calibrated=False`**
  (`default()` still returns the neutral placeholder): reproducing RiboNN faithfully
  is not calibration *for BT4's regime* (RiboNN's own ablation puts only ~31% of
  per-nt signal in the CDS), so promotion needs a passing `verify_expression_gate`
  on a CDS-variant panel — human-only, data-gated. New `bt4[expression-ribonn]`
  extra (torch + pandas, lazily imported so `import bt4` stays light). **First real
  end-to-end runs against the licensed weights have now happened** (on a
  maintainer's machine, human-run — the weights are non-commercial and never
  bundled/CI-run), validating the adapter and surfacing two integration bugs that
  only appear once the live forward pass executes, both since fixed: (1) RiboNN
  returns its ensemble as **one row per cross-validation model**, so the per-input
  realignment must **group by `tx_id` and average** (mean over cell types *and* the
  ensemble) — a plain `set_index` left duplicate labels and `float(Series)` raised;
  (2) scoring now **requires non-empty `utr5`/`utr3`** and refuses empty ones up
  front with a clear message (RiboNN's loader reads an all-empty UTR column as NaN
  and its `.str` preprocessing crashes — and the UTRs carry most of RiboNN's signal,
  so an empty-UTR score is not meaningful anyway). Both are property-tested against a
  synthetic RiboNN output table; `calibrated` stays `False` (unchanged).
- **tAI — landed (real data).** The deferred tAI item is now shipped honestly:
  `biomodels/codon/tai.py` builds relative adaptiveness from **real human tRNA
  gene copy numbers** (GtRNAdb hg38, 431 genes/47 anticodons, bundled with a
  content-SHA-256 provenance sidecar and independently re-counted) via a faithful
  port of the dos Reis (2004) `get.ws` wobble model (verified s-values; the
  bacterial lysidine path gated on `sking=1`). `TaiTerm` is the exact-DP
  `tai_logw` objective (additive, `delta==score`), wired through `tai_weight` in
  the pipeline/CLI/app and as a frontier axis, with the tRNA table's hash entering
  the provenance stamp when tAI is active. Real tRNA tables for **eight
  organisms** now ship (human, mouse, rat, zebrafish, *Drosophila*, *C. elegans*,
  *Arabidopsis*, and *S. cerevisiae* -- all from GtRNAdb, independently
  re-counted, content-hashed, and stamped citation-gated academic use, not
  CC/public-domain); organisms without bundled tRNA data raise -- no fabricated
  tables.
- **Phase 4 — Learned expression & polished app.** Expression predictor head,
  frontier reranking with calibration/uncertainty. Polish BT4 Studio (theming,
  accessibility) and ship **packaged installers** (PyInstaller/Briefcase) for
  macOS/Windows/Linux; optionally expose the `service/` HTTP API. External-
  validation report vs real gene distributions and published tools.
- **Phase 5 — Scale & ecosystem.** 🔶 **Opened.** **Library / degenerate-design
  mode has landed** (`optimize/sample.py` + `pipeline/library.py`, exposed as
  `api.library` and `bt4 library PROTEIN --n N`): instead of a single MFC
  optimum it draws a library of sequences by **sampling** each residue's
  synonymous-codon distribution (organism frequencies raised to `1/temperature`,
  keeping only codons that pass every LOCAL constraint's `ok_suffix`). It is an
  honest **stochastic sampler, not an optimizer** — every member round-trips
  (#1), carries metrics recomputed from its own DNA (#2), is fully deterministic
  from its seed (#7), and carries the new `OptimalityStatus.SAMPLED` certificate
  that makes **no** optimality or expression claim (§1/§10.6). GLOBAL rules
  (max-repeat, uORF) are not enforced during sampling but are validated and any
  residual violation reported honestly per member. Still ahead: more organisms
  with authoritative provenance, restriction-enzyme catalogs, tissue/condition-
  specific tables.

---

## 10. BT3 mistakes BT4 must NOT repeat

1. **No silent optimality loss.** BT3's `_MAX_CONTEXT=12` over-merged prefixes
   differing at bases 13–50 (Repeat window 18, GC window 50) and beam pruning
   dropped optima — neither flagged. → true per-constraint context; every
   cap/prune/relaxation in the certificate.
2. **No prose-only invariants.** The load-bearing `ok_suffix⇔validate` agreement
   was documented, not tested, and the cap violated it. → property-tested, CI
   fails otherwise.
3. **No non-additive term masquerading as additive.** → enforce `delta == score`
   (fix the GC re-scoring).
4. **No inert-constraint hacks** (`window=0, penalty=0` for splice). → declare
   scope; non-local terms live in the refinement layer with real deltas.
5. **No single magic-weighted scalar** (`splice_weight=5.0`). → multi-objective,
   returned frontier, calibrated/learned weights.
6. **No placeholder model presented as a feature.** BT3's splice model was
   consensus/synthetic — `--refine` optimized noise. → ship validated models or
   loudly refuse to claim calibration.
7. **No CAI-as-truth.** CAI is a weak proxy. → validate vs real expression /
   functional data; add tAI/CPB/folding/ramp/expression.
8. **No quadratic refinement.** → incremental O(context) per move.
9. **No private symbols across layers** (`_REGISTRY`/`_spec_to_dict`). → public,
   content-addressed registries and loaders.
10. **No provenance that lies.** BT3 hashed config field *names*, so a swapped
    TSV read as identical provenance. → hash actual table/model contents + git SHA.
11. **No missing CI, no stale docs, no legacy cruft.** BT3 had no CI, a BT2-era
    `DIRECTORY_STRUCTURE.md`, no-op legacy scripts, a "TBD" license, and no
    `--version`. → CI day one, docs synced in CI, no legacy scripts, licensed,
    single-sourced version, benchmark + determinism gates.
12. **No unbounded, unguarded service.** BT3's FastAPI was sync, unauth'd,
    unbounded. → async, queued, bounded, authenticated, streaming.
13. **No caller override bypassing feasibility.** → re-validate pinned/override
    stop codons through `ok_suffix`.
14. **No saturating noisy-OR splice aggregation.** → top-k / log-odds pooling on
    a real per-nucleotide model.
15. **No network service in the optimization inner loop, and no *unversioned,
    fallback-less* dependence on one.** BT3's sin was making the live ASSP scrape
    the *only* splice path, scored in-loop, with no offline model behind it. →
    the default splice path is a local, versioned, calibrated, offline-
    reproducible model; **ASSP stays as an opt-in, out-of-loop cross-check /
    validator** (see §6) — network-labeled, cached, non-blocking, never in the
    per-move scorer.
16. **No fragmentation across divergent branches.** → single-trunk development
    with CI gates; keep the UI and benchmark corpus in-tree, not stranded on
    side branches.

---

**Bottom line.** BT4 keeps BT3's one correct, load-bearing idea — constrained
combinatorial optimization over a codon trellis — and rebuilds everything around
it to be *honest* (auditable optimality, content-hashed provenance, enforced
invariants), *multi-objective* (a Pareto frontier over calibrated,
expression-relevant terms), and *real* (shipped, validated ML for splice /
folding / expression, benchmarked against published tools). BT3 optimized an
unvalidated proxy with a solver that quietly cheated and models that were never
biology. BT4 optimizes a validated objective vector with a solver that proves its
work and models that pass acceptance gates on real data — or refuses to claim
otherwise.
