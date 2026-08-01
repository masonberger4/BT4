# CLAUDE.md — BT4

Guidance for Claude Code (and humans) building **BT4**, the from-scratch
successor to BT3. This file is the constitution of the new repository: read it
before writing code, and keep it current as the architecture evolves.

> Status: **Phase 0 complete; Phase 1 done and early Phase 2 landed.** The pure
> `domain` layer, provenance manifest, packaging, layering contract, and CI are
> in place, and on top of them an **honest exact-DP core** now ships: a codon
> trellis with true per-constraint context and a real optimality certificate,
> the `CaiTerm`/`GcProximityTerm` objectives and `Homopolymer`/`ForbiddenMotif`
> constraints (with their `delta==score` and `ok_suffix⇔validate` property
> tests), a **CAI/GC Pareto frontier**, the stable `bt4.api`, a `bt4` CLI, and
> the first cut of **BT4 Studio** (the PySide6 desktop app). Phase 2 has since
> added codon-pair/ramp/CpG/%MinMax objectives, tandem/inverted-repeat
> constraints, and **two budget backends** — OR-Tools CP-SAT and an honest
> **Lagrangian relaxation** that keeps local constraints under a global GC budget.
> Still ahead: tAI, and the
> validated splice/folding/expression models — see §9. This document was written
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
  hash-pinned, calibrated models (SpliceAI/Pangolin-class Δsplicing, ViennaRNA
  folding, an optional learned expression head), gated on real held-out data.

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

Backends swap behind the contract; consumers never change. All models are
content-addressed, SHA-256-verified, pickle-free, and calibrated.

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

GC, homopolymer, tandem **and reverse-complement/inverted** repeats,
forbidden-motif, **restriction sites (REBASE-style, auto reverse-complement +
IUPAC)**, internal-ATG → proper **Kozak context** + **uORF pairing** (out-of-frame
ATG..in-frame stop), **CpG content** (elevate for vaccines / deplete for stealth),
**5′ folding**, splice. Hard constraints degrade gracefully via `relax()` and
report *which* constraints conflict rather than aborting with "no feasible codon".

### ML models (real or explicitly labeled baseline)

- **Splice — SpliceAI/Pangolin-class.** Wide-context (~10 kb) dilated-residual
  CNN, dense per-nucleotide donor/acceptor/neither. Objective reframed as
  **Δsplicing**: `P(site|designed) − P(site|reference)`. Validated against MFASS
  / Vex-seq / MPSA and cross-checked vs SpliceAI/Pangolin/SpliceBERT. Replace
  BT3's saturating noisy-OR aggregation with top-k / log-odds pooling.
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
- **Design controls:** toggle each constraint with real controls — GC window
  slider, max homopolymer, CpG mode (elevate/deplete), forbidden motifs,
  **restriction-enzyme picker** (REBASE catalog); choose objective **weights**
  *or* **ε-budgets** per term.
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

- **Stack:** Python 3.11+ for orchestration/API; a compiled **Rust core (PyO3 /
  maturin)** for the trellis inner loop and incremental scorers, with a
  pure-numpy fallback when the extension isn't built. Rust over Cython for
  memory safety, GIL-free DP, and rayon parallelism (Pareto sweep + batch).
  ILP/CP-SAT via OR-Tools; folding via ViennaRNA bindings; ML via torch (lazy).
- **Performance:** everything **incremental** — each refinement move is
  O(context), not O(L). *(BT3's SA recomputed full CAI + full validation per
  move — quadratic.)* Length-scaled iteration counts, block/segment moves,
  adaptive / parallel-tempering schedules. **Runtime and peak-memory scaling are
  regression-tested in CI** with an asserted curve and a wall-clock ceiling —
  BT3 had no performance test despite runtime being BT2's original weakness.
- **Testing:** Hypothesis property tests for every §5 invariant; **golden tests**
  (fixed protein panel → pinned DNA + metrics, regenerable by one command);
  **optimality tests** (beam/relaxed must match exact DP on small instances; the
  certificate must not lie); the committed **benchmark harness**
  (`scripts/benchmark.py`) comparing BT4 vs input vs GeneOptimizer/IDT/Twist on
  CAI/tAI/GC/CpG/repeats/splice/folding over a pinned panel (revive
  `almost-there`'s benchmark corpus).
- **CI from commit #1** (BT3 had none): GitHub Actions matrix (3.11/3.12/3.13)
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
  `max_homopolymer_run`, `reverse_complement`) available and CI-built; porting the
  full trellis inner loop to Rust and adding `build-table` remain. *This alone
  already beats BT3* (honest optimality, correct incremental GC).
- **Phase 2 — Multi-objective, richer biology & first app.** 🔶 **In progress.**
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
  Lagrangian backend whenever a local constraint or pairwise term is present. Also
  landed: an **internal-ATG strong-Kozak constraint**
  (`constraints/kozak.py`, LOCAL, `ok_suffix⇔validate`-tested) that forbids
  internal `ATG` in a strong Kozak context (purine at -3, G at +4) — with the
  genuinely non-local **uORF pairing** part *honestly deferred* to Phase 3, not
  faked with a padded window; an **enriched benchmark harness** reporting CpG,
  tandem-repeat, and %MinMax alongside GC/homopolymer/CAI (still naive-vs-BT4, no
  fabricated competitor numbers); and BT4's first **performance/scaling
  regression test** (`tests/test_performance.py`, §7) asserting sub-quadratic
  exact-DP runtime under a wall-clock ceiling. Remaining: tAI (needs authentic
  tRNA data); uORF pairing and CpG/whole-sequence *count* budgets (both non-local
  / not per-codon-decomposable, deferred); and the published comparison vs
  GeneOptimizer/IDT.
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
  `folding_calibrated=False` when the baseline is used (per §6/§10.6). Still
  ahead: the SpliceAI/Pangolin-class model trained on real GENCODE (Δsplicing
  objective, held-out-chromosome gate, hash-pinned artifact), parallel-tempering
  / block moves, the opt-in **ASSP** cross-check (§6) with offline fixtures, and
  per-site risk tracks in the UI.
- **tAI — landed (real data).** The deferred tAI item is now shipped honestly:
  `biomodels/codon/tai.py` builds relative adaptiveness from **real human tRNA
  gene copy numbers** (GtRNAdb hg38, 431 genes/47 anticodons, bundled with a
  content-SHA-256 provenance sidecar and independently re-counted) via a faithful
  port of the dos Reis (2004) `get.ws` wobble model (verified s-values; the
  bacterial lysidine path gated on `sking=1`). `TaiTerm` is the exact-DP
  `tai_logw` objective (additive, `delta==score`), wired through `tai_weight` in
  the pipeline/CLI/app and as a frontier axis, with the tRNA table's hash entering
  the provenance stamp when tAI is active. Other organisms raise until their
  tRNA data is added -- no fabricated tables.
- **Phase 4 — Learned expression & polished app.** Expression predictor head,
  frontier reranking with calibration/uncertainty. Polish BT4 Studio (theming,
  accessibility) and ship **packaged installers** (PyInstaller/Briefcase) for
  macOS/Windows/Linux; optionally expose the `service/` HTTP API. External-
  validation report vs real gene distributions and published tools.
- **Phase 5 — Scale & ecosystem.** More organisms with authoritative provenance,
  library/degenerate-design mode (sample the codon distribution, not just a
  single MFC target), restriction-enzyme catalogs, tissue/condition-specific
  tables.

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
