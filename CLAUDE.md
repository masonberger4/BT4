# CLAUDE.md — BT4

Guidance for Claude Code (and humans) building **BT4**, the from-scratch
successor to BT3. This file is the constitution of the new repository: read it
before writing code, and keep it current as the architecture evolves.

> **This file is the durable constitution** — the thesis (§1), architecture (§3),
> the three contracts (§4), the honesty invariants (§5), scientific scope (§6),
> engineering rules (§7), and the BT3 anti-patterns (§10). It changes rarely and
> deliberately, and it is **not** where per-PR status is tracked.
>
> **Live status and the next-task queue live in
> [`docs/NEXT_SESSION.md`](docs/NEXT_SESSION.md)** — a machine-scannable component
> board + an ordered, precondition-tagged queue with a single "start here". Shipped
> history lives in [`CHANGELOG.md`](CHANGELOG.md). Read NEXT_SESSION.md to learn
> *where we are and what to build next*; read this file to learn *the rules*. §9
> below is the durable **phase intent**, not the running status. No status fact
> should be written in both places.
>
> **Current phase (one line):** Phases 0–2 complete · Phase 3 groundwork landed ·
> Phase 4 **app polish + the expression promotion seam landed**, learned head still
> calibration-blocked (no attestation is bundled, because none has been earned) · Phase 5
> open, with **declared CAI reference sets** (highly-expressed by default) landed.
> The honest exact-DP + Pareto core, the Rust trellis port, the
> wrapped-but-**uncalibrated** RiboNN/SpliceAI/Pangolin models, the opt-in
> out-of-loop **ASSP** cross-check, and the full expression/splice design flow are
> all on `main` — and all of them now have a first-class BT4 Studio surface
> (*Design* with the consented ASSP cross-check and the reference-set picker,
> *Candidates & splice audit* with the opt-in RiboNN head, *Library (sampled)*,
> plus menus, shortcuts and runtime theming). **A measured 2026-08 audit
> ([`docs/REVIEW_2026-08_expression_and_context.md`](docs/REVIEW_2026-08_expression_and_context.md))
> re-pointed what remains**, and its four measured defects — three of which broke a
> §5 invariant — are now **fixed** (honest `RELAXED` certificates for unenforced
> GLOBAL rules on the frontier, GLOBAL-aware `validate`, a windowed `folding_dg`, and
> a new opt-in constraint `relax()` with culprit-named infeasibility; see
> [`docs/REVIEW_2026-08_sota_and_roadmap.md`](docs/REVIEW_2026-08_sota_and_roadmap.md)
> §3, which also benchmarks BT4 against 2023–2026 SOTA). **The defensible-default
> pass and the construct-context core have since landed too** (§4 of the same doc):
> a **windowed-GC** constraint routed honestly by tractability, the IUPAC
> extra-sites path, regime-tagged **application presets** (none applied by default),
> real Studio controls (objective weights, hard budgets, FASTA open, a validate
> panel, a splice track) — and then the standing architectural gap itself:
> **`ConstructContext`** now carries the 5′UTR / vector backbone, `SeededConstraint`
> makes every LOCAL rule junction-correct without touching the `Constraint`
> protocol, `UorfConstraint.cds_offset` catches a leader ATG reading into the CDS,
> one shared `junction_window()` folds the initiation region, and
> `api.audit_construct` audits the assembled construct including **restriction-site
> uniqueness**. What remains: real flanks for the wrapped splice CNNs (still
> `N`-padded), data and human-gated calibration, packaged installers, and breadth —
> see NEXT_SESSION.md.
>
> This document was written after a full review of the BT3 codebase and *every*
> BT3 branch (`master`, `almost-there`, `gemini`, `streamlit`, and the merged
> `claude/ultracode-app-redesign` line); the lessons are folded in below.
>
> **Keep the split clean.** When a phase lands, a contract changes, or the
> architecture evolves, update the durable rule *here in the same change* (§10.11)
> — but put the *status* in NEXT_SESSION.md, not in a running paragraph here. A
> stale constitution, and a fact duplicated across two docs that then drift, are
> both BT3 anti-patterns (§10.11, §2).

---

## 1. What BT4 is

BT4 back-translates a **protein** into a **coding DNA / mRNA** sequence that is
**optimized for expression-relevant objectives** in a target organism (default *Homo
sapiens*) **subject to biological constraints** (GC content, homopolymers,
tandem/inverted repeats, forbidden & restriction motifs, internal ATG / Kozak /
uORF, CpG budget, cryptic splice sites, 5′ mRNA folding).

**Scope, stated honestly (§10.6 applied to BT4's own framing).** BT4 optimizes a
**coding sequence**. It no longer optimizes it *in isolation*: `ConstructContext`
carries the 5′UTR and vector backbone, `SeededConstraint` makes every LOCAL rule
junction-correct, `junction_window()` folds the initiation region, and
`api.audit_construct` audits the assembled construct — so the initiator Kozak
context is reachable and a junction defect is visible. **What has not changed is
the claim BT4 is entitled to make.** Supplying context is optional and often
omitted, and on the bare-CDS path folding still sees only `CDS[0:48]` and the
splice CNNs still see the CDS padded with literal `N` — a measured *lower bound* on
those models' response, not their estimate (§6). Real flanks for the wrapped splice
CNNs remain unbuilt. So "expression-relevant objectives" stays the accurate claim —
a validated *expression outcome* is not one BT4 can make, and must not be written
as though it were. For the measured evidence behind the context work, see
[`docs/REVIEW_2026-08_expression_and_context.md`](docs/REVIEW_2026-08_expression_and_context.md)
for the measured evidence and
[`docs/NEXT_SESSION.md`](docs/NEXT_SESSION.md) for the queue.

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
  io/            fasta, json_out (versioned self-describing schema), run_manifest,
                 genbank (annotated reader/writer: residual violations ride out as
                 `misc_feature` spans; a vector map can be read back in as a
                 ConstructContext). A standalone manifest file is still INTENDED,
                 not built — do not cite it as shipped.
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
| CAI (Sharp & Li), against a **declared reference set** (§8) | additive per-site | exact DP |
| tAI (dos Reis) | additive per-site | exact DP |
| Codon-pair bias (CPS) | pairwise (prev codon) | DP, state extended by prev codon |
| 5′ ramp / %MinMax | positional, first ~30–50 codons | DP, position-aware weights |
| GC / GC-window | cumulative, bounded window | ILP linear constraint **or** incremental running-GC counter (never re-scored per window) |
| CpG / UpA budget | whole-sequence count | ILP budget or Lagrangian dual |
| 5′ folding ΔG | non-local, range-bounded near start | refinement layer (ViennaRNA) |
| Cryptic-splice Δrisk | non-local, whole-sequence | refinement layer (learned model) |
| Learned expression | non-local, learned | frontier reranking |

### Constraints (each a new file + registry entry)

**functional poly(A) signal** (`constraints/polya.py` — an `AATAAA`/`ATTAAA`
hexamer forbidden *only* when a downstream U/GU-rich element follows it, the
bipartite signal CPSF/CstF actually recognise; strictly more permissive than the
blunt `poly_a_signal` hexamer preset, which remains available as the stricter
option, and refinement-enforced because its ~45 nt footprint is far too wide for
the trellis — a structural rule making no calibrated cleavage claim),
GC, homopolymer, **max GC-run** (longest run of consecutive G/C, the "max GC
length" knob — LOCAL, exact), tandem **and reverse-complement/inverted** repeats,
**max repeat length** (any dispersed direct/inverted/palindromic repeat, anywhere,
RC-aware — genuinely GLOBAL, so refinement-enforced and honestly reported, *never*
merged into the trellis, §10.1), forbidden-motif (with named, documented
**forbidden-sequence presets**), **restriction sites (a 584-enzyme catalog derived from
a version-pinned REBASE release — every commercially available *restriction*
enzyme with a single fully-specified site — Type II, IIG, Type IIS, and the
modification-dependent IIM ones (DpnI's `GATC` avoidance is mainstream precisely
because a dam+ plasmid *is* Dam-methylated and *is* cut); methyltransferases and
homing endonucleases excluded; auto reverse-complement + IUPAC; content-hashed, entering the run
manifest when active, and re-derivable, never hand-typed; the recognition
*sequence* only, not cut position/star activity/methylation)**, **strong splice-consensus donor/acceptor motifs**
(`SpliceSiteMotifConstraint`/`avoid_splice_sites` — LOCAL, exact, IUPAC,
**sense-strand only** so *no* RC-banning; an honest structural heuristic that bans
only the strong consensus — donor `GTRAGT`, acceptor `YYYYYYNYAGG` — never the bare
`GT`/`AG`, and makes no calibrated splice claim, with the real audit deferred to
the CNNs, §10.6), internal-ATG → proper **Kozak context** + **uORF
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
  with **top-k / log-odds pooling** (never saturating noisy-OR, §10.14). **That
  pooling is hinged at an uncalibrated background, and the hinge is load-bearing:**
  `pool_log_odds` counts only positions above `DEFAULT_SITE_PROBABILITY = 0.5`, and
  **that one constant fails in opposite directions for the two backends that use it.**
  Measured on the 93-sequence designed-CDS panel: the opt-in **Pangolin** clears 0.5 on
  only **6 of 93** sequences (all six designs of one protein; no native clears it in any
  group), so for two of three proteins every pooled risk and every `delta_splicing` is
  identically zero while the underlying scores vary more than twofold — the signal is
  discarded. The **PWM baseline that `default()` actually returns** clears it on **93 of
  93**, peaks 0.981–1.000, so its hinge never binds and it flags a site on *every*
  designed CDS including the natural one — the signal is saturated. BT4's headline
  splice objective is therefore **mute on one path and indiscriminate on the other**
  until Part B derives a real operating point on data; one threshold is standing in for
  two different score scales. The rule that
  follows: **a pooled risk of zero is never reported alone.** `pooled_risk_detail`
  carries `below_background` / `max_score` so a floored zero is distinguishable from a
  measured one, `pool_top_k_logit` is the background-free ranking statistic that
  survives the hinge (**not a risk** — it goes negative and has no calibrated zero),
  and lowering the background to make the signal reappear is forbidden: it is the same
  uncalibrated knob pointed somewhere more flattering. **And the `N`-padding is not
  neutral either:** measured, replacing the adapters' 5,000 literal `N` with real human
  genomic flank raises the median peak score inside the CDS from **0.276 to 0.369**
  on the same 9-sequence set, enough to move designed sequences across the 0.5 cutoff — while three *different* real
  regions agree to three decimals and a shuffled control (same composition) inflates
  scores in 9 of 9, so the effect is distribution shift, not "any bases beat `N`". A
  splice number computed on the `N`-padded path is a **lower bound** on that model's
  response, not its estimate; supplying the real `ConstructContext` changes the answer
  rather than refining it. None of this licenses moving the cutoff — there are no labels,
  and a higher score is not a more correct one. **And the models are not inert here:**
  a textbook donor consensus planted into a designed CDS lifts the local peak from 0.052
  to 0.570 (~11x) at exactly the anchor base, while a composition-matched scramble and a
  `GT`->`CT` ablation keeping 7 of 9 bases both sit at host baseline — so the response is
  the splice signal, not a reaction to an edit. **But the floor is high:** a *weakened*
  real donor scores 0.357 and clears nothing, so BT4's 0.5 cutoff sits above the
  intermediate-strength sites cryptic splicing actually uses. This bounds inertness from
  below and **does not** license "a clean designed CDS has no cryptic site" — detecting a
  site BT4 planted is not evidence about sites nobody put there. See
  [`docs/REVIEW_splice_calibration.md`](docs/REVIEW_splice_calibration.md). *(License
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
  independently-trained CNNs reachable. **Pangolin has now passed its
  integration-fidelity gate** (2026-08, on a maintainer machine holding the GPL
  weights): 18 cases, **max abs deviation exactly 0.0** — the adapter reproduces
  upstream's per-position scores bit-for-bit. Its `FidelityAttestation` is
  committed (`biomodels/splice/data/pangolin.attestation.json`), and
  `promote_if_attested` honors it **only under an explicit opt-in**
  (`BT4_SPLICE_USE_ATTESTED=1` / `--use-attested-splice`), so `default()` still
  returns the PWM baseline. **SpliceAI has now passed its gate too** (2026-08, 18
  cases on the same panel, max abs deviation exactly 0.0), so both wrapped CNNs
  reproduce their published models bit-for-bit and both attestations ship. Promotion
  stays behind the same opt-in. The captured panel itself is never committed: it *is* the licensed
  model output, and only the attestation's eight license-clean scalars plus the
  public weight SHA-256s ship. **Honest scope:** these predict splice-*site
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
- **ASSP — optional online cross-check backend (landed).** (Alternative Splice
  Site Predictor.) A *kept, supported feature*, behind the same `SplicePredictor`
  contract as every other backend, with guardrails that make it honest rather than
  the BT3 liability it was. **Status: landed** — `AsspSplicePredictor`
  (`biomodels/splice/assp.py`), the graceful `run_splice_crosscheck`
  (`pipeline/splice_crosscheck.py`, exposed as `bt4.api.splice_crosscheck`), the
  `bt4 validate --splice-backend assp` / `bt4 optimize --check-splice assp` CLI
  wiring, and committed offline fixtures (`tests/fixtures/assp/`) all ship:
  - **Opt-in and out-of-the-inner-loop.** ASSP is a network service, so it is
    never used to score per-move inside the optimizer, and is **never** returned by
    `bt4.biomodels.splice.default` or `available_splice_backends` (it is requested
    explicitly by name, never auto-discovered). It runs only as a **final audit /
    validation pass** on a delivered sequence (`bt4 validate --splice-backend
    assp`, or `--check-splice assp` on `optimize`).
  - **Never silent, never blocking.** Gated behind an explicit flag *and* the
    `bt4[assp]` extra (httpx, lazily imported); polite rate-limiting + exponential
    backoff (`_throttle` / `_with_retries`); responses cached by sequence hash
    (`CachingAsspTransport`); if the service is unreachable or returns a garbled
    body the *raw* predictor raises an `AsspError`, but `run_splice_crosscheck`
    catches it and reports "unavailable" — it can never fail an optimization. CI
    uses stored **offline fixtures** (`FixtureAsspTransport`, selected via
    `$BT4_ASSP_FIXTURE_DIR`), no live calls. **Honest wire-format caveat:** the
    live transport targets ASSP's documented tabular site report but is *unverified
    against the live service* (unreachable during development), so the committed
    fixtures are *synthetic ASSP-format reports* (not real captures) — the same
    "no bundled panel ships" posture as the wrapped CNNs; the promotion path is a
    maintainer confirming the live transport.
  - **Labeled non-reproducible.** `AsspSplicePredictor.network_derived` is `True`
    and `calibrated` is `False`; any ASSP-derived number is stamped network-derived
    and **excluded from the reproducible-from-manifest guarantee** — the
    cross-check is reported as a separate advisory section (the CLI prints it to
    **stderr**, never into the stdout FASTA/JSON artifact or a `Result` manifest).
    The local baseline (and, when installed, the wrapped CNNs) remain the default,
    reproducible splice path.

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
  clearly labeled network-derived) on the delivered sequence. **Landed.** It is
  the *only* control that leaves the machine, so it is explicit in every
  direction: consent is asked before anything is sent (naming the service), the
  call runs on a background worker, and the report is led by its tags —
  network-derived, UNCALIBRATED, advisory, **not** part of the run manifest and
  never exported. An outage degrades to a labeled "unavailable" banner and can
  never fail a run (§10.15). The panel is cleared whenever the delivered sequence
  changes, so one sequence's sites can never be shown beside another's.
- **Expression head:** the wrapped **RiboNN** head is opt-in from the Candidates
  tab (toggle + species + the fixed 5′/3′ UTR context it requires), enabled only
  when `available_expression_backends()` reports the user's own checkout and
  weights actually resolve — never a dead control, and it says what is missing
  otherwise. It is `calibrated=False`, so the set stays **discovery order, not a
  ranking** and the solver's pick stays delivered (§10.6).
- **Library mode:** `api.library` has its own tab — members / temperature / seed,
  a per-member table, and a multi-record FASTA export — banner-led with
  **sampled, not optimized**, its badge coloured directly from the `SAMPLED`
  certificate so the label cannot drift from the engine's claim.
- **Export:** FASTA / JSON (the run manifest rides inside the JSON's `audit`) /
  **annotated GenBank**, so anything the app shows is reproducible from its stamp
  (except explicitly network-derived ASSP numbers, which never reach an export —
  regression-tested). The GenBank writer is the honest-residual surface: every
  residual violation is emitted as a `misc_feature` span at the base where it
  occurs, so a defect the optimizer could **not** remove reaches the map the user
  actually opens (SnapGene/Benchling/ApE) instead of living only in a JSON audit;
  overlapping findings merge into one readable span that names how many it covers,
  and the true count stays in the COMMENT block. The record carries **no
  timestamp** (invariant #7) and stamps the config hash + git commit instead. A
  standalone manifest file remains a design target, not shipped.

Accessibility, light/dark theming, and responsive layout are in scope — "nice"
is a requirement, not a nice-to-have. **Landed:** a File/Run/View/Help menu bar
with standard shortcuts makes every action keyboard-reachable, **View → System /
Light / Dark** switches theme at runtime (restyling stylesheet, plots, badges,
and the sequence viewers' violation bands from the still-live results), and every
control carries an accessible name, a buddy label, and an explanatory tooltip.
Only one engine flow runs at a time, gated from a single set of running-flags
rather than from thread references — so a missed reference clear cannot strand a
control.

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
- **Agent CI check-in cadence (convention).** When an automated agent opens a PR
  to this repo and watches its CI, the default is **auto-merge**: immediately after
  opening the PR (ready for review), enable GitHub **auto-merge with the squash
  method** (`enable_pr_auto_merge`) so GitHub itself merges the PR the moment its
  required status checks pass — fully event-driven, no polling and no timed merge.
  The rest stays **webhook-first**: the PR subscription wakes the session on CI
  *failures* and review comments, and a failing check simply blocks the pending
  auto-merge, so the drive-to-green posture is unchanged (diagnose + push a fix, or
  reply with the blocker; each new push keeps auto-merge armed). The merge happens
  only on all-green either way. *Why this is the default and not shell polling:* in
  the agent's execution environment GitHub is reachable **only** through the MCP
  tools (which require the agent to be invoked to call them) — there is no `gh` CLI
  or shell-level API access a background watcher could poll with, and the PR
  subscription does **not** reliably push CI *success* — so letting GitHub own the
  merge trigger is strictly better than any timer the agent could set.
  - **Repo prerequisites (one-time, owner-set).** Auto-merge needs two GitHub
    settings enabled: *Settings → General → Pull Requests → Allow auto-merge*, and
    branch protection on `main` with the CI checks marked **required**. The second
    is what makes auto-merge *wait for green* rather than merge a mergeable PR
    instantly — without a required check, a PR with no other block is immediately
    mergeable, so `enable_pr_auto_merge` returns a clean-status error and the agent
    uses the fallback below. With both on, `enable_pr_auto_merge` arms the PR and
    GitHub merges it the moment the required checks pass.
  - **Stale branch (`mergeable_state: "behind"`).** If branch protection also has
    *Require branches to be up to date before merging* on, an armed PR whose CI is
    already green can still sit unmerged at `mergeable_state: "behind"` once `main`
    advances under it — GitHub holds the merge rather than merging stale. Clear it
    with a **single `update_pull_request_branch`** call (the "Update branch" action:
    it merges the current base into the PR head, re-triggering CI). This is **not**
    a manual merge and does not defeat auto-merge — it satisfies the up-to-date gate
    so auto-merge can fire once CI is green again. Do not call `merge_pull_request`;
    let GitHub own the merge. (Same posture as a merge-conflict notice — drive the
    PR to mergeable, then let auto-merge complete it.)
  - **Fallback (only when auto-merge is unavailable).** If `enable_pr_auto_merge`
    errors — the repo has no required-status-check branch protection, or auto-merge
    is disabled at the repo level — fall back to the prior pattern: arm **at most a
    single self-check-in** (`send_later`) sized to when CI actually finishes (a few
    minutes for this repo, not a 60 s round) and **perform the merge from that
    firing**. This deliberately avoids the old linear-backoff loop, which armed a
    new trigger every round and left pending triggers to hand-delete. Rules that
    keep trigger churn at zero — **never call `delete_trigger` at all**; one-shots
    clean themselves up: (1) a `send_later` one-shot **self-disables after firing**
    (`ended_reason: run_once_fired`), so a fired trigger is already dead — nothing to
    delete; (2) do not stack a second check-in while one is pending — re-arm only
    *after* the current one fires and the PR still isn't resolved; (3) if the PR
    resolves (merges/closes) before a pending fallback fires, **do not delete it
    either** — let it fire once and no-op (see the PR is merged, do nothing). A
    single harmless no-op wake is cheaper than a permission-gated `delete_trigger`
    call, and auto-merge already removes the timer from the common path — so an
    agent should never need to delete a trigger.
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
- **A codon table must declare its reference set, and the default is
  highly-expressed.** `w = f/f_max` is meaningless without saying *which genes it
  was counted over*, so every table carries a `reference_set` label that travels
  with it into the result's audit, the CLI/app output, and the manifest — a
  number and the question it answers are never separated. BT4 ships two, both
  counted from public pinned sources by the same filtering rules so they differ
  only in *which genes*: `highly_expressed` (the top-N most abundant proteins by
  measured proteomics — CAI in Sharp & Li's original sense, `w = 1` marking the
  codon translation *prefers*) and `genome_wide` (every gene, marking the codon
  that is merely most *common*). Highly-expressed is the default wherever it
  exists, because codon optimization targets a highly-expressed protein and a
  genome-wide count is dominated by genes under no translational selection.
  Requesting a reference set an organism lacks **raises**; it never silently
  substitutes the other, and an organism whose abundance data cannot be joined to
  the pinned annotation gets **no** highly-expressed table rather than a guessed
  one. Neither table is a validated expression predictor — a highly-expressed
  reference makes CAI a better-founded proxy, not a true one (§10.7).
- Reproducible **table-build pipeline**: `scripts/build_organism_tables.py`
  (genome-wide) and `scripts/build_highly_expressed_tables.py` (highly-expressed)
  recount every bundled table from release-pinned public sources, each with a
  `--verify` mode that rebuilds and diffs against the committed bytes *and*
  sidecars; `bt4 build-table --cds genes.fasta --organism X` recomputes a table
  from a user's own declared CDS set. A test asserts a bundled table reproduces a
  *published* CAI for a benchmark gene set. *(BT3 asserted "representative" with
  no such check.)*
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

> This section is the durable **phase intent** and a narrative record of what each
> phase established. For **live status and the next task**, read
> [`docs/NEXT_SESSION.md`](docs/NEXT_SESSION.md) (the single source of truth for
> volatile state) — do not treat the ✅/🔶 markers below as the authoritative
> current board.

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
  upstream model's per-nucleotide scores **bit-for-bit**. **Pangolin's gate has
  been run and passed** (18 cases, deviation 0.0); its attestation is committed and
  honored under an explicit opt-in. **SpliceAI's gate has since run and passed as
  well** — same 18-case panel, deviation 0.0 — so both attestations ship and both are
  honored under that one opt-in. `default()` keeps returning the PWM baseline either way
  — it needs no per-user weight configuration, so it cannot presume a licensed
  model is installed or wanted. A pass is **regime-scoped**: it was captured on the
  bare-CDS, `N`-padded path, so `score_in_context` and `_FlankedPredictor` clear
  `calibrated` when real flanks are supplied. And it is *integration* fidelity, not
  statistical calibration — across eight predictors these models score median prAUC
  **0.419 on exonic** variants vs 0.773 intronic (Smith & Kitzman 2023), and BT4
  designs coding sequence, so its whole regime is the weaker half. The **two-backend
  agreement harness** (`backend_agreement` + `scripts/compare_splice_backends.py`)
  reports pairwise Spearman rank / sign agreement across whichever backends are
  available — with both CNNs installed, agreement between two real,
  independently-trained models — the first-class uncertainty signal of §6. On top
  of it the **localize-and-flag splice audit** (design-flow step 4,
  `biomodels/splice/audit.py::audit_splice`, `pipeline/splice_audit.py`,
  `bt4.api.splice_audit`) runs the available backends over a step-3 candidate set to
  **localize** residual cryptic sites (one flag per contiguous above-threshold run,
  at its peak) and attach that pooled agreement as the authoritative cross-backend
  signal. It is **out-of-loop and advisory — it never edits** (a targeted
  synonymous auto-edit is a deferred, per-backend-`calibrated`-gated future step);
  `all_calibrated` is `False` today, every `SpliceFlag` carries its emitting
  backend's `calibrated`, the site threshold is a heuristic display knob (the PWM
  baseline's score an arbitrary-units pseudo-score), per-flag `added_risk_vs_reference`
  is positive-worse and intra-backend (distinct from the larger-is-better panel
  `delta_splicing`), and cross-backend `also_flagged_by` is a raw ±window positional
  co-occurrence, never a kind-level agreement (Pangolin's combined track can't
  disagree on kind, so its flags are labelled `"splice"`). The
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
  Pangolin (GPL) and SpliceAI (CC BY-NC) are eligible to certify. **Pangolin has now
  passed**, and a committed attestation ships; promotion is opt-in
  (`BT4_SPLICE_USE_ATTESTED`), so `default()` still returns the PWM baseline.
  **A fidelity attestation is not a calibration claim, and the separate gate for that
  now exists too** (`biomodels/splice/gate.py` + `pipeline/splice_gate.py`,
  `api.splice_panel_gate` / `bt4 splice-gate`): a **per-stratum** verdict over
  PR-AUC / prevalence-normalized skill / ROC-AUC / top-k / MCC / Brier + Brier-skill /
  ECE — deliberately **not** Spearman, which on a binary label adds nothing over
  ROC-AUC and whose expression-gate threshold a *perfect* classifier fails at splice
  prevalence. It demands a declared `negative_construction` (average precision's floor
  is the prevalence, so a threshold without a pinned denominator is passable by
  thinning negatives), reads a strict panel format whose position convention is
  **verified against the sequence** rather than trusted (`api.read_splice_panel`
  refuses a mis-anchored panel and names the shift that would have worked), and runs
  four permanent baselines a backend must beat *in every stratum* — `permutation`,
  `gt_ag` (the canonical dinucleotide ~99% of introns follow), `pwm` (BT4's own shipped
  baseline, the free incumbent) and `constant`. Running the PWM backend as the head
  ties the `pwm` baseline exactly, so **BT4's default can never be evidence for
  itself**. Still ahead: the same fidelity gate for **SpliceAI**, and the regime-matched
  panels Part B needs (human-only, data-gated). The opt-in **ASSP** cross-check (§6) with
  offline fixtures has now **landed** — `AsspSplicePredictor` +
  `run_splice_crosscheck` + `bt4 validate --splice-backend assp` / `bt4 optimize
  --check-splice assp`, opt-in / out-of-loop / cached / rate-limited-with-backoff /
  network-derived-and-manifest-excluded / never-blocking, CI-driven from committed
  synthetic offline fixtures (block/segment moves and parallel tempering also
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
  **candidate-set assembler** (design-flow step 3,
  `pipeline/candidates.py::assemble_and_rank_candidates`, `bt4.api.candidates`) now
  builds the finalist set that rerank ranks: the **Pareto frontier** plus, when a
  GLOBAL rule is active *and* the delivered exact-DP seed violates it, a small
  **deterministic library of repeat-refined variants**; de-duplicated and
  batch-scored via the new `BatchExpressionPredictor` contract (`score_many`, e.g.
  RiboNN) when available. It obeys the same calibrated-gating rule - **discovery
  order + solver-delivered `chosen` when uncalibrated**, reorder + top pick when
  calibrated - and hardens it: the delivered (`chosen`) sequence is **invariant to
  `n`** (uncalibrated, the solver-delivered is pinned first; calibrated, the head's
  top pick tops the top-n keep; cap applied after scoring), variants are labelled
  `repeat_refined` (process, not a guaranteed fix) with residual GLOBAL violations
  disclosed per member, and de-dup/cap counts plus the predictor identity (in the
  manifest, invariant #9) are reported. **The design flow now has a first-class
  desktop surface**: BT4 Studio's new **Candidates & splice audit** tab
  (`app/worker.py::CandidatesWorker`, `app/studio.py`) runs `api.candidates` →
  `api.splice_audit` on a background `QThread` and renders the ranked, honestly-
  labeled candidate set (delivered pick starred, per-member source / CAI / GC /
  expression+units / calibration / hard-violation / distinct-splice-site counts)
  with two advisory banners — the same calibrated-gating honesty as the API: an
  **uncalibrated head is shown as discovery order, NOT a ranking** with the
  solver's pick starred and scores annotating only (a calibrated head switches to
  ranked-by-expression + top pick), and the splice banner leads with **UNCALIBRATED
  (advisory)** whenever `all_calibrated` is `False`, reports cross-backend rank/sign
  agreement, and states the flags localize sites heuristically and edit nothing
  (§10.5/§10.6/§6). The splice-flags column counts **distinct** sites (co-located
  cross-backend flags merged within the audit's match window), every table metric
  is recomputed per candidate from its own DNA (invariant #2), and an opt-in
  toggle routes the installed SpliceAI/Pangolin CNNs into the audit (PWM baseline
  only when off). The
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
  is not calibration *for BT4's regime* — **and the reason is not the one this
  document used to give.** RiboNN's ablation reports ~31% of per-*nucleotide*
  signal in the CDS, but its **length-integrated** attribution is 22/73/5
  (5′UTR/CDS/3′UTR), so the CDS is the *majority* of the total attributed signal;
  quoting only the per-nt figure is a misread (see
  [`docs/RESEARCH_codon_optimization_SOTA.md`](docs/RESEARCH_codon_optimization_SOTA.md)
  §0). The load-bearing gap is sharper: RiboNN has never been shown to
  **discriminate synonymous CDS variants of the same protein under a fixed UTR**,
  which is exactly BT4's regime. So promotion needs a passing
  `verify_expression_gate` on a CDS-variant panel — human-only, data-gated. New `bt4[expression-ribonn]`
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
  synthetic RiboNN output table; `calibrated` stays `False` (unchanged). **Batched
  scoring has now landed** (the expression/splice design flow's step 1,
  `docs/DESIGN_expression_splice_flow.md`): public `score_many(dnas)` /
  `delta_logte_many(designed, reference)` methods route a whole candidate set
  through the existing batched `_predict_te` (one TSV, one `predict` invocation over
  RiboNN's `top_k`-model ensemble), amortizing
  RiboNN's large fixed *per-invocation* overhead so scoring a frontier costs roughly
  one call rather than N; `delta_logte_many` scores the shared reference **once**;
  both keep the per-input validation and the `tx_id` realignment and return results
  in input order, and `score_sequence`/`delta_logte` now delegate to them.
  **Correction (2026-08, verified against upstream):** an earlier note here claimed
  RiboNN's predict entry point "exposes no worker-count parameter" and that a
  `num_workers=0` path was therefore left out. That was **wrong** --
  `predict_using_nested_cross_validation_models` takes both `batch_size` (default
  1024) and `num_workers` (default 4). The adapter now **forwards both**, defaulting
  to `batch_size=64` / `num_workers=0`, and neither can change a score: RiboNN pads
  every transcript to a *fixed* width (`max_utr5_len + max_cds_utr3_len` = 13318),
  not to a batch's longest member, and builds its predict dataloader with
  `shuffle=False`. The defaults are a **correctness requirement, not tuning**: the
  adapter scores from a mutated `sys.path` and a temporary working directory, which a
  *spawned* worker (Windows, macOS) does not inherit -- so `num_workers>0` hangs or
  fails there, and RiboNN rebuilds the dataloader once per ensemble member (up to 50
  times), paying the spawn cost every time; `batch_size=1024` allocates 1024
  fixed-width `(channels, 13318)` float32 tensors at once and OOMs an ordinary CPU
  box. `calibrated` remains `False` (no calibration claim -- a knob is not a gate).
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
- **Phase 4 — Learned expression & polished app.** 🔶 **App polish landed; the
  learned head stays calibration-blocked.** The BT4 Studio polish pass has
  shipped (§6.6): the two engine-ready backends that had no UI are surfaced — the
  wrapped **RiboNN** head opt-in on the Candidates tab (reached through the new
  public expression-backend registry, `available_backends` / `resolve_backend`,
  re-exported as `api.available_expression_backends` /
  `api.resolve_expression_backend`, so the app selects a head by name without
  importing `biomodels` across a layer, §3/§10.9) and the opt-in **ASSP**
  cross-check on the Design tab; **library mode** has its own tab; and the
  keyboard/menu/theming/accessibility work is done. Every one of these is
  calibrated-gating-preserving plumbing over the stable `bt4.api` — no engine
  change and no new claim: RiboNN and every splice backend remain
  `calibrated=False`, so they annotate and advise but never steer delivery.
  Remaining: the expression predictor head **earning** calibration (frontier
  reranking with calibration/uncertainty — blocked on a regime-matched
  CDS-variant panel), **packaged installers** (PyInstaller/Briefcase) for
  macOS/Windows/Linux, and the external-validation report vs real gene
  distributions and published tools.
  **The apparatus for that calibration has now landed, and only the data step
  remains.** The gate can finally judge the regime BT4 deploys in rather than the
  one RiboNN was trained on: `within_group` scores *inside* each protein (pooled
  scoring credits between-protein skill, which is exactly what a natural-gene-trained
  head has and exactly what BT4 cannot use — a regression test pins that a
  gene-identity-only head passes pooled and fails within-group); `recalibrate` fits
  the affine link on the calibration fold **only**, because a head reporting a CLR
  compositional residual cannot be compared to an assay's units by subtraction; the
  rank metric stays on **raw** predictions, so a head that ranks backwards cannot be
  rescued by a negative fitted slope while a deployed BT4 hands the user the worst
  candidate; `width_over_iqr` exposes a vacuous interval, since split conformal is
  valid for *any* score function and a **constant predictor passes coverage**; and a
  cluster bootstrap resamples whole proteins, because one protein's variants are a
  dependent cluster. Around it: a strict panel format that **refuses** a row RiboNN
  would silently drop (`api.read_panel`), five permanent baselines a head must beat —
  CAI above all, since BT4 already optimizes it in-loop (`api.expression_gate`,
  `bt4 expression-gate`), and `ExpressionAttestation`, the single scope-carrying seam
  that can flip `calibrated=True`, replacing a bare `dataclasses.replace`. The
  zero-data checks that decide whether a panel is worth acquiring at all live in
  `scripts/ribonn_sensitivity.py`. Evidence, corrections and protocol:
  [`docs/RESEARCH_ribonn_calibration.md`](docs/RESEARCH_ribonn_calibration.md).
  RiboNN remains **`calibrated=False`** — none of this is a claim, it is only the
  apparatus that could earn one.
  **The seam that would let such a claim reach a user has now landed too, and its
  scope is bound rather than declared.** `verified_predictor` was the sanctioned way
  to flip `calibrated`, but **nothing in `src/` called it** — so a passing gate and a
  committed record would have changed nothing for anyone, unlike the splice side,
  whose `promote_if_attested` has been wired into production behind
  `BT4_SPLICE_USE_ATTESTED` since #121. `biomodels/expression/attestations.py`
  mirrors it: `BT4_EXPRESSION_USE_ATTESTED` (or an explicit `use_attested=` on
  `api.resolve_expression_backend`) honours a resolvable attestation, `default()`
  still returns the neutral placeholder, and **nothing auto-promotes**. Because an
  expression attestation is earned against a *maintainer's own measured panel* rather
  than published weights, `$BT4_EXPRESSION_ATTESTATION` reads a local record, so using
  a result never requires committing the panel's identity; a mis-pointed path
  **refuses** rather than falling back, since a typo would otherwise be
  indistinguishable from "no attestation". **No attestation is bundled** — none has
  been earned, and shipping one that had not is exactly §10.6's fabricated
  placeholder. BT4 Studio carries the opt-in per run (never by mutating process env),
  shows the attestation's **scope on the page** (species / cell types / readout /
  `top_k` / UTR contexts / panel hash), **pins the head to that scope** while it is
  honoured, and flips the Candidates banner from *discovery order, NOT a ranking* to a
  ranking **with its scope named**. **And the scope is now the run's, not the
  caller's:** `attest_expression` derives species / cell types / readout from
  `GateComparison.scope` and treats a declared value as a cross-check that **refuses**
  on mismatch — a gate run averaging all 78 cell types can no longer be filed as a
  HEK293T result — while `verified_predictor` additionally binds `top_k` and the UTR
  context, and deliberately does **not** bind `batch_size` / `num_workers`, which
  provably cannot change a score. A promoted head carries the attestation's content
  hash into the manifest, so two runs steered by different claims cannot share a stamp
  (invariant #9).
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
  residual violation reported honestly per member. **Organism breadth has landed
  for the six stranded species:** mouse, rat, zebrafish, *Drosophila*, *C.
  elegans* and *Arabidopsis* now ship **genome-wide recounted** codon tables
  (`scripts/build_organism_tables.py`), taking BT4 from three selectable organisms
  to nine. This closed a real gap rather than padding a list — all six already had
  authentic GtRNAdb tRNA tables, but tAI is only offered for an organism you can
  *select*, and selection needs a codon table, so six of the eight tRNA tables were
  **unreachable** (a regression test now forbids that state). The builder is the
  §8 reproducible table-build pipeline made concrete: a **release-pinned** Ensembl
  CDS FASTA (never a moving `current` link), documented validity filtering, **one
  representative CDS per gene** so usage is not weighted by isoform-annotation
  depth, counting via BT4's own `count_codons`, and a provenance sidecar carrying
  the source URL, the **source file's own SHA-256**, assembly, release, total
  codons, and the per-filter drop tally — so a third party rebuilds the exact
  shipped bytes (`--verify` proves it). It **refuses rather than fabricates**:
  no valid CDS, or any of the 64 codons unobserved, aborts the write instead of
  smoothing in an invented number. The tables are checked against **external
  ground truth** (§8), not just self-consistency — GC3 ordering across species,
  mouse/rat independently landing within 1.5 points of human, and the known
  preferred Leu/stop codons per species. **The three legacy hand-curated tables
  are done:** human, *E. coli* and *S. cerevisiae* are now recounted through the
  same pinned-Ensembl pipeline as the other six, so **all** bundled tables
  are re-derivable counts and no organism — least of all the default — rests on
  undocumented numbers. That rebuild changed no delivered sequence (byte-identical
  across a four-protein × three-organism panel; CAI moved ≤ 0.0003), because CAI
  normalizes within each synonymous group and the most-preferred codon per amino
  acid was unchanged.
  **Declared CAI reference sets have now landed, and highly-expressed is the
  default** (§8) — the change that moves BT4 off codon *commonness* as its
  headline readout. `scripts/build_highly_expressed_tables.py` counts each
  organism's 300 most abundant proteins, ranked by **PaxDb v6.1** whole-organism
  integrated proteomics (CC BY 4.0) and joined to the *same* release-pinned
  Ensembl CDS the genome-wide tables use, through that release's own peptide
  FASTA — no third-party mapping layer, ambiguous identifiers dropped rather than
  guessed, organelle-encoded genes excluded (different genetic code and tRNA
  pool), and all three sources hash-pinned with a `--verify` rebuild. **N = 300 is
  evidence, not taste:** it is the smallest size on a tested grid at which every
  bundled organism observes all 64 codons, so no shipped table needed smoothing —
  and far above it the reference set dilutes back into the genome-wide answer.
  The tables reproduce the classic *E. coli* and yeast optimal codons and show
  **stronger codon bias than genome-wide in all eight** organisms, with the gap
  largest in yeast/fly and smallest in human/rat — the ordering dos Reis (2004)
  predicts. Eight organisms have one; *A. thaliana* does not, because
  PaxDb identifies its proteins by UniProt accession that the pinned Ensembl
  Plants annotation does not carry, so BT4 ships none for it rather than one built
  on a guess. What remains in this phase is enumerated in
  [`docs/NEXT_SESSION.md`](docs/NEXT_SESSION.md), not here — a "still ahead" list
  in two documents is exactly the pair of status facts §10.11 forbids.
  **The three industrial expression hosts have now landed too** — **CHO**
  (*Cricetulus griseus*, CHOK1GS_HDv1), ***B. subtilis*** 168 and ***K. phaffii***
  GS115 (*Pichia pastoris*) — taking BT4 to **twelve** selectable organisms, each
  with a recounted genome-wide table *and* a GtRNAdb tRNA table, because a codon
  table without tRNA data would make tAI silently unavailable exactly where a user
  asked for it (a shipped invariant, not a preference). Two honest limits ride with
  them, both recorded in their sidecars rather than smoothed over: CHO's tRNA set is
  GtRNAdb's CriGri_1.0 while its codon table is Ensembl's CHOK1GS_HDv1 — **the one
  organism whose two inputs are not assembly-matched**; and none of the three ships a
  highly-expressed reference set, for three *different* measured reasons (PaxDb has
  no data at all for *K. phaffii*; only a single study, not the whole-organism
  integrated set, for CHO; and for *B. subtilis* the integrated set exists and joins
  at 99.8% via a declared `BSU` → `BSU_` locus-tag rewrite the builder does not yet
  support). *B. subtilis* also loses **22.5%** of its CDS to the shared ATG-start
  filter — more than double *E. coli*'s 9.6%, since it genuinely uses TTG/GTG starts
  — and that was **measured for this organism** rather than inherited: counting the
  dropped genes back in moves **no** amino acid's top codon and shifts `w` by at most
  0.023, so it is a precision gap, not a wrong answer.

  **Tissue/cell-type-specific tables are deliberately out of scope** (maintainer
  decision): the work is large, the resulting number is hard to qualify honestly,
  and the upside over a whole-organism highly-expressed reference is small.

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
