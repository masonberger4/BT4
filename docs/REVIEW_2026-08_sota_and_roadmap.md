# BT4 review — SOTA benchmark, honesty fixes, and the forward roadmap

**Date:** 2026-08 · **Companion to:**
[`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md)
(the measured audit) · **Scope:** where BT4 sits against the state of the art, what
this change fixed, and what to build next.

> **What this document adds.** The measured audit
> ([`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md))
> already scores BT4 against its own requirements and pins the four honesty defects
> with reproduction commands (its §4, §10). This companion does three things that
> audit does not: (1) benchmarks BT4 against the **published state of the art** in
> codon optimization (2023–2026), which is the research the review was asked to do;
> (2) records that the **four Tier-0 honesty defects are now fixed** in this change,
> with before/after reproduction; and (3) states the **sequenced roadmap** the live
> queue in [`NEXT_SESSION.md`](NEXT_SESSION.md) is re-pointed from. It does not
> restate the audit's scorecard or its measured defect evidence — it links to them
> (§10.11: no fact lives in two docs).

---

## 1. What state-of-the-art codon optimization looks like (2023–2026)

The field has two schools. **BT4 sits in the stronger one — and, on honesty and
hard guarantees, ahead of everyone in it.**

### 1.1 Exact / constrained combinatorial optimization — BT4's school

The landmark is **LinearDesign** (Baidu, *Nature* 2023): it frames mRNA design as
*exact* optimization over the exponential synonymous space by coupling the MFE
folding DP to a **codon lattice/DFA**, optimizing a single objective balancing MFE
and CAI via a weight λ, with beam search as an explicit speed knob. Follow-ons
(**EnsembleDesign**, expected-partition-function design, 2024–25) minimize *ensemble*
free energy. This is exactly BT4's codon-trellis-DP thesis — and BT4 generalizes it
to a **vector** objective + a returned **Pareto frontier** + honest **optimality
certificates**, which is genuinely ahead on rigor.

One architectural difference to name plainly: LinearDesign folds RNA secondary
structure *into the DP itself*; BT4 treats folding as a **refinement post-pass**.
That is defensible (honest, incremental) but means BT4 does not *co-optimize*
structure the way the structural school does. Joint codon+structure design is
tracked separately (NEXT_SESSION Tier 5), deliberately demoted below construct
context because joint folding that does not know the real 5′UTR is optimizing the
wrong window.

### 1.2 Generative deep-learning / language-model optimizers — the research frontier

A wave of transformer/LM optimizers now generate host-specific CDSs:
**CodonTransformer** (*Nat Commun* 2025 — multispecies BERT, 164 species,
context-aware, explicitly *minimizes negative cis-regulatory elements*),
**CodonBERT** (Genome Research 2024 — cross-attention, human/vaccine),
**CodonTranslator** (bioRxiv 2025 — decoder-only, 2,100+ species), **DeepCodon**
(2025 — preserves rare-codon clusters), **CodonRL** (RL, multi-objective),
**ColiFormer** (2025 — balances CAI/GC/tAI/RNA-stability/cis-element minimization),
and **RiboDecode** / deep-generative optimizers trained directly on **ribosome
profiling** (*Nat Commun* 2025).

- **What they do that BT4 doesn't:** learn expression signal end-to-end from large
  data, and produce "natural-like" sequences.
- **What BT4 does that they can't:** give **hard guarantees**. An LM
  *probabilistically* reduces cis-regulatory elements; it does **not** guarantee
  zero BsaI sites, a bounded homopolymer/GC-run/repeat, or an enforced GC-count
  window — the exact things a synthesis / cloning / AAV-LVV workflow requires.
  **This is BT4's real moat**, and it is the same guarantee IDT and Twist sell that
  the academic LMs do not.

### 1.3 The multi-objective consensus (BT4 is aligned)

The field agrees **CAI alone is a weak predictor**. SOTA optimizers jointly weigh
CAI, tAI, GC/GC-window, mRNA-structure stability (MFE), codon-pair bias, a 5′ ramp,
and cis-regulatory-element avoidance (splice sites, poly-A, AU-rich elements, uORFs).
BT4's objective-vector + constraint-registry design already embodies this. One
caveat the audit and the literature agree on: BT4's **5′ ramp term encodes a
mechanism the field falsified** (Goodman/Church/Kosuri 2013 — reduced 5′ RNA
*structure*, not codon rarity, drives the 5′ effect). The term is honestly hedged in
code as a heuristic; relabeling its claim and routing the real 5′ lever through
folding is a Phase-1 item.

### 1.4 Construct context is now first-class — and it is BT4's biggest real gap

- **The 5′UTR is a dominant, heavily-modeled expression lever with its own SOTA
  subfield:** **UTR-Insight** (2025, explains ~89% of mean-ribosome-load variance,
  up to +319% protein vs α-globin 5′UTR), **UTRGAN**, **UTailoR**, **UTR-LM**,
  **Optimus 5-Prime**. Reviews now push toward **coordinated / construct-level
  design** ("from UTR and codon optimization to coordinated design"). A CDS-only
  optimizer leaves the single best-characterized expression lever on the table.
- **Translation-efficiency / ribosome-load models** are the learned-expression SOTA:
  **RiboNN** (which BT4 already wraps, uncalibrated), cell-type-specific codon
  optimization (*NAR* 2025). RiboNN's own length-integrated signal attribution is
  ~5′UTR/CDS/3′UTR = 22/73/5 — i.e. it *needs the UTRs*, reinforcing that scoring a
  bare CDS is under-powered.
- **Splice:** **SpliceAI** (Illumina) and **Pangolin** (usage-rate quantification,
  RNA-seq across 4 tissues/species) are the accepted SOTA — **exactly what BT4
  already wraps.** Correct posture; the gaps are calibration and *real flanking
  context* (both CNNs are currently fed the CDS padded with literal `N`). Cryptic
  splice sites in transgenes are a documented clinical failure mode in AAV/LVV (e.g.
  the HMGA2 truncation in an X-SCID lentiviral trial), so this is not academic.

### 1.5 Net verdict

BT4's **philosophy is current and its rigor is ahead**: exact optimization + Pareto
frontier + content-hashed provenance + enforced invariants + *hard* manufacturability
guarantees the academic LMs cannot make. BT4 is **behind on** exactly three things,
in priority order: (1) **construct context** (5′UTR + backbone) — now first-class in
the field; (2) a **calibrated learned expression head** (wraps RiboNN but ships it
uncalibrated for its regime); (3) optionally, co-optimizing mRNA structure inside the
solver. Full sourcing at the end of this document.

---

## 2. How BT4 measures on the requested functions

The measured scorecard lives in
[`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md)
§3 (six requirements), §6 (GC), and §7 (construct context). In one line per axis:

- **Homopolymer, restriction sites (584-enzyme REBASE catalog, IUPAC + RC,
  content-hashed into the manifest), GC-run, forbidden sequences + presets,
  tandem/inverted repeats, dispersed max-repeat** — all **implemented and honest**;
  LOCAL rules are exact in the trellis, the GLOBAL max-repeat is refinement-enforced
  with residuals disclosed. Small wiring gaps remain (the IUPAC `extra_sites` path
  and a few repeat knobs are unreachable from the UI/CLI) — Phase 1.
- **Core objectives** (CAI with a declared, default-highly-expressed reference set;
  tAI on real GtRNAdb data; codon-pair bias; %MinMax; CpG/dinucleotide) — **strong
  and honest**, `delta==score` and `ok_suffix⇔validate` property-tested for every
  term and constraint.
- **Windowed GC** — **missing** (only a global count budget, a per-codon proximity
  prior, and a GC-run cap exist); this is the single most important missing
  *manufacturability* constraint (IDT/Twist spec GC per 50-bp window). Phase 1.
- **5′UTR / vector backbone (construct context)** — **absent**: `OptimizeConfig`
  carries no sequence outside the CDS, folding sees only `CDS[0:48]`, the splice CNNs
  are `N`-padded, and the initiator Kozak context is unreachable. This is the
  headline architectural gap (audit §7). Phase 2.
- **UI** — honest surfaces are excellent (certificate badge from the engine's own
  status, inline violation spans, calibrated/uncalibrated labels, ASSP consent), but
  coverage is thin: no objective weights / ε-budgets, no FASTA upload, per-site
  tracks are GC+CpG only, no validate-an-existing-sequence panel, no construct-context
  input. Phase 1 + Phase 3.

---

## 3. Shipped in this change — the four Tier-0 honesty defects, fixed

The audit's §4 measured four defects; three broke a §5 honesty invariant on `main`.
All four are fixed here, each with a regression test that did not previously exist.
Before/after is reproducible with the audit's own commands (its §10).

| # | Defect (audit §4) | Invariant | Fix in this change |
|---|---|---|---|
| A.1 | `run_frontier` reported GLOBAL rules but never enforced them, shipping `PROVEN_OPTIMAL` over dozens of hard violations (the Studio Optimize path) | §5 #6 | A frontier point that violates a GLOBAL rule now gets a **`RELAXED`** certificate naming the unenforced rule, with the residual in the audit; a clean point keeps its certificate. The frontier stays a pure explorer — *repair* lives in `run_optimize` and the candidate assembler. |
| A.2 | `run_validate` / `POST /validate` silently dropped every GLOBAL constraint, returning `feasible=True` on a violating sequence | §5 #2 | `run_validate` now audits the GLOBAL rules too; `bt4 validate --max-repeat-length N` catches the repeat it always should have. |
| A.3 | `folding_dg` was computed whole-sequence but labelled `5' dG` (the SA optimized the 48-nt window) | §5 #2 | The audit now reports the **same 5′ window the SA optimized** (`model.score_sequence`), so reported == computed. |
| A.4 | `avoid_internal_start` infeasible on most proteins; `InfeasibleError` named every active constraint, not the culprit; `relax()`/`RELAXED` promised but absent | §4.2 contract | New opt-in `relax()` (`bt4.domain.relax.SoftConstraint`, on `InternalStartConstraint`): an infeasible-but-relaxable rule degrades to a **`RELAXED`** result naming the dropped rule, with residuals still audited against the original hard rule. `InfeasibleError` now names the **failing residue and only the culprit constraints**. A constraint that does not opt in (e.g. restriction sites) is never silently dropped. |

Reproduce (before this change each line lied; after, each is honest):

```bash
# A.2 — was feasible:True, now flags the 12-nt repeat
bt4 validate ATGAAGGTTTCCTTATGAAGGTTTCC --max-repeat-length 4

# A.4 — was InfeasibleError blaming "homopolymer, internal_start"; now a RELAXED result
python -c "from bt4 import api; from bt4.pipeline.optimize import OptimizeConfig; \
r=api.optimize('MAAMG', OptimizeConfig(avoid_internal_start=True)); \
print(r.certificate.status.value, r.certificate.relaxed_terms)"
```

```python
# A.1 — no frontier point may be proven_optimal while violating a GLOBAL rule
from bt4 import api
from bt4.pipeline.optimize import OptimizeConfig
fr = api.frontier("KKKKKKKKKKKK", OptimizeConfig(max_repeat_length=6, max_homopolymer=None), steps=5)
assert not any(r.certificate.is_proven_optimal and (r.audit.get("max_repeat_residual") or 0) > 0
               for r in fr.results)

# A.3 — reported folding dG equals the optimized 5' window, not the whole sequence
from bt4.biomodels.folding import default, DEFAULT_FIVE_PRIME_WINDOW
res = api.optimize("MAALKHETQWYCDEFGHIKLMNPQRS", OptimizeConfig(refine=True, max_homopolymer=None))
m = default()
assert res.audit["folding_dg"] == m.five_prime_dg(res.dna, DEFAULT_FIVE_PRIME_WINDOW)
```

New regression tests: `tests/test_pipeline_api.py` (validate-enforces-global,
frontier-never-certifies-optimal-over-violation, folding-window, relax + culprit
naming), `tests/test_relax.py`, `tests/test_certificate.py`. Also fixed here: the
`gc_target` silent-no-op is now honestly documented across CLI/Studio/API (it is a
soft objective, active in a single solve only when `gc_weight > 0`; always swept on
the frontier), and the stale "gap-bounded" Lagrangian-budget docstrings are corrected
to `PROVEN_OPTIMAL`.

---

## 4. Forward roadmap

The live, precondition-tagged queue is [`NEXT_SESSION.md`](NEXT_SESSION.md); this is
the durable shape. Two maintainer decisions gate Phase 2 (recorded below).

- **Phase 1 — make the default defensible** (decision-independent). Windowed-GC
  constraint (LOCAL, exact in the trellis — cheaper and stronger than a refinement
  rule); reach the IUPAC `extra_sites` path so "ban the degenerate site directly" is
  actionable; regime-tagged application presets (`pipeline/presets.py`, **no default
  selected** — see Decision 1); Studio objective weights + ε-budgets, FASTA upload, a
  validate-an-existing-sequence panel, splice + 5′-folding per-site tracks; relabel
  the ramp term's claim. Presets depend on this change's `relax()` so a preset can
  never ship an infeasible `avoid_internal_start`.
- **Phase 2 — construct context** (the headline capability, strictly serial). A pure
  `ConstructContext` domain object + an additive `OptimizeConfig.context` (byte-
  identical when absent); oORF pairing across the UTR-CDS junction (cheapest,
  best-evidenced, no ML); a `SeededConstraint` wrapper so every LOCAL rule is
  evaluated across the junction (and the initiator Kozak becomes reachable); a shared
  `junction_window()` folding function + a `SplicePredictor.score_region()` fed **real
  flanks** instead of `N`; a whole-construct audit + restriction-site uniqueness with
  a precomputed k-mer index (so a 6 kb backbone is not re-scanned per SA move).
- **Phase 3 — construct-context UI + provenance display.**

### Maintainer decisions on record

1. **Regime — user-defined, no default.** BT4 stays regime-agnostic; regime-tagged
   presets (AAV/LVV transgene, plasmid/production, IVT mRNA) ship as *data* with none
   selected by default. CpG direction, 5′-structure sign, and length ceilings are
   preset fields, never hardcoded engine branches.
2. **Backbone hash in the manifest — ask per run.** When a user supplies a
   5′UTR/backbone, whether a content hash enters the run manifest is a per-run choice
   (privacy-safe default). Always-on rider: the ASSP network cross-check is
   hard-blocked from transmitting any flank/backbone bytes.

---

## Sources

- **Exact-optimization school:** LinearDesign — *Nature* 2023 (s41586-023-06127-z);
  EnsembleDesign / expected partition function (PMC12261492, arXiv 2401.00037).
- **Generative / LM optimizers:** CodonTransformer — *Nat Commun* 2025
  (s41467-025-58588-7); CodonBERT — Genome Research 2024 (github FPPGroup/CodonBERT);
  CodonTranslator — bioRxiv 2025.11.24.690310; DeepCodon — ScienceDirect
  S2693125725000433; CodonRL — PMC12918928; ColiFormer — MDPI bioengineering 13/1/114.
- **5′UTR subfield:** UTR-Insight (BMC Genomics 2025, PMC11796101); UTRGAN (Bioinf
  Adv 2025); UTailoR (iScience 2025, PMC12506572); coordinated-design review
  (ScienceDirect S2090123226004959).
- **Expression / ribosome load:** RiboNN / cell-type codon optimization — *NAR* 2025
  gkaf233 (PMID 40156867); deep-generative — *Nat Commun* 2025 s41467-025-64894-x.
- **Splice:** Pangolin (PMID 35449021); SpliceAI Lookup (Broad); lentiviral cryptic-
  splice HMGA2 truncation (PMC9240040).
- **Industry & benchmarking:** IDT & Twist codon-tool pages; multi-criteria
  comparative framework (PMC12010093); mammalian glycoprotein codon-strategy eval
  (bioRxiv 2026.03.18.712111).
