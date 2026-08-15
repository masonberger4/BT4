# BT4 review — expression, construct context, and the honesty ledger

**Date:** 2026-08 · **Tree reviewed:** `3805c83` · **Scope:** the whole product,
audited against what BT4 is *supposed* to be.

> **What this document is.** A measured audit, not an opinion piece. Every number
> below was produced by *running the code* in this tree — the CLI, `bt4.api`, and
> the committed comparison harness — not read from a docstring. Where the code and
> the docs disagree, the code wins and the disagreement is recorded. Reproduction
> commands are in [§10](#10-reproducing-every-number).
>
> **What this document is not.** A changelog, and not a statement of shipped
> behaviour changes. Nothing here was fixed by writing it. The forward queue lives
> in [`NEXT_SESSION.md`](NEXT_SESSION.md); the durable rules live in
> [`../CLAUDE.md`](../CLAUDE.md). This review is the evidence those two documents
> are re-pointed from.

---

## 0. The requirement being audited

BT4 should back-translate a protein into a coding sequence that is:

1. **optimized for protein expression** (the primary objective),
2. subject to **user-specified sequence restrictions** — homopolymer length,
   restriction sites, GC runs, arbitrary forbidden sequences, repeats of a given
   length — the rules that decide whether a construct can be *synthesised*,
   *propagated as a plasmid*, *packaged into AAV or LVV*, and *expressed*,
3. aware of the **5′UTR** that sits in front of the CDS in the real transcript,
4. aware of the **vector backbone** the CDS is cloned into,
5. aware of **splice sites**,
6. and delivered through a **simple UI that reports every metric honestly** and
   hands back the final sequences.

---

## 1. Verdict

**BT4 is an unusually well-engineered, unusually honest CAI optimizer with
manufacturability constraints. It is not yet an expression optimizer.**

The reason is structural, and it is one sentence: **the optimizer only ever sees
the CDS.**

- Folding is scored on `CDS[0:48]` with no leader in front of it.
- The splice CNNs are scored on the CDS floating in 5,000 literal `N` bases.
- RiboNN is the one model that *can* accept UTRs — and it lives in a side tab
  where it is structurally forbidden from influencing the delivered sequence.

Nothing that knows anything about expression is allowed to touch the answer. What
ships by default is an exactly-solved, well-provenanced, honestly-certified
**CAI maximizer with a homopolymer cap**.

That is a narrower claim than `CLAUDE.md` §1 makes ("optimized for real expression
outcome"), and — to the project's credit — exactly the claim `README.md` makes
("maximize codon adaptation"). The README is right.

---

## 2. What is genuinely strong

This section is not politeness. These are the parts that should be protected while
the rest changes.

| Strength | Evidence |
|---|---|
| Exact codon-trellis DP with a real optimality certificate, over the **true** per-constraint context with no global cap | `optimize/exact_dp.py:144` |
| Content-hashed provenance — table bytes, model SHAs, git SHA, seed; reproducible from the stamp | `provenance/manifest.py` |
| Nine organism codon tables, all recounted from release-pinned public CDS sets, with **declared reference sets** (highly-expressed by default) | `biomodels/codon/data/`, `scripts/build_organism_tables.py` |
| Calibration honesty enforced *structurally*: every wrapped model ships `calibrated=False` and cannot steer delivery until it earns the flag | `pipeline/rerank.py:72`, `pipeline/candidates.py:354` |
| 584-enzyme REBASE-derived restriction catalog — IUPAC-aware, reverse-complement-aware, exact in the DP, content-hashed and re-derivable | `constraints/restriction.py`, `scripts/build_enzyme_catalog.py` |
| Every metric recomputed from the delivered DNA rather than trusted from an accumulator | invariant #2, throughout |
| **1,185 non-GUI tests pass**, including the `ok_suffix⇔validate` and `delta==score` property tests | `tests/` |

*(The 51 GUI tests in `tests/test_app_smoke.py` could not run in the review
sandbox — it lacks `libEGL`, so PySide6 will not import. They were skipped, not
failed. The Rust extension is likewise not built here, so all timings below are the
pure-Python fallback and would only improve.)*

No surveyed commercial tool exposes an optimality certificate, and most are
explicitly non-reproducible. On **method and rigor**, BT4 is ahead of the field.
The gap is entirely in *what it optimizes* and *what it can see*.

---

## 3. Scorecard against the six requirements

| Requirement | State | Anchor |
|---|---|---|
| Homopolymer limit | **Works** — exact in the DP, on by default (6) | `constraints/rules.py:30` |
| Restriction sites | **Works** — exact, both strands, IUPAC | `constraints/restriction.py:190` |
| GC-run limit ("max GC length") | **Works** — exact; off by default | `constraints/gc_run.py:56` |
| Forbidden sequences | **Works** — exact; 4 named presets | `constraints/rules.py:105`, `constraints/forbidden.py:47` |
| Tandem / inverted repeats | **Works** — exact | `constraints/repeats.py:43,139` |
| Dispersed repeat length | **Best-effort — and a silent no-op in two paths** | §4.1, §4.2 |
| GC content | **No GC-content constraint exists** | §6 |
| Windowed GC | **Missing** — promised in `CLAUDE.md` §6, never built | §6 |
| **Optimize for expression** | **No** — the default objective is CAI alone | §5 |
| **5′UTR** | **Absent** from the optimizer entirely | §7 |
| **Vector backbone** | **Absent everywhere** in `src/` | §7 |
| Splice sites | Opt-in motif heuristic + an uncalibrated CNN audit run on the wrong context | §7 |
| Simple, honest UI | Honest: **yes**. Complete: **no**. Simple: **no** | §8 |

Six of the six user-specified restrictions are genuinely implemented and exactly
enforced — that half of the requirement is in good shape. The failures are
concentrated in expression, context, and the surfaces.

---

## 4. Four defects, each measured

### 4.1 A constraint the GUI offers that the GUI's own code path does not enforce

BT4 Studio's **Optimize** button calls `api.frontier` (`app/worker.py:123`), and
`run_frontier` only *reports* the dispersed-repeat rule — it does not enforce it
(`pipeline/optimize.py:1026-1029`).

```
run_optimize (CLI path)     max_repeat_length=10 → longest repeat  9 nt
                                                   cert=heuristic, enforced=clean, residual=0

run_frontier (Studio path)  max_repeat_length=10 → longest repeat 58 nt
                                                   cert=proven_optimal, 96 violations
                                                   max_repeat_enforced=None
```

`app/theme.py:111` maps `proven_optimal` to green `#1e8b4d`. So a user who sets
*Max repeat length* in the GUI is shown a **green "proven optimal" badge** on a
sequence carrying a 58-nt perfect direct repeat they explicitly asked to remove,
beside a metrics row reading `Hard violations: 96`. `avoid_uorf` behaves
identically in that path — `uorf_enforced` returns `None`, meaning the rule was
never applied at all.

This is the failure mode `CLAUDE.md` §10.1 exists to prevent: *no silent optimality
loss*. The certificate is not lying about the DP it ran; it is lying by omission
about the problem the user posed.

### 4.2 `validate` silently drops every GLOBAL constraint

`run_validate` builds only the LOCAL constraint set (`pipeline/optimize.py:1105`
calls `_build_constraints`, never `_build_global_constraints`). Measured on a
sequence containing a 9-A homopolymer **and** a 31-nt exact direct repeat:

```
validate(max_homopolymer=6)    → 1 violation    (LOCAL — works)
validate(max_repeat_length=8)  → 0 violations   (GLOBAL — silent no-op)
validate(avoid_uorf=True)      → 0 violations   (GLOBAL — silent no-op)
```

`POST /validate` inherits the same behaviour. A user auditing a sequence for long
repeats gets a clean bill of health for a sequence with a 31-nt repeat in it, with
no warning that the check did not run. **No test covers this path.**

### 4.3 `folding_dg` is computed over the whole CDS and labelled `5' dG`

`pipeline/optimize.py:827` calls `model.five_prime_dg(result.dna)` with no
`window`, so `window=None` and `five_prime_window` returns the **entire sequence**
(`biomodels/folding/base.py:65-66`). Meanwhile the SA refinement optimized
`model.score_sequence(dna)`, which folds the **48-nt** window. The CLI then prints
the result as `folding 5' dG` (`cli/__main__.py:293`).

Measured on a 156-nt CDS with the baseline backend:

```
what the SA actually optimized (48-nt window) : -39.0
what the audit reports and the CLI prints     : -138.0   ← labelled "5' dG"
```

Reported ≠ computed. That is a direct violation of invariant #2, shipping today.
With the uncalibrated baseline it is additionally an O(n³) Nussinov over the full
length — which `biomodels/folding/baseline.py:25-27` explicitly warns against.

### 4.4 `avoid_internal_start` is infeasible on most real proteins, and the error names the wrong constraints

Methionine has exactly one codon, so an internal `ATG` cannot be moved. If the
**preceding** residue's codons all begin with a purine (G, A, V, E, D start with G;
I, M, T, N, K start with A) and the **following** residue's codons all begin with G
(G, A, V, E, D), then both strong-Kozak conditions are *forced* and no synonymous
assignment can satisfy the constraint.

Minimal reproducer: **`MAAMG`**. Ala is `GCN` (so −3 is a forced G) and Gly is
`GGN` (so +4 is a forced G).

Measured over random proteins:

```
length  100 aa  → infeasible on  35%
length  200 aa  → infeasible on  62%
length  400 aa  → infeasible on  82–100%
length  700 aa  → infeasible on 100%
```

The GUI exposes this as a plain checkbox labelled *Internal ATG*. Ticking it on a
700-aa protein raises:

```
InfeasibleError: no feasible codon under constraints: homopolymer, internal_start
```

`InfeasibleError` is constructed from **every active constraint**
(`optimize/exact_dp.py:183`), not the conflicting one, and carries no residue
position — so the user cannot tell which box to untick. `tandem_unit=3` fails on
~25% of 400-aa proteins the same opaque way.

`CLAUDE.md` §4.2 and §6 promise the opposite behaviour: *"Hard constraints degrade
gracefully via `relax()` and report which constraints conflict rather than aborting
with 'no feasible codon'."* **`relax()` does not exist anywhere in the codebase**
and is not on the `Constraint` protocol (`domain/contracts.py:87-110`).
`OptimalityStatus.RELAXED` is defined (`domain/certificate.py:33`) and never used.

---

## 5. What a default run actually optimizes

```
bt4 optimize <protein>   →   maximize  1.0 × Σ log w(codon)
                             subject to  max_homopolymer = 6
```

Every other term is weight `0.0` or off (`pipeline/optimize.py:218-256`). `CaiTerm`
is LOCAL with `context_len() == 0` (`objectives/terms.py:57-59`), so on a default
run the objective decomposes **per codon with no coupling at all** — the exact
trellis is solving a nearly trivial problem, and only earns its keep once a
PAIRWISE, POSITIONAL, or real LOCAL constraint is switched on.

CAI comes back at exactly **1.000**. On a protein containing an internal repeat,
back-translating deterministically to the top codon reproduces that repeat in DNA:
a 50-residue test protein with a duplicated 25-residue block yields a **58-nt
perfect direct repeat** — a synthesis-vendor flag and a recombination risk during
plasmid propagation.

### The whole industry deliberately does not do this

`python scripts/compare_tools.py` — BT4's own committed harness, on the KRas4B
panel from Ranaghan et al. 2021 (CC BY 4.0), every metric recomputed by BT4:

| name | optimizer | CAI | tAI | GC % | CpG | max homo |
|---|---|---|---|---|---|---|
| Native | native human KRAS CDS | 0.592 | 0.315 | 37.7 | 7 | 5 |
| IDT | Integrated DNA Technologies | 0.632 | 0.329 | 45.2 | 35 | 8 |
| OPTIMIZER | OPTIMIZER | 0.666 | 0.353 | 48.0 | 44 | 18 |
| GeneArt | Thermo GeneArt | 0.671 | 0.328 | 42.5 | 28 | 8 |
| JCAT | Java Codon Adaptation Tool | 0.702 | 0.378 | 46.4 | 35 | 18 |
| Twist | Twist Bioscience | 0.716 | 0.350 | 51.5 | 47 | 4 |
| COOL | Codon Optimization OnLine | 0.720 | 0.339 | 48.1 | 45 | 19 |
| GeneWiz | GENEWIZ (Azenta) | 0.748 | 0.357 | 54.0 | 54 | 4 |
| GenScript | GenScript (OptimumGene) | 0.831 | 0.386 | 54.0 | 38 | 4 |
| **BT4** | **BT4 (this tool)** | **1.000** | **0.413** | **62.4** | **47** | **4** |

**No shipping tool maximizes CAI.** Nine independent tools cluster at CAI 0.63–0.83
and GC 42–54%. BT4's default sits alone at 1.000 / 62.4%.

That is not nine tools being sloppy. It is nine tools respecting the literature
this repo already cites in
[`RESEARCH_codon_optimization_SOTA.md`](RESEARCH_codon_optimization_SOTA.md) §1 —
Kudla et al. (*Science* 2009), Welch et al. (*PLoS ONE* 2009), and the
over-optimization / co-translational-misfolding results. **The engine is not
wrong. The default operating point is.**

---

## 6. GC: no constraint exists, and neither control can do the job

There are ten constraint classes — Homopolymer, Forbidden, GcRun, Restriction,
Tandem, Inverted, MaxRepeat, Kozak, uORF, SpliceMotif. **None of them constrains GC
content**, windowed or otherwise (`grep gc_window src/` returns nothing).

GC has exactly two controls, and I measured both on KRas4B:

**The soft objective term saturates and never reaches its own target.**

```
gc_target = 0.50
gc_weight =  0   →  GC 62.4%   worst 50-nt window 76%   CAI 1.000
gc_weight =  1   →  GC 61.9%   worst 50-nt window 76%   CAI 0.996
gc_weight =  2   →  GC 56.4%   worst 50-nt window 66%   CAI 0.933
gc_weight = 50   →  GC 56.4%   worst 50-nt window 66%   CAI 0.933   ← saturated
```

`GcProximityTerm` is *separable* — `context_len() == 0`
(`objectives/terms.py:99-105`) — so every codon independently picks its own
closest-to-target synonym. It cannot trade a forced high-GC codon (Gly `GGN`, Pro
`CCN`, Ala `GCN` are all ≥ 66% GC) against a lower-GC one elsewhere. A "GC target"
that provably cannot reach its target is a modelling defect, not a tuning problem.
It is also **inert by default** (`gc_weight = 0.0`), so `--gc-target 0.70` alone
changes nothing.

**The hard count budget reaches the target but cannot control clustering.**

```
gc_max = 50%  →  GC 49.9%   50-nt window 36–74%   CAI 0.915   proven_optimal
gc_max = 45%  →  GC 45.0%   50-nt window 30–62%   CAI 0.858   proven_optimal
```

It hits the target exactly, with a proven-optimal certificate — genuinely good. But
a whole-sequence count says nothing about *where* the GC sits: at 50% total GC the
worst 50-nt window is still **74%**. It is also absent from the GUI, ignored by
`run_frontier`, and mutually exclusive with refinement.

**Consequence, measured** — default output over 30 random 300-aa proteins:

```
50-nt GC window above 75% GC   24/30    (worst window 82%)
longest repeat ≥ 20 nt          0/30
homopolymer > 6                 0/30
```

Every failure is on the high side: CAI-max against a human highly-expressed table
drives GC3 up, and nothing constrains where that GC clusters. Native KRAS sits at
37.7% GC with a 20–56% window; BT4's default at 62.4% with a 44–76% window.

A windowed-GC rule is the most common synthesis-vendor screen, is **already listed
in `CLAUDE.md` §6** ("GC / GC-window … bounded window"), and would be **LOCAL and
exact in the trellis** (`context_len = window − 1`). The windowed computation even
exists already on the reporting side — `bt4 tracks` prints a 50-nt GC window
(`pipeline/tracks.py`). Only the constraint is missing.

### The engine is strong; nothing leads the user to it

The same protein, three configurations:

```
native human KRAS CDS                    CpG   7
bt4 default                              CpG  47   CAI 1.000  GC 62.4%   proven_optimal
bt4 --cpg-weight 1 (deplete)             CpG   0   CAI 0.923  GC 54.9%   proven_optimal
bt4 --cpg-max 10 (hard count budget)     CpG  10   CAI 0.948  GC 56.6%   proven_optimal
```

Driving CpG from 47 to 0 costs 0.077 CAI and **keeps a proven-optimal
certificate** — better than anything the vendor tools expose. But for an AAV or
lentiviral transgene, where CpG content drives innate sensing and transgene
silencing, the default ships **6.7× the native CpG count**, and neither the CLI nor
the GUI hints that one flag would fix it.

Configured deliberately, BT4 is excellent. A plausible AAV/LVV profile
(`cai_weight=1, gc_weight=2, gc_target=0.50, cpg deplete, max_homopolymer=6,
max_gc_run=6, max_repeat_length=12, avoid_splice_sites, avoid_uorf`, poly-A + TATA
presets) on a **700-aa** protein:

```
20.7 s   CAI 0.882   GC 53.0%   50-nt window 42–64%   CpG 5
zero hard violations · repeats clean · uORFs clean
```

That is a defensible, vendor-safe, CpG-depleted transgene — produced entirely by
knobs that already ship. **No user will ever find that combination.** It is nine
flags with no guidance, three of which the GUI does not expose. That is the whole
argument for named application presets.

### Runtime is not the constraint

Default solve: 300 aa in 0.05 s, 700 aa in 0.12 s, 1,200 aa in 0.19 s. With heavy
rules on, a 400-aa protein takes ~2.4 s (`avoid_splice_sites`), ~2.7 s (six
enzymes), ~7.9 s (`max_repeat_length=8`, refinement-bound). Comfortably
interactive, before the Rust extension.

---

## 7. The architectural gap: the CDS is optimized in a vacuum

`OptimizeConfig` has **40 fields and not one of them is sequence outside the CDS**
(`pipeline/optimize.py:218-256`). The DP seeds from an empty prefix —
`optimize/exact_dp.py:161`: `layer = {"": (0.0, "")}`. A search of `src/` for
`genbank`, `plasmid`, `vector`, `backbone` returns nothing.

Seven verified consequences:

1. **Folding sees the wrong 48 nucleotides.** It folds `CDS[0:48]`
   (`biomodels/folding/base.py:64-69`). The duplex that actually occludes 40S
   scanning forms across the cap-proximal leader ↔ start-codon junction — sequence
   BT4 has no field to accept. Every SA move at codon 17 or beyond has *zero*
   folding effect, while moves are proposed uniformly across the sequence.
2. **Kozak is half-modelled.** `InternalStartConstraint` scores *internal* ATGs
   only. The real initiator's −3 base lives in the 5′UTR and is invisible; the +4
   base is the first base of codon 2, which BT4 already chooses. The most
   actionable initiation signal in the whole design is unreachable.
3. **SpliceAI and Pangolin are fed 5,000 `N` bases on each side**
   (`biomodels/splice/spliceai.py:387`, `pangolin.py:484`; `N` maps to an all-zero
   one-hot at `:141` / `:192`). In the user's construct those flanks are real
   vector sequence, often carrying real splice signals. This is SpliceAI's
   documented custom-sequence convention, but it is a nucleotide vacuum the
   networks essentially never saw in training. The two-backend **agreement** check
   — BT4's headline uncertainty signal — cannot detect a context artefact both
   backends share.
4. **Δsplicing's reference defaults to BT4's own output**
   (`pipeline/splice_audit.py:92`), so the delivered sequence scores exactly `0.0`
   added risk by construction. The `SplicePredictor` contract's framing ("the
   *added* splice risk a synonymous redesign introduces relative to a reference",
   `biomodels/splice/base.py:6-8`) implies a natural baseline the pipeline never
   supplies.
5. **Restriction-site logic is absolute, not contextual.** BT4 can only ban a site
   outright. It cannot answer the question a cloner actually asks: *is this site
   still unique in my plasmid?*
6. **Junction-spanning violations are invisible.** A forbidden motif or enzyme site
   formed across the leader ↔ codon-1 boundary is unreachable by `ok_suffix`,
   because the prefix starts empty.
7. **RiboNN's UTR fields reach no constraint.** They exist only in the Candidates
   tab (`app/studio.py:1000-1019`), feed an uncalibrated head, and are explicitly
   *"never exported"*. `resolve_backend` defaults both UTRs to `""`
   (`biomodels/expression/__init__.py:112-113`) and neither `rerank.py` nor
   `candidates.py` ever passes them — so on every path BT4 itself drives, RiboNN
   either is not reached or would raise.

`SplicePredictor.score_sequence(dna)` takes exactly one sequence
(`biomodels/splice/base.py:200`) — there is no place to pass flanking context at
all. Making the models context-aware is a **contract change**, not a parameter.

The project's own design-of-record already names this boundary:
[`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md) lists
*"Beyond-CDS design (5′UTR / initiation-region) — where the literature says most of
the expression signal lives"* under **Future / out of scope for v1**. This review
is the argument for moving it in.

---

## 8. The surfaces

### The frontier shows options the user cannot take

On KRas4B the default frontier is real and informative:

```
0: CAI 0.9331  GC 56.4%
1: CAI 0.9790  GC 60.1%
2: CAI 0.9959  GC 61.9%
3: CAI 1.0000  GC 62.6%   ← delivered
```

The delivered point is always the CAI corner, and BT4 Studio has **no
frontier-point picker**. The plot shows the user three better-balanced designs and
gives them no way to select one.

### BT4 Studio

- **27 undifferentiated controls** in one flat column, no grouping and no
  basic/advanced split. Exactly **two** are ones a bench scientist must set
  (Protein, Organism); a third (GC target) is usually set. The other 24 are
  constraint and solver knobs.
- **Protein input is paste-only.** `QFileDialog` appears three times in the app,
  all `getSaveFileName` — there is **no file-open dialog and no drag-and-drop
  anywhere**. Pasted FASTA *text* is tolerated and record 0 silently taken.
- **No GenBank reader or writer exists in the package.** `io/` is FASTA + JSON
  only. `CLAUDE.md` names GenBank three times — the branch-history table (`:105`),
  the `io/` architecture diagram (`:148`), and the export promise (`:463`) — while
  the README correctly does not. GenBank I/O is already in the *intended*
  architecture and was simply never built. Bench scientists clone from annotated
  vector maps.
- **Seven engine capabilities are CLI/service-only** and unreachable from the GUI:
  the CpG/UpA count budget (a headline README feature), the GC budget, the entire
  5′-ramp frontier axis, `--refine` and therefore all folding reporting, codon-pair
  *de*optimization (negative `cpb_weight`, the attenuated-vaccine use case),
  `tandem_copies` / `inverted_loop`, and `seed`.
- **The metrics table is a hard-coded 9-row list** (`app/studio.py:62-75`), so
  audit keys the CLI prints never appear in the GUI: `tai`, `cg_count`/`ta_count`,
  and — most importantly — `max_repeat_enforced` / `max_repeat_residual` and
  `uorf_enforced` / `uorf_residual`, the two "partial enforcement" notes. In the
  GUI a partially-enforced global rule surfaces only as an unexplained violation
  count.
- **The honesty machinery that *is* present is genuinely good**: the certificate
  badge colour is derived from the engine's own status and cannot drift; ASSP
  numbers are regression-tested never to reach an export; every uncalibrated score
  sits beside a `Calibrated` column under a banner stating the order is discovery
  order, not a ranking. I found **no uncalibrated number displayed without a
  caveat**. The honesty gap is the opposite one — caveats the GUI never shows
  because its metrics table is fixed.

### Smaller gaps

- **The stop codon is not selectable.** It is optimized through the trellis like
  any other codon, so you get the highest-CAI stop. Invariant #8 discusses a
  "caller-pinned/overridden stop" but no config field exists. Tandem stops
  (`TAATAA`) to suppress readthrough — routine in transgene design — cannot be
  expressed.
- **`bt3_synthesis_artifacts`** ships six motifs inherited from BT3 with no cited
  source, in a repo that otherwise hash-pins and cites every dataset.
- **`avoid_reverse_complement`** (default `True`) is exposed only by the FastAPI
  service, not the CLI or GUI. For a sense-strand rule like a poly(A) signal in a
  lentiviral genome, banning the reverse complement over-constrains.
- **`extra_sites`** exists on `RestrictionSiteConstraint` (`restriction.py:205`)
  but is wired to nothing, so a user cannot supply a recognition sequence the
  catalog lacks — and `forbidden_motifs` rejects IUPAC, so there is no workaround.
- **`run_frontier` ignores the GC budget** and says so only in a code comment
  (`pipeline/optimize.py:1013-1014`); an API caller gets no warning.
- **Doc/code drift:** the Lagrangian budget route is documented as "gap-bounded
  (not proven-optimal)" (`pipeline/optimize.py:673-675`) but returns
  `PROVEN_OPTIMAL` (`optimize/lagrangian.py:364`). Understated rather than
  overstated — but still wrong.

---

## 9. Roadmap

Ordered by expected value. Tier 0 fixes things that are wrong *today*; Tier 1 is
the cheapest large improvement; Tier 2 is the actual ask.

### Tier 0 — honesty defects (each violates a §5 invariant)

| # | Item |
|---|---|
| 0.1 | **`run_frontier` must enforce GLOBAL rules or refuse a config that sets one.** Suggested shape: refine the *delivered* point (degrading that point's certificate to `HEURISTIC`) while other frontier points stay exact and honestly labelled — the badge becomes per-point. **No frontier point may claim `proven_optimal` while violating a GLOBAL rule.** §4.1's reproducer becomes the regression test. |
| 0.2 | **`run_validate` / `POST /validate` must apply GLOBAL constraints.** Today they are silent no-ops. Add the missing test. |
| 0.3 | **`folding_dg` reported == optimized.** Fix structurally — one shared window function used by both the objective and the audit — not by making two call sites happen to agree. |
| 0.4 | **Diagnostic `InfeasibleError` + `relax()`.** On an empty layer, re-test the failing residue's codons against each constraint individually (O(codons × constraints), paid only on failure) and name the residue and the culprit. Then implement `relax()` using the already-defined but unused `OptimalityStatus.RELAXED`, so `avoid_internal_start` degrades to a reported soft violation where Met makes it unsatisfiable. `MAAMG` is the golden test. **Blocks 1.2.** |

### Tier 1 — make the default defensible

| # | Item |
|---|---|
| 1.1 | **Windowed-GC constraint** — LOCAL, `context_len = window − 1`, exact in the trellis. Promised in §6 of the constitution, absent from the code, the top vendor screen, and the windowed computation already exists in `pipeline/tracks.py`. |
| 1.2 | **Application presets** — `mammalian_plasmid`, `aav`, `lentiviral`, `mrna_ivt`, `ecoli`. A preset sets objective weights *and* constraints, so it produces an `OptimizeConfig` and belongs in `pipeline/presets.py`, mirroring the `constraints/forbidden.py` catalog pattern. Four honesty rules: the manifest carries the **resolved field values**, never the preset name alone; explicit user knobs win and the report says which were overridden; a test asserts every override key is a real config field; and no preset may enable `avoid_internal_start` until `relax()` exists. The `mrna_ivt` preset shows why presets earn their place — it turns `avoid_splice_sites` **off**, because a cytoplasmically delivered mRNA never meets a spliceosome and charging for cryptic splice sites there is a pure feasibility tax. Encoding that judgment *with its reason* is what a preset adds over raw knobs. |
| 1.3 | **Frontier-point picker** — the plot already shows designs the user cannot select. |

### Tier 2 — the construct-context model (architectural)

The framing: **BT4 models a CDS; it needs to model a CDS *sited in a construct*.**
One value object, then a per-rule decision — seeded prefix, whole-construct audit,
or both.

| # | Item |
|---|---|
| 2.1 | **`ConstructContext` in `domain`** — upstream, downstream, topology, `complete`, and `masked_spans`. Pure and dependency-free (reverse-complement stays in `io`). Circularity handled by linearising with the CDS in the *middle*, so both junctions are interior to a linear scan. Accepts `ACGTN` and truncates each flank at the nearest `N` — which keeps `N` out of every constraint without touching one. |
| 2.2 | **Wrap the constraints; do not seed the DP layer.** A `SeededConstraint` passes `upstream_tail(k) + prefix` to `ok_suffix` while the trellis accumulator stays CDS-only — so the emitted DNA, the layer key, the tie-break and the Rust path are all untouched, and with no context the path is byte-identical to today. Two shipped constraints are **not** translation-invariant (`kozak.py` skips index 0; `uorf.py` computes frame and a 5′ window from index 0) and must declare a `cds_offset`, enforced by a registry test. |
| 2.3 | **Objectives stay CDS-local, deliberately.** Seeding `delta` without seeding `score` would break invariant #4. The codon-pair term correctly returns 0 at position 0: a 5′UTR has no reading frame, so there *is* no codon pair at the start codon — refusing to score it is the honest answer, not a limitation. Report whole-construct dinucleotide content as a separate recomputed metric instead. |
| 2.4 | **3′ junction** — `solve_exact(..., tail_context)` filters the final layer against the downstream flank. Exact (the layer key is a sufficient statistic) and one pass. It does widen `ok_suffix`'s implicit "`next_codon` is 3 nt" assumption; every shipped implementation is already length-generic, but that must become a stated contract plus a property test. |
| 2.5 | **Whole-construct audit** — repeats and motifs *between insert and backbone*, where plasmid recombination actually comes from. `masked_spans` is load-bearing: AAV ITRs are a 145-nt palindrome and LVV LTRs are duplicated by construction, so without masking the repeat report is pure noise for exactly the two systems that need it most. ⚠️ Re-scanning a 6 kb backbone per SA move is the §10.8 quadratic trap — precompute the backbone k-mer index once, with a CI scaling assertion as a release blocker. |
| 2.6 | **Restriction-site uniqueness** — the genuinely novel capability; no surveyed tool offers it. It *compiles down*: if the backbone has exactly one site, "unique in the plasmid" reduces to "zero in the insert" — a plain seeded LOCAL constraint, still `proven_optimal`. If the backbone has zero and you want one, it is a count budget, which the shipped amount-bucketed DP already solves exactly. If the backbone has two, refuse honestly — the CDS cannot fix it. Ship a per-enzyme cloning report either way. |
| 2.7 | **Fold the real 5′UTR ↔ CDS junction** (`utr5[-45:] + cds[:48]`), through one shared window function used by both objective and audit, so invariant #2 is structurally enforced. |
| 2.8 | **`SplicePredictor.score_region(sequence, region)`** — not optional `upstream=`/`downstream=` kwargs, which would let a backend silently ignore them and return the wrong answer. `score_sequence` becomes `score_region(dna, (0, len(dna)))`, so no backend breaks. Then feed real flanks instead of 5,000 `N`s. ⚠️ **This collides with a queued calibration gate**: the bit-for-bit fidelity attestation rests on the N-padded path, so the flanked path must ship as a separately labelled mode and **cannot inherit that attestation**. |
| 2.9 | **Δsplicing reference policy** — native CDS → naive back-translation → self (labelled). The GenBank flow makes the first free: the feature you chose to replace *is* the wild-type. The reference must be scored in the same flanking context, or the Δ conflates a CDS change with a context change. |
| 2.10 | **Real initiator Kozak** from the supplied 5′UTR: −3 is `upstream[-3]` (fixed), +4 is the first base of codon 2 (BT4's to choose). Report always; prefer G at +4 only as a *soft* preference — residue 2 may have no G-initial codon (Lys, Phe). Split uORFs into fixable (ATG in the CDS) and fixed-context (ATG in the UTR); reporting violations the user cannot act on is how a report becomes unread. |

### Tier 3 — I/O and UI

| # | Item |
|---|---|
| 3.1 | **stdlib GenBank reader** in `io/`, keeping the layer's no-dependency policy. Parse `LOCUS` topology (load-bearing for circularity) and features with `complement`/`join`/`order`/partial locations; **refuse loudly** on remote accessions and `gap()` rather than mis-placing a CDS. Accept `ACGTN`, never coerce. |
| 3.2 | **GenBank writer** — the assembled construct with coordinate-shifted original features, a `CDS` feature carrying certificate and config hash, and **residual GLOBAL violations emitted as `misc_feature` annotations**, so honest residuals travel into the map the user opens in SnapGene instead of dying in terminal output. This is the highest value-per-line item in the roadmap. No wall-clock date (invariant #7). |
| 3.3 | **SnapGene `.dna` stays out of scope**, with a precise error telling the user to export GenBank from SnapGene. It is an undocumented, reverse-engineered binary format; a hand-rolled parser is a correctness liability this project's honesty posture cannot back. Revisit later as an optional extra wrapping a third-party reader — hash-pinned, never vendored, the posture already used for ViennaRNA / SpliceAI / RiboNN. |
| 3.4 | App: open a vector map, pick the feature to replace from a clickable feature table, export an **annotated construct** rather than a bare CDS. |
| 3.5 | App: basic/advanced split; protein file open and drag-drop; drive the metrics table from the audit dict instead of a fixed 9 rows; surface the seven hidden engine capabilities. |

### Tier 4 — biology the four target systems need

Independent of Tier 2, so it can land early for visible value.

| Item | Placement | Note |
|---|---|---|
| **Functional** poly(A) signal — hexamer **plus** a downstream GU/U-rich element within ~10–30 nt | LOCAL, exact (`context_len ≈ 36`) | The current preset bans the bare hexamer, which occurs by chance every ~4 kb — a feasibility tax with a poor hit rate. The full predicate still fits a bounded window. Catastrophic inside an LVV genome. |
| Sense-strand splice **donor × downstream-acceptor pairing** | LOCAL motif (exists) + GLOBAL reporting | The existing motif constraint is already sense-strand-only, exactly right for LVV. The new capability: enumerate acceptors in the *fixed* downstream flank once, then report each CDS donor as a paired hazard with the span it would splice out. Reporting-only while every backend is `calibrated=False`. |
| AAV packaging-size accounting (ITR-to-ITR vs ~4.7 kb; ~2.4 kb for scAAV) | **reporting only** | BT4 controls zero levers — CDS length is fixed by the protein. Wiring it as a pseudo-constraint would be a §10.4 inert-constraint hack. As a report it is genuinely actionable. |
| Uridine depletion for IVT mRNA | LOCAL additive term (`delta == score` trivially) | Offer total-U *and* third-position-only and make the user choose; both are defensible and third position is where synonymous freedom lives. Also available as an exact whole-sequence budget via the shipped bucketed DP. |
| m1Ψ slippery-sequence avoidance | LOCAL motif preset + a **per-base** homopolymer limit | `HomopolymerConstraint(6)` cannot express "U runs ≤ 3 but other runs ≤ 6". Applies **only** to m1Ψ-modified mRNA — meaningless for a plasmid, which is exactly the kind of context a preset makes legible. |
| Codon optimality / CSC (mRNA stability) | LOCAL additive term | Distinct from tAI; a validated axis BT4 does not carry. |

### Tier 5 — the flagship, separate track

**LinearDesign-class joint codon + secondary-structure optimization** in the
trellis. The field's strongest validated in-vivo result and architecturally native
to BT4 — but a research project, not a PR, and it deserves its own design doc
(already identified in
[`RESEARCH_codon_optimization_SOTA.md`](RESEARCH_codon_optimization_SOTA.md) §4).
Note the dependency runs the right way: joint folding that knows the real 5′UTR is
strictly better than joint folding that does not, so **context first**.

### Sequencing

Tier 0 is independent and lands first. Tier 1 is independent of the context work
and delivers most of the visible improvement for the least risk — with one hard
edge: **1.2 cannot ship before 0.4**, or the presets will enable a constraint that
fails on most proteins. Tier 4 is fully independent. Tier 2 is strictly serial
(2.1 → 2.2 → 2.4 → the rest), with 2.4 the riskiest single change and worth its own
PR. Tier 3's reader runs in parallel with Tier 2 and joins it at the
insertion-site step.

---

## 10. Reproducing every number

```bash
pip install -e '.[dev,app]'
QT_QPA_PLATFORM=offscreen python -m pytest tests/ -q       # baseline

python scripts/compare_tools.py                             # §5 the CAI/GC outlier table

# §4.1 — the frontier path does not enforce max_repeat_length
python -c "
from bt4 import api; from bt4.pipeline.optimize import OptimizeConfig
p='MAKQLEDKVEELLSKNYHLENEVARLKKLVKLLEDKVEELLSKNYHLENEV'
c=OptimizeConfig(max_repeat_length=10)
print('optimize:', api.optimize(p,c).audit.get('max_repeat_enforced'))
d=api.frontier(p,c,steps=5).delivered()
print('frontier:', d.audit.get('max_repeat_enforced'), d.certificate.status.value, len(d.violations))"

# §4.2 — validate drops GLOBAL constraints
python -c "
from bt4 import api; from bt4.pipeline.optimize import OptimizeConfig
d='ATG'+'A'*9+'GCTGGAGGACAAGGTGGAGGAGCTGCTGTCC'*2+'TAA'; d=d[:len(d)//3*3]
for lbl,c in [('homopolymer(LOCAL)',OptimizeConfig(max_homopolymer=6)),
              ('max_repeat(GLOBAL)',OptimizeConfig(max_homopolymer=None,max_repeat_length=8))]:
    print(lbl, len(api.validate(d,c).violations))"

# §4.3 — folding reported != optimized (the 156-nt CDS from the default run above)
python -c "
from bt4.biomodels.folding import default
m=default()
d=('ATGGCCAAGCAGCTGGAGGACAAGGTGGAGGAGCTGCTGTCCAAGAACTACCACCTGGAGAACGAGGTG'
   'GCCCGCCTGAAGAAGCTGGTGAAGCTGCTGGAGGACAAGGTGGAGGAGCTGCTGTCCAAGAACTACCAC'
   'CTGGAGAACGAGGTGTAA')
print('optimized (48-nt window):', m.score_sequence(d))
print('reported by the audit    :', m.five_prime_dg(d))
print('honest 5-prime number    :', m.five_prime_dg(d, 48))"

# §4.4 — avoid_internal_start infeasibility
python -c "
from bt4 import api; from bt4.pipeline.optimize import OptimizeConfig
api.optimize('MAAMG', OptimizeConfig(avoid_internal_start=True))"
```

---

## Cross-references

- [`../CLAUDE.md`](../CLAUDE.md) — the constitution: §1 thesis, §4 contracts, §5
  invariants, §6 scientific scope, §10 anti-patterns.
- [`NEXT_SESSION.md`](NEXT_SESSION.md) — live status and the forward queue, which
  this review re-points.
- [`RESEARCH_codon_optimization_SOTA.md`](RESEARCH_codon_optimization_SOTA.md) —
  the field survey; §7's "the CDS is a minority of the expression signal" is the
  scientific ceiling this review argues BT4 should stop accepting as fixed scope.
- [`COMPARISON.md`](COMPARISON.md) — honest positioning vs the vendor tools; its
  TL;DR already concedes the expression point this review measures.
- [`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md) — the
  design-of-record whose "Future / out of scope" section names beyond-CDS design.
