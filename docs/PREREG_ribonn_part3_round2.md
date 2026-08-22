# Pre-registration — RiboNN Part 3, round 2

**Status: DRAFT until merged. Nothing has been run against this design.** Merging this
file freezes it. After that, no threshold, endpoint, protein, UTR, seed or analysis
choice may change, and a failure is recorded as a failure.

## 0. Why there is a round 2

Round 1 (guide Steps 10–13, 2026-08-21, one UTR pair, HEK293T) returned a split verdict
against its own pre-registered gate:

| endpoint | value | threshold | verdict |
|---|---|---|---|
| `within_over_between` | 0.4062 | ≥ 0.20 | pass |
| `median_abs_gc3_spearman` | 0.7385 | ≤ 0.70 | **fail** |

The rule was "continue only if both hold", so **round 1 says stop, and that verdict
stands on its own terms whatever round 2 finds.** Both are reported together, always.

Round 2 is **not** round 1 with a looser bar. The failed threshold is not reused,
adjusted, or appealed. Round 2 exists because three *measurement* limitations leave the
finding under-determined, and all three are fixable with no new laboratory data:

1. **n = 3 proteins.** The failing statistic is a median over three values, so one
   protein (KRas4B, 0.7385) set it.
2. **One UTR pair.** A single 5′/3′ context cannot separate a CDS property from a
   CDS×UTR interaction.
3. **A statistic with no meaningful zero.** `ρ(GC3)` asks "is the response correlated
   with GC3?", which forces an arbitrary cut-off. The decision needs "does the response
   carry stable information *beyond* what BT4 already computes for free?" — a question
   whose null is genuinely zero.

## 1. What this round can and cannot establish

**It uses no measured expression data.** Every number is a RiboNN *prediction*. A change
in that prediction is a statement about the model's output and **not** evidence that
translation efficiency changed. Whether RiboNN's ordering is *correct* is unanswerable
here by construction; only a measured panel can answer it.

**"The score moved" is not a result.** RiboNN is deterministic and every variant is a
different input, so a non-zero response is guaranteed a priori. Round 1 already recorded
`responds_to_synonymous_change: true`, which is true whenever variation exceeds one part
in a billion. Any endpoint whose null is "no movement at all" is therefore vacuous, and
none is used below.

The question that *is* answerable without lab data, and that gates spending:

> For synonymous CDS variants of one protein under a fixed UTR — BT4's deployment
> regime — does RiboNN's score carry information that is **(a)** not already available
> from features BT4 optimizes directly and for free, and **(b)** a stable property of
> the CDS rather than of the UTR it happened to be scored in?

(a) fails ⇒ the head duplicates the trellis. (b) fails ⇒ a gate passed on a panel's UTRs
says nothing about a user's UTRs, so no achievable panel transports. A pass licenses
*acquiring or building* a panel. It is not evidence of skill, promotes nothing, and
`calibrated` stays `False` throughout.

## 2. Design

**Proteins — 16.** Fixed-seed pseudorandom draw (`seed=20260821`) from Ensembl **MANE
Select** human transcripts, protein length ∈ **[150, 1200] aa**, excluding round 1's
genes (BECN1, KRAS, PDE3A). The draw script and the resulting symbols + transcript IDs
are committed **with this file, before any scoring**. Proteins are never re-drawn.

**Variants — 40 per protein (primary).** `api.library(protein, n=40, temperature=1.0,
seed=<fixed per protein>)` — BT4's own synonymous sampler, default `OptimizeConfig` for
*Homo sapiens*.

**Variants — frontier (secondary).** `api.candidates(protein, steps=11, n=24)`, what BT4
actually delivers. Reported, never decisive: frontier points are extremal and
correlated, so this is a realism check, not a powered test.

**UTR contexts — 4.** Real human 5′/3′ pairs from MANE Select: **HBB** and **ACTB**
(in hand), plus two drawn by the same seeded rule, constrained so 5′UTR lengths span
≥ 2× across the set. Every variant is scored in all four.

**Cell type.** Primary `HEK293T` (matching round 1); the all-cell-type mean is recorded
as a secondary read, which costs nothing since RiboNN emits every cell type in one pass.

**Model configuration — bound, not varied.** `species=human`, `top_k=5`,
`batch_size=64`, `num_workers=0`. `top_k` is part of the scope an attestation binds.

**Free-feature set** — what BT4 gets without RiboNN: GC, GC3, CAI, tAI, CpG density,
CDS length.

Scale: 16 × 40 × 4 = **2,560 primary scorings** ≈ 3.5 h at the measured 12–13 seq/min;
secondary ≈ 2 h. Baselines are computed from sequence and cost nothing.

## 3. Endpoints

Every endpoint below is a **comparison against a null or a baseline computed on the same
data**, never an absolute cut-off on a quantity with no natural zero. Uncertainty is a
**cluster bootstrap resampling whole proteins** (10,000 draws), because one protein's
variants are a dependent cluster.

### E1 — Is there information beyond the free features?

Within each protein and UTR context, fit RiboNN's score on the free-feature set in rank
space; keep the **residual fraction** `1 − R²_adj`. Do the same for two baseline scorers
— **pure GC** and **pure CAI** — and for `NullExpressionModel` (constant).

> **E1 passes iff** the bootstrap 95% CI for
> `median(residual_RiboNN) − max(median(residual_GC), median(residual_CAI))`
> **excludes 0 from above.**

No arbitrary margin: the test is whether RiboNN reliably leaves more unexplained than the
free features do. `NullExpressionModel` must come back degenerate; if it does not, the
harness is broken and the run is **void**, not negative.

### E2 — Is that extra information a stable property of the CDS?

Take E1's residual — RiboNN's score with the free features removed. Within each protein,
correlate the residual across UTR contexts (Spearman, all 6 pairs), and take each
protein's median.

> **E2 passes iff** the bootstrap 95% CI for the across-protein median cross-UTR residual
> correlation **excludes 0 from above.**

**This is the endpoint with a genuinely meaningful zero.** If RiboNN's non-free component
is CDS×UTR interaction, residuals computed under different UTRs are uncorrelated and the
median sits at 0. If it is a real CDS property, they correlate. Unlike raw score
movement, "greater than zero" is here a substantive claim rather than an arithmetic
certainty.

### E3 — Does the response survive holding GC fixed?

Within each protein and UTR context, bin variants into GC quintiles; take the ratio of
pooled within-bin score SD to overall within-protein score SD. Null: the same statistic
under 1,000 within-protein permutations of the GC labels.

> **E3 passes iff** the bootstrap 95% CI for the across-protein median retention
> **excludes the permutation null.**

E3 is design-based where E1 is model-based. They are reported together **whether or not
they agree**; disagreement is itself a reportable finding, not something to resolve by
picking one.

### Sanity floor — not an endpoint

`within_over_between ≥ 0.20`, recomputed. Failing it voids the run (the harness stopped
working) rather than making it negative.

## 4. Detection is not sufficiency

With 16 proteins × 40 variants this design is well-powered, so a **very small** effect
will clear zero. Statistical detectability is not practical usefulness, and the two are
kept apart deliberately:

- E1–E3 decide **whether a non-free, stable signal exists at all.**
- The following **magnitudes are reported alongside, and are a human judgement, not an
  automatic gate**: the median residual fraction; the median cross-UTR residual ρ; the
  median GC-stratified retention; and each one's bootstrap CI.

A pass with a residual correlation of, say, 0.05 is a real but tiny signal, and the
honest reading is "detectable, probably not worth a five-figure panel." That call is made
by a person, on the record, with the numbers in front of them — it is not delegated to a
threshold.

## 5. Decision rule

- **All three pass** → proceed to Part 4, *subject to* the §4 magnitude judgement, and
  carrying forward the measured fact that no identified public panel reaches the size
  floor (round 1: PERSIST-seq ~77 usable rows / 4 groups; Mauger ~82 rows / 3 proteins).
  A pass licenses **building or hunting for** a panel — it does not make either rejected
  panel usable.
- **Any fails** → stop. Record the negative result. RiboNN stays `calibrated=False` and
  `NEXT_SESSION.md` item 11 closes with reasons.

## 6. What is not allowed

- Changing any threshold, protein, UTR, seed, or endpoint after any score is seen.
- Re-running with different UTRs, temperatures, cell types, or `top_k` and reporting the
  better outcome. The four UTR contexts are the design, not four attempts.
- Reporting E1/E2/E3 selectively. All three are reported with their values and CIs
  whatever the verdict.
- Treating a pass as evidence of predictive skill, or any score movement as evidence
  about translation efficiency.

## 7. Analysis code is frozen with this document

Scoring and analysis scripts are committed **in this PR, before the first run**, and
re-analysis with different code is a new pre-registration rather than a correction.
Content hashes of the scripts and of the panel manifest (16 proteins, 4 UTR contexts,
seeds) are recorded here at freeze time.
