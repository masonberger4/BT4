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

**Neither half can be gated on its own**, and §3.2 records the known-answer regimes that
prove it: (a) alone passes per-context noise, which is exactly what no fixed-UTR panel
could calibrate; (b) alone passes a pure GC/CAI blend, because once the fit removes the
free features such a model leaves deterministic dust that is identical in every context
and so correlates at 1.0. The gate is therefore their **product**, measured against a
family of free-feature blends run through the same pipeline.

A pass is a **necessary condition**, not permission to spend. It is not evidence of
skill, promotes nothing, and `calibrated` stays `False` throughout.

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

**Cell type — `HEK293T` only, matching round 1.** An earlier draft of this document also
promised the all-cell-type mean as a free secondary read, on the reasoning that RiboNN
emits every cell type in one pass and `cell_types` merely selects columns. That is true
of the model and **false of BT4's adapter**: `score_many` reduces internally with
`self.cell_types`, and no public method returns two reductions from one pass, so a second
readout means a second forward pass — 7 h rather than 3.5 — or reaching into the private
`_predict_te`, which is the cross-layer private-symbol pattern CLAUDE.md §10.9 forbids.
Cell-type robustness is therefore **out of scope for this round** and would need its own
pre-registration. Recorded here rather than quietly dropped, because the claim was wrong
and the correction is the kind this document exists to force.

**Model configuration — bound, not varied.** `species=human`, `top_k=5`,
`batch_size=64`, `num_workers=0`. `top_k` is part of the scope an attestation binds.

**Free-feature set** — what BT4 gets without RiboNN: GC, GC3, CAI, tAI, CpG density,
CDS length.

Scale: 16 × 40 × 4 = **2,560 primary scorings** ≈ 3.5 h at the measured 12–13 seq/min;
secondary ≈ 2 h. Baselines are computed from sequence and cost nothing.

## 3. The gate

The gate is a comparison against a null **computed on the same data**, never an absolute
cut-off on a quantity with no natural zero. Uncertainty throughout is a **cluster
bootstrap resampling whole proteins** (10,000 draws), because one protein's variants are
a dependent cluster.

### 3.1 The free-feature fit

Within each protein and UTR context, RiboNN's score is regressed on the free features
(GC, GC3, CAI, tAI, CpG density, CDS length) **rank-transformed, plus their squares**.
What is left is "the residual".

**Void condition.** If the across-protein median |Spearman(residual, GC)| exceeds
**0.10**, the fit did not remove GC, and the run is **void** -- neither pass nor fail.

### 3.2 The gate: stable non-free signal, against a free-feature null

For each protein, combine **size** and **stability** of the residual:

    stable_non_free = residual_fraction x max(0, median cross-UTR residual rho)

and compare it against the same statistic computed for a family of **free-feature
baseline scorers** pushed through the identical pipeline: pure GC, pure CAI, and 24
seeded random linear blends of all six features. The floor is the **worst case** over
that family.

> **PASS iff** the cluster-bootstrap 95% CI for the across-protein median of
> `stable_non_free(RiboNN) − max over baselines` **excludes 0 from above.**

Both halves are load-bearing, and **neither alone works** -- established by running
known-answer regimes through the analysis (`scripts/prereg_round2_selftest.py`), not by
reasoning about the definitions:

- **Stability alone passes a pure GC/CAI blend.** After the fit removes the free
  features, such a model leaves only numerical dust -- but the dust is deterministic, so
  it is *identical* in every UTR context and correlates at exactly 1.0.
- **Size alone passes per-context noise**, which is precisely the thing that cannot be
  calibrated.
- **Size x stability against *single-feature* baselines still passes a blend.** The fit
  is in rank space, and a monotone transform of a blend is not linear in the ranked
  features, so a blend leaves a small stable residual that a `score = GC` baseline does
  not. The null has to be *any* blend, which is why the baseline family exists.

The residual limitation, stated rather than hidden: this design cannot separate
"information beyond the free features" from "structure the fit failed to capture". The
blend family bounds that misspecification empirically; it does not eliminate it.

### 3.3 Reported alongside, never gating

`residual_fraction`; median cross-UTR residual rho; GC-stratified retention (0 = pure GC
detector, 1 = GC-independent -- an earlier draft had this direction backwards and gated
on it); |rho(residual, GC)| as the fit diagnostic.

### Sanity floor -- not an endpoint

`within_over_between >= 0.20`, recomputed on the new panel. Failing it voids the run --
the harness stopped working -- rather than making it negative.

## 4. What this round can rule out, and what it cannot

**It can only establish a necessary condition.** Without measured expression there is no
way to show RiboNN is *useful*; there is only a way to show it cannot be. So:

- **Gate fails** -> decisive. The non-free component does not survive a change of UTR, so
  it is not a CDS property, and no panel built on one set of UTRs would transport to a
  user's. Stop.
- **Gate passes** -> a necessary condition holds. **This is not permission to spend.** At
  16 proteins x 40 variants the design is well powered, so a *very small* effect clears
  zero, and statistical detectability is not practical usefulness.

The spending decision is therefore an explicit **human judgement on magnitudes, recorded
on the record**, not an automatic consequence of the gate. Reported with cluster-bootstrap
CIs:

| magnitude | reads as |
|---|---|
| median cross-UTR residual rho | how much of the non-free component is a stable CDS property |
| median residual fraction `1 - R2_adj` | how much of RiboNN's within-protein signal the free features miss |
| median GC-stratified retention | 0 = pure GC detector, 1 = GC-independent |
| median \|Spearman(residual, GC)\| | fit diagnostic; must be <= 0.10 or the run is void |

A pass at a residual correlation of, say, 0.05 is real and almost certainly not worth a
five-figure panel. That call is made by a person looking at these numbers.

## 5. Decision rule

- **Gate fails** -> stop. Record the negative result. RiboNN stays `calibrated=False` and
  `NEXT_SESSION.md` item 11 closes with reasons.
- **Gate passes** -> record the magnitudes and make the spending judgement explicitly,
  carrying forward the measured fact that no identified public panel reaches the size
  floor (round 1: PERSIST-seq ~77 usable rows / 4 groups; Mauger ~82 rows / 3 proteins).
  A pass licenses **building or hunting for** a panel; it does not make either rejected
  panel usable, and it is not evidence of predictive skill.

## 6. What is not allowed

- Changing any threshold, protein, UTR, seed, or endpoint after any score is seen.
- Re-running with different UTRs, temperatures, cell types, or `top_k` and reporting the
  better outcome. The four UTR contexts are the design, not four attempts.
- Reporting the gate or the §4 magnitudes selectively. Every one is reported with its
  value and CI whatever the verdict, including the ones demoted in §3.3.
- Treating a pass as evidence of predictive skill, or any score movement as evidence
  about translation efficiency.

## 7. Analysis code is frozen with this document

Every script and the drawn panel are committed **in this PR, before the first scoring
run**. Re-analysing with different code is a new pre-registration, not a correction.
SHA-256 at freeze:

| file | sha256 |
|---|---|
| `scripts/prereg_round2_draw_panel.py` | `18a46434c99c906fc9dacbbb60e23741c87e6e62d3a4c362c57994e48a694181` |
| `scripts/prereg_round2_score.py` | `56655ed6c8a6442f91ff0e0454016aa071a1d08759c3731974c404e218aa5a10` |
| `scripts/prereg_round2_analyze.py` | `d7b4d7929f650ea80e32c367a11089aa9d9c96c85772d8577e0ae27168480b78` |
| `scripts/prereg_round2_selftest.py` | `508d428eb9eb0d3bb5c300d10beb519d4d8383811d548f12148082f402af6009` |
| `scripts/data/prereg_round2_panel.json` | `77a2d3721346515fe3dc4c51d81678df7400e23a2fb9dcbe326aca45ea48c0e5` |

The panel additionally carries its own `content_hash`
(`b7da49ea1bd241994cd875b13101215259a7c81f055444461801e8f6916a5359`) over the drawn
proteins, UTR contexts and seeds — deliberately *excluding* the incidental skip tallies,
so `--verify` proves the **draw** reproduces rather than that the file is byte-identical.
The source it was drawn from is pinned too: NCBI MANE release 1.5,
`MANE.GRCh38.v1.5.summary.txt.gz`, sha256
`d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c`, 19,363 MANE Select
rows.

**The panel as drawn, for the record** (16 proteins, 156–765 aa): ZNF286A, ACADL, SNRPA,
B3GNT9, GGA2, WRNIP1, IL1RN, TMEM150A, GKN1, MRGPRE, DBX1, IER3, C17orf99, PRSS57, PAK5,
NOX5. UTR contexts: HBB (5′ 50 / 3′ 134), ACTB (84 / 600), DEFB118 (50 / 753), LY6E
(115 / 620).

**A known weakness, recorded rather than re-drawn.** The four 5′UTRs are 50, 84, 50 and
115 nt — a 2.3× span, which satisfies the pre-registered ≥ 2× rule but clusters three of
the four between 50 and 84. The 3′UTRs spread more (134–753). Re-drawing to obtain a
prettier spread would be selecting the design for its appearance, so the draw stands; the
consequence is that the cross-context half of the gate is tested over a narrower range of
5′ contexts than the rule's spirit intends, and a null result should be read with that in
mind.
