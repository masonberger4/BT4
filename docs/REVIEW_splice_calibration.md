# Splice calibration — measured results

**Status: both halves of Part B have now been run** — the variant-effect half against a
published benchmark's own scores, and the site-prediction half against BT4's own wrapped
Pangolin on real genomic sequence. Nothing here changes any backend's `calibrated` flag,
and nothing here is a calibration claim. This is a record of what was measured, by whom,
on which bytes.

Runbook: [`DESIGN_splice_cnn_calibration.md`](DESIGN_splice_cnn_calibration.md) Part B.
Run on a maintainer's Windows machine, 2026-08-19, with BT4's own
`bt4 variant-gate`. The scores are `splicebench2023`'s **own pre-computed** SpliceAI and
Pangolin columns, so no model weights were involved: this measures BT4's gate against a
published benchmark, not BT4's adapters against a model.

## Panels

| file | variants | genes | held out | content hash |
|---|---|---|---|---|
| `variants.tsv` | 3,616 | BRCA1, FAS, MST1R, POU1F1, WT1 | **no** | `da59449b27dec3a4…` |
| `all6.tsv` | 3,912 | + MLH1 | **no** | `fcd1b87094922943…` |
| `heldout.tsv` | 1,539 | MST1R, POU1F1 | **yes** | `3718d96361f32acb…` |

Source: Smith & Kitzman, *Genome Biol* 24:294 (2023), Zenodo record 8351879
(`splicebench_data.tar.gz`, md5 `e628ca38209064be73d28d5bddf1ae80`, verified on
download). MIT licensed. Negative construction: assayed variants the assay called
non-disruptive. Label `sdv_fc2`; stratifier `exon`.

**Over half of this benchmark is not held out.** BRCA1 (chr17), FAS (chr10) and WT1
(chr11) are on chromosomes both SpliceAI and Pangolin trained on — 2,077 of 3,616
variants. Only the chr3 genes are held out.

## Results

### Not held out — all five genes (`variants.tsv`)

| score | exonic AP | exonic skill | intronic AP | intronic skill |
|---|---|---|---|---|
| `spliceai_masked` | 0.535 | 0.350 | 0.777 | 0.715 |
| `pangolin_masked` | 0.545 | 0.365 | 0.784 | 0.724 |
| `spliceai_unmasked` | 0.523 | 0.334 | 0.763 | 0.698 |
| `pangolin_unmasked` | 0.609 | 0.454 | 0.796 | 0.739 |

Pangolin beats SpliceAI in every comparison, matching Zeng & Li's published ordering.

### Held out — chr3 only (`heldout.tsv`), `pangolin_masked`

| stratum | n | prevalence | AP | skill | ROC | ECE |
|---|---|---|---|---|---|---|
| exonic | 796 | 0.423 | 0.665 | **0.419** | 0.668 | **0.345** |
| intronic | 743 | 0.227 | 0.598 | **0.480** | 0.770 | **0.181** |

## What the numbers say

**1. On the held-out subset, intronic skill is much lower and exonic is not.** Same
backend, same score column:

| | exonic skill | intronic skill |
|---|---|---|
| all five genes | 0.365 | 0.724 |
| held out only | **0.419** | **0.480** |

Intronic skill falls by a third (−33.7%); exonic *rises* slightly (+0.054). Whatever the
full-panel intronic figure reflects, a third of it does not survive restricting to the
chromosomes these models did **not** train on — while the exonic figure, already the
weaker one, does not fall at all.

> **This comparison is confounded, and the confound cannot be removed from this panel.**
> Holding out chr1/3/5/7/9 here does not remove training genes from a fixed set — it
> *replaces* the gene set. The five-gene row is BRCA1 + FAS + MST1R + POU1F1 + WT1; the
> held-out row is MST1R + POU1F1 alone. So the drop is training-chromosome overlap
> **and** gene difficulty, inseparably, and no split of this benchmark can separate them:
> the two are perfectly collinear because each gene sits on exactly one chromosome. The
> direction is worth recording, and it is the direction leakage would predict — but
> "a third of intronic skill is training leakage" is **not** what these numbers show, and
> must not be written that way. Separating the two needs a panel carrying held-out and
> non-held-out variants in the *same* genes, which this benchmark does not contain.

**2. The exonic/intronic gap nearly closes on held-out data.** On skill the gap is 0.359
across all five genes and **0.061** on chr3 alone. The published pair is computed over a
panel that likewise includes training-chromosome genes, so if the gap there is inflated
by the same mechanism, it is inflated too — but per the confound above, *this panel
cannot tell you the mechanism*. What it does show is that a widely-cited gap between two
strata is not stable across subsets of the very benchmark it came from, which is reason
enough not to design an operating point around its size.

**3. Calibration is poor, and worst exactly where BT4 operates.** Expected calibration
error, held-out exonic: **0.345**. That is the average gap between a score and the
observed frequency of disruption at that score.

> **Read this carefully.** `DS_maxm` and `pang_max_abs` are **delta** scores — differences
> between probabilities — not probabilities of disruption. An ECE computed on them asks
> "does a delta of 0.3 mean 30% of such variants are disruptive?", which is a fair
> question but not the same as "is the model's output a calibrated probability". What it
> does establish is that **these numbers cannot be read as probabilities of anything**,
> which matters because BT4's own operating point (`DEFAULT_SITE_PROBABILITY = 0.5`)
> treats scores as if they were. It is direct evidence against importing published
> cutoffs, and consistent with OpenSpliceAI's finding that SpliceAI-architecture models
> are overconfident.

**4. Two comparisons that turned out not to be comparisons.** The published
0.419 / 0.773 pair is a **median across tools** pooled over all six datasets. Matching
the composition with `--include-mlh1` moved both figures *further* from it (exonic
0.535 → 0.575, intronic 0.777 → 0.841), because MLH1's variants are easier than average.
Neither the levels nor the gap are comparable between a single tool and a median over
tools; only the ordering is.

And on the held-out panel raw **AP inverted** (exonic 0.665 above intronic 0.598) purely
because the exonic stratum has nearly double the prevalence — 42.3% against 22.7% — and
average precision's floor *is* the prevalence. On skill and on ROC-AUC the expected
ordering held. The panel was correct; three successive versions of BT4's own reporting
text were not (fixed in #103, #104, #105).

## What this does NOT establish

- **Nothing about BT4's adapters.** These are the benchmark's own pre-computed scores.
  BT4's wrapped Pangolin has passed an *integration-fidelity* gate (bit-for-bit
  reproduction, 18 cases, deviation 0.0); that is a different question and remains the
  only thing its `calibrated` flag has ever meant.
- **Nothing about BT4's regime.** Every variant here is a natural single-nucleotide
  change in a natural gene. BT4 designs **synonymous variants of a coding sequence in a
  vector** — a regime no panel in this benchmark covers. This is the same gap CLAUDE.md
  documents for RiboNN, for the same reason.
- **No promotion.** Every run reported `PROMOTABLE` unreachable or `False`, and no
  threshold was declared. `bt4.biomodels.splice.default()` still returns the PWM baseline.

---

# Site prediction — measured (2026-08-19)

Run on the same maintainer machine, with **BT4's own wrapped Pangolin** against the
hash-verified GPL weights. Unlike the variant half, this exercises BT4's adapter end to
end: the one-hot encoding, the 12-model ensemble, the anchor convention, and the gate.

Panel built by `scripts/make_gencode_splice_panel.py` from **GENCODE v44 / GRCh38
primary assembly**, both downloaded from the same release directory and md5-verified
(`7450ef42cf9cb3d29625320b22d4bb45`, `9c3fc2ca260a767530dddb0f26721a6b`). MANE Select
transcripts only, ±5,000 nt of real flank, antisense-overlapping windows skipped.

| panel | windows | positions | sites | prevalence | content hash |
|---|---|---|---|---|---|
| `trial.tsv` | 5 | 372,634 | 129 | 0.000346 | `cb9c8f519e279268…` |
| `panel20.tsv` | 20 | 861,096 | 333 | 0.000387 | `d30ba2cc7dbe1d13…` |

Both are **fully held out**: groups are chr1/3/5/7, none of which either model trained
on. Motif consistency **100%** on both — every annotated site carries its canonical
GT/AG, and none of these genes has a minor-spliceosome intron.

## Results — `panel20.tsv`, `--cnn-anchors`

### Pangolin (one combined track)

| | AP | skill | ROC | top-k | ECE |
|---|---|---|---|---|---|
| **Pangolin** | **0.983** | **0.983** | 1.000 | **0.940** | 0.050 |
| `pwm` (best baseline) | 0.096 | 0.096 | 0.994 | 0.159 | 0.080 |
| `gt_ag` | 0.003 | 0.003 | 0.940 | 0.003 | 0.120 |
| `permutation` | 0.000 | 0.000 | 0.509 | 0.003 | 0.050 |
| `constant` | 0.000 | 0.000 | 0.500 | 0.000 | **0.000** |

Bar declared **before** the run: `--min-pr-auc-skill 0.75`. All four conditions held, so
the run reports `PROMOTABLE on this panel: True` — **the first time any BT4 splice
backend has reached that state.** It is not a promotion, and §"What this does not
establish" below is the reason.

### SpliceAI, and the comparison that needed making comparable

SpliceAI emits separate donor and acceptor tracks, so its default run scores **two**
strata against Pangolin's one — a different, harder task, because in the donor stratum
an acceptor site is a *negative*. Putting those figures side by side is the comparison
two consecutive gate runs invite and do not support; `--combined-track on` (added for
this) scores it on the shared task.

| run | stratum | skill | top-k | ECE |
|---|---|---|---|---|
| SpliceAI, default | acceptor | 0.969 | 0.910 | 0.000 |
| SpliceAI, default | donor | 0.961 | 0.904 | 0.000 |
| **SpliceAI, `--combined-track on`** | splice | **0.965** | **0.907** | 0.000 |
| **Pangolin** | splice | **0.983** | **0.940** | 0.050 |

**Separating cost SpliceAI nothing.** The combined figures are the exact mean of the
separated pair to three decimals — skill `(0.969+0.961)/2 = 0.965`, top-k
`(0.910+0.904)/2 = 0.907` — so its kind discrimination is effectively perfect: its donor
track is already near zero at acceptor sites. The comparability caveat is still right in
principle; its measured magnitude here is nil, and an earlier draft of this document
implied otherwise.

**The between-model gap reproduces the published one where the metric has room.**

| | observed | published (Zeng & Li 2022) |
|---|---|---|
| top-k gap | **0.033** (0.940 − 0.907) | 0.040 (79% − 75%) |
| AP gap | 0.018 (0.983 − 0.965) | 0.080 (0.85 − 0.77) |

Both absolute levels sit far above published, as §"What this does not establish" records
— but the *ordering* and, on top-k, the *magnitude* are faithful. The AP gap is
compressed because both models are near ceiling there and average precision has nowhere
left to go, while top-k at 0.94/0.91 still does.

An internal check that had to hold and did: `pwm` scores **0.096 / 0.159 / 0.080** and
`gt_ag` **0.003** in *both* combined runs, identical, because the baselines are
sequence-derived and backend-independent. `permutation` differs (ROC 0.482 vs 0.509),
correctly — it shuffles the *head's* scores, so a different head gives a different null.

### Two-backend agreement — do they point at the same bases?

`bt4 splice-agreement`, both CNNs over the same panel. Neither gate report above can
answer this: two backends can each score ~0.97 while being confident about different
positions.

| | |
|---|---|
| Jaccard of called positions | **0.855** (307 shared of 359 union) |
| Spearman over called positions | 0.820 |

| annotated site recovered by | count | of 333 |
|---|---|---|
| both backends | 300 | 90.1% |
| only Pangolin | 15 | 4.5% |
| only SpliceAI | 6 | 1.8% |
| **neither** | **12** | **3.6%** |

Recall implied by the 2×2 — Pangolin 315/333 = **94.6%**, SpliceAI 306/333 = **91.9%** —
tracks the gate's top-k (0.940, 0.907) as it should, the small differences being the two
metrics' different denominators.

**Running both is not redundant.** On **21 sites (6.3%)** exactly one model finds the
site. Those are precisely the positions an audit should surface as uncertain, and no
single-backend run can identify them.

**But agreement is not correctness, and the 12 is the number that shows it.** Two
independently-trained models miss the *same* 12 sites. Agreement on a miss is a
**correlated blind spot**, not reassurance — these architectures are similar, their
training corpora overlap, and both learned from the same style of annotation. That is
the standing limit on reading cross-backend agreement as an uncertainty signal: it
bounds *independent* error, not shared error.

**Seven positions both models call are not annotated sites** (307 shared calls, 300 of
them real). On a MANE-Select-only panel those are as likely to be genuine sites of
non-MANE isoforms — which this panel scores as negatives — as they are to be shared false
positives. Worth checking before either reading is adopted; it is a property of the panel
construction, not a measured model error.

## What the numbers say

**1. The per-kind anchors are confirmed on real data.** Donors peak at **−1** for 100% of
sites and acceptors at **+1** for 99%. Those offsets were derived in #102 from SpliceAI's
training-label construction and Pangolin's CLI — from reading source, not measuring — and
this is the first check against real genomic sequence with real weights. It was the
largest correctness risk in this half: a wrong scalar anchor previously drove the PWM
control from 0.853 to 0.0001.

**2. The result is stable across a 2.3× larger panel.** Between the two panels head skill
moved −0.005 (0.988 → 0.983), top-k −0.006, and the `pwm` baseline −0.011. Five windows
were not a lucky draw.

**3. Pangolin beats every baseline by roughly 10×** on skill (0.983 against `pwm` 0.096).
That margin is the point of the permanent baselines: BT4 ships the PWM for free and with
no licence, so a wrapped CNN that could not clear it would not have earned a PyTorch
dependency, a hash-pinned weight set, and a GPL term.

**4. `gt_ag` is the clean illustration of why this gate leads on skill.** The canonical
dinucleotide scores **ROC 0.940** and **AP 0.003**: it ranks acceptably and predicts
appallingly, because at 1-in-2,600 prevalence there are thousands of GT/AGs that are not
splice sites. ROC is near-saturated for everything here (`pwm` 0.994, head 1.000) and
carries almost no information.

**5. The ECE column is not evidence, and the run now says so.** `constant` and
`permutation` both **match or beat** Pangolin's ECE of 0.050 — a base-rate predictor and
a shuffled null are better calibrated than the model. That is what ECE measures at this
prevalence, and it is why an ECE ceiling was removed from what counts as a declared bar
(#109). The note naming the offending baselines is generated from the run's own numbers.

## What this does NOT establish

- **Not a promotion, and the panel is easier than the published benchmark.** Zeng & Li
  report Pangolin at AUPRC **0.85** and top-1 **79%**; this run reads 0.983 and 0.940,
  **+0.133 and +0.150 above published**, consistently across both panel sizes. A
  systematic gap in the flattering direction is a statement about the panel, not the
  model. Theirs is genome-wide across the test chromosomes with a far larger and harder
  negative pool; this is 20 MANE gene bodies at 100% canonical motifs. The honest reading
  is *"Pangolin locates splice sites within MANE gene bodies"*, which is narrower than
  what the published figure measures.
- **No exonic/intronic split is available here.** Pangolin emits one combined P(splice)
  track, so the panel scores a single `splice` stratum. The 0.419-exonic / 0.773-intronic
  penalty that matters most to BT4 is a *variant-effect* finding and is not checkable on
  this panel shape — so this result must not be read as relieving it.
- **Nothing about BT4's regime.** Every site here is a natural splice site in a natural
  gene. BT4 designs **synonymous variants of a coding sequence in a vector**, and a model
  that finds real splice sites has not thereby been shown to correctly stay silent on a
  designed CDS, or to flag one that creates a cryptic site. This is the same gap CLAUDE.md
  documents for RiboNN, for the same reason, and it is the one this benchmark cannot close.
- **`calibrated` is unchanged.** Pangolin's flag still reports *integration fidelity*
  only, `default()` still returns the PWM baseline, and promotion remains behind
  `BT4_SPLICE_USE_ATTESTED`. Both gates now have a passing result on their own terms;
  whether that warrants changing the default is a deliberate human decision, not a
  consequence of this run.

---

# Designed synonymous CDS — measured (2026-08-19), and the defect it exposed

The panel that reaches BT4's own regime: `scripts/make_designed_cds_panel.py` over
Ranaghan et al. 2021 Table 4 (CC BY 4.0) — three proteins, each the native human CDS
plus 30 designs from three anonymized commercial optimizers. 93 members, content hash
`fa0df04f46e53dd9…`. It carries **no splice labels and cannot**: designed coding
sequence has no splice ground truth. Everything below is label-free.

## The first run was reported wrong, and the error was BT4's, not the model's

The probe's headline field is the **Δsplicing spread** across a group's designs. On the
first run, with the hash-verified Pangolin weights:

| protein | Pangolin Δ spread | reported as |
|---|---|---|
| Beclin1 | `0.0000` | "flat across all 30 designs" |
| KRas4B | `0.0000` | "cannot rank these candidates at all" |
| PDE3A | `1.0885` | the one responsive group |

That reading was wrong. Running Pangolin directly on the same sequences showed its raw
per-position scores are **not** flat and **not** zero:

| sequence | peak score | positions > 0.5 | pooled risk @ 0.5 | pooled @ 0.01 |
|---|---|---|---|---|
| KRas4B native | 0.128 | **0** | 0.0000 | 7.89 |
| KRas4B design0 | 0.276 | **0** | 0.0000 | 9.38 |
| Beclin1 design0 | 0.435 | **0** | 0.0000 | 10.96 |
| PDE3A design1 | 0.445 | **0** | 0.0000 | 9.97 |

Every position was nonzero and the native differed from its designs by more than
twofold. **No position on any sequence reached 0.5** — and `pool_log_odds` sums
`max(0, logit(p) − logit(background))` with `background = DEFAULT_SITE_PROBABILITY =
0.5`, so the hinge floored every score to zero before the Δ was taken. The zeros were a
property of BT4's pooling, reported as a property of the CNN.

The same hinge explains the rest of that run: the `+0.000` rank agreements were Spearman
correlations of constants, and the `0.000` sign agreement on KRas4B was the sign of zero.

**Why it was invisible.** `DEFAULT_SITE_PROBABILITY` is documented as *"a display /
localization knob, not a calibrated cutoff"* — but it was wired in as a hard gate inside
risk pooling, where instead of shifting a display it silently zeroed the output. A pooled
risk of `0.0` meant either "no risk" or "nothing cleared an admittedly-uncalibrated
cutoff", with nothing distinguishing them. In BT4's regime the second is the universal
case.

Note this is *not* suppression of a long tail: with `top_k = 3` only three positions
ever contribute, so top-k already excludes the tail. The hinge's only effect on those
three is to floor them — buying non-negativity at the price of blindness below 0.5.

## The fix

**Lowering the background was rejected.** It is the same uncalibrated knob pointed
somewhere more flattering, and the constant's own docstring says the number is a
convention. Deriving a real operating point is Part B's job, on data.

What landed instead:

- `pool_top_k_logit` — the same top-k log-odds with the hinge **and the background
  removed**. It takes no background parameter at all, so it introduces no operating
  point; it is monotone in the model's scores everywhere, so it still separates
  sequences the risk has flattened. It is **not a risk**: it goes negative, and it has
  no calibrated zero. Where the hinge does not bind it equals the pooled risk exactly
  (`logit(0.5) == 0`), so it adds a statistic without moving one.
- `PooledRisk` / `pooled_risk_detail` — the same number plus what makes it attributable:
  `n_above_background`, `max_score`, and `below_background`, which is true exactly when
  the risk is zero *by construction*.
- Every consumer that reports a risk now reports which zero it is: `bt4 designed-probe`,
  `bt4 validate --splice-backend`, the Studio ASSP banner, `BackendCandidateAudit`, and
  `AgreementReport.degenerate`.
- The audit and the cross-check now **pool against the threshold they localize at**. They
  passed the caller's `threshold` to localization and let pooling keep the default 0.5,
  so `--threshold 0.2` flagged sites at 0.35 and pooled them as zero. At the default
  threshold nothing changes — the two constants are the same number.

No shipped number moved: `pool_log_odds` and `pooled_risk` are byte-identical, pinned by
`test_default_background_reproduces_the_legacy_expression`.

## Results, re-measured with the fix

`bt4 designed-probe`, Pangolin (4 tissues) and the PWM baseline, 30 designs vs native:

| protein | backend | Δ spread (risk) | response spread | rank agree (risk) | rank agree (response) |
|---|---|---|---|---|---|
| Beclin1 | Pangolin | `0` **floored**, peak 0.435 | 4.2301 | `+0.000` | `+0.614` |
| Beclin1 | PWM | 7.4284 | 7.4284 | — | — |
| KRas4B | Pangolin | `0` **floored**, peak 0.323 | 3.8850 | `+0.000` | `+0.195` |
| KRas4B | PWM | 5.5951 | 5.5951 | — | — |
| PDE3A | Pangolin | 1.0885 | 5.9209 | `+0.381` | `+0.162` |
| PDE3A | PWM | 5.3478 | 5.3478 | — | — |

Sign agreement, risk vs response: Beclin1 `0.000` / `0.400`, KRas4B `0.000` / `0.600`,
PDE3A `0.100` / `0.667`.

**What the numbers say.**

- **Pangolin does respond to synonymous change.** Response spreads of 3.9–5.9 log-odds
  across designs of the same protein, on the axis BT4 varies and nothing else. The
  earlier "cannot rank these candidates" conclusion is withdrawn.
- **BT4's shipped Δsplicing is mute in BT4's regime.** Two of three groups floored
  entirely; the third (PDE3A, spread 1.0885) only because a few designs happened to
  clear 0.5. Routing `delta_splicing` into candidate selection today would contribute
  nothing on most proteins — not because the model is silent, but because the pooling is.
- **The PWM baseline's risk and response are identical to four decimals in all three
  groups.** That is the arithmetic working as expected, not a coincidence: its top-3
  positions always clear 0.5, so the hinge never binds and the two poolings coincide. It
  is also why this defect never showed up against the default backend.
- **Cross-backend agreement is low where it is now measurable** — `+0.614`, `+0.195`,
  `+0.162` response rank agreement between Pangolin and the PWM baseline. Two backends
  that would often pick different candidates. That is an uncertainty signal, and it
  argues against routing either into selection before Part B.

## What this does NOT establish

- **No labels, so no accuracy claim of any kind.** Nothing here was assayed and none of
  it is annotated. "Responds to synonymous change" is not "responds correctly".
- **The response statistic is not calibrated and is not a risk.** It is a ranking
  quantity with no meaningful zero. It must never be quoted as how spliceogenic a
  sequence is, and it does not move any `calibrated` flag.
- **Low agreement does not say which backend is right**, and high agreement would not
  either — both pad with 5,000 literal `N`, so a shared artifact is invisible to it.
- **This is not the specificity panel.** Whether a model stays correctly silent on a
  clean designed CDS and correctly flags one carrying a cryptic site still needs labels.

## Still to run

- **A specificity panel in BT4's regime**: designed synonymous CDS with *known* splice
  outcomes — does a model stay correctly silent on a clean designed CDS, and correctly
  flag one that creates a cryptic site? Every labeled measurement above is recall on
  natural sites in natural genes, and the designed-CDS panel has no labels at all.
- **Derive BT4's own operating point** (Part B). The measurements above are the concrete
  argument for it: the shipped 0.5 makes the whole splice objective inert on designed
  coding sequence, and no defensible replacement can be picked without labeled data.
- **Resolve the 7 shared calls that are not annotated sites** — non-MANE isoform sites
  the panel scores as negatives, or shared false positives. The two readings have
  opposite implications and the panel as built cannot separate them.
- Deciding whether the attested-promotion opt-in should become the default, and a Studio
  checkbox for it.
