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
twofold. None of these four reached 0.5 — and `pool_log_odds` sums
`max(0, logit(p) − logit(background))` with `background = DEFAULT_SITE_PROBABILITY =
0.5`, so the hinge floored every score to zero before the Δ was taken. The zeros were a
property of BT4's pooling, reported as a property of the CNN.

> **Correction (2026-08-20).** The sentence that stood here — *"No position on any
> sequence reached 0.5"* — was **wrong**, and it propagated into `CLAUDE.md` §6,
> `NEXT_SESSION.md`, the `CHANGELOG`, and three docstrings in `base.py`. It was
> generalized from the four sequences in the table above to all 93, and this document
> already contradicted it two sections below: PDE3A's Δ spread of `1.0885` is
> arithmetically impossible under a hinge at 0.5 unless something cleared 0.5. The
> full-panel count is in the next section. The finding it was offered as evidence for
> — that the hinge discards the CNN's signal — survives; the universal quantifier did
> not, and generalizing a spot check is exactly the §5 failure this repo is built to
> catch.

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

## The full-panel count, and the defect it exposes in the *other* direction

Every sequence in the panel, both backends, counting positions above the 0.5 pooling
background:

| backend | Beclin1 | KRas4B | PDE3A | total | peak range |
|---|---|---|---|---|---|
| **Pangolin** (opt-in CNN) | 0/31 | 0/31 | **6/31** | **6/93** | 0.088 – 0.748 |
| **PWM baseline** (`default()`) | 31/31 | 31/31 | 31/31 | **93/93** | 0.981 – 1.000 |

Two findings, and the second was not previously recorded anywhere.

**Pangolin is floored, but not uniformly.** All six clearing sequences are *designs*,
all of one protein; **no native CDS clears 0.5 in any group**. So the hinge silences the
CNN on two proteins out of three and passes a sparse, erratic signal on the third —
which is why PDE3A alone showed a nonzero Δ spread (1.0885) and why that spread should
never have been read as "the responsive group". Six sequences with one or two clearing
positions each is not a measurement of anything; it is the tail of a distribution
poking over an arbitrary line. *(With no labels, the tempting reading — designs
introducing sites their native lacks — is exactly what this panel cannot support.)*

**The shipped default has the inverse defect, and it is arguably worse.** The PWM
baseline clears 0.5 on **every sequence**, native included, at peaks of 0.981–1.000. Its
top-`k` contributions are therefore never hinged, its risk and response coincide to four
decimals, and it flags a site on **100% of designed coding sequences**. A detector that
fires on everything is uninformative in precisely the way a detector that fires on
nothing is — and this one is the backend `bt4.biomodels.splice.default()` returns, so it
is the path most users are on, and it drives Studio's "distinct splice sites" column.

**One constant is standing in for two different score scales.** 0.5 is simultaneously
too high for the CNNs and too low for the PWM baseline. No single value fixes both,
which is the concrete reason deriving an operating point needs data per backend rather
than one better guess.

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

## The N-padding is not neutral — measured (2026-08-20)

The open question after the correction above was whether Pangolin's low band on designed
CDS is the model's view of the sequence or an artifact of the **5,000 literal `N`** both
adapters pad with (`pangolin.py`, `spliceai.py` — upstream's own documented convention).
A low structureless band is what a *broken input* produces too, and nothing had
separated the two.

Four arms, same designed CDSs, differing **only** in what surrounds them. Scores are
sliced back to the CDS by `score_in_context`, so every number below describes coding
positions:

| arm | flank | median peak in CDS | range |
|---|---|---|---|
| **A** | `N` × 5,000 (shipped default) | 0.2757 | 0.113 – 0.445 |
| **B** | random uniform ACGT | **0.1731** | 0.059 – 0.391 |
| **C** | real human chr1, site-free | **0.3691** | 0.125 – 0.749 |
| **D** | shuffled chr1 (same composition) | 0.4944 | 0.204 – 0.759 |

*(9 sequences: native + 2 designs per protein. Real flank from the 160 kb chr1 window in
`big_panel.tsv`, from a stretch carrying no annotated site.)*

**Replication**, because the surprise (D > C) rested on one shuffle draw — 3 designs ×
{3 *different* real regions, 3 *independent* shuffles of one region}:

| protein | real #1 / #2 / #3 | shuffled #1 / #2 / #3 |
|---|---|---|
| Beclin1 | 0.4622 / 0.4622 / **0.4638** | 0.5652 / 0.4913 / 0.5070 |
| KRas4B | 0.3691 / 0.3691 / **0.3644** | 0.5155 / 0.4349 / 0.4718 |
| PDE3A | 0.5840 / 0.5840 / **0.5855** | 0.6347 / 0.5820 / 0.6821 |

### What this establishes

- **`N`-padding systematically deflates the model's CDS scores.** On the **same
  9-sequence set**, median peak **0.2757 padded vs 0.3691** in real genomic context
  (+0.093) — and under real flanks several designed CDSs reach or cross 0.5 that never
  did padded (arm C's range extends to 0.749, arm A's stops at 0.445).

  > **Correction (2026-08-20).** This bullet, and the copies in `CLAUDE.md` §6 and
  > `NEXT_SESSION.md`, previously read **"0.276 → 0.462"**. That paired arm A's median
  > over the **9-sequence** main set with **0.4622**, which is *one protein's* value from
  > the **3-sequence replication** below — a set mismatch, and in the flattering
  > direction, roughly doubling the stated effect. The same-set figure is **+0.093**. The
  > finding survives; its magnitude was overstated. **So the "everything floors below 0.5"
  picture is partly an input artifact, not purely a misplaced threshold.** The operating
  point is less absurdly located than it looked; the input was wrong.
- **Real genomic context behaves like a stable background.** Three *different* 10 kb
  regions give peak scores agreeing to three decimal places on every protein. That is
  the behaviour wanted from a model whose score should describe the CDS: the flank
  supplies context without dictating the answer.
- **Shuffled sequence inflates scores, in 9 of 9 comparisons.** A permutation of real
  genomic sequence is *less* biologically realistic than the original, not more, so a
  higher score there is a distribution-shift response, not better detection. It is the
  control that rules out "any real bases beat `N`" — random ACGT (arm B) scores *below*
  `N`, so the effect tracks composition and structure, not mere non-`N`-ness.

### What it does NOT establish

- **Not that real flanks make the model *right*.** There are no labels here. A higher
  score is not a more correct score, and "real context raises the peak" is equally
  consistent with it becoming more false-positive-prone. Direction of *correctness* is
  unmeasurable on this panel — that is the same wall every label-free result here hits.
- **Not a threshold.** Nothing here licenses moving `DEFAULT_SITE_PROBABILITY`. It says
  the number that threshold is applied to depends materially on the flank, which is an
  argument for supplying the **real construct context** (`ConstructContext`), not for
  moving the cutoff.
- **Scope:** one locus (three regions of one chr1 window), 9 sequences in the main arm
  and 3 in the replication, peak score only. The flank/CDS **seam** is a genuine
  confound — a real construct has such a junction, but the one used here is arbitrary.

### The consequence for BT4

BT4 already has the seam to fix this: `ConstructContext` carries the 5′UTR and backbone,
and `score_in_context` routes them to the CNNs. What this measurement adds is that using
it is **not a refinement — it changes the answer**, by ~0.19 of median peak score, in the
direction that matters at a 0.5 cutoff. A splice number computed on the `N`-padded path
should be read as a lower bound on that model's response, not as its estimate.

## The detection floor — measured (2026-08-20)

#123 showed the `N`-padding suppresses scores but left the central question open: is
Pangolin **blind** inside designed coding sequence, or **correctly silent** because a
clean designed CDS contains no strong splice site? Those produce identical output on a
label-free panel. The way to separate them is to plant a site whose location you know
and see whether the model finds it.

Substitution (not insertion) of a 9-mer spanning the junction, so CDS length and reading
frame are preserved: exon `MAG` | `GTRAGT` intron. Two third-party designed hosts
(KRas4B, Beclin1 — Ranaghan designs, so **not** BT4 output and therefore not pre-depleted
by BT4's own `avoid_splice_sites` rule), three plant positions each, real chr1 flank.

| rung | 9-mer | median peak | cleared 0.5 |
|---|---|---|---|
| **L0** host, unmodified | — | 0.0524 | 0/6 |
| **L5** scrambled | `GAAGCTGAT` (composition-matched) | 0.0525 | 0/6 |
| **L4** GT→CT ablation | `CAGCTAAGT` (7/9 bases = L1) | 0.0543 | 0/6 |
| **L3** weaker | `AAAGTATCT` | 0.0547 | 0/6 |
| **L2** weakened consensus | `CAGGTATGT` | **0.3567** | 0/6 |
| **L1** full consensus | `CAGGTAAGT` | **0.5700** | **5/6** |
| L1 on the **`N`-padded shipped path** | `CAGGTAAGT` | 0.5328 | 4/6 |

### What this establishes

- **Pangolin is not inert in designed coding sequence.** A planted textbook donor lifts
  the local peak from 0.052 to 0.570 — about **11×** — and the peak lands on **exactly**
  the base `CNN_ANCHOR_OFFSETS` predicts (donor = −1) in 5 of 6 plants. The anchor
  convention, previously established on natural genomic panels, is confirmed to hold in
  a designed-CDS host.
- **The response is signal-specific, not change-driven.** Three independent controls all
  sit at host baseline: a composition-matched scramble (0.0525), a *weaker* motif
  (0.0547), and — the decisive one — **GT→CT ablation at 0.0543**, which keeps 7 of the
  9 bases and destroys only the invariant `GT`. Changing one dinucleotide removes the
  entire signal. The model is reading the splice signal, not reacting to an edit.
- **There is a real dose-response**: L1 0.570 → L2 0.357 → everything else ≈ 0.053.

### The floor is high, and that is the finding that matters

**L2 — a genuinely weakened but real donor — scores 0.357 and clears nothing.** The
detection floor at a 0.5 cutoff sits *above* the intermediate-strength sites that cryptic
splicing actually uses. Only a textbook-perfect consensus is reported. Combined with
#123's ~0.19 median deflation on the `N`-padded path, the shipped configuration misses
one plant in six even at full consensus (4/6 vs 5/6).

So the operating point is not merely uncalibrated in the abstract: **it is demonstrably
above the model's own response to a real site of less-than-textbook strength.**

### What it does NOT establish

- **Not correct silence.** A planted textbook consensus is the *easiest possible*
  stimulus, and its ground truth is assumed rather than assayed. Detecting it says the
  model is not inert; it does **not** license "a clean designed CDS is free of cryptic
  sites". That inference is affirming the consequent, and no label-free panel can supply
  it.
- **Not sensitivity.** BT4 can generate unlimited designed sequence; it cannot generate a
  designed sequence *known* to carry a real cryptic site. This ladder measures response
  to a site BT4 chose and placed — which bounds inertness from below and nothing else.
- **Scope:** 2 hosts × 3 positions, donors only, one motif per rung, Pangolin only.
  Acceptors and SpliceAI are untested here. `L5` had to be chosen from the permutations
  of `L1` that introduce no in-frame stop **and do not recreate `GT` at the junction** —
  a first pass used an unconstrained shuffle that hit a stop codon at every plant site
  and was skipped everywhere, i.e. a control that never ran. That constraint is now
  enforced in code (`build_scramble`) rather than left to a lucky draw.

The ladder is committed as `scripts/probe_splice_detection_floor.py`, so this is
re-runnable by anyone holding the weights rather than a one-off. It prints **derived
scalars only** — peak, offset, threshold-cleared — never per-position arrays, which are
the licensed model's output. A probe that would introduce an in-frame stop is reported as
skipped, never dropped silently.

### Where this leaves the "should we train?" question

It substantially supports the **"correct silence"** branch over the **"blind"** branch:
the models demonstrably detect a strong planted site in this exact regime, so their
silence on clean designed CDS is not an inability to see. Training was only the right
instrument under the "blind" branch. What the ladder points at instead is the **operating
point** — a cutoff sitting above the model's response to a real weakened donor — and
that is derived on labelled data, not trained.

## Is the candidate ranking signal or ensemble noise? — measured (2026-08-20)

Step 4. BT4 hands the user **one** candidate. Before asking whether a splice model ranks
*correctly* — which needs labels nobody has — there is a prior question that needs none:
**is the ranking stable at all, or is it ensemble noise?** A ranking that changes with the
training seed cannot be right even in principle, and that is refutable without ground truth.

Pangolin's prediction is the mean of **12 members = 3 CV folds × 4 tissues**. Those tracks
were retained separately (the adapter averages them at `pangolin.py:511-514` and exposes no
per-member seam) and scored under real chr1 flank, per #123.

**A naive split-half was the wrong design and was discarded.** Folds are re-training
replicates — their disagreement is *noise*. Tissues are different biological targets —
their disagreement is *signal about heterogeneity*. A random 6/6 split averages the two
into an uninterpretable number. The correct instrument is a **two-facet generalizability
study** that separates them.

### Floor census — why the ranking statistic had to change first

| protein | (candidate, member) cells floored by the 0.5 hinge | candidates floored on **all 12** members |
|---|---|---|
| KRas4B | 242/360 (67.2%) | 7/30 |
| Beclin1 | 204/360 (56.7%) | 5/30 |
| PDE3A | 138/360 (38.3%) | 6/30 |

Between a third and two thirds of every measurement is destroyed by the hinge, and 5–7
candidates per protein are identically zero on every member. Nothing can be ranked on that.
All statistics below use the background-free `pool_top_k_logit`.

### Variance components and the generalizability coefficient

`y_ift = μ + a_i + φ_f + τ_t + (aφ)_if + (aτ)_it + (φτ)_ft + e_ift`, n=30 designs, 3 folds,
4 tissues, one observation per cell. `Eρ² = σ²_a / (σ²_a + σ²_δ)`.

| protein | σ²_a (candidate) | σ²_a×fold | σ²_a×tissue | σ²_res | **Eρ² U1** | **Eρ² U2** |
|---|---|---|---|---|---|---|
| KRas4B | 2.4984 | 0.0811 | 0.1737 | 0.4300 | 0.975 | **0.959** |
| Beclin1 | 1.6203 | 0.2125 | 0.3466 | 0.2363 | 0.947 | **0.901** |
| PDE3A | 1.7866 | 0.1020 | 0.2267 | 0.2249 | 0.971 | **0.942** |

*U1 = retraining universe (tissue fixed); **U2 = tissue-general**, the headline, because BT4
never asks the user for a tissue and so implicitly claims tissue-generality.*

**The ranking is not noise.** True candidate variance exceeds every error term by 5–10×, and
Eρ² = 0.90–0.96 under the stricter universe. Corroborated by the structured splits, which
are reported instead of a random one: **fold-vs-fold Spearman +0.861 to +0.970** (three
pairwise, each half holding all four tissues), **tissue-vs-tissue median +0.827 to +0.891**
(range down to +0.685). Tissue disagreement exceeds fold disagreement in all three proteins
— σ²_a×tissue is about **2× σ²_a×fold** — which is heterogeneity behaving like heterogeneity,
not a defect.

*(A Jensen check: pooling is convex, so average-then-pool and pool-then-average need not
agree. Here they do — Spearman +0.989 to +0.996 — so the ordering does not depend on which
was used.)*

### But the delivered **pick** is not stable, and that is what BT4 ships

A reliable ordering does not imply a stable winner, because the top candidates are near-ties:

| protein | pick changes across 3 folds? | across 4 tissues? | where those picks sit in the full ensemble's own ranking |
|---|---|---|---|
| KRas4B | **yes**, 2 distinct | **yes**, 2 distinct | ranks 0–1 of 30 |
| Beclin1 | **yes**, 3 distinct | **yes**, 3 distinct | ranks 0, 2, **6, 7** of 30 |
| PDE3A | no — stable | no — stable | rank 0 |

In 2 of 3 proteins the candidate BT4 would deliver **depends on Pangolin's fold and tissue
configuration**. Beclin1 is the worst case and also the lowest Eρ² (0.901): under one tissue
the winner is a sequence the full ensemble ranks **7th of 30**.

### What this establishes, and what it does not

- **Establishes:** the instrument is not noise-limited. One failure mode — "the ranking is
  ensemble noise" — is **excluded**, and no assay was needed to exclude it.
- **Sharpens an earlier result:** the low cross-backend agreement with the PWM baseline
  (+0.614 / +0.195 / +0.162) cannot be explained by Pangolin being unstable, since Pangolin's
  own ranking is highly reliable. The two backends genuinely disagree, and **one of them is
  wrong** — with the prior against the baseline that flags a site every ~14 nt and rates the
  natural gene worst.
- **Does NOT establish that the ranking is correct.** Reliability is not validity: a ranking
  can be perfectly reproducible and perfectly wrong. Every member shares an architecture and
  most of its training data, so a shared blind spot is invisible to all of this.
- **Scope:** Pangolin only (SpliceAI is not installed on this machine), 3 proteins × 30
  designs, one flank locus, donors and acceptors pooled. Aggregating Eρ² across only three
  protein clusters is not meaningful, so per-protein values are reported and no pooled CI is
  offered.

### Consequence

If splice Δ were ever routed into candidate *selection*, the delivered sequence would change
with a retrained Pangolin or a different tissue in 2 of 3 proteins here. That is an argument
for reporting the ranking with its near-ties visible — not for treating the argmax as a
decision.

## How much flank do these models actually need? — measured (2026-08-20)

Prompted by a design question: gene-therapy payloads cap at ~4.7 kb (AAV) to ~10 kb (LVV),
so if a model with a 10,000 nt receptive field *requires* 5,000 nt each side, it is the
wrong tool for BT4. **Receptive field is not requirement**, and the difference is testable.

Detection of a planted textbook donor (the #124 ladder) against available real flank:

| flank/side | host baseline | planted donor | lift |
|---|---|---|---|
| **0** (`N`-padded, shipped) | 0.0532 | 0.5447 | 10.2× |
| 100 | 0.0539 | 0.6042 | 11.2× |
| 250 | 0.0538 | 0.6224 | 11.6× |
| 1,000 | 0.0532 | 0.6129 | 11.5× |
| **1,500** (AAV-scale) | 0.0528 | 0.6216 | 11.8× |
| 5,000 | 0.0522 | 0.6561 | 12.6× |

**Detection saturates around 100–250 nt.** Twenty times more context (250 → 5,000) buys
about 5%. A 1.5 kb CDS inside a 4.7 kb AAV payload has ~1.6 kb of real promoter/UTR/polyA
on each side — comfortably past saturation.

**Independently corroborated by published ablations**, which BT4 had not previously cited:
Jaganathan et al. 2019 trained SpliceAI at 40 / 200 / 1,000 / 5,000 nt per side for top-k
**0.57 / 0.90 / 0.93 / 0.95**; OpenSpliceAI (Chao et al., eLife 2025) retrained at all four
lengths and reports 80→400 nt as the large gain (+62% donors / +74% acceptors) with
400 nt→2 kb and 2 kb→10 kb each worth only ~3–4%. *(Cited from a literature sweep, not
verified against the primary sources here — check before relying on the exact figures.)*

### And the earlier flank effect was an extreme-value artifact, not a change in response

Two statistics on the *same unmodified sequences*:

| flank/side | **LOCAL** (max over a 13 nt window) | **GLOBAL** (max over the whole CDS) |
|---|---|---|
| 0 | 0.0536 | 0.2757 |
| 250 | 0.0535 | 0.3666 |
| 1,500 | 0.0537 | 0.4333 |
| 5,000 | 0.0536 | 0.4622 |

**The model's response at any given position is flat to four decimal places across the
entire flank range.** What moves is the maximum over hundreds of positions — an
extreme-value statistic, which a small distribution shift relocates substantially. BT4's
`pooled_risk` is top-3 over the whole CDS, so **it is BT4's aggregation that is
flank-sensitive, not the model's detection.** Same shape as the pooling-hinge and
saturating-baseline findings above: the defect is in what BT4 computes from the scores.

### Consequence for tool fit

The 10 kb receptive field exists to model **gene architecture** — pairing a donor with an
acceptor across a long intron, using exon/intron context to reject decoys. BT4 does not ask
that question. It asks *"did my redesign create something that looks like a splice site?"*,
which is **local** and saturates well inside any construct BT4 will ever design.

So the construct-size objection does not disqualify these models for **site localization**.
It does bear on the whole-CDS pooled risk, which was already the part with no defensible
operating point.

**Caveat:** n=4 per flank length (2 hosts × 2 positions), planted positives only. The
saturation point is approximate; the qualitative conclusion — well below 5 kb — is robust
across every host and position tried.

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

---

# The operating point — what the literature can and cannot settle (2026-08-21)

Prompted by a direct design question: set the site threshold "low enough that we remove
any real risk of splicing occurring while preserving protein expression." Two halves —
where the threshold should sit, and what it costs. A six-facet primary-source sweep plus
measurements run here.

## The decisive negative result

**No published source thresholds the quantity BT4 thresholds.** BT4 compares a *raw
per-position* pseudo-probability `P(site)` on a *single* designed CDS against a fixed
constant. Neither model's authors ever do this:

| published number | what it actually thresholds | transfers to BT4? |
|---|---|---|
| SpliceAI 0.2 / 0.5 / 0.8 | `max(DS_AG, DS_AL, DS_DG, DS_DL)` — a **delta** between two sequences | **no** |
| Walker 2023 LR bands (≤0.1, 0.1–0.2, ≥0.2) | SpliceAI **max delta** score | **no** |
| Pangolin `-s SCORE_CUTOFF` | `|change in score|` between ref and alt | **no** |
| Pangolin 0.14 (~5% false-sign rate) | **delta** in splice-site *usage*, human vs chimp | **no** |
| SpliceAI top-k 0.95, Pangolin top-1 79% | raw per-position `P(site)` — but with a **floating** threshold set so #predicted = #true | **no** (not a fixed cutoff) |

The one family of published results that *is* about the raw per-position score
characterizes it with **top-k accuracy**, where the threshold floats with the number of
true sites. Neither paper ever says "call a site when raw `P(site) > c`."

So `DEFAULT_SITE_PROBABILITY = 0.5` does not inherit SpliceAI's authority. It is the same
*number* as SpliceAI's recommended cutoff attached to a **different quantity** — and that
coincidence is very likely where it came from.

**Consequence:** the literature cannot hand BT4 this threshold. It has to be derived.

## What the literature *does* settle: direction, and the cost of being wrong

Two opposing pressures, both primary-verified.

**Toward a lower threshold — the cost of a miss is large.**

- A single strong cryptic site removed up to **78%** of protein expression (PrP study,
  PMC11342779); canonical-site mutant fell to 33% of parent.
- **13 of 17** randomly chosen transgenes spliced aberrantly at exon–exon junctions;
  **17 of 17** at the V5 tag (PMC9379414). Cryptic splicing on transgene insertion is the
  *common* case, not an edge case.
- Empirically-optimal thresholds were consistently **lower** than the tool authors'
  recommended cutoffs (Smith & Kitzman 2023).
- Walker 2023's ClinGen calibration concluded 0.5 "may be calibrated too high."

**Toward a higher threshold — the risk is concentrated, and decoys are everywhere.**

- In a lentiviral vector, **4 cryptic sites accounted for ~82%** of aberrant fusion
  transcripts (Cesana et al., *JCI* 2012). The burden is concentrated in a few *strong*
  sites, which a conservative cutoff still catches.
- Consensus-matching pseudo-sites outnumber real donors roughly **21:1**; a CDS is
  saturated with `GT[R]AG`-like matches that never splice.
- At SpliceAI deltaMax ≥ 0.06, **10% of a background set** is called splice-disrupting.

**And the regime penalty applies to every call BT4 makes:** median prAUC **0.419 exonic**
vs **0.773 intronic** (Smith & Kitzman 2023). BT4 designs coding sequence, so it operates
entirely in the weaker half.

## Deriving the threshold — measured in-container against the real weights

The Pangolin weights are present in this container and **all 12 pins match**
`PINNED_WEIGHT_SHA256`, so the following are first-hand measurements, not quotations.

### The donor ladder reproduces, and 0.5 catches nothing

Third-party designed host, 300 nt, shipped `N`-padded path, three plant positions:

| rung | median peak | cleared 0.5 |
|---|---|---|
| L0 host unmodified | 0.0569 | 0/3 |
| L5 composition scramble | 0.0562 | 0/3 |
| L4 `GT`→`CT` ablation | 0.0552 | 0/3 |
| L3 weaker motif | 0.0582 | 0/3 |
| **L2 weakened real donor** | **0.1890** | 0/3 |
| **L1 full consensus** | **0.3206** | **0/3** |

The ordering and the ~6× separation reproduce exactly. But on this host **not even a
textbook consensus donor clears 0.5.**

### The absolute score is flank-dependent; the ratio is not

Same motif, same position, varying only how much designed CDS surrounds it
(non-repetitive flank built from distinct BT4-optimized proteins):

| flank/side | L0 host | L2 weak | L1 full | L1/L0 |
|---|---|---|---|---|
| 0 | 0.0514 | 0.4357 | **0.6554** | 12.7× |
| 100 | 0.0493 | 0.3393 | 0.5484 | 11.1× |
| 250 | 0.0486 | 0.3040 | 0.5084 | 10.5× |
| 500 | 0.0493 | 0.2521 | **0.4487** | 9.1× |
| 1,000 | 0.0473 | 0.1994 | **0.3742** | 7.9× |

**The identical textbook donor scores 0.655 bare and 0.374 inside 1 kb of designed CDS —
crossing 0.5 in one direction purely because of surrounding sequence that contains no
splice signal.** The absolute peak falls 1.75× while the ratio to local background stays
in a tight 7.9–12.7× band.

*This runs opposite to the earlier real-genomic-flank result in this document (0.545 →
0.656 as flank grew). Both are consistent: real genomic flank supplies intron/exon
architecture that makes a site more plausible, while more exon-like designed CDS makes it
less so. The direction of the flank effect depends on what the flank is — which is
precisely why a fixed absolute cutoff cannot be regime-independent.*

### The acceptor arm, built to the published architecture

The prior ladder was donor-only. A realistic acceptor needs branch point (YNYURAC,
18–40 nt upstream) + polypyrimidine tract + `YAG`|`G`, assembled frame-safe (33 nt,
11 codons, zero in-frame stops):

| rung | peak | vs host |
|---|---|---|
| A0 host unmodified | 0.0536 | 1.0× |
| **A1 full (BP+PPT+YAG)** | **0.4879** | 9.1× |
| A2 branch point killed | 0.3486 | 6.5× |
| A3 polypyrimidine tract broken | 0.0856 | 1.6× |
| A4 `AG`→`AC` ablation | 0.0513 | 1.0× |
| A5 composition scramble | 0.2328 | 4.3× |

Biologically correct throughout: the invariant `AG` is essential (A4 collapses to
baseline), the pyrimidine tract is nearly as essential (A3), and the branch point
contributes without being required (A2 partial). **A1 — a fully-formed acceptor — scores
0.4879 and misses the 0.5 cutoff.**

**A5 is a failed control and the failure is instructive:** permuting a pyrimidine-rich
cassette tends to *recreate* a tract followed by `AG`, so the scramble is accidentally a
weak acceptor. An acceptor scramble must be constrained to avoid a terminal `AG` the way
the donor scramble is constrained to avoid a junction `GT`.

### Pooling all 21 measurements

| statistic | lowest true site | highest negative | separation | usable interval |
|---|---|---|---|---|
| **absolute peak** | 0.1890 | 0.0856 | 2.21× | (0.086, 0.189) |
| **peak / local background** | 3.32× | 1.60× | 2.08× | (1.60, 3.32) |

**At the shipped 0.5, 11 of 14 true sites are missed — 79%**, including full-consensus
donors at three of five flank lengths and the fully-formed acceptor. A ratio cutoff of
**3.0×** produces **zero false positives and zero false negatives** across all 21.

An absolute cutoff near **0.13** (the geometric midpoint of the pooled interval) also
separates this set — but its usable window is narrow and, as the flank sweep shows, it
drifts with construct context in a way the ratio does not.

**This converges with the literature's negative result from the other direction.** Neither
model's authors threshold a fixed absolute either: they use top-k (a *floating* threshold)
or a delta against a reference. A background-relative statistic is the same family of
answer.

## Both of BT4's splice defenses are tuned to catch only the strongest sites

Measured here, and the two findings compound.

**1. The in-loop constraint is inert.** `avoid_splice_sites=True` changed **0 of 12**
random proteins' delivered sequence, and cost **0.0000** CAI on KRas4B (1.0000 → 1.0000)
while removing **zero** flagged sites. Across 20 default designs (9,000 nt) it fired
**once**. Both halves are *correct* — positive controls fire on `GTAAGT` (donor) and
`TTTTTTTTTTCTTCTAGG` (acceptor) — but CAI-max drives GC3 to **99.5%** (native KRAS:
31.7%), and the strong consensus motifs are AT-rich, so they essentially never arise.

**2. The out-of-loop CNN threshold misses intermediate sites**, as the ladder shows.

So the constraint catches only full consensus (which never appears) and the threshold
catches only near-consensus (0.57), while the sites cryptic splicing actually uses are
intermediate-strength (0.357). **Neither defense covers the middle of the distribution.**

## A correction to an earlier claim in this document's own framing

Prior sessions asserted that BT4's high GC3 mechanically strips AT-rich splice motifs, so
apparent splice safety was a GC3 artifact. **Directly tested, that claim is not
supported** as stated:

| | native KRAS | bt4 default |
|---|---|---|
| GC3 | 31.7% | 99.5% |
| bare `GT` | 32 | 29 |
| bare `AG` | 57 | 43 |
| `GTAAGT` (textbook donor) | 0 | 0 |
| pyrimidine tract ≥8 + `AG` | 0 | 0 |
| PWM donor sites > 0.5 | 16 | 16 |
| PWM acceptor sites > 0.5 | 32 | 28 |

Motif counts barely move and the PWM flags the *same* number of donors in both. The
accurate statement is narrower: **CAI-max never generates the AT-rich consensus in the
first place**, which is why the constraint is inert — not that GC3 removes motifs from a
sequence that had them. Wobble-base composition does shift drastically (A 34.9% → 0.5%,
T 33.3% → 0.0%), but that did not translate into fewer flagged sites on this instrument.

*Caveat: the PWM is a saturated instrument (below), so this weakens the GC3 claim rather
than settling it. The discriminating test needs the CNN weights.*

## The PWM baseline is barely more selective than the dinucleotide

| | native KRAS | bt4 design |
|---|---|---|
| acceptor flags / all `AG` | 32/57 = **56.1%** | 28/43 = **65.1%** |
| donor flags / all `GT` | 16/32 = **50.0%** | 16/29 = **55.2%** |
| one flagged site per | 11.8 nt | 12.9 nt |

At the shipped 0.5 it calls roughly **half to two-thirds of every `GT` and `AG`** a splice
site. On the acceptor cassette ladder it scored the *host baseline* (0.7916) **higher**
than the `AG`-ablated cassette (0.7479) — it cannot grade the controls at all. This is
what `default()` returns today.

## What this does NOT establish

- **No threshold is calibrated by any of this.** Every positive is a motif planted here,
  on one backend (Pangolin), with **ground truth assumed rather than assayed**. It bounds
  an interval; it does not calibrate a probability. The 3.0x ratio cutoff has zero errors
  on 21 points that were *chosen to be separable* -- that is a sanity floor, not a
  measured error rate.
- **The pooled interval mixes regimes.** Absolute peaks were pooled across flank lengths
  and across donor and acceptor cassettes, which is what narrows the absolute window to
  (0.086, 0.189). A single-regime interval is wider; the narrow one is the honest number
  precisely because BT4 applies one constant across all regimes.
- **It says nothing about specificity.** Every positive is a motif BT4 planted. Detecting
  a planted site is not evidence about sites nobody put there, and lowering the threshold
  to 0.14–0.21 would flag Beclin1 and KRas4B designs **with no way to know whether those
  flags are true**.
- **The two pressures were not jointly optimized.** Cesana's concentration finding (4
  sites → 82%) and the ladder's floor finding point in opposite directions, and no
  measurement here adjudicates between them. The direct test is available and unrun:
  **score Cesana's four named LVV cryptic sites (SA1/SA3/SA4/SD5) and see where they
  fall.** If they clear 0.5, the shipped constant is defensible for LVV; if they sit near
  0.3, it is not.
- **The GC3 correction is measured on one protein with a saturated instrument**, and does
  not license the opposite claim either.
