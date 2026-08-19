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

## Results — `panel20.tsv`, Pangolin, `--cnn-anchors`

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

## Still to run

- The **SpliceAI** integration-fidelity gate (needs the CC BY-NC weights and a TF 2.15
  environment), and the same site-prediction panel scored by SpliceAI — which would give
  a genuine **two-backend agreement** figure, the first-class uncertainty signal of §6.
- A panel that reaches BT4's actual regime: designed synonymous CDS variants, where the
  question is specificity rather than recall.
