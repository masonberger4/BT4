# Splice calibration — measured results (variant-effect half)

**Status: the variant-effect half of Part B has been run. The site-prediction half has
not.** Nothing here changes any backend's `calibrated` flag, and nothing here is a
calibration claim. This is a record of what was measured, by whom, on which bytes.

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

## Still to run

The **site-prediction** half: build a panel from a pinned GENCODE release with
`scripts/make_gencode_splice_panel.py` and gate it with `bt4 splice-gate --cnn-anchors`.
That half needs the ~3 GB genome download and the licensed weights, and it is the half
that exercises BT4's own adapters rather than the benchmark's numbers.
