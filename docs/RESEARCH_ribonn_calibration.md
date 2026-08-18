# Calibrating RiboNN — what the evidence says, and the protocol it implies

*Status: research record + protocol. **No claim here promotes anything.** RiboNN ships
`calibrated=False` and stays that way until a gate passes on a real panel, recorded as an
`ExpressionAttestation` (see `bt4.biomodels.expression.attestation`).*

`docs/NEXT_SESSION.md` item 11 asks whether the wrapped RiboNN head can be promoted to
`calibrated=True`. This document is the evidence behind that question and the protocol
that would answer it. It was assembled from primary sources — the paper's full text, the
upstream repository, and the dataset records themselves — with an adversarial
verification pass over every dataset and statistical claim. Claims an adversarial review
**overturned** are marked **CORRECTED** — several of which corrected statements BT4's
own documentation had been carrying.

The operational runbook that goes with it — environment setup, the free sanity checks,
the panel hunt, pre-registration, the gate command, and the decision tree — is
[`DESIGN_ribonn_calibration.md`](DESIGN_ribonn_calibration.md). This document is the
*why*; that one is the *how*.

## 1. What RiboNN actually is

Zheng D, Persyn L, Wang J, Liu Y, Ulloa-Montoya F, Cenik C, Agarwal V. **"Predicting
the translation efficiency of messenger RNA in mammalian cells."** *Nature
Biotechnology* 44(5):783–796 (2026). DOI `10.1038/s41587-025-02712-x`, PMCID
`PMC12323635`. ("RiboNN" is the repo/Zenodo name, not the paper's title.)

| | |
|---|---|
| **What it predicts** | TE defined as *the residual of a compositional linear regression per transcript per sample*, on centred-log-ratio (CLR) normalised counts |
| **Units** | Dimensionless, log-scale, mean-centred. Bare `nn.Linear` head, MSE loss, **no output activation — never exponentiate** |
| **Outputs** | **78 human / 68 mouse per-cell-type columns**, no aggregate column |
| **Input limits** | 5′UTR ≤ 1381 nt; CDS+3′UTR ≤ 11937 nt; over-length rows are **silently dropped** |
| **Ensemble** | 90 models/species (10 outer × 9 inner folds); `predict` returns **one row per outer fold (10)**, each already a 5-model mean — **CORRECTED**, it is *not* one row per model |
| **Reported accuracy** | r² = 0.62 human / 0.61 mouse for mean TE. **CORRECTED:** the often-quoted r = 0.78/0.74 belongs to their **LGBM baseline**, not RiboNN |
| **Train/test splits** | Random ten-fold over transcripts |

**On the 5′UTR/CDS/3′UTR attribution.** Both numbers in BT4's docs are real and
verbatim: per-nucleotide density **67/31/2**, total global contribution
**22/73/5**. Length-integrated, the CDS carries the **majority**. BT4's existing
framing is correct. Two caveats worth recording: it is Captum saliency×input, not a
retrain-without-region ablation; and the best *measured* decomposition (Höllerer &
Jeschek, 1.2 M 5′UTR–CDS pairs in *E. coli*) puts the CDS at 20%.

## 2. Why "it reproduces the published model" is not calibration

- **The paper never tests synonymous CDS variants of one protein under a fixed UTR.**
  Full-text search: "synonymous" appears only in the introduction and reference
  titles; "codon optimization", "shuffle", "recode", "variant library" appear in no
  experiment. Its codon analysis is *insertional* — codons inserted into the reading
  frame, which **changes the protein**.
- **The one CDS-attributable datapoint is weak:** swapping ORFs under shared 5′UTRs
  gives **r² = 0.11**.
- **Zero-shot on designed reporter mRNAs: r² = 0.17–0.19**, rising to 0.49–0.50 only
  after fine-tuning — against r² = 0.62 on natural genes.
- The authors themselves concede the reverse transfer problem: reporter-trained
  models "offer limited insights into the translation of endogenous mRNAs."
- **No independent replication or adversarial benchmark of RiboNN exists** (22 citing
  works, none critical).

**CORRECTED — drop the leakage argument.** If BT4 ever argues "random CV folds ⇒
homology leakage", that does not survive: the paper reports a homology control
("removing highly homologous test mRNAs led to highly similar r²"). It is a
*qualitative* control with no threshold or numbers, so it is weak — but the charge as
stated is refuted. **The argument that survives is simply that the synonymous-variant
regime was never evaluated.**

**A subtler problem that is real.** RiboNN saw essentially every natural human/mouse
transcript; BT4 deploys it on novel designed CDSs. That breaks the *exchangeability*
that conformal coverage assumes, and because you cannot audit what a frozen
third-party model saw, you cannot bound the penalty. This is the one place a frozen
model genuinely hurts.

## 3. What "calibrated" has to mean, precisely

Two separate claims that must never be conflated:

**Discrimination (ranking).** Spearman ρ is the honest primary metric — invariant to
any strictly **increasing** transform. **CORRECTED:** not "any monotone" — a
*decreasing* map flips the sign, and the sign convention of a CLR residual against an
assay readout is not guaranteed, so the sign must be checked, not assumed.

**Calibration proper (intervals).** Split conformal:
`q̂ = ⌈(n+1)(1−α)⌉/n` empirical quantile of calibration residuals, giving
`1−α ≤ coverage ≤ 1−α + 1/(n+1)`.

Two facts that discipline the whole design:

1. **CORRECTED — coverage stays *valid* under arbitrary unit mismatch.** The theorem
   holds for any score function, however uninformative. Unit mismatch destroys
   **sharpness**, not validity. So "the units differ, therefore coverage is
   meaningless" would be wrong. The honest statement: unrecalibrated intervals are
   valid and **useless**. The metric that exposes this is **median interval width ÷
   label IQR**, which must be reported next to coverage.
2. **A constant predictor passes the coverage test exactly.** For `f(x) = c`, the
   interval `[c−q̂, c+q̂]` has exactly valid coverage. So the **constant predictor must
   be a permanent baseline**, precisely so that its coverage "pass" is visible in the
   report.

**CORRECTED — a frozen model needs only two splits, not three.** There is no fitting
step, so the whole panel is calibration data: split 1 fits the link `measured ≈ a ×
pred + b`, split 2 conformalizes the *recalibrated* residuals. This halves the data
cost. (Fitting the link on the conformal split would break the independence the
theorem requires.)

**CORRECTED — do not report "calibration slope = 1".** Once you fit `a`, `b`, the
slope is 1 *by construction*. Report instead the **stability of (a, b) across held-out
protein groups** and the residual spread.

## 4. The data problem — no dataset fully qualifies

There is **no** mammalian panel of synonymous CDS variants of one protein under a
fixed UTR, with ribo-seq-derived log-TE and downloadable per-variant sequences +
measurements, at usable n. mRNABench, the field's curated benchmark, states it has no
synonymous-variant task.

| Dataset | System | Measured | n variants | n proteins | UTRs fixed | Download | Licence | Verdict |
|---|---|---|---|---|---|---|---|---|
| **PERSIST-seq** (Leppek 2022, *Nat Commun* 13:1536) | HEK293T | **ribosome load** + in-cell half-life | 121 CDS variants (~64 strictly synonymous Nluc) | 1 per arm | ✅ hHBB 5′/3′ | ✅ HuggingFace `morrislab/mrl-hl-lbkwk`, 203 rows, Parquet | CC BY 4.0 | **PARTIAL — best mammalian option.** Right readout, small n, one protein |
| **Mauger 2019** (*PNAS* 116:24075, Moderna) | HeLa + mouse | luminescence / fluorescence / ELISA | 43 Luc + 30 eGFP-degron + 9 hEpo + 4 eGFP | **4** | ✅ "identical 5′ and 3′ UTRs" | Sequences in Dataset S2; **whether Dataset S1 holds per-variant measurements is unresolved** | CC BY | **PARTIAL / OPEN — highest-value check remaining** |
| **iCodon library** (Diez 2022, GSE207584) | Zebrafish embryos | mRNA decay | ~1,395 measured | **100** | ✅ one shared UTR pair | ✅ GEO FASTA + CSVs | CC BY 4.0 | **PARTIAL.** Ideal *design* (100 groups!), wrong species, wrong readout |
| **Mordstein 2020** (*Cell Syst* 10:351) | HeLa / HEK293 | flow-seq protein + polysome | 217 | 1 (GFP) | ✅ | ❌ SRA raw reads only (PRJNA596086) | CC BY | NOT USABLE off-the-shelf; largest n, needs reprocessing |
| **CodonBERT MLOS** (Li 2024) | HeLa | ELISA protein | 167 released (543 in paper) | 2 | claimed, **UTRs never disclosed** | `CDS,Value` only | Sanofi non-commercial | **NOT USABLE** — a full-length model needs UTRs |
| RiboDecode 2025 | HEK293T + mice | luminescence/ELISA | ~7/protein | 4 | ✅ | no accession | CC BY-NC-ND | NOT USABLE — n=7 gave **p = 0.077** |
| Nieuwkoop 2023 / Kudla 2009 / Schmitz 2021 / Cambray 2018 | *E. coli* | GFP/RFP | 154 – 228k | 1 | ✅ | mostly ✅ | mixed | NOT USABLE (species) — but the **best design templates** |
| **RiboNN's own labels** (`TE_classic_ML`, CenikLab) | natural genes | CLR-TE + a `fold` column | ~11k | many | native | ✅ GitHub | GPL-3.0 | **Adapter validation only** — not a calibration panel |

**Also verified: the two panels already in this repo carry no measurements.**
Ranaghan et al. 2021 measured expression for exactly **one** sequence (GeneArt
KRas4B, 23 ± 4 mg/L soluble yield, ***E. coli***). No per-variant expression data
exists for Beclin 1 or PDE3A. So `scripts/data/ranaghan2021_tab4.fasta` is a
**sequence-only** resource — perfect for the free sensitivity checks in Stage 1,
useless as a validation panel.

**Generating one yourself is out of scope but worth the number:** ~300 full-length
variants ≈ **$15–19k in synthesis alone**, before IVT/transfection/sequencing;
6–12 months; mid-five to low-six figures.

## 5. What a defensible protocol looks like

- **Grouping unit = the protein.** Variants of one protein share length, amino-acid
  composition, UTR and assay batch. The effective sample size for any cross-protein
  claim is **G (number of proteins)**, not N (rows). Report both.
- **Two group-disjoint splits** (fit-the-link, then conformalize).
- **Sample sizes — two different answers.**
  - *Coverage:* hard floor `n ≥ 1/α − 1` (**9** at 90%) — matches exactly what
    BT4's `conformal_quantile` does (returns `+inf` below it). For slack ±0.10 you
    need n≈22; **±0.05 needs n≈102**; ±0.01 needs n≈2491.
  - *Ranking:* for 80% power at α=0.05, ρ=0.5→n≈33, ρ=0.4→n≈51, **ρ=0.3→n≈89**,
    ρ=0.2→n≈198. Assay noise attenuates by √reliability, so a true 0.4 at ribo-seq
    replicate reliability ~0.5 reads as 0.28 → n≈100.
  - **A novel-*protein* coverage claim needs ~100 held-out proteins — unreachable. A
    novel-*variant-of-a-known-protein* claim needs ~100 variants — reachable.** This
    is the single most important sizing fact in the whole plan.
- **Primary metric:** Spearman computed *within* protein, aggregated across proteins
  (ProteinGym's aggregation stage). Secondary: Kendall τ-b, precision@k at BT4's real
  candidate-set size, and coverage with a **Clopper-Pearson CI at ≥3 α levels**, plus
  median interval width ÷ label IQR.
- **Permanent mandatory baselines** (all must be beaten, and kept in the report
  forever): group-level label permutation, **CAI-only**, GC3-only, CDS-length-only,
  and the **constant predictor**.
- **CORRECTED — there is no community-standard absolute Spearman cutoff.** Any bare
  "0.4" or "0.5" would be invented. The threshold must be phrased as *"the
  cluster-bootstrap CI lower bound exceeds every baseline"*, with a pre-registered
  absolute floor stated as a pre-commitment, not a standard.
- **Homology grouping:** MMseqs2 `--min-seq-id 0.5 -c 0.8`. **CD-HIT cannot do 30% —
  its floor is ~40%.**
- **Readout choice matters enormously.** Mauger: protein output correlates r = 0.90
  with mRNA half-life but only **r = 0.45** with translation rate. Comparing raw
  protein output to a CLR-residual TE re-introduces exactly the term TE divides out.
  **Use log(protein / mRNA), or mean ribosome load — never raw protein output.**
- **Confounders to regress out:** GC/GC3, CpG count *and clustering*, CDS length,
  start-region ΔG, cryptic splice / premature polyA creation.
- **Holdout reuse is a real cost:** adaptive reuse degrades error from
  O(√(log k/n)) to O(√(k log n/n)), so a few-hundred-row panel supports only a
  handful of honest runs. One locked run is the discipline.

## 6. The honest prior

r² = 0.62 on natural genes, **0.17–0.19 zero-shot on designed reporters**, **0.11**
on the only CDS-attributable test. Plan for the gate to fail, and make failing cheap.

---

---

## What this implies for BT4, concretely

The machinery this document motivated has landed and is described in `CHANGELOG.md`:

| Need identified above | Where it lives now |
|---|---|
| Pooled scoring credits between-protein skill | `verify_expression_gate(within_group=True)` |
| Arbitrary units make residuals meaningless | `verify_expression_gate(recalibrate=True)` |
| A constant predictor passes coverage | `ExpressionGateReport.width_over_iqr` |
| Variants of one protein are a dependent cluster | cluster bootstrap over whole groups |
| A head must beat CAI to have added anything | `bt4.pipeline.expression_gate.BASELINES` |
| A result must bind to exact bytes | `ExpressionPanel.content_hash()` |
| A claim is scoped, not a bare boolean | `ExpressionAttestation` (species/cell types/readout) |
| Is it even worth buying a panel? | `scripts/ribonn_sensitivity.py` |

**The honest prior remains that this gate fails.** RiboNN scores r² = 0.62 on natural
genes, 0.17–0.19 zero-shot on designed reporters, and 0.11 on the only CDS-attributable
test its own paper reports. The exercise is designed so that a failure is a cheap,
recorded, publishable outcome rather than a sunk cost — and so that a *pass* means
something when it happens.
