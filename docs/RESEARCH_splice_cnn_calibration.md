# Reproducing SpliceAI + Pangolin locally, and building a validation panel

Background research for the splice-CNN calibration work
([`DESIGN_splice_cnn_calibration.md`](DESIGN_splice_cnn_calibration.md) is the
procedure; this file is the evidence behind it). Relates to CLAUDE.md §6 (wrapped
splice backends) and §4.3 (the calibration gate).

Everything below was verified against primary sources on **2026-08-17** unless
marked **[UNVERIFIED]**.

---

## 0. Headline findings

1. **`Illumina/SpliceAI` was archived by its owner on 2026-04-20 and is now
   read-only.** No upstream fix will ever land for the TF/Keras 3 breakage. The
   maintained path is the Broad's fork `bw2/SpliceAI` (v1.3.4) plus the Broad's
   published Docker images.
2. **All 17 weight hashes BT4 currently pins are correct.** The upstream bytes
   were independently downloaded and re-hashed: 5/5 SpliceAI `.h5` and 12/12
   Pangolin `.3.v2` match `PINNED_WEIGHT_SHA256` exactly.
3. **Pangolin ships two disjoint model sets and the repo's own two entry points
   disagree about which to use.** `pangolin/pangolin.py` (the CLI) loads
   `final.{1,2,3}.{0,2,4,6}.3.v2` (12 fine-tuned models); `scripts/custom_usage.py`
   loads `final.{1..5}.{i}.3` (base models). Same repo, different numbers. BT4 pins
   the CLI's v2 set — correct, but the divergence must be documented.
4. **`pip install pangolin` installs an unrelated package** (a probabilistic
   programming language). Pangolin-the-splice-model is GitHub-only.
5. **The regime BT4 designs in — exonic / synonymous variants — is where both
   models are weakest.** Smith & Kitzman 2023 measured median prAUC 0.773 for
   intronic vs **0.419 for exonic** variants. This must be stated wherever BT4
   surfaces splice numbers on a CDS.

---

## 1. SpliceAI — install and run

### 1.1 Environment (the fragile part)

`spliceai` 1.3.1 (2020-03-07) declares only `keras>=2.0.5`, `tensorflow>=1.2.0`,
`pyfaidx>=0.5.0`, `pysam>=0.10.0`, `numpy>=1.14.0`, `pandas>=0.24.2` — no upper
bounds and no `requires_python`. Those unbounded pins are the whole problem.

Known-good pin (recommended starting point):

```bash
conda create -n spliceai python=3.10 -y && conda activate spliceai
pip install "tensorflow==2.15.*"     # last TF whose bundled Keras is 2.x
pip install "numpy<2" "pandas<2.2" "setuptools<81"
pip install spliceai==1.3.1
```

Failure modes to expect:

- **TF >= 2.16 defaults to Keras 3**, which will not load these legacy Keras-2
  `.h5` graphs. Either stay on TF 2.15, or install `tf_keras` and export
  `TF_USE_LEGACY_KERAS=1`. (BT4's `_import_keras` already falls back to
  `tf_keras`; TF 2.15 is the verified-safe path.)
- `WARNING:absl:No training configuration found in the save file...` on load is
  **benign** for inference — the Broad's own container prints it.
- `spliceai/utils.py::one_hot_encode` calls `np.fromstring(seq, np.int8)`. Still
  present in NumPy 2.x but long-deprecated; pin `numpy<2` or patch to
  `np.frombuffer(seq.encode(), np.int8)`.
- `pkg_resources` (used by `utils.py` and `__init__.py`) is removed in
  setuptools >= 81 → pin `setuptools<81`.

Alternative that sidesteps all of the above:

```bash
docker run -p 8080:8080 docker.io/weisburd/spliceai-38:latest
curl "http://localhost:8080/spliceai/?hg=38&variant=chr8-140300616-T-G"
```

### 1.2 Licensing

- **Code: PolyForm Strict License 1.0.0** — noncommercial only; forbids
  distributing the software or making derivative works.
- **Weights (`spliceai/models/`): CC BY-NC 4.0**, academic / non-commercial;
  commercial use needs an Illumina licence.
- `setup.py` and the PyPI metadata say `license='GPLv3'`. **This is stale and is
  contradicted by the authoritative LICENSE file.** CLAUDE.md §6 already flags
  this; the flag is confirmed correct.

### 1.3 Weights and reference data

Weights ship **inside the pip package** — no separate download.
`spliceai-1.3.1.tar.gz` is 16,669,756 B, SHA-256
`65c76b012ffd2ca97ca96d7f4c0897c78b9aba4d4ca4068331f7fb5cd5c3b7e1`.

| File | Bytes | SHA-256 |
|---|---|---|
| `spliceai/models/spliceai1.h5` | 3,131,720 | `e1fd5adcef7489d604b10e79c40078ef790d51ef048c4ce3869c9119ac5de42b` |
| `spliceai/models/spliceai2.h5` | 3,131,720 | `6ab042b82ab966b6d3582cb31b96f0859ea08a864f168d69e83aa14450a3b66e` |
| `spliceai/models/spliceai3.h5` | 3,131,720 | `e2e790bde53dfdf410c6dc434a86122a7d12f3f38dc2ef45d85986e9ecf22fad` |
| `spliceai/models/spliceai4.h5` | 3,131,720 | `ca88ac9e58e69ba6fdeed319b72f063f164c9abf7392eaccef903e94c1d99dd6` |
| `spliceai/models/spliceai5.h5` | 3,131,720 | `791cd22c62a80a08d2ca674615a93ce8159d7b55bd157cfef2983b1bd6b41391` |
| `spliceai/annotations/grch37.txt` | 4,418,184 | `16e495da074e965f46f60ced4e73292843e949051d2dc31dd521da3fd7615dce` |
| `spliceai/annotations/grch38.txt` | 4,250,439 | `565c0ea83f66b2210182364c3d626fef526229e3d02f444bab247329b6695f2c` |

Annotation format (TSV, 20,275 / 19,306 lines incl. header):
`#NAME  CHROM  STRAND  TX_START  TX_END  EXON_START  EXON_END`.
`-A grch37` = GENCODE **V24lift37** canonical; `-A grch38` = GENCODE **V24**
canonical (per `__main__.py` help text — more precise than the README's "V24").

Genome FASTA (use `hgdownload.soe.ucsc.edu`, not `.cse.`, which has a certificate
mismatch):

- `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz` (948,731,419 B)
- `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` (983,659,424 B)

### 1.4 Programmatic use (matches BT4's adapter convention)

```python
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np

context = 10000                       # 5000 nt each side
models = [load_model(resource_filename('spliceai', f'models/spliceai{i}.h5'))
          for i in range(1, 6)]
x = one_hot_encode('N'*(context//2) + seq + 'N'*(context//2))[None, :]
y = np.mean([m.predict(x) for m in models], axis=0)   # (1, len(seq), 3)
acceptor_prob, donor_prob = y[0, :, 1], y[0, :, 2]    # channel 0 = null
```

Ensemble = unweighted mean of 5 models. Encoding is **position-major `[L][4]`**
(Keras channels-last). Output is a per-position 3-way softmax
`[null, acceptor, donor]`, so BT4 populating **both** `donor` and `acceptor` is
correct.

### 1.5 A free, committed fidelity oracle

Upstream ships `examples/input.vcf` (10 records) **and** `examples/output.vcf`
with exact 2-decimal delta scores — e.g.
`19:38958362 C>T` → `T|RYR1|0.00|0.00|0.91|0.08|-28|-46|-2|-31`;
`2:179415988 C>CA` → `CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29`;
`21:47406854 CCA>C` → `C|COL6A1|0.04|0.98|0.00|0.00|-38|4|38|4`.

```bash
spliceai -I examples/input.vcf -O out.vcf -R hg19.fa -A grch37   # -D 50, -M 0 defaults
diff <(grep -v '^##' out.vcf) <(grep -v '^##' examples/output.vcf)
```

**Caveat:** this exercises the *variant* path (`get_delta_scores`, which does
ref/alt one-hot, strand flip and indel realignment), **not** `score_sequence`. A
raw per-position panel still needs local capture. A clean diff does prove the
install is sound.

---

## 2. Pangolin — install and run

### 2.1 Environment

```bash
conda create -n pangolin python=3.10 -y && conda activate pangolin
# PyTorch per https://pytorch.org/get-started/locally/ — the repo pins no version
conda install -c conda-forge pyvcf     # only for VCF CLI input; see caveat
pip install gffutils biopython pandas pyfastx
git clone https://github.com/tkzeng/Pangolin.git && cd Pangolin && pip install .
```

- **Do not `pip install pangolin`** — a different project entirely.
- License: **GPL-3.0** (the LICENSE file is verbatim GNU GPLv3). Version 1.0.2,
  repo actively maintained (not archived).
- `setup.py` does not list `torch` as a dependency at all, and README says only
  "Python 3.6 or higher and conda".
- **PyVCF caveat:** PyPI `PyVCF` 0.6.8 last shipped 2016-03-18 with classifiers
  topping out at Python 3.4. On modern Python use conda-forge's `pyvcf`, or the
  maintained fork `PyVCF3` (1.0.4). CSV input avoids `vcf` entirely.
- Checkpoints are `torch.save` zip archives (254 entries), so
  `torch.load(..., weights_only=True)` works and is the safe default on torch >= 2.6.
- Colab (GPU, no install):
  `https://colab.research.google.com/github/tkzeng/Pangolin/blob/main/PangolinColab.ipynb`

### 2.2 Model files

`pangolin/models/` holds **64** files, every one 2,877,321 B:

- `final.{1..5}.{0..7}.3` — 40 base models (5 training runs × 8 heads)
- `final.{1..3}.{0..7}.3.v2` — 24 fine-tuned models (3 runs × 8 heads)

Head index → `0/1` Heart P(splice)/usage, `2/3` Liver, `4/5` Brain, `6/7` Testis.
The `.v2` files are the paper's fine-tuned models (label smoothing, human-only,
3 runs per tissue), used for its variant-effect figures — which is why the CLI
ships them.

The **CLI** uses only the 12 fine-tuned P(splice) models — verified SHA-256:

| File | SHA-256 |
|---|---|
| `final.1.0.3.v2` | `f0478fab173b75f7f7e9fe96688bad6c50fa4a46d70557f423b110caaf565501` |
| `final.2.0.3.v2` | `c4c6bb4880fa6fb28b14182ae3ea0600edb07056158f55325b5e6e6e48fc9f26` |
| `final.3.0.3.v2` | `ec685a6e7105a4486c1f89a005458a13deb3fe7171f13d434f4877e386d10676` |
| `final.1.2.3.v2` | `559c05de3e1ce65c2515ca3e92ef85edb0ec2e47686ca58060e25891ce06eb3a` |
| `final.2.2.3.v2` | `48758ba8b95eee9aa9feea52672ef06ca1b34111299c27f8a710f734d8b9aae5` |
| `final.3.2.3.v2` | `7cb576c2b24db4fdd6970c4ca4fb7c20ae1b1d8ae80645ebbe689848b5743129` |
| `final.1.4.3.v2` | `c50b12e0c0af776d5674ca5e346493f8265783494d4df383364de9c1136657f6` |
| `final.2.4.3.v2` | `e03303bed4fd6f135ec0f6c1b192cce954ea42d0646f44d17b4a6fbb2b1f610e` |
| `final.3.4.3.v2` | `9476d2e25520d7ff15bece0cd5d3b657e3b1dd3cc5fcab1d9c3b62bea7a0c5b6` |
| `final.1.6.3.v2` | `2aae563fa18a8a9b6699c6c96e0d32b8ec7543f8f805fb3bc9de77302cc9f66e` |
| `final.2.6.3.v2` | `7d3c0b1b2a60067b940dec315567874fbc8bcd322f1b7c76bf969f51f0f53f7f` |
| `final.3.6.3.v2` | `756e7721a382cace24e9bfea5b543af5623f2487d9a3efe7385e9c76367005fd` |

### 2.3 Architecture and output

`L=32`, `W=[11×8, 21×4, 41×4]`, `AR=[1×4, 4×4, 10×4, 25×4]`, 16 residual blocks,
skip connections added before blocks 1, 5, 9, 13. `CL = 2·Σ AR·(W−1) = 10,000` —
the same receptive field as SpliceAI. Input must be `N >= 10,001`; output covers
the middle `N − 10,000` bases.

`forward()` returns **12 channels** = 4 tissues × (2-way softmax + 1 sigmoid):
P(splice) at indices `[1, 4, 7, 10]`, usage at `[2, 5, 8, 11]`. P(splice) is a
*binary* softmax — donor and acceptor are **not** separated, which is why BT4 puts
it in `SpliceResult.donor` with `acceptor` all-zero. Encoding is **channel-major
`[4][L]`** (note the transpose), the opposite of SpliceAI.

### 2.4 Reference data + run

```bash
python scripts/create_db.py gencode.v38.annotation.gtf.gz   # filters Ensembl_canonical
pangolin examples/brca.vcf GRCh37.primary_assembly.genome.fa.gz \
         gencode.v38lift37.annotation.db brca_pangolin
```

Output format `gene|pos:largest_increase|pos:largest_decrease|Warnings:`. Options:
`-c`, `-m {True,False}` (default True), `-s CUTOFF`, `-d 50`, and an
undocumented-in-README `--score_exons`. Prebuilt GENCODE v38 databases are linked
from the README (Dropbox).

**The repo ships `examples/brca.vcf` and `examples/brca.csv` as inputs only — there
is no committed expected-output file** (`examples/brca_pangolin.vcf` is a 404). A
Pangolin fidelity panel must be captured locally, not diffed against the repo.

---

## 3. Published validation numbers

### SpliceAI (Jaganathan et al., *Cell* 2019, doi:10.1016/j.cell.2018.12.015)

- **Train:** genes on chr **2, 4, 6, 8, 10–22, X, Y**. **Test:** chr **1, 3, 5, 7,
  9**, paralog-filtered. Confirmed from three independent open-access sources
  (Pangolin 2022; CI-SpliceAI 2022; OpenSpliceAI 2025), since the paper itself is
  paywalled and not in PMC.
- Trained on GENCODE v24 / GRCh37, one primary isoform per gene, enriched with
  novel junctions common in GTEx.
- **Top-k accuracy 0.95** — the paper's headline claim.
- **PR-AUC 0.98 — [UNVERIFIED]** against the primary source. Widely quoted, but
  cell.com is Cloudflare-blocked and the paper is not in PMC/Europe PMC.
- Delta-score thresholds (README, authoritative): **0.2 = high recall, 0.5 =
  recommended, 0.8 = high precision**, where delta = max(DS_AG, DS_AL, DS_DG, DS_DL).
- OpenSpliceAI's paralog audit found **0.71% of MANE transcripts on chr 1,3,5,7,9
  are in fact paralogous to training sequences** — the split is not perfectly clean.

### Pangolin (Zeng & Li, *Genome Biology* 23:103, 2022, CC BY — full text read)

- Same train/test chromosome split as SpliceAI; test genes with mean TPM < 2.5
  excluded; non-human training genes filtered to exclude orthologs/paralogs of test
  human genes.
- Trained on **4 species** (human, rhesus macaque, rat, mouse) × **4 tissues**
  (heart, liver, brain, testis); RNA-seq from ArrayExpress E-MTAB-6798 / -6811 /
  -6813 / -6814.

| Benchmark | Pangolin | SpliceAI |
|---|---|---|
| Splice-site prediction, top-1 | **79%** | **75%** |
| Splice-site prediction, top-0.5 | **94%** | **87%** |
| Splice-site AUPRC | **0.85** | **0.77** |
| MFASS variant effect, AUPRC | **0.56** | **0.47** |
| MaPSy, Pearson r | **0.61** | **0.50** |

MMSplice / HAL / MaxEntScan all scored < 37% top-1 and < 0.30 AUPRC. On MFASS at
80% precision Pangolin's recall is 29%; near splice sites (0–9 nt) AUPRC is 0.75,
beyond 9 nt it falls below 0.35. FAS exon 6: Spearman 0.79 (189 single
substitutions), 0.80 (3,059 combinations). Tissue-specific usage: Spearman
0.35–0.50, median 0.43.

### The number that matters most for BT4

Smith & Kitzman, *Genome Biol* 24:294 (2023), PMC10734170, CC BY: 8 predictors on
3,616 variants across 5 genes. **Median prAUC 0.773 for intronic vs 0.419 for
exonic variants.** SpliceAI and Pangolin were the top two overall, but the paper
concludes improvements are needed *"especially within exons"*, and that concordance
is *"lower for exonic than intronic variants, underscoring the difficulty of
identifying missense or **synonymous** SDVs."*

Also from the same paper: **gene-model annotation choice altered SpliceAI's
predictions for > 10% of variants** in POU1F1 and WT1, projected to affect roughly
1 in 5 human genes — so which annotation a panel used belongs in its provenance.

---

## 4. Statistical calibration (distinct from integration fidelity)

- **Walker et al., *AJHG* 2023** (ClinGen SVI Splicing Subgroup), PMC9980257 — the
  only likelihood-ratio calibration of SpliceAI found. On 2,736
  non-canonical-splice-site variants with in-vitro splicing assays across BRCA1,
  BRCA2, MLH1, MSH2, MSH6, PMS2, NF1, POU1F1:

  | SpliceAI Δ | Likelihood ratio | Interpretation |
  |---|---|---|
  | ≤ 0.1 | 0.17 [0.14–0.21] | moderate evidence **against** spliceogenicity |
  | 0.1 – 0.2 | 1.00 [0.71–1.39] | **uninformative** |
  | ≥ 0.2 | 15.99 [13.23–19.32] | moderate evidence **for** |

  On BRCA1 SGE, SpliceAI reached AUC 0.959, but the authors conclude **"SpliceAI
  score ≥ 0.5 may be calibrated too high,"** excluding many genuinely spliceogenic
  variants; they recommend 0.2 / 0.1 operating points.

- **Chao et al., OpenSpliceAI, *eLife* 2025**, PMC12575001 — the only true
  reliability/ECE analysis. Applies **class-wise temperature scaling** (a
  one-parameter Platt variant), reporting reliability diagrams, ECE and NLL
  before/after on a 10% held-out calibration split. Finding: SpliceAI-architecture
  models are **slightly overconfident** (temperature > 1); after calibration,
  donor/acceptor scores move away from the saturating 0/1 extremes. Ships a
  `calibrate` subcommand — a ready-made recipe to mirror. Also the source of the
  N-padding artifact finding already cited in
  [`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md).

- **Sullivan et al., *AJHG* 2025**, PMC12081236 — argues the Δ score *"does not
  consistently mirror the actual outcome"*: both tools correctly flagged a PKD1 3′SS
  disruption but could not say which alternative site would be used.

- **Smith & Kitzman 2023** — optimal Youden-J thresholds **varied widely across
  exons and variant classes** for most tools (SpliceAI, Pangolin and SQUIRLS were
  the least variable). Recommends normalizing tools by genome-wide call rate — the
  score threshold at which each calls 5/10/20% of a 500k background set — rather
  than by nominal score.

**Net:** there is no published evidence that raw SpliceAI/Pangolin scores are
calibrated probabilities. Both are ranking scores whose thresholds must be
re-derived per use case. **[UNVERIFIED]** whether any published work applies
isotonic regression to SpliceAI deltas; temperature scaling (OpenSpliceAI) and LR
binning (Walker) are the two approaches found.

---

## 5. Panels for a sites / non-sites validation set

| Panel | Content | Licence | Location |
|---|---|---|---|
| **Smith & Kitzman 2023** | 3,616 variants, 5 genes, **8 tools pre-scored** (incl. SpliceAI + Pangolin), plus 500,000 background exonic SNVs | **MIT** ("All code and data are provided under an MIT license") | `github.com/kitzmanlab/splicebench2023` + Zenodo + GB additional files 2–4 |
| **SpliceVarDB** | 50,715 experimentally assessed variants, 8,000+ genes; 13,673 splice-altering, 55% outside canonical sites | paper CC BY 4.0; DB redistribution terms **[UNVERIFIED]** | `splicevardb.org` |
| **MFASS** | 27,733 ExAC variants across 2,198 exons; 1,050 SDVs (3.8%); 83% outside canonical sites | **no LICENSE file** — check before redistributing | `github.com/KosuriLab/MFASS`, `processed_data/snv/snv_data_clean.txt` |
| **Vex-seq** | 2,059 variants, 110 alternative exons | GEO terms | GEO `GSE113163` |
| **MaPSy** | 4,964 **exonic** disease-causing (HGMD) variants; ~10% affect splicing | publisher | Soemedi 2017, *Nat Genet*, supplementary |
| **BRCA1 SGE** | 13 exons saturation genome editing | publisher | Findlay 2018, suppl. table 1 |
| **ClinVar** | splice-altering variants | public domain | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` |
| **GENCODE** | annotation + genome, for building labels directly | open | `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/` |

Verified GENCODE URLs (all HTTP 200):

```
.../release_44/gencode.v44.annotation.gtf.gz          49,721,965 B
.../release_44/gencode.v44.basic.annotation.gtf.gz    29,570,410 B
.../release_44/GRCh38.primary_assembly.genome.fa.gz  844,691,642 B
.../release_38/gencode.v38.annotation.gtf.gz          46,556,621 B   # Pangolin's release
.../release_24/gencode.v24.annotation.gtf.gz          38,708,036 B   # SpliceAI's release
```

**Most CDS/synonymous-relevant for BT4:** the BRCA1-SGE synonymous subset inside
the MIT-licensed Kitzman panel; MaPSy (exonic, protein-code-altering); and the
exonic slice of MFASS. All three sit in exactly the regime where prAUC drops to
~0.42.

Donor/acceptor extraction from GENCODE: take `exon` features of a chosen
transcript set (filter `tag "Ensembl_canonical"`, or use MANE as OpenSpliceAI
does); for a `+`-strand transcript the donor is the base after each exon end and
the acceptor the base before each exon start (reverse for `−`); drop the
transcript's first acceptor and last donor. Negatives are the other positions in
the same gene bodies — exactly the top-k denominator both papers use. SpliceAI's
own bundled `annotations/grch38.txt` is already a flat `EXON_START`/`EXON_END`
table if you want to skip GTF parsing.

---

## 6. Follow-ups this research implies for BT4

1. **Correct the archival status** in CLAUDE.md §6 and `biomodels/splice/spliceai.py`:
   `Illumina/SpliceAI` archived 2026-04-20; name `bw2/SpliceAI` + the Broad Docker
   images as the live install path.
2. **Document the Pangolin `.3` vs `.3.v2` split** in `biomodels/splice/pangolin.py`'s
   module docstring — the pins are right, but the reason is subtle and a maintainer
   following `custom_usage.py` would silently produce different numbers.
3. **Tighten `bt4[splice-spliceai]`**, which declares `tensorflow>=2.6` with no
   upper bound — a fresh install pulls a TF that cannot load the weights.
4. **Keep "integration fidelity" and "statistical calibration" as separate words**
   in BT4's docs. Walker 2023 and OpenSpliceAI 2025 show the second is a genuinely
   different, still-unmet gate.
5. **State the exonic penalty** (prAUC 0.419 vs 0.773) wherever BT4 reports splice
   risk on a CDS — it is the honest scope limit for a codon-optimization tool.

---

## Sources

[tkzeng/Pangolin](https://github.com/tkzeng/Pangolin) ·
[Zeng & Li 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02664-4) ·
[Illumina/SpliceAI](https://github.com/Illumina/SpliceAI) ·
[Jaganathan et al. 2019](https://doi.org/10.1016/j.cell.2018.12.015) ·
[Smith & Kitzman 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10734170/) ·
[kitzmanlab/splicebench2023](https://github.com/kitzmanlab/splicebench2023) ·
[Walker et al. 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9980257/) ·
[OpenSpliceAI 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12575001/) ·
[Sullivan et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081236/) ·
[CI-SpliceAI 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9165884/) ·
[MFASS](https://github.com/KosuriLab/MFASS) ·
[SpliceVarDB](https://pmc.ncbi.nlm.nih.gov/articles/PMC11480807/) ·
[Vex-seq GSE113163](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113163) ·
[SpliceAI-lookup](https://github.com/broadinstitute/SpliceAI-lookup) ·
[GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)
