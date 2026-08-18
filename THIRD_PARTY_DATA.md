# Third-party data & model provenance

BT4's own code is MIT-licensed (see [`LICENSE`](LICENSE)). This file records the
provenance and licensing status of the **data bundled in the repository** and the
**external models BT4 wraps but does not bundle**, so the licensing picture is
visible without opening each sidecar. Every bundled dataset ships a machine-readable
`*.provenance.json` / `*.LICENSE.md` sidecar next to it; this file summarizes them.

## Bundled data (in the tree)

| Data | Location | Source | License / status |
|---|---|---|---|
| tRNA gene copy numbers (8 organisms) | `src/bt4/biomodels/codon/data/*.trna.tsv` (+ `*.trna.provenance.json`) | GtRNAdb (Chan & Lowe 2016; tRNAscan-SE 2.0), independently re-counted from a SHA-256-pinned FASTA | **Citation-gated academic use — NOT a CC/public-domain grant.** GtRNAdb states no explicit data license. These are raw factual gene counts by anticodon (generally not copyrightable), honestly provenanced. |
| Codon-usage tables — **genome-wide** (9 organisms) | `src/bt4/biomodels/codon/data/<organism>.tsv` (+ `*.provenance.json`) | Real codon counts over one representative CDS per gene, recounted from **release-pinned Ensembl** CDS FASTAs by `scripts/build_organism_tables.py` | **Ensembl** data, made freely available by EMBL-EBI ([terms](https://www.ebi.ac.uk/about/terms-of-use)); cite Harrison et al., *NAR* 2024, doi:10.1093/nar/gkad1049. *A. thaliana* additionally requires the TAIR10 assembly citation (Lamesch et al. 2012) and its Araport11 annotation (Cheng et al. 2017). Each sidecar carries the source URL, the source file's own SHA-256, and the per-filter drop tally; `--verify` rebuilds and diffs the shipped bytes. |
| Codon-usage tables — **highly-expressed** (8 organisms) | `src/bt4/biomodels/codon/data/<organism>.highly_expressed.tsv` (+ `*.provenance.json`) | Codon counts over each organism's 300 most abundant proteins, ranked by **PaxDb v6.1** whole-organism integrated proteomics and joined to the same pinned Ensembl CDS through that release's peptide FASTA, by `scripts/build_highly_expressed_tables.py` | **PaxDb: CC BY 4.0** (von Mering Lab, SIB / University of Zurich) — cite Huang et al., *Mol Cell Proteomics* 2023, doi:10.1016/j.mcpro.2023.100640. Sequences and the protein→gene mapping are **Ensembl** (same terms and citation as the row above). All three sources are SHA-256-pinned in every sidecar. |
| KRas4B benchmark FASTA | `scripts/data/kras_ranaghan2021.fasta` (+ `.LICENSE.md`) | Ranaghan et al. 2021 | **CC BY 4.0** (cited; changes documented; ATUM/DNA2.0 truncation flagged). |
| Reproducibility panel FASTA | `scripts/data/ranaghan2021_tab4.fasta` (+ `.LICENSE.md`) | Ranaghan et al. 2021, Table 4 | **CC BY 4.0** (cited). |
| RiboNN weight manifest | `src/bt4/biomodels/expression/data/ribonn_sha256.json` | Public SHA-256 **content hashes only** (90 human + 90 mouse) | Hashes, not weights — no model bytes bundled. |
| RiboNN gate attestations | (none committed yet) | Derived scalars + public content hashes only, never raw model outputs | An `ExpressionAttestation` records a passing gate license-cleanly; RiboNN's per-sequence scores are non-commercial outputs and are structurally excluded (`_ALLOWED_FIELDS`). |
| ASSP cross-check fixtures | `tests/fixtures/assp/*.txt` | **Synthetic** ASSP-format reports (not real service captures) | BT4-authored test fixtures. |

## External models — wrapped, never bundled

BT4 lazily imports the user's **own** installed package and weights (pointed at via
env vars such as `$BT4_RIBONN_DIR`), hash-verifies the weights against a pinned
SHA-256 **before** loading, and never vendors or reimplements them. No model weights
are in this repository or its git history.

| Model | License (upstream) | How BT4 uses it |
|---|---|---|
| RiboNN (expression) | Sanofi **non-commercial** | Optional `bt4[expression-ribonn]`; user-supplied checkout + weights, hash-pinned. |
| SpliceAI | code PolyForm Strict 1.0.0 · weights **CC BY-NC 4.0** (noncommercial) | Optional `bt4[splice-spliceai]`; user-supplied install. |
| Pangolin | **GPL-3.0** | Optional `bt4[splice-pangolin]`; user-supplied install. |
| ViennaRNA | GPL | Optional `bt4[fold]`; user-supplied install. |
| ASSP (online) | Web service (Alternative Splice Site Predictor) | Opt-in `bt4[assp]` network cross-check; results labeled network-derived, never bundled or reproducible-from-manifest. |

Because BT4 is open-source and non-commercial in scope, the non-commercial SpliceAI
and Sanofi RiboNN terms are compatible with the wrapped, user-supplied-weights
approach; the GPL components (Pangolin, ViennaRNA) are handled the license-clean way
BT4 handles all of them — imported from the user's install, never linked or vendored.
