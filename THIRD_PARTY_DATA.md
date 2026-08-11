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
| Codon-usage tables | `src/bt4/biomodels/codon/data/*.tsv` (+ `*.provenance.json`) | Kazusa-style **representative published** per-thousand frequencies | Representative published values, not an authoritative per-genome recount; use ratios within a synonymous group only. |
| KRas4B benchmark FASTA | `scripts/data/kras_ranaghan2021.fasta` (+ `.LICENSE.md`) | Ranaghan et al. 2021 | **CC BY 4.0** (cited; changes documented; ATUM/DNA2.0 truncation flagged). |
| Reproducibility panel FASTA | `scripts/data/ranaghan2021_tab4.fasta` (+ `.LICENSE.md`) | Ranaghan et al. 2021, Table 4 | **CC BY 4.0** (cited). |
| RiboNN weight manifest | `src/bt4/biomodels/expression/data/ribonn_sha256.json` | Public SHA-256 **content hashes only** (90 human + 90 mouse) | Hashes, not weights — no model bytes bundled. |
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
