# Ranaghan et al. 2021 — Table 4 reproducibility panel (`ranaghan2021_tab4.fasta`)

## Source

Ranaghan, M. J., Li, J. J., Laprise, D. M., & Garvie, C. W. (2021).
**Assessing optimal: inequalities in codon optimization algorithms.**
*BMC Biology*, **19**:36. DOI: [10.1186/s12915-021-00968-8](https://doi.org/10.1186/s12915-021-00968-8)

Sequences are from the paper's **Additional file 1 (MOESM1), Table 4** — the
repeatability experiment in which three codon-optimization algorithms each
optimized the same three human genes ten independent times.

## License

Published under the **Creative Commons Attribution 4.0 International License
(CC BY 4.0)** — <https://creativecommons.org/licenses/by/4.0/>. Reuse is
permitted with attribution. Attribution: *Ranaghan, Li, Laprise & Garvie, BMC
Biology 19:36 (2021), CC BY 4.0.* **Changes made:** none to sequence content;
records were reformatted into FASTA and the header fields normalized.

## Contents (93 records)

Three human proteins — **KRas4B** (P01116, 188 aa), **Beclin 1** (Q14457,
450 aa), **PDE3A** (Q14432, 1141 aa) — each with:

* one **Native** coding sequence (the reference), and
* **Algorithm1 / Algorithm2 / Algorithm3**, each run **ten** times
  (`run1`..`run10`).

Header format: `Protein|Source|[runN]|acc=...|uniprot=...|len=...`.

## Honesty notes

* **The tools are anonymized.** The paper reports these three algorithms as
  *Algorithm 1/2/3*, not named vendors. Do **not** map them to specific
  commercial tools — that mapping is not published. This panel is therefore a
  **run-to-run reproducibility / variability** comparison, a fundamentally
  different thing from the named-tool head-to-head in `kras_ranaghan2021.fasta`
  (which does attribute GeneArt/IDT/Twist/GenScript). Keep the two separate.
* The **ten runs of one algorithm are repeat runs of a single stochastic tool**,
  so their spread measures that tool's determinism, not a difference between
  tools. Ten runs ≠ ten tools.
* Every metric BT4 reports over this panel is **recomputed by BT4's own
  functions** from each delivered nucleotide sequence — nothing is copied from
  the paper, and BT4 is never claimed "better."
