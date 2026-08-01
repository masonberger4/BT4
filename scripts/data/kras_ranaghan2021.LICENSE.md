# Provenance and license: `kras_ranaghan2021.fasta`

The ten KRas4B coding sequences in `kras_ranaghan2021.fasta` are reproduced from
the published, openly-licensed dataset of Ranaghan et al. (2021). They are used
here purely as a fixed, real-world panel of commercial and academic
codon-optimized sequences to compare against BT4's own output.

## Source

> Matthew J. Ranaghan, Jeffrey J. Li, Dylan M. Laprise, and Colin W. Garvie.
> **"Assessing optimal: inequalities in codon optimization algorithms."**
> *BMC Biology* **19**, 36 (2021).
> DOI: [10.1186/s12915-021-00968-8](https://doi.org/10.1186/s12915-021-00968-8)
> PMC: [PMC7893858](https://pmc.ncbi.nlm.nih.gov/articles/PMC7893858/)

The sequences are the KRas4B case-study coding sequences from that paper: the
native human KRas4B CDS plus one codon-optimized variant per tool the authors
evaluated.

## License

This dataset is licensed under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

> This article is licensed under a Creative Commons Attribution 4.0 International
> License, which permits use, sharing, adaptation, distribution and reproduction
> in any medium or format, as long as appropriate credit is given to the
> original author(s) and source, with a link to the license and notation of any
> changes made.

**Changes made:** none to the sequence content. The nucleotide sequences were
extracted verbatim and reformatted into a plain multi-record FASTA file with the
tool names as record headers. No bases were added, removed, or edited.

## Tool-name mapping (FASTA header -> optimizer)

The FASTA headers use the short labels from the paper. Several correspond to
named commercial algorithms; the mapping to the vendor/product (including the
CLAUDE.md section 7 targets GeneOptimizer / IDT / Twist) is:

| FASTA header | Optimizer / vendor                                   | Kind       |
| ------------ | ---------------------------------------------------- | ---------- |
| `Native`     | Native human *KRAS* CDS (reference, not optimized)   | reference  |
| `GeneArt`    | Thermo Fisher GeneArt - **GeneOptimizer** algorithm  | commercial |
| `GeneWiz`    | GENEWIZ (Azenta Life Sciences)                       | commercial |
| `DNA2.0`     | **ATUM** (formerly DNA2.0)                            | commercial |
| `IDT`        | Integrated DNA Technologies (Codon Optimization Tool) | commercial |
| `Genscript`  | GenScript (OptimumGene)                              | commercial |
| `Twist`      | Twist Bioscience                                     | commercial |
| `JCAT`       | Java Codon Adaptation Tool                           | academic   |
| `OPTIMIZER`  | OPTIMIZER                                            | academic   |
| `COOL`       | Codon Optimization OnLine (COOL)                     | academic   |

`GeneArt == Thermo GeneOptimizer` and `DNA2.0 == ATUM` in particular: when
CLAUDE.md section 7 names "GeneOptimizer", that is the `GeneArt` row here.

## Important caveat: the `DNA2.0` sequence encodes a different protein

The `DNA2.0` (ATUM) sequence is **not** a synonymous optimization of the same
protein as the other rows. It is a genuine **C-terminal truncation**: it encodes
only 169 amino acids, whereas the native KRas4B and every other panel sequence
encode the full 188-amino-acid protein. It is therefore a *different protein*,
not a codon variant of the reference.

Any comparison tooling that uses this panel **must flag `DNA2.0` as a length /
protein mismatch** and must not treat its recomputed metrics as directly
comparable to the full-length rows. Per-codon geometric-mean indices (CAI, tAI)
and compositional metrics (GC%, CpG, homopolymer run) are still well-defined on
the truncated sequence, but they describe a shorter, different coding sequence -
so they are reported with the mismatch flagged, never silently pooled with the
full-length variants. See `scripts/compare_tools.py`.
