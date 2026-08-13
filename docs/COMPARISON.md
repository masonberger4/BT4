# How BT4 compares to other codon-optimization tools

An honest positioning of BT4 against the commercial and open-source codon-
optimization landscape. In keeping with BT4's governing principle — *never
present an unenforced constraint, an unvalidated number, or a heuristic result as
if it were real* ([`../CLAUDE.md`](../CLAUDE.md) §1) — this document states where
BT4 genuinely differs, **and where it does not**, and it foregrounds a scientific
caveat that applies to BT4 as much as to every other tool.

> **TL;DR.** BT4's real, near-unique advantage is one of **method and rigor** —
> exact multi-objective optimization, an optimality certificate, byte-level
> reproducibility, content-hashed provenance, and *honestly-labeled* validated
> models. It is **not** yet differentiated on the axis a bench scientist cares
> about most — "will this express better?" — because (a) the CDS is a minority of
> the expression signal, and (b) BT4's splice and expression models still ship
> `calibrated=False`. Both statements are true at once.

---

## 1. These aren't quite the same kind of product

Most "competitors" are **gene-synthesis vendors** whose optimizer is a design
front-end that funnels orders and guarantees the sequence is manufacturable *by
them*. BT4 is a **standalone research optimizer** with no synthesis business
behind it. That difference shapes what "added benefit" even means: a vendor tool
is judged by "can I order this and will it ship," BT4 by "can I trust, defend, and
reproduce this design." They optimize for different things, honestly.

---

## 2. Feature comparison

| Capability | IDT / Twist / GeneArt / GenScript | DNA Chisel (open source) | **BT4** |
|---|---|---|---|
| Optimization method | Greedy or stochastic single-shot | Local search + constraints | **Exact DP + ILP/CP-SAT + Lagrangian** |
| Multiple objectives | Collapsed into one hidden scalar | Weighted scalar | **True Pareto frontier** |
| Optimality guarantee | None | None | **Certificate: proven / gap-bounded / relaxed** |
| Reproducible output | **No — explicitly random** | Partial | **Byte-identical from a seed** |
| Provenance | Undisclosed tables, opaque model | Open code | **Content-hashed manifest (tables + model SHAs + seed + git SHA)** |
| Real folding thermodynamics | Heuristic ΔG windows | Some | **ViennaRNA (calibrated)** |
| Cryptic-splice model | None found in any shipping tool | No | **Wraps SpliceAI / Pangolin (`calibrated=False`)** |
| Learned expression model | ATUM only (proprietary, opaque) | No | **Wraps RiboNN (`calibrated=False`)** |
| Open source | Mostly no | Yes | Yes |

What the vendor tools actually are, from their own documentation:

- **IDT Codon Optimization Tool** — stochastic weighted-random codon sampling
  toward the organism's natural bias, dropping codons used <10% of the time.
  Single, **non-reproducible** output by design ("multiple optimization attempts
  of the same sequence will generate different results"). No learned expression
  model, no thermodynamic folding, no splice/CpG/uORF handling; screening is
  manufacturability-oriented. A 2025 peer-reviewed benchmark found it
  underperformed several free academic tools on codon metrics.
- **Twist codon optimization** — currently marketed as "LLM-based," though Twist's
  own older FAQ describes it as "randomly assigns synonymous codons… does not use
  machine learning," a tension Twist does not reconcile. **Non-deterministic**,
  single output, explicitly disclaims optimality ("optimal expression cannot be
  guaranteed"). Handles GC, homopolymers, restriction sites, 5′ hairpin ΔG,
  internal RBS/promoters; no splice/CpG/uORF modeling exposed. Its only public
  validation is a gated, non-peer-reviewed 32-antibody white paper.
- **GeneArt GeneOptimizer** (Raab et al. 2010) — a greedy **sliding-window**
  heuristic: at each window it scores synonymous combinations with a weighted
  multiparameter quality function, fixes one codon, and slides on. Explicitly not
  a global optimum; single output; proprietary; tied to synthesis.
- **ATUM / DNA2.0 GeneGPS** (Welch et al. 2009) — the one incumbent built on
  **real measured expression** of dozens of synonymous variants. Proprietary,
  opaque, single output, no calibration transparency — but empirically grounded in
  wet-lab data in a way BT4's CDS objectives are not (see §4).
- **DNA Chisel** (Zulkower & Rosser 2020, open source) — constraint-based design
  with weighted objectives via local stochastic/exhaustive search. The closest in
  spirit to BT4, but scalarizes objectives (no Pareto frontier) and gives no
  optimality certificate.

---

## 3. Where BT4 genuinely differs (real and near-unique)

1. **Exact multi-objective optimization with an honesty certificate.** No
   surveyed tool provides this. Everyone else is greedy or stochastic and
   collapses trade-offs into a single hidden number. BT4 returns a **Pareto
   frontier** and states *how optimal* each point is and *what was relaxed*.
2. **Reproducibility and provenance.** Vendor tools are explicitly
   non-deterministic with undisclosed codon tables. BT4 is byte-reproducible from
   a seed and stamps a content-hashed manifest (table contents, model weights,
   seed, git SHA). This matters most for **regulated / therapeutic** work, where a
   sequence must be justified and re-derived, not just ordered.
3. **Validated ML integrated with honest calibration.** The literature says
   folding and splicing matter (§5), yet **no shipping codon tool integrates a
   validated splice model**. BT4 wraps published SpliceAI / Pangolin / ViennaRNA /
   RiboNN — and labels each by calibration status rather than asserting it.

---

## 4. Where BT4 does *not* have an edge (the honest part)

1. **Empirical grounding.** ATUM/GeneGPS fit their models to **measured
   expression** of synonymous variants. BT4's CAI / tAI / codon-pair objectives
   are sequence statistics, not fits to wet-lab yield. BT4's answer is the RiboNN
   wrap — but it is `calibrated=False` and CDS-only, so *today* ATUM can claim
   something BT4 cannot: expression predictions fit to real data.
2. **Manufacturability / synthesis integration.** IDT and Twist guarantee a
   sequence they can actually synthesize, with vendor-grade complexity screening.
   BT4 has partial constraints (max-GC-run, homopolymer, repeats) but no
   synthesis-vendor complexity model.
3. **The paradigm itself isn't novel.** DNA Chisel already does open-source
   constraint-based multi-objective design. BT4's edge over it is **exactness,
   the Pareto frontier, certificates, and reproducibility** — not constraint-based
   optimization as such.

---

## 5. The scientific caveat that applies to *every* codon tool (including BT4)

The peer-reviewed consensus on codon optimization is blunt, and BT4's honesty
ethos requires stating it plainly:

- **CAI barely predicts expression.** Kudla et al. 2009 (*Science*) found codon
  bias **did not correlate** with expression across 154 synonymous GFP variants;
  **5′ mRNA folding near the start dominated**. Welch et al. 2009 (*PLoS ONE*)
  state flatly that "CAI has no value in predicting gene expression" in their data.
  Note that this verdict lands on CAI computed against a *highly-expressed*
  reference set — the strong form of the metric, and the one BT4 now ships as its
  default. Moving from genome-wide counts to a highly-expressed reference makes
  CAI a **better-founded proxy**; it does not make it predictive — in Welch's
  own data the *E. coli* variant built by maximizing exactly it expressed at a
  fraction of the alternatives they tested.
- **Expression is multifactorial.** tRNA adaptation (tAI), codon-pair bias, the 5′
  translational ramp, GC, CpG/UpA content, cryptic splicing, and — above all — 5′
  UTR and initiation-region folding all contribute (Tuller 2010; Coleman 2008;
  Hanson & Coller 2018).
- **The CDS is a minority of the signal.** RiboNN's own analysis attributes
  ~**67% of per-nucleotide signal to the 5′UTR and only ~31% to the CDS**.
- **Optimization can be harmful.** Maximizing yield can disrupt co-translational
  folding, causing misfolding, altered function, or cryptic immunogenic epitopes
  (Mauro & Chappell 2014). Harmonization is often preferable to maximization
  (Mignon et al. 2018).

**BT4 optimizes the CDS.** So even a *perfectly exact* CDS optimizer is optimizing
a minority of the biological signal. This is not unique to BT4 — every codon tool
shares it — but it has two honest consequences:

- BT4's **honesty framing is vindicated** by this literature: it refuses to
  present CAI or an un-gated model as predictive, which is precisely the field's
  documented failure mode.
- BT4's *validated* advantage today is therefore narrower than the feature list
  suggests: **honest constraint enforcement, exact multi-objective trade-offs, and
  real folding (ViennaRNA)**. The splice and expression models — the pieces that
  would let BT4 optimize what actually drives expression — are the ones still
  `calibrated=False` (their promotion path is the §9 acceptance/fidelity gates).

---

## 6. Bottom line

BT4's added benefit is real but specific: it is the only tool that does **exact,
multi-objective, reproducible, certificate-backed** optimization with
**honestly-labeled** validated models — a *rigor-and-trustworthiness* proposition
that the synthesis-vendor incumbents structurally cannot occupy (non-determinism
is fine for them; it is disqualifying for a design you must defend).

BT4 is **not** yet differentiated on "will this express better?" Earning that
claim is the roadmap: **calibrate** the splice and expression models on real,
regime-matched data (§9 gates), and **expand beyond the CDS** toward 5′UTR /
initiation design, where the literature says the signal actually lives.

---

## Sources

Vendor documentation (IDT, Twist, GeneArt, GenScript, ATUM) and the following
peer-reviewed / primary sources, current as of this writing:

- Kudla et al. 2009, *Science* — 5′ folding dominates; codon bias does not
  correlate with per-gene expression. <https://www.science.org/doi/10.1126/science.1170160>
- Welch et al. 2009, *PLoS ONE* — "CAI has no value in predicting gene
  expression" in their data. <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007002>
- Tuller et al. 2010, *PNAS* — 5′ translational ramp. <https://www.pnas.org/doi/10.1073/pnas.0909910107>
- Coleman et al. 2008, *Science* — codon-pair bias. (Cell Reports follow-up: <https://www.cell.com/cell-reports/fulltext/S2211-1247(20)30535-0>)
- Hanson & Coller 2018, *Nat Rev Mol Cell Biol* — codon optimality via elongation
  and mRNA stability. <https://www.nature.com/articles/nrm.2017.91>
- Mauro & Chappell 2014, *Trends Mol Med* — a critical analysis of codon
  optimization in human therapeutics. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4253638/>
- Mignon et al. 2018, *FEBS Letters* — codon harmonization vs maximization. <https://febs.onlinelibrary.wiley.com/doi/10.1002/1873-3468.13046>
- Raab et al. 2010 — the GeneOptimizer sliding-window algorithm. <https://link.springer.com/article/10.1007/s11693-010-9062-3>
- Zulkower & Rosser 2020, *Bioinformatics* — DNA Chisel. <https://academic.oup.com/bioinformatics/article/36/16/4508/5869515>
- RiboNN (Zheng et al. 2024) — translation-efficiency CNN; ~31% of per-nt signal
  in the CDS. <https://pubmed.ncbi.nlm.nih.gov/39149337/>
- Independent multi-tool benchmark (2025). <https://pmc.ncbi.nlm.nih.gov/articles/PMC12010093/>

> This document is a positioning and research-integrity statement, not a
> marketing claim. If a comparison here drifts from what BT4 actually enforces or
> from the cited literature, the code and [`../CLAUDE.md`](../CLAUDE.md) are
> authoritative — fix this file in the same change (§10.11).
