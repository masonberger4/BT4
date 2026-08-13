# State-of-the-art codon & mRNA sequence optimization — a grounded survey

**Purpose.** A cited map of where codon / mRNA-coding-sequence design actually is
in the biotech industry and literature as of 2026, written to steer BT4's
roadmap. It is deliberately honest about what is *validated* versus *promising*,
in keeping with the project's governing principle ([`../CLAUDE.md`](../CLAUDE.md)
§1): never present an unvalidated number as if it were real. Where a claim is a
headline result from one group, it is attributed as such, not stated as settled
fact.

This document is **research and gap analysis, not a changelog** — nothing here is
claimed to ship in BT4 unless the "BT4 today" column says so. Decisions that
follow from it belong in [`NEXT_SESSION.md`](NEXT_SESSION.md); durable rules that
follow from it belong in `../CLAUDE.md`.

---

## 0. The one-paragraph takeaway

The field has moved decisively off **CAI-maximization as the objective**. Two
things replaced it, and they are complementary rather than competing:

1. **Multi-objective, structure-aware exact optimization.** The landmark is
   **LinearDesign** (Zhang et al., *Nature* 2023), which jointly optimizes codon
   usage *and* whole-mRNA secondary-structure stability (minimum free energy) in
   a single exact dynamic program over a codon lattice, and reported **up to 128×
   antibody titre, ~5× mRNA half-life and ~3× protein** for a COVID vaccine
   antigen in mice versus a codon-optimized benchmark. This is the single most
   important result for BT4, because BT4's exact codon trellis is *the same class
   of algorithm* — BT4 already has the two ingredients (a codon trellis and a
   ViennaRNA model) but keeps folding in a post-hoc refinement layer rather than
   jointly in the DP.
2. **Learned, expression-validated models.** Transformer / language-model designers
   (CodonTransformer, CodonBERT, DeepCodon, mRNA-GPT, mRNABERT, codonGPT, RiboNN)
   trained on large sequence and ribosome-profiling corpora. These are the
   research frontier for *predicting* expression, but the honest caveat the field
   itself repeats is that **the CDS explains a minority of expression variance**
   (RiboNN's own ablation: ~31% CDS vs ~67% 5′/3′UTR) and that **public codon
   tools do not reliably raise expression and can lower it** (misfolding, over-
   optimization).

BT4's existing thesis — *a Pareto frontier over calibrated, locality-partitioned
terms, honest about optimality* — is exactly where the multi-criteria literature
says design should go (§6). The biggest concrete gap is **joint codon + structure
optimization** (§4), followed by an **mRNA-therapeutic design mode** (§5).

---

## 1. The paradigm shift: from CAI to a validated objective vector

Every recent comparative study reaches the same verdict: single-metric CAI
optimization is inadequate, and often actively harmful.

- **CAI barely predicts expression.** Kudla et al. (*Science* 2009) found codon
  bias did not correlate with GFP expression across 154 synonymous variants; 5′
  mRNA folding dominated. Welch et al. (*PLoS ONE* 2009) stated flatly that "CAI
  has no value in predicting gene expression" in their data — the same *E. coli*
  variant built by maximizing highly-expressed CAI expressed at a fraction of the
  alternatives. (BT4 already foregrounds this in [`COMPARISON.md`](COMPARISON.md).)
- **Public tools can reduce expression.** The 2025 multi-criteria comparative
  analysis (Sci Rep / PMC12010093) tested ten tools (JCat, OPTIMIZER, ATGme,
  TISIGNER, GenSmart, ExpOptimizer, IDT, Genewiz, GeneOptimizer, VectorBuilder)
  and found they "frequently produced divergent results," that CAI-only tools
  "failed to account for mRNA secondary structure stability and codon-pair bias,"
  and that **there is no universally best tool** — the right choice depends on
  host and target. Its recommended best practice is "a multi-criteria framework
  that integrates CAI, GC content, mRNA folding energy, and codon-pair
  considerations," host-tailored.
- **Over-optimization harms folding.** Rare-codon clusters can be functionally
  required for co-translational folding; erasing them causes misfolding
  (§7). This is why DeepCodon (2025) is explicitly built to *preserve* conserved
  rare-codon clusters rather than maximize CAI.

**Implication for BT4.** BT4's multi-objective Pareto core is the correct
architecture for this world; the work is in *which* axes it carries and which are
*calibrated*, not in the paradigm.

---

## 2. The design axes the field actually uses

Each row: what it is, the key evidence, and BT4's current status. "BT4 today"
reflects `main` at the time of writing.

| Axis | What it captures | Key evidence | BT4 today |
|---|---|---|---|
| **Codon adaptation (CAI)** | Match to a highly-expressed reference set | Sharp & Li 1987; weak expression proxy (Kudla, Welch) | **Have** — highly-expressed reference sets now default (§8 of CLAUDE.md) |
| **tRNA adaptation (tAI)** | Supply/demand of charged tRNAs; elongation speed | dos Reis 2004 | **Have** — real GtRNAdb tables, 8 organisms |
| **Codon optimality / CSC** | Codon identity → **mRNA stability** via elongation-coupled decay | Presnyak & Coller, *Cell* 2015 ("Codon optimality is a major determinant of mRNA stability"); CSC metric (NAR 2019) | **Not directly** — tAI is a related but distinct proxy; no CSC term |
| **Codon-pair bias (CPB)** | Over/under-represented adjacent codon pairs; used for viral attenuation (SAVE) | Coleman et al. 2008 | **Have** — `CpbTerm`, PAIRWISE, exact in trellis |
| **mRNA secondary structure (MFE)** | Whole-transcript folding stability → half-life | **LinearDesign, *Nature* 2023** (joint CAI+MFE, 128× antibody); more structure ⇒ longer half-life | **Partial** — ViennaRNA in a *refinement* layer, not jointly in the DP; **the flagship gap (§4)** |
| **5′ start-region structure** | Structure near the cap/start **inhibits** 40S scanning & initiation | Kudla 2009; 40S 1-D scanning work (2024–25) | **Partial** — 5′-folding refinement targets this |
| **GC content / GC3** | Folds into stability, manufacturability, host bias | Host-specific (E. coli high-GC, yeast A/T-rich) | **Have** — GC target + window + max-GC-run |
| **CpG content / ZAP** | CpG-rich mRNA degraded by zinc-finger antiviral protein; depletion evades, elevation is immunostimulatory | Takata et al. *Nature* 2017; CpG/ZAP HIV work | **Have** — CpG deplete/elevate term + budget |
| **Uridine content (U-depletion)** | Lowers reactogenicity / innate sensing in mRNA therapeutics | mRNArchitect, mRNAid design principles | **Missing** — no uridine-depletion term |
| **m1Ψ slippery sequences** | N1-methylpseudouridine causes **+1 ribosomal frameshifting** at slippery motifs, producing off-target products immunogenic in humans | Mulroney et al. *Nature* 2024; m1Ψ decoding (Nat Commun 2024) | **Missing** — no slippery-sequence avoidance |
| **Cryptic splicing** | Introduced donor/acceptor sites mis-splice the transgene | SpliceAI (Illumina), Pangolin | **Have (uncalibrated)** — wrapped CNNs + PWM baseline + ASSP cross-check |
| **Co-translational folding / harmonization** | Match host *translation-rate profile* (not maximize) to preserve folding | Angov 2008/2011 (harmonization); Pechmann & Frydman | **Partial** — 5′ ramp + %MinMax terms, but no harmonization mode |
| **Learned expression / TE** | Directly predict translation efficiency from sequence | RiboNN; cell-type-specific model (NAR gkaf233, 2025) | **Have (uncalibrated)** — RiboNN wrapper, `calibrated=False` |

---

## 3. The two methodological camps

BT4 sits squarely in the first and should own it while wrapping the second.

### 3a. Exact / combinatorial optimization — **BT4's lineage**

- **LinearDesign** (Zhang et al., *Nature* 2023, Baidu / Oregon State / StemiRNA /
  Rochester). Couples a minimum-free-energy folding DP with a codon-level lattice
  graph — each layer is a codon position, each path a synonymous mRNA — and finds
  the *exact* optimum trading CAI against MFE, in ~11 minutes for the spike
  protein. **This is architecturally the same idea as BT4's codon trellis**, and
  its result (§0) is the strongest single argument in the field for structure-
  aware exact design.
- **2026 successors** show the lineage is active: the *Montparnasse Algorithm*
  (Cazenave, arXiv 2606.07562, 2026) and *tensor-based secondary-structure*
  RNA-design models (arXiv 2604.19718, 2026) push exact/near-exact joint codon +
  structure design further, and there is even a quantum-computing formulation of
  the co-optimization (arXiv 2507.18817). BT4's honest-optimality-certificate
  framing is a natural fit here.

**Why this matters for BT4:** BT4 already returns a Pareto frontier with an
optimality certificate. Adding a folding-aware trellis state (the LinearDesign
lattice) would let BT4 do *jointly and exactly* what it currently does in two
stages (exact codon DP → SA refinement over folding), and would put BT4 on the
same axis as the strongest validated result in the field.

### 3b. Learned / generative — **wrap, don't reinvent (BT4's existing stance)**

A fast-moving 2024–2026 wave of transformer / language-model designers:

- **CodonTransformer** — multispecies, 164 organisms, natural-like codon
  distributions.
- **CodonBERT** — BERT trained on highly-expressed human transcripts for vaccine
  codon choice.
- **DeepCodon** (2025) — preserves conserved rare-codon clusters; E. coli.
- **mRNA-GPT / mRNABERT / codonGPT** (2025) — GPT/BERT-scale generative mRNA
  designers; codonGPT uses reinforcement learning to satisfy biological
  constraints. GEMORNA reported up to **15× EPO** expression.
- **RiboNN** — ribosome-load / translation-efficiency regression; the model BT4
  already wraps (`calibrated=False`).
- **Cell-type-specific codon preferences** — a deep model on tissue-resolved
  expressed transcripts (NAR gkaf233, 2025) finds genuine cell-type codon-
  optimization signal. *(Relevant to the maintainer's decision to scope tissue-
  specific tables out — the evidence exists; the decision was cost/qualifiability,
  not absence of signal. Recorded so the decision stays informed.)*

**Honest read.** These predict/generate well in-distribution but (a) are usually
not validated for the *specific* CDS-variant regime a codon optimizer operates
in, and (b) inherit the CDS-is-a-minority-of-signal ceiling. BT4's posture —
wrap the user's own installed model, hash-pin weights, ship `calibrated=False`
until a regime-matched gate passes — is exactly right and matches how the field's
more careful practitioners talk about their own models.

---

## 4. The flagship gap: joint codon + structure optimization

**Recommendation (highest impact, architecturally native).** Add a
**folding-aware objective that the trellis optimizes jointly**, LinearDesign-style,
rather than only as SA refinement over an exact-DP seed.

- **Why it's the top priority.** It is the single design change with the strongest
  *validated in-vivo* evidence behind it (§0), and BT4 is unusually well-placed to
  do it honestly: it already has the codon trellis, a ViennaRNA `FoldingModel`, and
  an optimality-certificate framework the exact-design camp lacks a clean version of.
- **Structure is context-dependent — this is the subtlety BT4 can get right.** The
  literature is not "more structure is better" uniformly:
  - **Body of the transcript:** more secondary structure ⇒ longer half-life ⇒ more
    total protein (LinearDesign's core finding).
  - **5′ start region:** structure near the cap/start **inhibits** 40S scanning and
    initiation (Kudla 2009; scanning studies). The right objective *minimizes*
    start-proximal structure while *allowing/maximizing* structure in the body.
  BT4's existing split (a 5′-folding refinement term separate from whole-sequence
  concerns) already gestures at this; a joint formulation should make the
  position-dependence explicit.
- **Honesty constraints that must survive.** ViennaRNA MFE is calibrated
  thermodynamics (keep it labeled so). A joint objective must still emit an honest
  certificate — LinearDesign is exact for its scoring, but if BT4 approximates
  (beam, windowed folding) it must say so, per §5/§10.1.

This is Phase-3's "non-local models & refinement done right" taken to its
conclusion, and it is the recommendation most worth a design doc of its own.

---

## 5. A second gap: an mRNA-therapeutic design mode

The mRNA-medicine tools (**mRNArchitect**, Garvan 2024; **mRNAid**, 2024) converge
on a design profile distinct from "express a transgene in a cell line":

- **Uridine depletion** — count U at codon third positions and minimize it, to cut
  innate reactogenicity. **BT4 lacks this**; it is a clean additive per-codon term
  (like CpG), exact in the trellis.
- **Avoid m1Ψ slippery sequences** — for N1-methylpseudouridine-modified
  therapeutic mRNA, slippery motifs cause +1 frameshifting and off-target,
  immunogenic products; synonymous disruption of slippery sequences mitigates it
  (Mulroney et al. *Nature* 2024). **BT4 lacks this**; it is a LOCAL motif
  constraint, exactly the shape BT4's constraint engine already handles.
- **Structure context (§4)** and **CpG/U handling** are application-switched:
  vaccines may *want* structure and CpG (adjuvant effect); frequently-dosed
  therapeutics avoid both. This argues for a small number of **named design
  presets** (vaccine / therapeutic / cell-line expression) layered over the
  existing knobs, not new machinery.

None of these require abandoning BT4's exact-DP core; they are additional
LOCAL/additive terms and constraints plus a preset layer.

---

## 6. What the research says BT4 already gets right

Worth stating, because it means the direction is sound and the work is additive:

- **Multi-objective over single-metric.** The exact recommendation of the 2025
  multi-criteria framework.
- **Honest optimality.** No surveyed commercial tool (GenScript GenSmart's "200+
  factors / population immune algorithm," IDT, GeneArt) exposes an optimality
  certificate or is reproducible from a stamp; most are explicitly stochastic.
- **Calibrated-vs-labeled ML.** BT4's `calibrated=False`-until-gated posture is the
  responsible version of what the ML-designer papers should (but often don't) say.
- **CpG as deplete/elevate with a biological rationale** (ZAP) rather than a magic
  knob.
- **Codon-pair bias and %MinMax / ramp** are already the "context and profile"
  terms the harmonization literature cares about.

---

## 7. The honest ceiling (applies to BT4 and every competitor)

- **The CDS is a minority of the expression signal.** RiboNN attributes ~67% of
  per-nucleotide signal to UTRs, ~31% to the CDS. A *perfect* CDS optimizer still
  optimizes a minority of the biology. UTR design, cap, poly(A), and delivery
  dominate and are out of BT4's current scope.
- **Optimization can harm.** Maximizing yield can disrupt co-translational folding
  (misfolding, altered function, cryptic epitopes). Harmonization is often
  preferable to maximization (Mignon et al. 2018). This is the argument for a
  harmonization *mode*, not just a maximization one.
- **No tool is universally best** (the multi-criteria study's headline). Host and
  target decide. BT4's per-organism, per-reference-set, per-objective transparency
  is the right response to this — it makes the choice the user's, explicitly.
- **Validation is the currency.** The field's credible claims (LinearDesign's 128×,
  GEMORNA's 15×) are *wet-lab in-vivo* results. BT4's differentiator is method and
  honesty, not a validated expression claim — and it must keep saying so until it
  has its own gated data.

---

## 8. Ranked recommendations for BT4

1. **Joint codon + secondary-structure optimization (LinearDesign-class), position-
   aware** (§4). Highest validated impact; native to the trellis. Deserves its own
   design doc.
2. **mRNA-therapeutic terms** (§5): uridine-depletion objective + m1Ψ slippery-
   sequence constraint, both exact/LOCAL and cheap.
3. **Codon optimality / CSC** as a first-class term distinct from tAI (§2) —
   mRNA-stability-through-codon-choice is a validated axis BT4 doesn't yet carry.
4. **Named application presets** (vaccine / therapeutic / expression) that set the
   structure/CpG/U trade-offs coherently, rather than leaving them as unguided
   knobs (§5).
5. **Codon harmonization mode** — a "match the host profile" objective alongside
   maximization, for co-translational folding (§2, §7).
6. **Keep wrapping learned models, keep them `calibrated=False`** until a
   regime-matched gate passes (§3b). The E. coli tRNA table (already queued) also
   unlocks tAI where translational selection is strongest.

Items 2–5 are additive terms/constraints/presets over the existing exact core.
Item 1 is the one architectural investment, and the one the evidence most
strongly rewards.

---

## Sources

Primary results:
- Zhang et al., "Algorithm for optimized mRNA design improves stability and
  immunogenicity," *Nature* 621:396 (2023). doi:10.1038/s41586-023-06127-z
  (LinearDesign; up to 128× antibody, ~5× half-life, ~3× protein).
- Presnyak et al., "Codon optimality is a major determinant of mRNA stability,"
  *Cell* 160:1111 (2015).
- Kudla et al., "Coding-sequence determinants of gene expression in
  *Escherichia coli*," *Science* 324:255 (2009).
- Welch et al., "Design parameters to control synthetic gene expression in
  *Escherichia coli*," *PLoS ONE* 4:e7002 (2009). doi:10.1371/journal.pone.0007002
- Mulroney et al., "N1-methylpseudouridylation of mRNA causes +1 ribosomal
  frameshifting," *Nature* 625:189 (2024). doi:10.1038/s41586-023-06800-3
- Takata et al., "CG dinucleotide suppression enables antiviral defence targeting
  non-self RNA," *Nature* 550:124 (2017).
- Coleman et al., "Virus attenuation by genome-scale changes in codon pair bias,"
  *Science* 320:1784 (2008).
- dos Reis et al., "Solving the riddle of codon usage preferences," *NAR* 32:5036
  (2004).

Tools, models, and reviews:
- "Comparative Analysis of Codon Optimization Tools: Advancing toward a
  Multi-Criteria Framework," *Sci Rep* (2025). PMC12010093.
- "mRNArchitect: sequence design of mRNA medicines," bioRxiv 2024 (Garvan
  Institute). doi:10.1101/2024.12.03.626696 *(first author not verified here;
  cited by tool + venue)*
- "mRNAid, an open-source platform for therapeutic mRNA design," *NAR Genom
  Bioinform* 6:lqae028 (2024).
- Learned designers (cited by tool + venue; first-author surnames not
  independently verified): **CodonBERT**; **CodonTransformer** (multispecies,
  164 organisms); **DeepCodon** (rare-codon-cluster preserving, 2025,
  PMC13109293); **mRNA-GPT**; **mRNABERT** (*Nat Commun* 2025,
  s41467-025-65340-8); **codonGPT** (reinforcement learning, bioRxiv 2025).
- "Deep generative optimization of mRNA codon sequences," *Nat Commun* (2025).
  s41467-025-64894-x (PubMed 40875799).
- Cell-type codon-optimization deep model, *NAR* 53:gkaf233 (2025).
- Exact-design frontier: Cazenave, "The Montparnasse Algorithm for RNA Design,"
  arXiv 2606.07562 (2026); tensor-based RNA design under codon constraints,
  arXiv 2604.19718 (2026); quantum co-optimization, arXiv 2507.18817 (2025).
- Angov et al., codon harmonization, *PLoS ONE* 3:e2189 (2008) and *Methods Mol
  Biol* 2011.
- GenScript GenSmart product documentation (200+ factors, "population immune
  algorithm"); accessed 2026.

*Compiled 2026-08. Claims attributed to a single group are marked as such; treat
in-vivo fold-improvements as that group's reported results, not settled
constants.*
