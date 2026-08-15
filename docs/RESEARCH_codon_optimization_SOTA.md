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
   a single exact dynamic program over a codon lattice. It is architecturally the
   same class of algorithm as BT4's trellis. **Its headline number is routinely
   misquoted, including by an earlier draft of this document** — the 57–128× is
   *anti-spike IgG in n = 6 mice* against a **vendor codon optimizer's** output;
   the HEK293 protein figure was **2.9×**, and the durable claim is **stability**
   (~5–8× in-solution half-life), not expression. It has not been independently
   replicated, and the counter-evidence is now substantial (§4).
2. **Learned, expression-validated models.** Transformer / language-model designers
   (CodonTransformer, CodonBERT, DeepCodon, mRNA-GPT, mRNABERT, codonGPT, RiboNN)
   trained on large sequence and ribosome-profiling corpora. These are the
   research frontier for *predicting* expression, and the honest caveat the field
   repeats is that **public codon tools do not reliably raise expression and can
   lower it** (misfolding, over-optimization).

**A third finding, added after a 2026 verification sweep, outranks both:
context.** Every documented failure mode that has reached patients — cryptic
splicing in a transgene cassette, an internal poly(A) that destroys a lentiviral
genome, a uORF that overlaps into the CDS — is a *junction* phenomenon. The CDS
supplies one half and the vector supplies the other. The 2026 statement of
direction is Shi et al., *J Adv Res*: the field is moving "from isolated
optimization of sub-regions to **coordinated design strategies that account for
cross-regional dependencies**." A CDS-only optimizer, however exact, is
structurally blind to that class.

BT4's existing thesis — *a Pareto frontier over calibrated, locality-partitioned
terms, honest about optimality* — is exactly where the multi-criteria literature
says design should go (§6). The biggest concrete gap is **construct context**
(§3c), not joint codon+structure optimization; see §4 for why that ranking
changed.

> ### Two corrections this document previously got wrong
>
> **1. "The CDS is a minority of the expression signal" is a misreading of
> RiboNN — quote both numbers or neither.** RiboNN reports two different
> attributions and they point opposite ways:
>
> | | 5′UTR | CDS | 3′UTR |
> |---|---|---|---|
> | **Per-nucleotide** information density | **67%** | 31% | 2% |
> | **Length-integrated total** attribution (human) | 22% | **73%** | 5% |
>
> *(mouse: 23 / 73 / 4.)* Per nucleotide the 5′UTR is far denser; integrated over
> length **the CDS carries ~73% of the total attributed translation-efficiency
> signal**. Quoting 67/31/2 as a variance decomposition is the single most common
> misreading of the paper — and it would argue against BT4's own existence. The
> real ceiling on a CDS optimizer is a *different* fact: **mRNA abundance, not
> translation rate, is the majority channel for protein abundance** (Li, Bickel &
> Biggin, *PeerJ* 2:e270 (2014): mRNA levels explain ≥56%, and measured
> ribosome-profiling translation-rate variance is only 12% of what had been
> inferred), and **integration site alone spans ~1,000×** (Akhtar et al., *Cell*
> 154:914 (2013), >27,000 random reporter integrations). State *that* as the
> ceiling, not a mis-halved RiboNN number.
>
> **2. LinearDesign's headline number is not an expression result, and it has not
> replicated.** See §4.

---

## 1. The paradigm shift: from CAI to a validated objective vector

Every recent comparative study reaches the same verdict: single-metric CAI
optimization is inadequate, and often actively harmful.

- **CAI barely predicts expression, with numbers.** Kudla et al. (*Science*
  324:255, 2009): 154 synonymous GFP variants spanning a **250-fold** protein
  range gave **CAI r = 0.14, p = 0.09 (not significant)**, while folding ΔG over
  nucleotides **−4 to +37** explained **44%** of the variance (r = 0.66) under T7
  and **59%** (r = 0.77) with a bacterial promoter. Note that window: it **spans
  the UTR↔CDS junction**, and a CDS-only tool cannot compute it.
- **The highly-expressed reference premise is falsified in *E. coli*.** Welch et
  al. (*PLoS ONE* 4:e7002, 2009), 21 Φ29 polymerase + 24 scFv variants over a
  >40-fold range: *"CAI has no value in predicting gene expression for either gene
  set."* The favourable codons turned out to be those read by the tRNAs **most
  highly charged under amino-acid starvation — explicitly *not* the codons
  abundant in highly expressed *E. coli* genes.** This matters to BT4 directly:
  declaring a reference set fixes *which genes were counted*; it does not fix *the
  fact that frequency is the wrong signal*. Highly-expressed remains the better
  default (it is better-founded than genome-wide), but no tooltip may imply it
  makes CAI predictive.
- **Human synonymous-site selection is not translational selection at all.**
  Radrizzani, Kudla, Izsvák & Hurst, *Nat Rev Genet* 25:431 (2024) — the
  "unwanted transcript hypothesis": selection for abundant-tRNA codons *"is not
  readily observed in humans."* What **is** selected is avoidance of spurious
  transcripts, mis-splicing, cryptic splice sites and TE/virus-like RNAs, plus a
  high-GC / low-CpG "self" signature. **This is the single most consequential
  paper for a CDS optimizer's objective design** — and it is good news for BT4's
  direction rather than bad, because it demotes CAI while promoting exactly the
  things BT4 already ships to first class: cryptic-splice avoidance, CpG budgets,
  and GC control.
- **The 5′ "rare codon ramp" is a structural effect, not a codon-rarity one.**
  Goodman, Church & Kosuri (*Science* 342:475, 2013), >14,000 reporters, ~14×
  maximum / 4× median effect: *"reduced RNA structure and not codon rarity itself
  is responsible."* BT4's `RampTerm` rewards **lower codon adaptiveness** over the
  first ~35 codons — i.e. it implements the mechanism this study falsified. The
  effect is real; the lever is wrong.
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
| **Codon optimality / CSC** | Codon identity → **mRNA stability** via elongation-coupled decay | Presnyak & Coller, *Cell* 2015; Wu et al., *eLife* 8:e45396 (2019): *"the regulatory information affecting mRNA stability is encoded in codons and not in nucleotides"*, with stabilizing codons tracking **charged/total tRNA ratios**. Independently recovered by RiboNN, whose learned codon influence correlates **strongly with CSC**, moderately with tRNA abundance, and **negatively** with A-site occupancy (optimal GCU/GGU/GAU/AAC; suboptimal AGG/AGA/UCA/UUA) | **Not directly** — tAI is related but distinct; no CSC term. **This is the best-founded codon axis in the 2024–26 literature and the one BT4 lacks.** Honest caveat: no published head-to-head quantifies how much protein a CSC-optimized CDS gains over a CAI-optimized one — swapping CAI for CSC is better *reasoning*, not a validated gain |
| **Start-codon (Kozak) context** | The 6 nt preceding the AUG | Shukla et al., *NAR* 54(14):gkag728 (2026): randomizing **only those 6 nt** spans ~**100-fold** steady-state protein across 4,042 variants, on an isogenic landing pad that removes integration confounding | **Cannot see it** — −3 lives in the 5′UTR. For a fixed protein the CDS owns a narrow but real slice: +4/+5 are codon-2 bases fixed by the amino acid except for the 6-fold residues (Leu/Ser/Arg), while **+6, the wobble base of codon 2, is always synonymously free** |
| **oORF (uORF overlapping into the CDS)** | A 5′UTR AUG whose stop lies **inside** the CDS | Johnstone, Bazzini & Giraldez, *EMBO J* 35:706 (2016): oORFs repress significantly more than non-overlapping uORFs (Wilcoxon P = 1.23e−3) because ribosomes cannot reinitiate downstream of one | **Cannot compute it** — the AUG is in the UTR and the **stop position is set by the CDS's synonymous choices**. Neither half owns the answer. Deterministic, no ML required, and **no shipping optimizer performs it** |
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

**Honest read.** These predict/generate well in-distribution but are usually not
validated for the *specific* CDS-variant regime a codon optimizer operates in.
The sharpest statement of that gap: **RiboNN has never been shown to discriminate
synonymous CDS variants of the same protein under a fixed UTR** — which is
precisely BT4's regime. BT4's posture — wrap the user's own installed model,
hash-pin weights, ship `calibrated=False` until a regime-matched gate passes — is
exactly right and matches how the field's more careful practitioners talk about
their own models.

Two later entrants have better wet-lab evidence but are not shipping tools:
**RiboDecode** (*Nat Commun* 2025) — generative CDS design trained directly on
ribosome profiling, reporting **>2× the best LinearDesign sequence** in vitro and
~10× stronger neutralizing antibody for influenza HA in BALB/c mice (single-lab,
not replicated); and **GEMORNA** (*Science* 390:6773, 2025) — up to 41× firefly
luciferase and 15× human EPO (wet-lab numbers confirmed; whether its CDS design
conditions on the chosen UTR is **UNVERIFIED**).

---

## 3c. The gap that outranks both camps: construct context

Added after the 2026 verification sweep, because it changes the ranking in §4.

Every failure mode in this literature that has **reached patients** is a junction
phenomenon — the CDS supplies one half and the vector supplies the other:

- **Cheng et al., *Int J Biol Sci* 18:4914 (2022).** 17 genes × 4 common vectors.
  **13 of 17 spliced at their own retained exon–exon junction** (`AN|GTRAG`, e.g.
  hFACI `CAG|GTAAG`); **17 of 17 spliced at the V5 tag** (`G|GTAAG`). **Every
  acceptor was supplied by the vector** — the mPGK–PuroR linker, SV40, Neo/KanR.
  Fixed by synonymous recoding of the donor. *(BT4's `GTRAGT` donor consensus is
  well-aimed at exactly these offenders.)*
- **Kowarz et al., *eLife* 11:e74974 (2022).** Codon optimization *created* donors:
  wild-type Spike had 6 donors / 5 acceptors; Ad5.S had 12/10, ChAdOx1-S 13/5,
  Ad26.COV2.S 9/10. The acceptor was the adenoviral **pIX** site, reached by
  read-through **past the cassette poly(A)**. *(Cite the splicing mechanism only —
  the paper's downstream clinical inference is contested.)*
- **De Ravin et al., *Nat Commun* (2022).** A cryptic acceptor in the **cHS4-400
  insulator** drove HMGA2 truncation and clonal expansion in a clinical trial —
  the clone reached 15–22% of CD34+/myeloid cells. Fixed by a **2-bp `AG`→`TG`
  change**.

The corollary for BT4's splice audit is mechanical and severe: **splice events
need two halves, and the CDS supplies only one.** A CDS-only optimizer, however
good its model, is structurally blind to the failure mode that has reached
patients. One fortunate asymmetry: LVV work (*BMC Biotechnol* 9:86, 2009) shows
that **mutating cryptic acceptors merely shifts usage to other acceptors** —
donors are the tractable target, and donors are what a CDS optimizer owns.

**And BT4's current splice input is worse than no context.** OpenSpliceAI (Chao et
al., *eLife* reviewed preprint 107454) documents that SpliceAI *"exhibits an
inherent bias near the starts and ends of transcripts which are padded with
flanking N's … predicting donor and acceptor sites in these boundaries with an
extremely high signal **that disappears when the sequence is padded with the
actual genomic sequence**."* That N-padding is SpliceAI's own reference default
(`'N'*5000 + seq + 'N'*5000`) and is what BT4 does today — so BT4's audit is most
likely to hallucinate sites at the two positions a designer cares about most: the
start of the CDS and the region near the stop. The same re-benchmark shows
accuracy rises steeply only through ~400 nt of context (**+62%/+74%** donor/acceptor
from 80→400 nt, then marginal gains to 2 kb and 10 kb), so **±400 nt of real flank
beats ±5 kb of N**. That is a strictly-better fix costing one input field.

**Prior art BT4 is behind on.** *TIsigner* (*NAR* 49:W654, 2021) takes `-u/--UTR`
and optimizes accessibility over windows with **negative coordinates** (*E. coli*
−24:+24, yeast −7:+89, **mouse −8:+11**). *DNA Chisel* (MIT) optimizes a
sub-region via `CodonOptimize(location=…)` while `EnforceGCContent(window=50)` and
`AvoidPattern` evaluate over the **entire record** — a 10 kb plasmid included. The
mRNA-medicine tools (mRNArchitect, mRNAid, mRNAdesigner) accept UTR fields but
treat them as **assembly slots**; mRNAid's own paper says transcripts are "fused
to already optimized UTRs, obviating the need for optimization of these regions."
So *context-as-optimization-substrate* is genuinely unclaimed — but
*context-as-constraint-scope* is not, and BT4 should not claim novelty it does not
have.

---

## 4. Joint codon + structure optimization — demoted, with the evidence

**An earlier draft of this document ranked this the top priority on the strength
of LinearDesign's headline number. The 2026 verification sweep does not support
that ranking.** The recommendation is retained, at lower priority and with the
counter-evidence attached.

### 4a. The counter-evidence, which must travel with any structure claim

1. **The titre result has never been replicated, and the expression gain is
   ~2–3×, not 128×.** The 57–128× is anti-spike IgG in n = 6 mice against a
   *vendor codon optimizer's* output; HEK293 protein was **2.9×**.
2. **In an independent 2026 mammalian bake-off, LinearDesign was the worst
   performer.** Yang/Zhu/Meijers (Institute for Protein Innovation, bioRxiv
   10.64898/2026.03.18.712111), 18 human and murine Wnt-pathway glycoproteins,
   five strategies: *"strategies prioritizing RNA stability consistently reduced
   expression"* and *"codon optimization did not provide a general advantage over
   native coding sequences."* The counterweight that must ship with this citation:
   a **skewed** scheme using the most abundant codons matched native and
   "occasionally enhanced protein output." *(Full text 403s; the abstract-level
   claims above are confirmed, finer details are UNVERIFIED.)*
3. **MFE-optimized mRNA gives monosome-dominated polysome profiles.** Leppek et
   al., *Nat Commun* 13:1536 (2022), PMC8940940: a LinearDesign construct showed
   an "unusual monosome-concentrated polysome profile" — low ribosome load — yet a
   long in-cell half-life and ~2× luciferase at 24 h. **There is a sweet spot;
   over-stabilization causes stalling → collisions → decay.**
4. **Folding free energy correlates only weakly with in-cell lifetime and protein
   expression** (Jin et al., *JBC* 301:108015, 2025, verbatim). The same paper
   notes a methodological asymmetry: mRNA-1273 and BNT162b2 use modified
   nucleotides, but their MFEs in the LinearDesign study were computed with the
   **unmodified** energy model.
5. **MFE minimization is substantially a GC-maximization proxy.** So **folding ΔG
   and GC are not independent Pareto axes** — a frontier that treats them as such
   will mislead. BT4 already ships GC target, GC-window and max-GC-run, which makes
   this a live risk rather than a hypothetical one.

### 4b. What survives, and is still worth building

- BT4 is unusually well-placed to do it honestly: it already has the codon
  trellis, a ViennaRNA `FoldingModel`, and an optimality-certificate framework the
  exact-design camp lacks a clean version of.
- **The position-dependence is the part with real, uncontested evidence**, and it
  is buildable without a joint solver — see §4c.
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
  (beam, windowed folding) it must say so, per §5/§10.1. **Never claim an MFE-
  minimized design improves *expression*** — the surviving claim is stability.

### 4c. The cheap version, which is where the uncontested evidence actually is

Position-dependent structure has hard numbers, and none of them need a joint
solver — they need a **defined cap distance**, which needs a 5′UTR:

- **Cap-proximal thresholds.** Babendure et al., *RNA* 12:851 (2006): translation
  efficiency drops **abruptly as hairpin ΔG goes −25 → −35 kcal/mol**, and
  shifting a hairpin **9 nt** relative to the cap modulates translation
  **>50-fold**. Kozak (1989): a −30 kcal/mol hairpin **12 nt** from the cap
  inhibits; the same hairpin at **50–60 nt** does not.
- **The two regions.** Mauger et al., *PNAS* 116:24075 (2019): the inhibitory
  region is the 5′UTR **plus the first ~30 nt of the CDS**; the facilitatory
  region is the remainder of the CDS plus the 3′UTR. Functional half-life ↔ total
  output **r = 0.90**; translation rate ↔ output **r = 0.45**.
- **The Kozak-1990 conditional, which no codon optimizer models.** A hairpin
  **inside the CDS at ~14 nt** downstream of the AUG *improves* recognition of an
  AUG in a **suboptimal** context, and is neutral or harmful with a strong one
  (Kozak, *PNAS* 87:8301, 1990). With the UTR supplied, BT4 can read the actual
  Kozak strength and **set the sign of the downstream-structure term accordingly**.

**Therefore: do not build a single global MFE/ΔG objective.** Structure's sign
flips somewhere around nt +15 to +40. LinearDesign itself excises 5 codons from
its own DP because one global objective cannot express this. A two-region folding
model, anchored to a real cap distance, is both more correct and far cheaper than
a joint solver.

**Complexity note, stated correctly.** MFE-based mRNA design is **polynomial**
(O(n³) with the standard ≤30-nt two-loop cap). The partition-function /
ensemble-free-energy version is **not proven NP-hard** — Dai et al.
(arXiv:2401.00037) call it *"wide open, and likely NP-hard."* The correct
certificate justification is *"no exact algorithm is known"*, never *"the problem
is proven intractable."* Relatedly, adding a whole-sequence budget (GC count, CpG
count) to a CFG∩DFA intersection is exactly what breaks LinearDesign-style
polynomiality — **which is why BT4's bucketed-DP / Lagrangian machinery is a
complexity-level differentiator, not a feature-list one.**

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

- **The ceiling is upstream of the CDS, but not where this document used to say.**
  Not "the CDS is 31% of the signal" (§0 — that is a misread of a per-nucleotide
  density). The real ceiling: **mRNA abundance, not translation rate, is the
  majority channel** for protein abundance (Li, Bickel & Biggin 2014: mRNA levels
  explain ≥56%), and **integration site alone spans ~1,000×** (Akhtar 2013), while
  **Kozak context alone spans ~100×** (Shukla 2026). Promoter, integration site and
  UTR design will dominate any CDS result. Instructive decomposition of a *real*
  optimization gain: Fath et al. 2011, stably integrated *mip-1α* in CHO —
  transcription **+30%**, half-life **+14%**, steady-state mRNA **+80%**, TE
  **+20%** → **~3× protein**. The gain came mostly through mRNA level.
- **Honest expected effect sizes**, which BT4 should show next to every result:

  | Starting point | Honest expectation | Basis |
  |---|---|---|
  | Native human gene, GC-balanced, good vector | **1.0–1.5×; often indistinguishable from noise** | Fath 2011 (<2× for ~43% of 50 genes); Radrizzani 2024; IPI 2026 |
  | Native gene that is AT-rich, CpG-rich, or carries cryptic splice sites / uORFs / internal poly-A | **2–10×** | Fath's right tail (JNK3 14.7×, AQP5 9.1×, JAK2 6.0×) |
  | Non-mammalian ORF in mammalian cells | **5–100×** | the historical big wins are all cross-kingdom |
  | Any gene in *E. coli* | **wide (10–250×) across synonymous space — but not predicted by CAI** | Kudla, Welch, Boël |
  | **Risk of making it worse** | **~4–15%** | Fath (SLC39A1 **0.29×**, VKORC1 0.76×); Konkle 2021 (clinical expression loss) |

  The Fath numbers have **vendor authorship** (GeneArt/Life Technologies), so treat
  the 45%-of-genes-≥2× success rate as an **optimistic bound** and say so.
- **Optimization can harm.** Maximizing yield can disrupt co-translational folding.
  Harmonization is *a* response, but state its trade honestly: Moss, Chamness &
  Clark (*Annu Rev Biophys* 53:87, 2024) find harmonization "can improve protein
  **folding yield**, despite **lower overall protein accumulation**," and **there
  is no large-N head-to-head showing it beats optimization on titer.**
- **No tool is universally best.** Host and target decide. BT4's per-organism,
  per-reference-set, per-objective transparency is the right response.
- **Validation is the currency.** BT4's differentiator is method and honesty, not a
  validated expression claim — and it must keep saying so until it has its own
  gated data.

---

## 8. Ranked recommendations for BT4

**Re-ranked by the 2026 verification sweep.** The previous #1 (joint
codon+structure) has moved down; the evidence behind it did not survive scrutiny
(§4a), while the context gap gained three patient-level citations (§3c).

1. **Accept an optional 5′UTR, and evaluate every existing constraint across the
   UTR⊕CDS junction** (§3c). Mechanical, and it is the precondition for items 2–4.
   Seeding the trellis makes a junction-spanning restriction site *avoidable*
   rather than merely *reportable*.
2. **uORF/oORF pairing across the junction** — *the highest value-per-cost item in
   this document.* Deterministic, no ML, no calibration gate: scan the user's
   5′UTR for uAUGs, walk each out-of-frame one into the CDS to its in-frame stop,
   and classify uORF (stops in the UTR) vs **oORF** (stops in the CDS). Because the
   oORF stop position is a function of synonymous choices, BT4 can **move** it. No
   shipping optimizer does this.
3. **Two-region folding anchored to a real cap distance** (§4c), replacing any
   single global ΔG. This is where the uncontested positional evidence lives.
4. **Feed the splice CNNs real flanks and never N-pad** (§3c). ±400 nt of real
   sequence beats ±5 kb of N, and N-padding is a documented artifact generator at
   exactly the positions that matter. Ship as Δ vs a reference cassette.
5. **Accept a vector backbone as audit-only context** — never as an optimization
   substrate, and never transmitted anywhere (it is the user's IP).
6. **Codon optimality / CSC** as a first-class axis distinct from tAI (§2), with
   the caveat that no head-to-head quantifies the gain over CAI.
7. **Vendor manufacturability rules BT4 lacks, exactly as published:** windowed GC
   **range** (max−min across 50 bp windows ≤ 50 points, Twist DOC-001081 REV4);
   a **Tm-based repeat trigger** (any repeat with Tm ≥ 60 °C regardless of length,
   Twist — length-only rules miss GC-rich short repeats); **per-base homopolymer
   limits** (IDT: A/T ≥10, G/C ≥6 — the only published asymmetry, and
   thermodynamically motivated). Report a *profile*, not pass/fail: Twist's own
   documentation says "there is no single reason for the rejection of a sequence."
8. **An expected-effect-size disclosure next to every result** (§7). Costs almost
   nothing, every number is cited, and it is the honesty failure that matters most
   — reporting "CAI 0.42 → 0.91" without it is precisely the unvalidated-number-
   presented-as-real that §1 of the constitution forbids.
9. **mRNA-therapeutic terms** (§5): uridine depletion + m1Ψ slippery-sequence
   avoidance, both exact/LOCAL and cheap.
10. **Named application presets** (§5), including vector-type budgets.
11. **Joint codon + secondary-structure optimization** (§4) — retained, demoted,
    and to be built only with §4a's counter-evidence attached.
12. **Keep wrapping learned models, keep them `calibrated=False`** (§3b).

Items 1–2 and 6–10 are additive over the existing exact core. Item 11 is the one
large architectural investment, and it is no longer the one the evidence most
strongly rewards.

### What NOT to do

A short list, because each of these is a plausible-sounding move the evidence
refuses:

- **Do not claim "the 5′UTR is 67% of expression."** Quote both attributions (§0).
- **Do not build a single global MFE objective** — structure's sign flips (§4c).
- **Do not put folding ΔG and GC on a frontier as independent axes** — MFE
  minimization is substantially a GC-maximization proxy (§4a.5).
- **Do not claim MFE-minimized designs improve expression** — the surviving claim
  is stability, and four independent results say over-stabilization *reduces*
  expression (§4a).
- **Do not assert ensemble-free-energy codon design is NP-hard** — it is explicitly
  open (§4c).
- **Do not flip any model to `calibrated=True` because context improved.** Feeding
  SpliceAI real flanks fixes an *input* defect. Feeding RiboNN a real UTR makes it
  *runnable*. Neither is a gate.
- **Do not treat "highly-expressed reference set" as a fix for CAI** — it fixes
  which genes were counted, not that frequency is the wrong signal (§1).
- **Do not ship generic internal Shine–Dalgarno avoidance** — Li/Oh/Weissman 2012
  is contested by Mohammad 2016 (ribosome-profiling artifact). Ship the σ70 + SD +
  in-frame-ATG **co-occurrence** toxicity rule instead, which has a documented case.
- **Do not ship the "~150 bp palindrome lethality" or "≥50 bp spacer" rules** —
  refuted as attributed. What survives: a **246 bp** interrupted palindrome is an
  SbcCD target; **85 bp** is not.
- **Do not hard-code vendor thresholds BT4 cannot cite**, and do not merge Twist's
  *global* 25–65% GC (Express Genes) with GenScript's *per-100 bp-window* 25–65% —
  different rule types.
- **Do not present vendor complexity scores as outcome predictors** — Nguyen 2024
  found no correlation with assembly fidelity, and IDT moved its own ceiling from
  10 → 100 → 150 within two months of 2026.
- **Do not send the vector backbone anywhere.** It is IP. The ASSP cross-check is
  BT4's only outbound control and must be hard-blocked from transmitting backbone
  bytes — the simplest honest answer is to disable it whenever a backbone is loaded.
- **Do not vendor LinearDesign** — redistribution is prohibited and a Baidu patent
  is filed. Lazy-drive the user's own build, exactly as BT4 does for SpliceAI.
- **Do not build an m6A objective** — deposition is 3′UTR/stop-proximal and the
  `DRACH` consensus is not a usable CDS lever.

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

Added by the 2026 verification sweep (each fact-checked against a primary source):
- Zheng et al. (RiboNN), *Nat Biotechnol* 44(5):783 — per-nt 67/31/2 **and**
  length-integrated 22/73/5; human r = 0.78; early codons ~2× attribution weight.
- Radrizzani, Kudla, Izsvák & Hurst, "Selection on synonymous sites: the unwanted
  transcript hypothesis," *Nat Rev Genet* 25:431 (2024).
- Goodman, Church & Kosuri, *Science* 342:475 (2013) — the 5′ ramp is structural.
- Wu et al., *eLife* 8:e45396 (2019) — CSC; "encoded in codons and not in
  nucleotides."
- Shukla et al., *NAR* 54(14):gkag728 (2026) — ~100× from the 6 nt before the AUG.
- Johnstone, Bazzini & Giraldez, *EMBO J* 35:706 (2016) — oORF > uORF repression.
- Akhtar et al., *Cell* 154:914 (2013) — TRIP; ~1,000× integration-site range.
- Li, Bickel & Biggin, *PeerJ* 2:e270 (2014) — mRNA level ≥56% of protein abundance.
- Mauger et al., *PNAS* 116:24075 (2019); Babendure et al., *RNA* 12:851 (2006);
  Kozak, *PNAS* 87:8301 (1990) — positional structure.
- Leppek et al., *Nat Commun* 13:1536 (2022), PMC8940940; Jin et al., *JBC*
  301:108015 (2025) — the over-stabilization counter-evidence.
- Yang/Zhu/Meijers (Institute for Protein Innovation), bioRxiv
  10.64898/2026.03.18.712111 — independent mammalian bake-off.
- Cheng et al., *Int J Biol Sci* 18:4914 (2022); Kowarz et al., *eLife* 11:e74974
  (2022); De Ravin et al., *Nat Commun* (2022) — vector-context cryptic splicing.
- Chao et al. (OpenSpliceAI), *eLife* reviewed preprint 107454 — the N-padding
  artifact and the context re-benchmark.
- Wu, Yang & Colosi, *Mol Ther* 18:80 (2010) — AAV never packages >~5.2 kb
  (supersedes Grieger & Samulski, *J Virol* 79:9933 (2005) on capsid content;
  Grieger remains authoritative on the post-entry penalty).
- *NAR* 53(2):gkae1170 (2025) — AAV ITR instability at CG direct repeats; ΔsbcC
  at 42 °C; Stbl3 does *not* protect the ITRs.
- Bzymek & Lovett (2001) — the 7 bp / <100 bp RecA-independent repeat regime.
- Moss, Chamness & Clark, *Annu Rev Biophys* 53:87 (2024) — harmonization's trade.
- Twist tech note DOC-001081 REV4 (windowed GC range, Tm-based repeat trigger);
  IDT per-base homopolymer guidance; Salis 2020 (synthesis-failure features).
- Shi et al., *J Adv Res* (2026) — "coordinated design strategies that account for
  cross-regional dependencies."

*Compiled 2026-08, revised 2026-08 after a six-lens verification sweep in which
every load-bearing claim was independently re-checked and anything unconfirmed was
marked or removed. Two claims the earlier draft carried did not survive and are
corrected in place: the RiboNN attribution (§0) and LinearDesign's headline (§4).
Claims attributed to a single group are marked as such; treat in-vivo
fold-improvements as that group's reported results, not settled constants.*
