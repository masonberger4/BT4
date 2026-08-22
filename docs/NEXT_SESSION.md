# BT4 — live status & next-task queue

**This file is the single source of truth for *volatile* state**: where BT4 is
right now, and what to build next. The durable rules — architecture, the three
contracts, the honesty invariants, the anti-patterns — live in
[`../CLAUDE.md`](../CLAUDE.md) (the constitution) and are **not** restated here;
this file links into it by section. Shipped history lives in
[`../CHANGELOG.md`](../CHANGELOG.md). Update this file **in the same PR** as the
code it describes (§10.11): flip a status row, re-point the queue — a small
bounded edit, not a prose rewrite.

Every `DONE` below is stamped so a claim can be checked against the tree.

---

## Status board

**Phases:** 0–2 complete · Phase 3 groundwork landed · Phase 4 **app polish
landed**, learned-expression head still calibration-blocked · Phase 5 open, with
**declared CAI reference sets** (highly-expressed by default) landed. All merged
and green on `main`.

**Status vocabulary:** `DONE` · `GROUNDWORK` (contract + baseline shipped,
calibration pending) · `BLOCKED-data` (needs a matched-regime panel) ·
`BLOCKED-human` (needs licensed weights / a maintainer machine) · `NOT-STARTED`.

> ✅ **The four Tier-0 defects (measured 2026-08) are FIXED.** See
> [`REVIEW_2026-08_sota_and_roadmap.md`](REVIEW_2026-08_sota_and_roadmap.md) §3 for the
> fix record and before/after reproduction, and
> [`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md)
> §4 for the original measured evidence. `run_frontier` now downgrades a point that
> violates a GLOBAL rule to a `RELAXED` certificate (never `proven_optimal`);
> `run_validate` / `POST /validate` now audit GLOBAL constraints; `folding_dg` now
> reports the optimized 5′ window; and `avoid_internal_start` degrades gracefully via a
> new opt-in `relax()` with a culprit-named `InfeasibleError`.
>
> ✅ **Tier 1 (defensible default) and the Tier 2 construct-context core have also
> landed** — windowed GC, IUPAC extra sites, application presets, Studio weights /
> budgets / FASTA open / validate panel / splice track, and then
> `ConstructContext` + junction-correct constraints + junction folding + the
> whole-construct audit with restriction-site uniqueness. See
> [`REVIEW_2026-08_sota_and_roadmap.md`](REVIEW_2026-08_sota_and_roadmap.md) §4.
> ⚠️ **BT4's splice risk pooling was structurally mute in BT4's own regime, and is now
> honest about it (2026-08-19).** `pool_log_odds` counts only positions above
> `DEFAULT_SITE_PROBABILITY = 0.5`; measured against the hash-verified Pangolin weights
> on the designed-CDS panel, **only 6 of 93 sequences carry any position above 0.5**
> (all six designs of one protein), so `delta_splicing` was identically zero for two of
> three proteins while the raw scores varied more than twofold — and the **PWM baseline
> `default()` returns clears 0.5 on 93 of 93**, peaks 0.981–1.000, so the same constant
> saturates the shipped path while flooring the opt-in one.
> The background was **deliberately not lowered** (same uncalibrated knob, better view).
> Landed instead: `pool_top_k_logit` (background-free ranking statistic, **not a risk**),
> `PooledRisk` / `pooled_risk_detail` (a zero is now attributable), pooling coupled to
> the localization threshold in the audit and cross-check, and every risk-reporting
> surface saying *which* zero it is. Re-measured, Pangolin's background-free response
> spread across synonymous designs of one protein is 3.9–5.9 log-odds — it responds, and
> BT4 was discarding it. Evidence:
> [`REVIEW_splice_calibration.md`](REVIEW_splice_calibration.md). **The operating point
> itself is still underived** — that is Part B, and it is `BLOCKED-data`.
>
> ⚠️ **And the `N`-padding is not neutral (measured 2026-08-20).** Replacing the
> adapters' 5,000 literal `N` with real human genomic flank raises the median peak score
> *inside the CDS* from **0.276 → 0.369** (same 9-sequence set), moving designed
> sequences across the 0.5 cutoff. Controls make it interpretable: three **different** real regions agree to three
> decimals (real context is a stable background), random uniform ACGT scores *below* `N`
> (so it is not "any bases beat `N`"), and a composition-matched **shuffle inflates
> scores 9/9** (so the lift is distribution shift, not restored function). Consequence:
> a splice number from the `N`-padded path is a **lower bound**, and passing a real
> `ConstructContext` changes the answer rather than refining it. It licenses **no**
> threshold change — there are no labels, and a higher score is not a more correct one.
>
> ✅ **The models are NOT blind in BT4's regime (measured 2026-08-20).** A planted
> textbook donor lifts the local peak **0.052 → 0.570 (~11×)** at exactly the anchor base
> `CNN_ANCHOR_OFFSETS` predicts, in a third-party designed CDS host. Signal-specific, not
> change-driven: composition-matched scramble **0.0525**, `GT`→`CT` ablation keeping 7 of
> 9 bases **0.0543**, both at host baseline (0.0524). **But the floor is high** — a
> *weakened* real donor scores **0.357** and clears nothing, so the 0.5 cutoff sits above
> the intermediate-strength sites cryptic splicing actually uses. This **retires the
> "train the models" branch**: they detect a strong site in this exact regime, so silence
> on clean designed CDS is not inability to see. What remains is the **operating point**,
> which is derived on labelled data, not trained. It does **not** establish correct
> silence — detecting a site BT4 planted says nothing about sites nobody put there.
>
> ✅ **The candidate ranking is reliable, but the delivered PICK is not (measured
> 2026-08-20).** A two-facet generalizability study over Pangolin's 12 members (3 folds ×
> 4 tissues) — the fold/tissue split matters, folds are re-training noise while tissues are
> biology, so a random split-half is uninterpretable — gives **Eρ² = 0.959 / 0.901 / 0.942**
> under the tissue-general universe. Candidate variance beats every error term 5–10×, so
> **"the ranking is ensemble noise" is excluded without any labels**. But the top candidates
> are near-ties: in **2 of 3 proteins the argmax changes** with fold or tissue, and Beclin1's
> worst case delivers a sequence the full ensemble ranks **7th of 30**. Also sharpens the
> low PWM agreement: Pangolin is *stable*, so the disagreement is substantive and one of the
> two is wrong. **Reliability is not validity** — a ranking can be reproducible and wrong.
>
> ✅ **Tier 3 (GenBank I/O) has landed too** — an annotated reader/writer whose
> `misc_feature` spans put residual violations on the map the user opens, and whose
> reader turns an existing vector map into a `ConstructContext`. **Tier 4 has now
> partly landed too** — the functional (bipartite) poly(A) constraint and AAV/LVV
> packaging accounting, and **tAI now reaches every bundled organism** (E. coli was the
> last gap; exercising its bacterial path surfaced and fixed a latent tAI
> fidelity bug).
>
> ✅ **The three industrial expression hosts have landed (2026-08-19)** — **CHO**
> (CHOK1GS_HDv1), ***B. subtilis*** 168 and ***K. phaffii*** GS115 (*Pichia
> pastoris*), taking BT4 from nine selectable organisms to **twelve**, each with a
> recounted genome-wide codon table *and* a GtRNAdb tRNA table (the shipped
> invariant: a codon table without tRNA data makes tAI silently unavailable exactly
> where a user asked for it). Two limits ride with them and are recorded, not
> smoothed: CHO is **the one organism whose two inputs are not assembly-matched**
> (tRNA on CriGri_1.0, codons on CHOK1GS_HDv1), and none of the three ships a
> highly-expressed set — for three *different* measured reasons. **Start here is now
> the *B. subtilis* highly-expressed table**, which is a real, evidence-backed next
> step rather than a wish: its PaxDb integrated dataset exists and joins at 99.8%
> via a declared `BSU` → `BSU_` rewrite (queue item 6).

| Component | State | Calibrated? | Primary file(s) |
|---|---|---|---|
| Exact-DP codon trellis + certificate | DONE | n/a | `optimize/exact_dp.py` |
| Rust trellis port (`trellis_solve`, regime-gated) | DONE | n/a | `rust/bt4_core`, `bt4_native` |
| Objectives: CAI, tAI, GC, ramp, CpG, %MinMax, codon-pair | DONE | n/a | `objectives/` |
| tAI (real GtRNAdb, **all 12 organisms**; prokaryotic `sking` from provenance) | DONE | n/a | `biomodels/codon/tai.py`, `scripts/build_trna_tables.py` |
| Codon tables: **all 12** recounted from pinned Ensembl CDS (`genome_wide`) | DONE | n/a | `biomodels/codon/data/`, `scripts/build_organism_tables.py` |
| **Highly-expressed reference sets** (PaxDb top-300, 8 of 12 organisms, the **default**) | DONE | n/a (a declared reference set, not a model) | `scripts/build_highly_expressed_tables.py`, `biomodels/codon/tables.py` |
| Constraints: homopolymer, GC-run, max-repeat, tandem/inverted, forbidden+presets, restriction, Kozak-ATG, uORF, splice-motif | DONE | n/a | `constraints/` |
| Budget backends: CP-SAT, Lagrangian, dinucleotide-count | DONE | n/a | `optimize/{cpsat,lagrangian}.py` |
| SA refinement + block moves + parallel tempering | DONE | n/a | `optimize/anneal_refine.py` |
| Folding (ViennaRNA + labeled baseline) | GROUNDWORK | ViennaRNA=yes, baseline=no | `biomodels/folding/` |
| Splice PWM baseline | GROUNDWORK | no (baseline) | `biomodels/splice/` |
| Splice CNNs: Pangolin (GPL) / SpliceAI (CC BY-NC) | **BOTH fidelity gates PASSED** (exact, 18 cases each, same panel `f3589fd1…`) | **yes, opt-in** for both (`BT4_SPLICE_USE_ATTESTED=1`) — integration fidelity only, NOT statistical calibration | `biomodels/splice/{pangolin,spliceai}.py`, `biomodels/splice/attestations.py` |
| Splice audit (localize-and-flag) + backend agreement | DONE | advisory (`all_calibrated=False`) | `biomodels/splice/audit.py` |
| **Splice risk pooling** (`pool_log_odds` + the background-free `pool_top_k_logit`) | DONE, and **honest about being degenerate on designed CDS** — a zero risk now reports whether it was floored | operating point **underived** (0.5 is a convention, not evidence) | `biomodels/splice/base.py` |
| Splice fidelity-attestation layer | **DONE and in use** — a committed Pangolin attestation ships | n/a | `biomodels/splice/{attestation,attestations}.py`, `biomodels/splice/data/` |
| **ASSP cross-check (opt-in, out-of-loop network validator)** | DONE | `network_derived`, not calibrated | `biomodels/splice/assp.py`, `pipeline/splice_crosscheck.py` |
| Expression: `ExpressionPredictor` + `NullExpressionModel` + rerank hook | GROUNDWORK | placeholder=no | `biomodels/expression/`, `pipeline/rerank.py` |
| Expression: wrapped RiboNN (Sanofi non-commercial) | GROUNDWORK | **no** (acceptance gate pending) | `biomodels/expression/ribonn.py` |
| **Expression promotion seam** (`promote_if_attested`, `BT4_EXPRESSION_USE_ATTESTED`, `$BT4_EXPRESSION_ATTESTATION`, Studio toggle + scope display) | **DONE and wired** — `verified_predictor` now has production callers | n/a (carries a claim, cannot make one). **No attestation is bundled**: none has been earned | `biomodels/expression/attestations.py`, `app/studio.py` |
| Expression attestation **scope binding** (species/cell types/readout derived from the run; `top_k` + UTR context bound; gate record reconstructable) | DONE | n/a | `biomodels/expression/attestation.py`, `pipeline/expression_gate.py::GateScope` |
| Candidate-set assembly + expression rerank | DONE | calibrated-gated | `pipeline/candidates.py`, `bt4.api.candidates` |
| Library / degenerate-design (SAMPLED) mode | DONE | n/a (sampler, not optimizer) | `optimize/sample.py`, `pipeline/library.py` |
| Restriction catalog (584 enzymes, REBASE-derived + content-hashed) | DONE | n/a | `constraints/restriction.py`, `constraints/data/` |
| Surfaces: `bt4.api`, `bt4` CLI, FastAPI service, provenance | DONE | n/a | `api/`, `cli/`, `service/`, `provenance/` |
| **GenBank I/O** — annotated writer (residual violations as `misc_feature`) + reader (vector map -> `ConstructContext`) | DONE | n/a | `io/genbank.py` |
| **Functional poly(A) constraint** (bipartite hexamer + downstream U/GU element) | DONE | n/a (structural, not a cleavage predictor) | `constraints/polya.py` |
| **AAV/LVV packaging accounting** (report only — BT4 controls no lever) | DONE | n/a | `pipeline/packaging.py` |
| **BT4 Studio** — Design / Candidates+splice-audit / Library tabs, RiboNN + ASSP surfaced, menus + runtime theming | DONE | n/a | `app/studio.py`, `app/worker.py`, `app/theme.py` |
| Expression backend registry (`available_expression_backends` / `resolve_expression_backend`) | DONE | n/a | `biomodels/expression/__init__.py`, `api/` |
| Packaged app bundle (PyInstaller one-file / `.app`) | **DONE for Linux, verified by launching it** (2026-08-20); macOS/Windows built by CI, unverified here | n/a | `packaging/bt4-studio.spec`, `tests/test_bundle_spec.py` |
| Signed, double-clickable *installers* (`.dmg` niceties, MSI, code signing) | NOT-STARTED (signing is a deliberate non-goal, `packaging/README.md`) | n/a | `packaging/` |

The **expression/splice design flow**
([`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md)) is
**complete through step 5** — batched RiboNN scoring, the splice-consensus motif
constraint, candidate assembly + rerank, the localize-and-flag splice audit, and
the BT4 Studio UI that surfaces them — including the opt-in RiboNN head, which
now reaches the Candidates tab through the public expression-backend registry.
Only step 6 (auto-edit / auto-select) remains, and it is calibration-gated (see
the queue). The opt-in **ASSP** cross-check (the last named Phase-3 "still ahead"
item) has landed as a network validator *and* as BT4 Studio's one explicitly
consented, clearly-labeled network control.

---

## Strategic direction (read before picking large new work)

Two documents set direction, and they now agree on the same conclusion from
different ends.

[`RESEARCH_codon_optimization_SOTA.md`](RESEARCH_codon_optimization_SOTA.md) is a
grounded 2026 survey of the field. Its top finding: the strongest *validated
in-vivo* result (**LinearDesign**, *Nature* 2023 — joint codon + mRNA-structure
exact optimization) is the **same class of algorithm as BT4's codon trellis**, and
BT4 already has both ingredients but keeps folding in a refinement layer rather
than jointly in the DP.

[`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md)
is a measured audit of this tree against what BT4 is *supposed* to be. Its verdict:
BT4 is an unusually well-engineered, unusually honest **CAI optimizer with
manufacturability constraints**, not yet an expression optimizer — because at the
time of that audit **the optimizer only ever saw the CDS**. *(Partly addressed
since: the CDS is no longer designed in isolation — a construct context now feeds
junction-correct constraints, junction folding, leader-aware uORF pairing, and a
whole-construct audit. `score_in_context` now feeds the splice CNNs real flanking
sequence when the user supplies it -- literal `N` padding remains only for the
residual width beyond that, which a ~10 kb-context CNN unavoidably needs. RiboNN
still sits where it cannot influence delivery.)* The audit also
found **four defects measured by running the code**, three of which broke a §5
invariant; all four are now fixed.

The review's roadmap is the queue below. Note the dependency the two documents
share: joint codon+structure optimization that knows the real 5′UTR is strictly
better than joint folding that does not — so **construct context comes first**, and
LinearDesign-class joint design stays a separate track with its own design doc.

## Next-task queue

> ⚠️ **NEW, unqueued defect (found 2026-08-20): the splice path is organism-blind.**
> Neither `api.splice_audit` nor `api.splice_crosscheck` takes an organism, so BT4 will
> run **Pangolin's human heart/liver/brain/testis heads** on an *E. coli* or *B. subtilis*
> design and print a risk. Bacteria have no spliceosome; the number is not weak evidence,
> it is **inapplicable**. This got worse in #120, which added *B. subtilis* and
> *K. phaffii* as selectable organisms. Same shape applies to IVT-mRNA/saRNA delivery,
> which never enters the nucleus and so never meets the spliceosome. Fix by labelling
> **inapplicable** rather than deleting — one CDS routinely moves plasmid → mRNA → AAV —
> but do not print a risk. §5-shaped: an inapplicable model output presented as a result.


Ordered. Each item is tagged by precondition. **Pick the first `self-contained`
item unless you have a reason not to.** Full evidence and file:line anchors for
items 1–3 are in
[`REVIEW_2026-08_expression_and_context.md`](REVIEW_2026-08_expression_and_context.md)
§4 and §9.

1. **[DONE 2026-08] Tier 0 — the four measured defects.** All four fixed in the
   `codon-optimization-review` change; fix record + before/after reproduction in
   [`REVIEW_2026-08_sota_and_roadmap.md`](REVIEW_2026-08_sota_and_roadmap.md) §3.
   - **`run_frontier` now downgrades a GLOBAL-violating point to `RELAXED`** (never
     `proven_optimal`) naming the unenforced rule, with the residual in the audit
     (`pipeline/optimize.py::run_frontier`). Repair of a violating seed stays in
     `run_optimize` and the candidate assembler; the frontier is a pure explorer.
   - **`run_validate` / `POST /validate` now audit the GLOBAL constraints**
     (`pipeline/optimize.py::run_validate`); covered by
     `tests/test_pipeline_api.py::test_validate_enforces_global_max_repeat`.
   - **`folding_dg` now reports the optimized 5′ window** (`model.score_sequence`),
     so reported == computed; covered by
     `test_refine_folding_audit_reports_optimized_window`.
   - **`avoid_internal_start` degrades gracefully** via a new opt-in
     `relax()` (`bt4.domain.relax.SoftConstraint`, `OptimalityStatus.RELAXED`
     now used) and `InfeasibleError` names the failing residue + culprit
     constraints (`optimize/exact_dp.py`). A non-opt-in rule (e.g. restriction) is
     never silently dropped. **Tier 1 presets are now unblocked.**
2. **[DONE 2026-08] Tier 1 — make the default defensible.** Landed: the
   **windowed-GC constraint** (`constraints/gc_window.py`) with honest
   tractability routing (<=12 nt exact in the trellis, wider refinement-enforced --
   the exact search grows ~exponentially in the window, measured); the IUPAC
   **`restriction_extra_sites`** path; **application presets**
   (`pipeline/presets.py`, none applied by default); Studio **objective weights**,
   **GC/dinucleotide budgets**, FASTA open, a **validate** panel, a **splice
   track**, and the previously hardcoded repeat/RC knobs; and the **ramp term
   relabelled** (its mechanism was falsified). See
   [`REVIEW_2026-08_sota_and_roadmap.md`](REVIEW_2026-08_sota_and_roadmap.md) §4.
   *(Superseded queue text for this item follows, kept for its evidence.)*

   **[superseded detail] Tier 1 — make the default
   defensible.** BT4's own `scripts/compare_tools.py` shows nine shipping tools
   clustered at CAI 0.63–0.83 / GC 42–54% while BT4's default sits alone at
   **CAI 1.000 / GC 62.4%**. The engine is not wrong; the operating point is.
   - **Windowed-GC constraint** — LOCAL, `context_len = window − 1`, exact in the
     trellis. It is listed in CLAUDE.md §6 and was never built; there is **no
     GC-content constraint of any kind** today. The soft term saturates at weight 2
     without reaching its target (it is separable, `context_len() == 0`) and the
     hard count budget cannot control clustering (74% window at 50% total GC). The
     windowed computation already exists in `pipeline/tracks.py`. Build the rule
     vendors actually publish: **GC *range*** — max − min across 50 bp windows
     ≤ 50 percentage points (Twist DOC-001081 REV4). Two cheap siblings, also
     published and also absent: a **Tm-based repeat trigger** (flag any repeat with
     Tm ≥ 60 °C regardless of length — length-only rules miss GC-rich short
     repeats) and **per-base homopolymer limits** (IDT: A/T ≥ 10, G/C ≥ 6, the only
     published asymmetry). Report a *profile*, not pass/fail — Twist's own docs say
     "there is no single reason for the rejection of a sequence" — and do not
     hard-code a threshold BT4 cannot cite.
   - **Application presets** (`mammalian_plasmid`, `aav`, `lentiviral`, `mrna_ivt`,
     `ecoli`) in a new `pipeline/presets.py`, mirroring the
     `constraints/forbidden.py` catalog pattern. The manifest must carry the
     **resolved field values**, never the preset name alone, and explicit user
     knobs must win with the override reported. A preset is a design profile, never
     a validated expression claim.
   - **Frontier-point picker** in Studio — the plot already shows better-balanced
     designs the user has no way to select.
3. **[self-contained, then architectural] Tier 2 — construct context.** The actual
   product gap: `OptimizeConfig` has 40 fields and not one is sequence outside the
   CDS; the DP seeds from an empty prefix (`optimize/exact_dp.py:161`). Strictly
   serial. See the review §9 Tier 2 for the full design, including why the right
   move is to **wrap the constraints (`SeededConstraint`) rather than seed the DP
   layer**, why objectives must stay CDS-local (invariant #4), the `cds_offset`
   trap in `kozak.py`/`uorf.py`, the `masked_spans` requirement for AAV ITRs and
   LVV LTRs, the whole-construct repeat **performance gate** (§10.8), and the fact
   that giving the splice CNNs real flanks **invalidates the queued N-padded
   fidelity attestation** rather than inheriting it.

   **Do the sub-items in this order — the cheapest one is also the best-evidenced:**
   - **oORF pairing across the junction** *(no ML, no gate, nothing else does it)*.
     An out-of-frame AUG in the user's 5′UTR whose in-frame stop lands **inside the
     CDS** is an oORF, which represses significantly more than a non-overlapping
     uORF (Johnstone et al., *EMBO J* 35:706 (2016), P = 1.23e−3). The stop position
     is a function of BT4's synonymous choices, so BT4 can **move** it. This is a
     prefix-seed on the existing GLOBAL `avoid_uorf` constraint.
   - **Junction-correct evaluation of every existing LOCAL constraint** (GC window,
     GC-run, homopolymer, restriction sites, forbidden motifs) — mechanical once the
     prefix is seeded.
   - **Real flanks for the splice CNNs, never N-padding.** ±400 nt of real sequence
     beats ±5 kb of N (OpenSpliceAI re-benchmark: +62%/+74% donor/acceptor from
     80→400 nt, marginal thereafter), and N-padding is a *documented artifact
     generator* at transcript boundaries. Fall back to the PWM baseline when real
     flanks are unavailable, and say why.
   - **Cap-distance-aware two-region folding** (§4c of the survey), which is where
     the uncontested positional evidence lives.
   - **Whole-construct audit + restriction-site uniqueness.** Note the prior art:
     DNA Chisel (MIT) already evaluates `EnforceGCContent(window=50)` and
     `AvoidPattern` over an entire plasmid record, and TIsigner already takes
     `-u/--UTR`. Context-as-*constraint-scope* is not novel; context-as-
     *optimization-substrate* is. Do not claim more than that.
4. **[DONE 2026-08] Tier 3 — GenBank I/O.** Landed as `io/genbank.py` (stdlib
   only). The writer emits **residual violations as `misc_feature` spans** at the
   base where they occur, so a defect the optimizer could not remove reaches the
   map the user opens rather than living only in a JSON audit; overlapping
   findings merge into one readable span labelled with how many it covers, and the
   true count stays in the COMMENT block alongside the certificate and the
   config-hash/git-commit stamp. The record carries **no timestamp** so it stays
   byte-reproducible (invariant #7). The reader (`parse_genbank` /
   `context_from_genbank`) closes the loop with Tier 2: **the vector map a user
   already has becomes the `ConstructContext` their CDS is designed inside**.
   Surfaced as `bt4 optimize --genbank`, `api.write_genbank`, and a Studio
   *Export GenBank* button (Ctrl+G). Verified against Biopython as the reference
   parser (`tests/test_genbank.py`, skipped when Biopython is absent). SnapGene
   `.dna` remains out of scope.
5. **[self-contained] Tier 4 — per-system biology.** *Partly landed 2026-08:* the
   **functional poly(A) constraint** (`constraints/polya.py`, `avoid_polya`) now
   forbids an `AATAAA`/`ATTAAA` hexamer **only when a downstream U/GU-rich element
   follows it** — the bipartite architecture CPSF/CstF actually recognise — so it is
   strictly more permissive than the blunt `poly_a_signal` hexamer preset, which
   stays available as the stricter option; its footprint (~45 nt) is far too wide
   for the trellis, so it is refinement-enforced and its residuals reported. Also
   landed: **AAV/LVV packaging accounting** (`pipeline/packaging.py`,
   `api.packaging_report`) — reporting only, since BT4 controls no lever over
   cassette size, and it names what it could **not** see rather than implying a
   partial count is a verdict. Both are wired into the AAV/LVV presets, CLI and
   Studio. **Still ahead in this tier:** donor×acceptor splice pairing as a report,
   uridine depletion, m1Ψ slippery motifs, and a codon-optimality/CSC term (the
   last is data-gated — do not ship a fabricated CSC table). Independent of Tier 2, so it
   can land early: a *functional* poly(A) constraint (hexamer **plus** downstream
   GU/U-rich element — still LOCAL at `context_len ≈ 36`), donor×acceptor splice
   pairing as a report, AAV packaging-size accounting (reporting only — BT4
   controls no lever), uridine depletion, m1Ψ slippery motifs, and a codon
   optimality/CSC term.
6. **[START HERE · self-contained] Phase-5 breadth, continued.** Nine organisms ship recounted
   genome-wide tables and eight also ship a highly-expressed reference set. What
   remains:
   - **Add further organisms — the three industrial hosts are DONE 2026-08.** CHO
     (`cricetulus_griseus_chok1gshd`, CHOK1GS_HDv1), *B. subtilis* 168 and
     *K. phaffii* GS115 all ship recounted genome-wide tables **and** GtRNAdb tRNA
     tables, taking BT4 from nine selectable organisms to twelve. None ships a
     highly-expressed set, for three different measured reasons — see the
     `WITHOUT_ABUNDANCE` note in `tests/test_highly_expressed_tables.py`, which
     records the evidence per organism rather than a blanket "not done".
     **The concrete follow-up this surfaced:** *B. subtilis* **can** have one. PaxDb
     v6.1 has its whole-organism integrated dataset (taxon 224308) and the join is
     available — PaxDb writes `BSU35360` where Ensembl Bacteria writes `BSU_35360`,
     and that `^BSU` → `BSU_` rewrite joins **4,042/4,052 = 99.8%**. It is a
     locus-tag punctuation difference derivable from the two pinned files alone, not
     a third-party mapping, so it is admissible under the *A. thaliana* rule; it
     needs `AbundanceSpec` to carry a **declared per-spec identifier rewrite** (plus
     a minimum-join-rate check so a bad rewrite fails loudly rather than thinning the
     reference set). That builder change is the next task here.
     Further organisms beyond these still follow the same rules: never fabricate a
     table, and never join abundance IDs through an unpinned mapping.
   - **Bacterial alternative start codons.** *(Now measured for a second organism:
     the filter costs *B. subtilis* **954 of 4,237 CDS = 22.5%** — TTG 553, GTG 387,
     ATT 8, CTG 5, ATC 1 — more than double *E. coli*'s share, because it genuinely
     uses alternative starts. Counting them back in moves **no** amino acid's top
     codon and shifts `w` by at most 0.023, so the same conclusion holds: a
     precision gap, not a wrong answer. This raises the value of fixing it.)*
     Both builders' shared validity filter
     requires an `ATG` start, dropping 409 of 4,239 *E. coli* CDS (9.6%) —
     including *tufA* and *hupB*. Measured impact on the shipped table: 16 of the
     300 selected genes change and **no** amino acid's top codon moves, so this is
     a quality gap, not a wrong answer. Fixing it means relaxing the filter in
     *both* builders together (they must stay identical or the two reference sets
     stop being comparable) and deciding whether an initiator codon is a codon
     *choice* at all, since the ribosome uses fMet-tRNA regardless.
   - **Tissue / cell-type-specific tables: dropped** (maintainer decision) — large
     effort, hard to qualify honestly, small upside over a whole-organism
     highly-expressed reference. Do not re-open without a new reason.
   - **tAI reach — DONE 2026-08.** *E. coli* now ships a GtRNAdb tRNA table (86
     genes, real data, source-SHA-pinned), so **all nine bundled organisms offer
     tAI**. Exercising the bacterial `sking=1` path surfaced a **latent fidelity
     bug**: BT4 computed `W[ATA] = p[8] * t[ATA]`, but the reference `get.ws` is
     `if(sking == 1) W[35] = p[9]` — a constant with *no* tRNA factor, because the
     AUA reader is tRNA-Ile2 (anticodon CAU, lysidine-modified), not a UAU tRNA.
     E. coli has zero TAT-anticodon genes, so the old form drove `W[ATA]` to 0 and
     the zero-filling step silently replaced it with the geometric mean. Fixed, and
     a test that had encoded the wrong model is corrected. `super_kingdom` now comes
     from the table's provenance instead of a hardcoded `sking=0`, and
     `scripts/build_trna_tables.py` makes tRNA tables re-derivable (`--verify`)
     for the first time.
7. **[self-contained] Remaining BT4 Studio work.** *(The control column no longer
   clips its widgets -- it is sized from the form's own size hint, which also
   un-hid the windowed-GC `max` spinbox that was previously off-screen.)*
   Beyond the frontier-point picker in item 2: a **basic/advanced split** (27 undifferentiated controls, of which
   only 2 are ones a bench scientist must set), **protein file open + drag-drop**
   (there is no `getOpenFileName` anywhere in the app today), *(**DONE 2026-08:**
   the metrics table is now driven by the audit dict instead of a hard-coded 9
   rows, so enforcement status, residual counts, tAI, CpG/UpA counts, relaxed
   rules and the folding read-out are no longer CLI-only. It keys off any
   `<rule>_enforced` audit key rather than a fixed list, so a rule added to the
   engine later surfaces without editing the GUI — which is how the poly(A) and
   windowed-GC rules would otherwise have been invisible here.)*
   surfacing the **seven engine capabilities the GUI hides** (CpG/UpA budget, GC
   budget, ramp axis, `--refine`, negative `cpb_weight`, `tandem_copies` /
   `inverted_loop`, `seed`), saving/restoring the control panel, richer per-site
   risk tracks, and a screenshot refresh.
8. **[self-contained] External-validation report** — compare BT4 output
   codon/GC/CpG distributions against real highly-expressed gene panels (§8), using
   public data and BT4's own recompute functions.
9. **[Linux DONE 2026-08-20 · rest human] Packaged installers.** The Linux
   one-file bundle now **builds and runs here** — the sandbox can install the Qt
   runtime and drive a real X display (`Xvfb`), so the old "can't self-test in this
   sandbox" limitation is retired: the packaged app was launched, given a protein,
   and produced a `PROVEN_OPTIMAL` design, a ranked candidate set and a sampled
   library. Doing that for the first time found a **release-blocking defect no
   from-source gate could see** — the spec collected `**/*.provenance.json`, so
   `ribonn_sha256.json` (read at import) was absent and the frozen app died before
   its first window. Fixed, pinned by `tests/test_bundle_spec.py`, and
   `--self-test` now runs a real design rather than only building the window.
   **Build the bundle and open it before every tag**; a green suite does not cover
   the artifact users download. What remains is human: macOS/Windows verification
   on real hardware, signing (a deliberate non-goal today), and the tag push /
   release cut (HTTP 403 in the sandbox).
10. **[Part A DONE 2026-08 — both CNNs passed · Part B `BLOCKED-data`] Splice CNN calibration.**
    **Pangolin passed its integration-fidelity gate**, on a maintainer machine
    holding the GPL weights: 18 cases, tolerance 1e-3, **max abs deviation exactly
    0.0** — BT4's adapter reproduces upstream's per-position scores bit-for-bit.
    The attestation is committed at
    `src/bt4/biomodels/splice/data/pangolin.attestation.json`
    (content hash `5176032c…`, eight license-clean scalars plus the public weight
    SHA-256s — never a raw score).

    **The pass is not vacuous.** All three donor probes peaked at position 302 with
    the consensus planted at 300–309, across three different random flank sets, so
    the model tracked sequence rather than the `N`-padding boundary. Designed CDSs
    averaged 0.085; consensus probes 0.61–0.71.

    **Promotion is opt-in** (`BT4_SPLICE_USE_ATTESTED=1`, or
    `--use-attested-splice`). `default()` still returns the PWM baseline, and a
    real-flank score still reports uncalibrated — the gate was captured N-padded,
    and regime scoping survives promotion.

    **What it does NOT establish:** that Pangolin's scores are calibrated
    *probabilities* for designed coding sequence. That is Part B of
    [`DESIGN_splice_cnn_calibration.md`](DESIGN_splice_cnn_calibration.md), still
    unmet, and the regime where these models are measured weakest (median prAUC
    **0.419 exonic** vs 0.773 intronic, Smith & Kitzman 2023) — which is BT4's
    entire regime. State that wherever BT4 reports splice risk on a CDS.

    **Part B now has a concrete, measured stake beyond "the score is uncalibrated".**
    The shipped operating point does not merely risk being *wrong* — at `0.5` it makes
    BT4's entire splice objective **inert in BT4's own regime**. Measured with the
    hash-verified Pangolin weights on the designed-CDS panel, only 6 of 93 sequences
    carry a position above 0.5, so `delta_splicing` was exactly zero for two of the
    three proteins — while the PWM baseline clears 0.5 on all 93, so the one constant
    floors the CNN and saturates the default. That is now
    *reported* honestly (`pool_top_k_logit`, `PooledRisk.below_background`; see
    [`REVIEW_splice_calibration.md`](REVIEW_splice_calibration.md)) but not *fixed* — a
    defensible replacement cannot be picked without labeled data, which is what makes
    this item `BLOCKED-data` rather than a knob to turn. Deriving the operating point is
    therefore part of this item's deliverable, not a follow-on.

    **Part B's machinery is now built; the data step is what is left.** Landed: the
    classification estimators (`pr_auc` as tie-grouped average precision, `roc_auc`,
    `mcc`, Brier + Brier skill, ECE, reliability bins, `top_k_accuracy`,
    `pr_auc_skill`), the acceptance gate itself
    (`biomodels/splice/gate.py` — two case types never mixed, no Spearman, a
    **per-stratum** verdict, and a mandatory `negative_construction`), the **panel
    format + strict reader** (`api.read_splice_panel`), and the baseline comparison
    (`api.splice_panel_gate`, `bt4 splice-gate`) over four permanent controls —
    `permutation` / `gt_ag` / `pwm` / `constant`. Run BT4's own PWM backend as the head
    and it ties the `pwm` baseline exactly, so `beats_every_baseline` is `False`: the
    shipped default can never be evidence for itself.

    The reader **refuses a mis-anchored panel** rather than scoring it — BT4 anchors a
    donor on the `G` of `GT` and an acceptor on the `G` of `AG`, and a panel built to
    the exonic-boundary convention (what the GENCODE recipe produces) is rejected with
    the exact shift that would have worked. The runner reports the matching backend-side
    trap: where each score actually peaked around a declared site, so a one-base anchor
    disagreement cannot masquerade as a hopeless model.

    **Panel acquisition is now scripted, and both sources are characterized.**
    `scripts/make_gencode_splice_panel.py` builds the site-prediction panel from a pinned
    GENCODE v44 + GRCh38 (URLs and md5s in the runbook's B1). Its arithmetic was
    *executed* against real GRCh38 -- 99.42% canonical GT/AG over 1,206 chr1 MANE sites,
    versus 0.08% and 44.2% for the two plausible wrong conventions -- and it handles the
    two traps that silently relabel true positives as negatives: a +/-5,000 nt window
    holds a **median of 8** annotated sites (only 2.8% hold just the centre one), and
    **27%** of gene-body windows contain opposite-strand sites.

    For variant effect, `kitzmanlab/splicebench2023` carries **no data in the repo** --
    it is one Zenodo archive (record 8351879), whose top directory must be renamed
    `for_zenodo` -> `data`. Label column `sdv_fc2`, stratifier `exon`, SpliceAI `DS_maxm`,
    Pangolin `pang_max_abs`; the paper's 0.419/0.773 exonic/intronic split reproduces
    from those columns to six decimals. **But 53% of it (BRCA1, FAS, WT1) is on
    chromosomes both models trained on** -- only the chr3 genes (POU1F1, MST1R, MLH1) are
    held out, and BT4's gate now correctly refuses to call the rest held out.

    *Provenance caveat:* Pangolin's held-out split is from its own paper; **SpliceAI's is
    not** -- the Cell paper is paywalled, so chr1/3/5/7/9 comes from OpenSpliceAI (eLife
    2025), which rebuilt its pipeline. Well-sourced second-hand, not primary.

    **The variant half has now been RUN** (2026-08-19, maintainer machine) and the
    measured results are recorded in
    [`REVIEW_splice_calibration.md`](REVIEW_splice_calibration.md). Headline findings, on
    `splicebench2023`'s own pre-computed scores so no weights were involved: removing the
    training-chromosome genes costs **a third of intronic skill** (0.724 → 0.480) while
    exonic *rises* slightly (0.365 → 0.419), collapsing the exonic/intronic skill gap from
    0.359 to **0.061**. The write-up is explicit that this is *confounded* — holding out
    chr1/3/5/7/9 replaces the gene set rather than filtering a fixed one, and the two are
    collinear — so it does **not** establish that the published gap is leakage. And held-out exonic **ECE is 0.345**, direct
    evidence these scores cannot be read as probabilities, which is what BT4's own 0.5
    operating point assumes. Three defects in BT4's *reporting* layer were found by those
    runs and fixed (#103, #104, #105) — each a claim about numbers that no test could
    reach.

    **The site-prediction half has now been RUN too** (2026-08-19, same machine), against
    BT4's own wrapped Pangolin on a GENCODE v44 / GRCh38 panel — so this half exercises the
    adapter end to end rather than a benchmark's pre-computed scores. On 20 held-out MANE
    windows (861,096 positions, 333 sites) Pangolin scores **skill 0.983 / top-k 0.940**
    against the `pwm` baseline's 0.096, and with a bar declared beforehand
    (`--min-pr-auc-skill 0.75`) the run reports **`PROMOTABLE on this panel: True`** — a
    first for any BT4 splice backend. The **per-kind anchors are confirmed on real data**
    (donor −1 for 100% of sites, acceptor +1 for 99%), which was this half's largest
    correctness risk. It is **not** a promotion: the figures sit +0.13/+0.15 *above*
    published, which says the panel (20 MANE gene bodies, 100% canonical motifs) is easier
    than the genome-wide benchmark, and every site is natural sequence rather than the
    designed synonymous CDS BT4 actually emits. `calibrated` is unchanged and `default()`
    still returns the PWM baseline.

    **SpliceAI's fidelity gate has now been run and PASSED too** (2026-08-19): 18 cases
    on the *same* panel (`content_hash f3589fd1…`), **max abs deviation exactly 0.0**,
    with the panel's peak scores spanning 0.029–0.925 so the pass is not the vacuous
    kind a flat panel gives. Its adapter had never executed against real weights before
    that run and was correct first time. Both attestations now ship, both honored under
    the one opt-in. **Part A is complete.**

    **The two-backend agreement has now been measured** (`bt4 splice-agreement`, both
    CNNs over panel20): Jaccard **0.855**, and of 333 annotated sites both find 300, only
    Pangolin 15, only SpliceAI 6, **neither 12**. Running both is *not* redundant — 21
    sites (6.3%) are found by exactly one model, which is what an audit should surface —
    but the 12 both miss are a **correlated blind spot**, and that is the standing limit
    on reading agreement as an uncertainty signal: it bounds independent error, not shared
    error. Like-for-like on the shared task (`--combined-track on`), Pangolin **0.983** vs
    SpliceAI **0.965**; the top-k gap (0.033) reproduces the published one (0.040) even
    though both absolute levels are inflated by an easy panel. All in
    [`REVIEW_splice_calibration.md`](REVIEW_splice_calibration.md).

    **Still ahead here:** a panel in BT4's own regime (designed synonymous variants, where
    the question is **specificity** not recall — every measurement so far is recall on
    natural sites); resolving the **7 positions both models call that the panel does not
    annotate** (non-MANE isoform sites, or shared false positives — opposite implications,
    and this panel cannot separate them); deciding whether the attested-promotion opt-in
    should become the default; and a Studio checkbox for it. *(The splice side already
    has one — `splice_attested_check`; what is missing is the CLI/promotion-default
    decision, not the control. The expression side's equivalent landed 2026-08.)*

    **The SpliceAI tooling proved out on first contact with real weights.**
    `scripts/capture_spliceai_panel.py` ships alongside the Pangolin one (same
    independence guard, statically enforced), and `run_splice_fidelity_gate.py`
    dispatches on the capture payload's own `backend` field, so a Pangolin capture can
    never be checked against the SpliceAI adapter. The SpliceAI capture **imports
    upstream's own `one_hot_encode`** instead of re-deriving it — SpliceAI ships it as a
    reusable function where Pangolin's CLI encodes inline — which makes a transposed
    layout or wrong base order in BT4's `_one_hot_rows` a gate *failure* rather than an
    artifact reproduced on both sides. It has **no fallback encoder** on purpose, and a
    test asserts it defines none: the obvious "fix" for the NumPy 2 `np.fromstring`
    breakage would silently destroy the independence — and on the real run the refusal
    fired for a *different* missing dependency and still refused correctly. The workflow
    is `make_splice_panel.py` → `capture_{pangolin,spliceai}_panel.py` (neither imports
    `bt4`) → `run_splice_fidelity_gate.py`; see the Splice CNN environment gotchas below
    for the `--no-deps` install that makes it run without pysam.

11. **[BLOCKED-data · human] Promote RiboNN to `calibrated=True`** — **the machinery
    is now built; only the data step is left.** Landed: `within_group` (the strict bar)
    and `recalibrate` on the gate, the cluster-bootstrap CI, the `width_over_iqr`
    vacuity check, the panel format + strict reader (`api.read_panel`), the baseline
    comparison (`api.expression_gate`, `bt4 expression-gate`,
    `scripts/run_expression_gate.py`), the `ExpressionAttestation` promotion seam, and
    the zero-data checks (`scripts/ribonn_sensitivity.py`). **The step-by-step runbook
    is [`DESIGN_ribonn_calibration.md`](DESIGN_ribonn_calibration.md)** (Windows/WSL
    install, the free Stage-1 checks, the panel hunt, pre-registration, the gate
    command, and the decision tree), backed by the evidence and corrections in
    [`RESEARCH_ribonn_calibration.md`](RESEARCH_ribonn_calibration.md). **Start instead
    from [`GUIDE_ribonn_calibration.md`](GUIDE_ribonn_calibration.md)** — the same
    procedure in plain language with every command verified against the code, the
    free/weights-free steps pulled to the front, and an Appendix B of the runbook's
    measured defects (the `bt4 expression-gate` shortcut cannot emit `gate_result.json`;
    the "~90 rows" panel floor fails a good head **56%** of the time, all of it on the
    coverage band, so size for **~200**; `$BT4_RIBONN_WEIGHTS` is undocumented anywhere
    else; and the job-1 "bit-for-bit" ✅ was an overclaim, now corrected).
    **Licence status: the maintainer reports the grant was received in writing (2026-08).**
    Recorded as reported, not as verified — no artifact in this repo can confirm it, and
    saying otherwise would be the kind of unbacked claim §5 exists to prevent. What *is*
    verifiable from upstream is the grant's shape: `LICENSE` and `MODEL WEIGHTS
    LICENSE.txt` grant use "to any person from academic research or non-profit
    organizations", so it is an **affiliation** grant, not merely a non-commercial one,
    and the guide's Step 1 has it resolved with `patent.gos@sanofi.com` before any
    download. **Guide Steps 5–9 (install + verify) are now DONE** on the maintainer's
    CPU-only Windows machine (2026-08): the weights are downloaded, the adapter scores a
    real CDS in `RiboNN CLR-residual TE` units with `calibrated = False`, the same sequence
    scores byte-identically twice, and `max_shift` is `0` in both weight sets. Getting
    there needed four corrections to this guide's own CPU recipe — see the gotchas section
    below, and guide Step 5 for the recipe. **The next step is therefore the panel, not the
    install.**
    **Remaining
    (human):** run the free sanity checks against the licensed weights, obtain a
    licence-clean regime-matched **CDS-variant** panel — no public dataset fully
    qualifies today, so read §4 of the research doc before spending anything — then
    pre-register thresholds and run the gate **once**. **The wiring caveat is GONE
    (2026-08):** `promote_if_attested` now has production callers behind
    `BT4_EXPRESSION_USE_ATTESTED` / `$BT4_EXPRESSION_ATTESTATION` and a BT4 Studio
    toggle, so a passing gate plus an attestation **does** change what a user gets —
    a ranked candidate set with the delivered pick chosen by the head, banner-led with
    the scope it was earned in. Nothing is bundled, so today it promotes nothing. Two
    things the run itself must now satisfy, both enforced rather than documented: the
    gate **refuses** a head configured to average every cell type against a panel that
    declares one (the old silent `--cell-type` trap, now caught *before* the scoring
    pass), and `attest_expression` **derives** the scope from the run instead of
    accepting it as free text, so a declared `cell_types` that disagrees is a refusal.
    Reproducing RiboNN faithfully is **not** calibration for BT4's
    CDS-variant regime — RiboNN **has never been shown to discriminate synonymous
    CDS variants of the same protein under a fixed UTR**, which is exactly the
    regime BT4 operates in. (Do *not* justify this with "only ~31% of signal is in
    the CDS" — that is the per-nucleotide density; the length-integrated
    attribution is 22/73/5, so the CDS is the majority. See the survey's §0.) Do
    not relabel a hand-weighted composite as "calibrated".
12. **[BLOCKED until #10/#11] Design-flow step 6** — targeted synonymous splice
    **auto-edit** and RiboNN **auto-select**, each unlocked only once its backend
    passes its gate. *(The delivery half of RiboNN auto-select is in fact already
    built and tested: a promoted head reorders the candidate set and re-picks
    `chosen`. What is still blocked is the only part that matters — a head that has
    **earned** the promotion. The remaining work here is therefore the splice
    auto-edit, plus whatever auto-select needs beyond re-picking.)*
13. **[long-horizon, demoted] Tier 5 — LinearDesign-class joint codon +
    secondary-structure optimization.** Previously ranked first; the 2026
    verification sweep does not support that. Its headline number is antibody titre
    in n = 6 mice against a vendor optimizer (HEK293 protein was **2.9×**), it has
    never replicated, an independent 2026 mammalian bake-off found
    stability-prioritizing strategies **reduced** expression, and MFE minimization
    is substantially a **GC-maximization proxy** — so folding ΔG and GC are not
    independent frontier axes. The part with uncontested evidence (position-
    dependent structure) is cheap and belongs in item 3. Keep this as an option,
    with the counter-evidence attached; never claim an MFE-minimized design
    improves *expression*, only stability. See the survey §4.

**Two open questions that block sequencing and are the maintainer's to answer:**

- **Does the supplied 5′UTR / backbone enter the run manifest?** Hashing it makes a
  context-dependent result reproducible from its stamp (invariant #9). But the
  backbone is the user's IP, a hash is a fingerprint, and manifests are meant to be
  shareable. This is a genuine conflict between two of BT4's own principles —
  provenance completeness vs. nothing leaves the machine — and needs a deliberate
  call, not a default. (Related and non-negotiable: the ASSP cross-check must be
  hard-blocked from ever transmitting backbone bytes; the simplest honest answer is
  to disable it whenever a backbone is loaded.)
- **What is BT4's regime — IVT mRNA therapeutic, vector transgene, or both?** It is
  upstream of almost everything: the 5′-structure sign, the CpG direction (deplete
  for durable AAV, elevate for a vaccine), the length ceilings, and whether
  "expression" means titre, half-life or immunogenicity all differ. A tool that
  silently serves both will give one of them wrong advice.

---

## Working agreements (operational; the durable rationale is in CLAUDE.md)

- **Honesty is structural** (§1, §5, §10.6). Never present an unenforced
  constraint, an unvalidated number, or a heuristic result as if it were real. A
  new model ships `calibrated=False` until it passes its gate. Never fabricate a
  data table — refuse and say why (as tAI / codon-pair bias do).
- **Adding a constraint/objective/model = a new file + a registry/export entry +
  its honesty property test** (`ok_suffix⇔validate` / `delta==score` / a
  calibration gate), never an engine edit (§4). import-linter enforces the layering
  (§3): `domain` imports nothing; `cli`/`service`/`app` import only `api`; heavy
  deps stay lazy behind contracts + optional extras.
- **Single-trunk + CI, auto-merge cadence (§7).** Branch off trunk (the harness
  assigns your branch name — don't reuse a prior one), open a PR ready for review,
  then **enable GitHub auto-merge (squash)** immediately (`enable_pr_auto_merge`) so
  GitHub merges the moment required checks pass — event-driven, no polling, no timed
  merge. Stay webhook-first: the subscription wakes you on CI *failures* and review
  comments; a failing check just blocks the pending auto-merge, so keep the
  drive-to-green posture (diagnose + push a fix, or reply with the blocker; each
  push keeps auto-merge armed).
- **Keep this file and CLAUDE.md current in the same change** (§10.11); update
  `README.md`/`CHANGELOG.md` when a user-facing surface or shipped behavior changes.
- **Sandbox limits (human-only, don't work around):** deleting remote branches,
  pushing git tags, cutting releases, and anything needing the licensed
  splice/expression weights.

## Runbook

Local gate (matches CI; run before pushing):

```bash
python -m ruff check src tests scripts
python -m mypy                      # whole-package; CI's dep-free quality job is authoritative
lint-imports                        # layering contract
QT_QPA_PLATFORM=offscreen python -m pytest tests/ -p no:cacheprovider
```

- **BT4 Studio (headless self-test):** `QT_QPA_PLATFORM=offscreen python -m bt4.app --self-test`
- **App deps not in the base env** — install `pip install -e '.[app]'` and the Qt
  system libs `libegl1 libgl1 libglib2.0-0 libxkbcommon0 libdbus-1-3`.
- **CLI:** `bt4 optimize`, `frontier`, `validate`, `tracks`, `library`, `presets`,
  `build-table` (+ `--check-splice` / `--splice-backend`, incl. `assp`).
- **Optional extras:** `[ilp]` `[fold]` `[ml]` `[service]` `[app]` `[assp]`
  `[splice-pangolin]` `[splice-spliceai]` `[expression-ribonn]` `[dev]`.
- **Rust ext:** built by CI; a byte-identical pure-Python fallback runs when absent.

## RiboNN environment gotchas (learned on real hardware)

Point `$BT4_RIBONN_DIR` at the user's own RiboNN checkout. Its stack needs
`numpy<2` (torch 1.13.1 ABI — note that is BT4's derived *ceiling*; upstream's
`environment.yml` pins numpy to 1.22.4 exactly, so this shorthand understates it and must
not be used as evidence about what else could share the environment) and
`setuptools<82` (its older `pytorch_lightning`
calls `pkg_resources`; upstream issue #10 / PR #11), and the Zenodo `weights.zip`
(record 17258709) extracted **into** a directory literally named `models/` under
`$BT4_RIBONN_DIR` — the zip's root holds `human/` and `mouse/`, so `-C models` is
the correct target — so the hard-coded
`models/<species>/<run_id>/state_dict.pth` path resolves without a symlink.
Scoring requires **non-empty** `utr5`/`utr3` (empty → refused; the UTRs carry most
of RiboNN's signal). Weights are non-commercial — never bundled or CI-run.

Three more, learned from reading upstream rather than from a crash:

- **The licence is an *affiliation* grant, not just a non-commercial one.** Both
  `LICENSE` and `MODEL WEIGHTS LICENSE.txt` grant use "to any person from academic
  research or non-profit organizations". Non-commercial *intent* may not qualify an
  unaffiliated maintainer. Resolve in writing (`patent.gos@sanofi.com`) before
  downloading weights.
- **`max_shift` is a determinism hazard.** RiboNN's `_stochastic_shift` is not gated
  on `self.training` and uses an unseeded `torch.randint`, so a nonzero `max_shift`
  in the shipped `models/<species>/runs.csv` MLflow params makes *inference*
  non-deterministic and breaks invariant #7. Check it before trusting any number.
  (Predict-time config comes from `runs.csv`, **not** `config/conf.yml`.)
- **Windows works but is not upstream-supported.** You never need `make` — its
  `install` target is only `mamba env create -f environment.yml -y`, and BT4 imports
  `src.predict` directly. But the weights folder *must* be named `models` (the
  symlink fallback needs Developer Mode), and `num_workers` must be `0` (now the
  adapter default). RiboNN's one tested environment is Ubuntu 20.04 + one NVIDIA GPU.
  Note `src/predict.py` calls `torch.load` with **no `map_location`**, so if the
  released state dicts carry CUDA tensors they will refuse to load on a CPU-only
  box — that is an upstream property, not a BT4 bug.

Four more, learned while walking a maintainer through a real Windows install (2026-08):

- **Clone into `%USERPROFILE%`, not `C:\`.** The guide's Step 9 `max_shift` check defaults
  to `~/RiboNN/models`, so a `C:\RiboNN` clone makes the one check written to fail loudly
  instead report that it cannot find `runs.csv`.
- **RiboNN and Pangolin should not share an environment** — installing `[splice-pangolin]`
  into RiboNN's env would upgrade torch off RiboNN's pin. The versions, and the fact that
  the `torch>=2.2` floor is **BT4's own** rather than an upstream Pangolin requirement,
  are stated in the guide's Step 5 and deliberately **not** repeated here: an earlier
  draft of this bullet restated them, and the copies drifted apart inside a single commit.
  The repo-side facts worth holding here are that nothing checks `torch.__version__` at
  runtime, and that the SpliceAI half is a *separate, untested* question (its extra is
  TensorFlow-only). A maintainer with a working Pangolin environment builds a second one
  and `pip install -e` BT4 into both — safe, because BT4's core has **zero** dependencies and
  the `expression-ribonn` extra is floors (`torch>=1.13`, `pandas>=1.5`), not pins, so it
  upgrades nothing already present. The consequence to state rather than let someone
  discover: inside the RiboNN environment the wrapped splice CNNs are unavailable and the
  audit falls back to the PWM baseline, so **one run cannot use both**.
- **Miniforge is not a prerequisite;** Miniconda/Anaconda work, and `mamba` is a speed
  convenience. The one real difference is that `conda env create -f environment.yml` names
  no channel, so `defaults` is consulted unless `nodefaults` is added to RiboNN's file —
  an earlier draft claimed "every command names its channel", which is false for the one
  command that matters most.
- **`conda activate` *and* `$BT4_RIBONN_DIR` are both shell state, and agent-driven shells
  lose both.** `conda run` replaces only the first — it inherits the parent environment —
  so a half-fix runs the right Python against an unset variable and fails with "RiboNN
  clone not found", which looks nothing like the cause. Guide Step 7 has the full form.

**The install has now been completed end-to-end on that machine (2026-08), and the guide's
CPU-only recipe did not survive contact with it.** Four measured findings; the working
recipe lives in guide Step 5 **only**, and is deliberately not restated here:

- **`pytorch=1.13.1` does not exist for win-64 on conda-forge** — its earliest build for
  that platform is **2.5.1**, so the guide's CPU-only `mamba create` line was
  `PackagesNotFoundInChannelsError` on Windows from the day it was written. It is fine on
  linux-64.
- **The official `pytorch` channel's win-64 `1.13.1` installs but will not import.** All
  10,672 files verify against their SHA-256, then `import torch` raises
  `OSError: [WinError 182]`. The traceback names **`shm.dll`**, which is a red herring —
  `shm.dll` fails only because it imports `torch_cpu.dll`, and `torch_cpu.dll` is what
  will not load, under every loader search-flag combination —
  with a 200-module import graph in which everything resolves to an x86-64 PE, all 722
  symbols imported from its four non-system dependencies present in their export tables,
  no delay-load or bound-import directory, and an image that maps cleanly as a datafile.
  **The cause was not established**; what is established is that the **pip** build
  (`torch==1.13.1` from `download.pytorch.org/whl/cpu`) works on the same machine.
- **The pip wheel and `numpy=1.22.4` are ABI-incompatible**, and it fails *quietly*: torch
  prints `Failed to initialize NumPy: module compiled against API version 0x10 but this
  version of numpy is 0xf` as a **warning** and runs on with `from_numpy` broken. This is
  the one gotcha above whose shorthand — "needs `numpy<2`" — is actively misleading, since
  the real constraint on the pip path is a *floor*, not only a ceiling. `numpy=1.23.5` is
  the version measured here — NumPy's `0x10` ABI starts at 1.23, and nothing between
  1.22.4 and 1.23.5 was tried.
- **A fresh *conda* solve pairs protobuf 5 with tensorboardX 2.5 and breaks `import
  pytorch_lightning`.** `environment.yml` pins no protobuf, and tensorboardX 2.5's
  generated `_pb2.py` raise `TypeError: Descriptors cannot be created directly` against
  protobuf ≥ 4. RiboNN's own `src/predict.py` imports `pytorch_lightning`, so **nothing**
  runs. conda cannot solve `protobuf<4` in that environment either (libmamba aborts with
  `RuntimeError: bad variant access`); `pip install "protobuf<4"` fixes it. The
  `<=2.5.1` ceiling is **conda-forge's recipe**, not upstream — pytorch-lightning 1.8.5's
  own metadata declares `tensorboardX (>=2.2)` unbounded — so the pip-based recipe needs
  no protobuf pin, verified both ways.

And one that is about tooling rather than RiboNN: **`conda run` refuses a multi-line
`python -c` argument** (`NotImplementedError: Support for scripts where arguments contain
newlines not implemented`, conda 26.5.3), so the stateless form recommended for agents
does not compose with the multi-line snippets the guide is written in. Write them to a
`.py` file.

**What the completed install produced** (guide Steps 8–9, human `HBB` UTRs and CDS from
Ensembl `ENST00000335295`): a score in `RiboNN CLR-residual TE` units with
**`calibrated = False`**, byte-identical across two runs, and `params.max_shift = [0]` for
both the human and mouse weight sets — so the determinism hazard above is **absent from
the released weights**, and `torch.load`'s missing `map_location` did **not** bite on a
CPU-only box. None of this is calibration: it is the adapter proven to run.

The replacement recipe was then rebuilt from scratch in a second environment and scored the
same sequence: **bit-for-bit the same number** as the environment that had been repaired
step by step. So the guide's recipe is the one that was measured, not one reconstructed
from the findings afterwards — and the score does not depend on which of the two ways the
environment was assembled.

## Splice CNN environment gotchas (learned on real hardware)

Point `$BT4_PANGOLIN_MODEL_DIR` / `$BT4_SPLICEAI_MODEL_DIR` at the weights. Neither
variable appeared in any doc before this note.

- **`pip install spliceai` fails on Windows, and it does not matter — but `--no-deps`
  alone is not enough.** It depends on **pysam**, which has no Windows wheels and cannot
  build (its `setup.py` shells out to a configure script). pysam serves SpliceAI's own
  VCF command line and BT4 never uses it: the adapter resolves the weights with
  `importlib.util.find_spec` — which locates the module *without executing it* — and
  loads the `.h5` files with Keras directly. Install both lines:

  ```
  pip install spliceai==1.3.1 --no-deps
  pip install "pandas<2.2" pyfaidx "setuptools<81"
  ```

  The second is **required for the capture step**, not polish. `spliceai/utils.py`
  imports `pandas`, `numpy`, `pyfaidx` and `keras`, and `--no-deps` skips all of them.
  Miss it and the failure is *late and confusing*: `available()` still reports `True`
  (it needs only Keras and the weight files) while `capture_spliceai_panel.py` dies at
  `from spliceai.utils import one_hot_encode`. Measured with those three installed and
  **pysam still absent**: the capture and the gate both run to completion.
- **`pip install pangolin` is a different package entirely** (a probabilistic
  programming language). Pangolin is GitHub-only: `pip install <checkout>`.
- **TensorFlow 2.16+ cannot load SpliceAI's weights**, because `tensorflow.keras` is
  Keras 3 from that release and the weights are 2019 Keras-2 `.h5` graphs. The extra
  pins `tensorflow>=2.6,<2.16`; to use a newer TensorFlow, install `tf_keras` alongside
  it and BT4 prefers that shim automatically (`_ambient_keras_is_v3`).
- **Set `TF_ENABLE_ONEDNN_OPTS=0` for both the capture and the gate.** TensorFlow warns
  on import that oneDNN "may see slightly different numerical results ... from different
  computation orders". Whatever it is set to, it must be the **same on both sides** — a
  capture with it on gated with it off produces a deviation caused by TensorFlow rather
  than by BT4's adapter, which is the one thing the fidelity gate must not confuse.
- **The two model stacks share one environment more easily than expected.** The plan
  assumed they could not coexist (Pangolin wants torch, SpliceAI an old TensorFlow).
  Measured: TF 2.15 forces `numpy<2`, and a CPU torch build runs fine on numpy 1.26.4, so
  a single venv serves both. A plain `python -m venv` is enough on any platform — nothing
  in the splice path needs a conda-only package — and the activated prompt is the check
  that a `pip install` is landing where it is meant to.
  *(Two corrections to this note, 2026-08. It recorded the measured build as
  `torch 2.13.0+cpu`, and **no such release exists** — PyTorch runs 1.13.x, then 2.0
  onward — so the version number cannot be recovered and has been dropped rather than
  guessed; the numpy finding it supports is unaffected, and the `splice-pangolin` extra's
  own floor is `torch>=2.2`. It also claimed "on Windows there is no `conda`", which is
  simply false and now contradicts the RiboNN Windows guidance above in this same file.)*

---

*History (what each session shipped) lives in
[`../CHANGELOG.md`](../CHANGELOG.md) — this file carries only current state and the
forward queue.*
