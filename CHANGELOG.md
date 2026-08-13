# Changelog

All notable changes to BT4 are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it cuts
its first tagged release.

## [Unreleased]

### Added
- **Highly-expressed CAI reference sets — and they are now the default.** A codon
  table's `w = f/f_max` only means something relative to a set of genes, and BT4's
  tables answered the wrong question: they were counted genome-wide, marking the
  codon that is most *common* (a quantity set largely by mutation and GC bias)
  rather than the codon translation *prefers*. `w` as Sharp & Li (1987) defined it
  comes from a reference set of **highly expressed** genes. BT4 now ships both,
  and every table, result, CLI line, app row, and manifest says which one it is.
  - **`scripts/build_highly_expressed_tables.py`** counts each organism's **300
    most abundant proteins**, ranked by **PaxDb v6.1** whole-organism *integrated*
    proteomics (CC BY 4.0 — a weighted consensus over many published studies, in
    ppm) and joined to the **same release-pinned Ensembl CDS** the genome-wide
    tables use, under the **same** filtering rules — so the two tables differ only
    in *which genes* they count, never in how a gene is read.
  - **The join uses no third-party mapping layer.** PaxDb protein IDs resolve
    against the pinned release's own peptide FASTA; an identifier that resolves to
    two genes is dropped as ambiguous and *counted separately* from one the
    annotation simply lacks. All three sources (abundances, peptide FASTA, CDS)
    are SHA-256-pinned, and `--verify` rebuilds and diffs the committed bytes and
    sidecars. Each sidecar also carries a digest of the ranked 300-gene roster, so
    a third party can prove they reproduced the same reference set.
  - **Organelle-encoded genes are excluded.** Mitochondria and plastids translate
    with a different genetic code and their own tRNA pool, and are never a BT4
    design target; they are negligible genome-wide but abundant enough to reach a
    top-300 list.
  - **N = 300 is evidence, not taste.** It is the smallest size on a tested grid
    (50…2000) at which *every* bundled organism observes all 64 codons, so no
    shipped table needed smoothing — an invented number in a reference table is
    exactly what BT4 refuses to ship. Below it, yeast alone leaves `CGA`/`CGG`
    unobserved; far above it the reference set dilutes back into the genome-wide
    answer (at N=2000 yeast and mouse agree with genome-wide at every amino acid).
  - **What changed in the output.** The most-used codon moves for 11 amino acids
    in *C. elegans*, 8 in *E. coli* (`TTT`→`TTC`, `CGC`→`CGT`, `GGC`→`GGT`,
    `ATT`→`ATC`, `CAT`→`CAC`, `GTG`→`GTT`, `AGC`→`TCT`, `TAT`→`TAC` — the classic
    *E. coli* optimal codons), 7 in zebrafish, 5 in yeast, and 2 in human
    (`AGA`→`CGC` for Arg, `AGC`→`TCC` for Ser, plus the preferred stop moving
    `TGA`→`TAA`). Two golden snapshots moved accordingly; both are regenerated and
    the reason is recorded beside them.
  - **External ground truth, not just self-consistency** (§8): the tables
    reproduce the classic *E. coli* and *S. cerevisiae* optimal codons, and codon
    bias is **stronger** in the highly-expressed set than genome-wide in all eight
    organisms — largest in yeast (+0.18) and *Drosophila* (+0.17), smallest in rat
    (+0.03) and human (+0.05), the ordering dos Reis et al. (2004) predict from
    translational selection being weak in large vertebrate genomes.
  - ***A. thaliana* deliberately has none.** PaxDb identifies its proteins by
    UniProt accession, which the pinned Ensembl Plants annotation does not carry,
    so a join would need an unpinned external mapping. BT4 ships no table rather
    than one built on a guess; the genome-wide table stays its only, honestly
    labeled, option — and asking for `highly_expressed` there **raises** instead of
    silently substituting the other table.
  - **Surfaces:** `--reference-set` on `bt4 optimize`/`frontier`/`validate`/
    `tracks`; `bt4 organisms` now prints each organism's default and available
    sets; a **Reference set** picker in BT4 Studio that repopulates from the engine
    per organism (and disables itself, with a reason, when only one exists);
    `reference_set` on the service's optimize request and in `/organisms`;
    `OptimizeConfig.reference_set`; and `api.available_reference_sets` /
    `api.default_reference_set`.
  - **Honest scope, unchanged.** A highly-expressed reference set makes CAI a
    better-founded proxy, not a validated expression predictor — Welch et al.
    (PLoS ONE 2009) found an *E. coli* variant built by maximizing exactly this
    quantity expressed at a fraction of alternatives. CAI stays one axis of the
    objective vector (§10.7).

### Changed
- **Results now report their CAI reference set.** `result.audit` carries
  `codon_reference_set`, the CLI prints it beside the CAI, and BT4 Studio shows it
  as its own metrics row — a CAI of 1.0 against highly-expressed counts and one
  against genome-wide counts are different claims about the same sequence, and the
  label travels with the number rather than living in a control the user may have
  changed since the run. The genome-wide sidecars gained a matching
  `reference_set` stamp, and a mis-filed sidecar (one claiming a reference set
  other than the table it sits beside) is now a load error rather than a silent
  mislabel.
- **`available_organisms()` now recognizes an organism table by shape, not by a
  suffix blocklist.** It accepts exactly `<organism>.tsv`; every other TSV in the
  data directory belongs to another axis of the same organism. The old rule
  excluded `.trna.tsv` by name and so failed *open* — the new reference-set tables
  would have appeared as organisms called `homo_sapiens.highly_expressed`.
- **All nine organism tables are now recounted from release-pinned public CDS
  sets.** Human, *E. coli* and *S. cerevisiae* were the last hand-typed
  "representative published values" with `cds_count: null` — which made BT4's
  **default organism** its least checkable table. They now go through the same
  `scripts/build_organism_tables.py` pipeline as the other six, so every bundled
  number is a real codon count with a source URL, the source file's own SHA-256,
  assembly, genebuild, and a per-filter drop tally (CLAUDE.md §8).
  - **The delivered sequences did not change.** Across a four-protein panel × the
    three rebuilt organisms, the optimized DNA is **byte-identical** and CAI moves
    by at most +0.0003. CAI normalizes within each synonymous group, and the
    most-preferred codon per amino acid is unchanged in all three — so the
    published values were qualitatively right, and this is a provenance upgrade
    rather than a change in behavior. Textbook biases still hold (*E. coli* CTG
    for Leu, yeast AGA for Arg), now as counted facts rather than asserted ones.
  - **Alternate haplotypes and patch scaffolds are excluded.** Ensembl ships
    alternate haplotypes and patch scaffolds alongside the primary assembly, each
    with its own gene IDs — so per-gene de-duplication does not collapse them and
    they enter a table as duplicate copies of genes already counted. Two species
    are affected: **human** (11,513 records — seven alternate MHC/HLA haplotypes,
    the chr19 KIR and LRC haplotypes, and the patch scaffolds) and **zebrafish**
    (6,029 records on `ALT_CTG*` contigs, which were **15.6% of that species'
    genes**; 98% of the symbolled ones duplicate a primary-chromosome gene).
    Both counts are stamped. Genuine unplaced contigs are kept, and the region
    label match covers `primary_assembly:`, which rat and *Drosophila* use — a
    naive "chromosomes only" rule would have discarded every record for those two.
  - **A blocklist alone was not enough, so the build now audits itself.** An
    earlier pattern that looked complete missed human's `HG*_NOVEL_TEST` patches
    (12 genes, 9 of them second copies of chr11 olfactory receptors) and every
    zebrafish `ALT_CTG*` contig. The builder now separately counts anything it
    *kept* whose region name still looks alternate/patch-like and stamps it as
    `kept_suspicious_region`; a test requires it to be zero for every organism.
    The next unknown naming variant fails in CI instead of quietly inflating a
    shipped table. Removing the duplicates changed **no** amino acid's preferred
    codon in any organism, so delivered sequences are unaffected.
  - `test_every_bundled_organism_is_recounted` now fails if any future organism
    reintroduces an undocumented table.

### Fixed
- **Unknown-enzyme suggestions no longer read as equivalent substitutes.** The
  near-miss list is a fuzzy match on the *name* with no notion of recognition
  sequence, so a suggestion usually cuts something entirely different (`NotI` is
  `GCGGCCGC`, `NcoI` is `CCATGG`). A bare list invited a user to accept a
  substitute that does not ban the site they care about — and the run would then
  report proven-optimal with zero hard violations while their real site sat in
  the delivered sequence, which is precisely the §1 failure BT4 exists to
  prevent. Every suggestion now carries **its own site** and is labelled a
  spelling match, so non-equivalence is *visible* rather than asserted, and the
  message points at banning the sequence directly instead. New public
  `enzyme_suggestions()` / `unknown_enzyme_message()` (re-exported through
  `bt4.api`) keep the CLI, API and BT4 Studio from drifting into telling the user
  different things about the same miss.
- **The catalog build's selection tally now closes.** It called itself "auditable
  rather than a black box" while length / unknown-site / multi-site rejections
  were folded silently into the gap between two counters — which is how an
  arbitrary length cap dropped SfiI without anything in the sidecar showing a
  loss. Each rejection reason is counted and stamped, and a test asserts the
  numbers add up. `rejected_site_length` is now `0`, so the stamp *proves* the
  sanity bound drops no real enzyme rather than asking to be trusted.
- **Enzyme catalog: corrected selection, and its hash now enters the manifest.**
  Found by an adversarial review of the catalog change below, and fixed before
  the numbers were relied on:
  - **The documented selection rule now matches the code.** The build described
    itself as "Type II only" while selecting any REBASE type starting with `R`.
    The composition is now stated precisely — Type II (`R2`), bifunctional/IIG
    (`RM2`), modification-dependent/IIM (`R2*`), and one Type III — in both the
    script and the sidecar. The modification-dependent entries are **kept
    deliberately**, with the reason recorded: an initial attempt to exclude them
    as "cannot cut unmethylated synthetic DNA" would have deleted **DpnI**, whose
    `GATC` avoidance is mainstream precisely *because* a plasmid from a dam+
    strain is Dam-methylated and is cut by it.
  - **SfiI and eleven others are no longer silently dropped.** An arbitrary
    12-base cap discarded `SfiI` (`GGCCNNNNNGGCC`) and `XcmI`
    (`CCANNNNNNNNNTGG`) with no signal. The bound is now a parse-sanity guard at
    20, justified by the longest real site in the source.
  - **Type IIB duplicates resolve to one site.** Enzymes like `AjuI`/`AloI`/`PsrI`
    appear as two REBASE records leading with opposite strands. The builder now
    verifies the two are reverse complements and keeps the first, instead of
    last-wins (which made the shipped site depend on file order).
  - **`enzyme_catalog_sha256` enters the run manifest** when restriction enzymes
    are active (invariant #9). The sites moved from Python source — covered by
    the manifest's `git_commit` — into a data file that no hash covered, so a
    swapped catalog would have changed the constraint while leaving byte-identical
    provenance: exactly the BT3 anti-pattern §10.10.
  - **Re-derivability stated honestly.** REBASE publishes only a *moving*
    current-release URL with no versioned permalink, so `source_url` pins nothing
    on its own; the digest is the pin, and the sidecar now says so rather than
    implying a stable link.

### Added
- **REBASE-derived restriction-enzyme catalog (17 → 584 enzymes)**
  (`bt4.constraints.restriction`, `scripts/build_enzyme_catalog.py`) — the
  catalog was seventeen hand-typed pairs described only as "textbook-correct",
  with no source, version, or way to check them. It is now **derived from a
  version-pinned REBASE release** and held as content-hashed package data, the
  same discipline the recounted codon tables get (CLAUDE.md §8):
  - **Selection is documented and auditable:** commercially available Type II
    enzymes (REBASE `ET`/`CR`) with a single fully-specified IUPAC recognition
    site of 4–12 bases. The provenance sidecar records the REBASE version, URL,
    the source file's own SHA-256, the stage-by-stage selection tally, and the
    shipped TSV's digest — so a third party re-derives and re-verifies it
    (`python scripts/build_enzyme_catalog.py --verify`).
  - **Type IIS enzymes included** — BsaI, BsmBI, BbsI, SapI, Esp3I, AarI, the
    Golden Gate workhorses. REBASE lists an asymmetric site once per strand; the
    builder *verifies* the second entry is the reverse complement of the first
    rather than assuming it, and takes one (BT4 bans both strands anyway).
  - **Verified against the values it replaces:** all 17 previously shipped
    enzymes resolve to byte-identical sites, cross-validating old and new.
  - **Isoschizomers kept** (`KpnI` and `Acc65I` both → `GGTACC`) so a user can
    name the enzyme they actually own.
  - New public `resolve_enzyme()` (case-insensitive, and on a miss offers the
    closest names instead of dumping the catalog) and `enzyme_provenance()`, re-exported
    through `bt4.api`. `ENZYMES` is now read-only shipped data.
  - **BT4 Studio** gains a searchable enzyme field: a completer that matches the
    *last* comma-separated token and substitutes it back, leaving earlier entries
    intact — a stock completer matches the whole line and breaks after the first
    enzyme.
  - **Honest scope, stated in the data:** BT4 models the recognition *sequence*
    only — not cut position, star activity, methylation sensitivity, or buffer.
    Some real entries are highly degenerate (`MspJI` is `CNNR`); banning one in a
    CDS can be genuinely unsatisfiable, and BT4 raises `InfeasibleError` naming
    `restriction_site` rather than returning a sequence that still contains it.
    A regression test pins that either/or across degenerate and ordinary sites.
- **Six new organisms, recounted from release-pinned public CDS sets** — mouse,
  rat, zebrafish, *Drosophila*, *C. elegans*, and *Arabidopsis*, taking BT4 from
  three selectable organisms to nine (Phase 5 organism breadth, CLAUDE.md §8/§9).
  This closes a real gap rather than padding a list: all six already shipped
  authentic GtRNAdb **tRNA** tables, but tAI is only offered for an organism you
  can *select*, and selection needs a codon-usage table — so six of the eight
  bundled tRNA tables were **unreachable**. A regression test
  (`test_every_trna_table_has_a_selectable_organism`) keeps that from recurring.
  - **Every number is a real count, never a curated summary.** New
    `scripts/build_organism_tables.py` downloads a **release-pinned** Ensembl /
    Ensembl Plants CDS FASTA (release 116 / plants 63 — pinned, not `current`,
    which moves), filters to complete unambiguous in-frame coding sequences
    (ACGT-only, length 3N, ATG start, terminal stop, no internal stop), takes
    **one representative CDS per gene** (the longest; ties broken by transcript id
    so the pick is deterministic) so codon usage is not weighted by how finely a
    gene happens to be annotated, and counts codons with BT4's own
    `count_codons`. The terminal stop is counted, since BT4 chooses the stop it
    appends.
  - **Re-derivable by a third party.** Each provenance sidecar now carries the
    source URL, the **downloaded file's own SHA-256**, assembly, **genebuild**
    (the gene annotation the CDS models come from — *not* the same thing as the
    assembly: Arabidopsis CDS are Araport11 models on the TAIR10 assembly, and the
    fly/worm models are FlyBase/WormBase; recording only the assembly would
    misattribute the very sequences that were counted), database release,
    total codons counted, the full per-filter drop tally, and the rebuild command
    — alongside the existing content hash of the TSV itself. `--verify` rebuilds
    into a temp directory and diffs against the committed bytes; all six verify
    byte-identically.
  - **Refuses rather than fabricates.** The builder aborts if a CDS set yields no
    valid sequences, or if any of the 64 codons goes unobserved — it will not
    smooth an invented number into a shipped table.
  - **Checked against external ground truth (§8), not just self-consistency.**
    The new tables reproduce independently-published facts: GC3 orders
    *Drosophila* (0.63) > zebrafish (0.54) > *Arabidopsis* (0.42) > *C. elegans*
    (0.40); mouse (0.573) and rat (0.578), counted from separate CDS sets, land
    within 1.5 points of human (0.587); preferred Leu is CTG in the GC3-rich
    genomes and CTT in the AT-rich ones; preferred stop is TGA in the GC-richer
    genomes and TAA in the AT-rich ones. Gene counts match published
    protein-coding counts (e.g. *C. elegans* 19,928; mouse 21,571).
- `write_table` gained `build` / `note` / `extra` parameters so a recount can
  describe what it actually did and attach a re-derivation trail. Reserved
  provenance keys cannot be shadowed by `extra`, so a sidecar can never disagree
  with itself.
- **BT4 Studio surfaces the engine-ready backends, gains library mode, and gets
  its Phase-4 polish** (`bt4.app`) — the two models that already existed behind
  `bt4.api` but had no UI are now wired in, plus the sampler and the accessibility
  work called for in CLAUDE.md §6.6. All of it is pure plumbing over the stable
  API (no engine change, no calibration claim):
  - **RiboNN in the Candidates tab.** An opt-in *Expression head* group (toggle,
    species, and the fixed 5'/3' UTR context the model requires) routes a
    `RiboNNExpressionModel` into `api.candidates`. The toggle is enabled **only**
    when `available_expression_backends()` reports the user's own checkout and
    weights actually resolve, so it is never a dead control, and it explains what
    is missing otherwise. Missing/non-DNA UTRs are refused *before* the run starts
    rather than raising mid-flight. RiboNN stays `calibrated=False`, so the banner
    still reads **discovery order, not a ranking** and the solver's pick stays
    delivered (§10.6).
  - **Validate with ASSP.** The one control that leaves the machine. It asks for
    consent first (naming the service and what is sent), runs
    `api.splice_crosscheck` on a background thread, and renders the report led by
    its tags — *network-derived, UNCALIBRATED, advisory, **not** part of the run
    manifest and never exported* — with the localized sites in a table. An outage
    degrades to a labeled "unavailable" banner and never fails a run (§10.15). The
    panel is cleared whenever the delivered sequence changes, so one sequence's
    splice sites can never be shown beside another's, and an export is
    byte-identical whether or not a cross-check ran (regression-tested).
  - **Library (sampled) tab.** `api.library` with members / temperature / seed
    controls, a per-member table, the selected member's sequence with its
    violation highlights, and a multi-record FASTA export whose every record is
    named `sampled`. The banner leads with **sampled, not optimized** — the
    `SAMPLED` certificate colours the badge directly, so it cannot drift from the
    claim the engine made — and reports measured diversity (distinct count, mean
    pairwise difference).
  - **Phase-4 polish.** A File/Run/View/Help menu bar with standard shortcuts
    makes every action keyboard-reachable; **View → System / Light / Dark**
    switches theme at runtime (restyling the stylesheet, both plots, the badges,
    and the sequence viewers' violation bands from the still-live results, via a
    new `SequenceViewer.set_dark`); tab order covers every new control and each
    carries an accessible name plus an explanatory tooltip.
  - **One source of truth for run gating.** All four flows (optimize, rank+audit,
    cross-check, library) share a `_wire_thread` helper and a `_update_run_buttons`
    gate driven by explicit running-flags rather than thread references — so a
    missed reference clear can no longer strand a control (the previous
    optimize-then-rank stuck-button class of bug is now structurally impossible).
  - New shared `_EngineWorker` base in `bt4.app.worker` (signal trio + the
    never-raise contract) with `CrossCheckWorker` and `LibraryWorker` alongside
    the existing two.

  Found by an adversarial review of the above and fixed in the same change (each
  with a regression test that fails without its fix):
  - **A late cross-check could be attributed to the wrong sequence.** A report
    describes exactly one sequence and carries it, so `_on_crosscheck_finished`
    now compares `report.dna` to the live delivered DNA and *discards* a report
    whose design changed while it ran, instead of rendering it. The panel-clearing
    rule covered only the other ordering.
  - **Menu shortcuts bypassed the single-flow gate.** `Ctrl+R` during an in-flight
    cross-check started a second engine flow, because only the buttons were gated.
    The Run actions are now gated alongside them, and each `_start_*` refuses via a
    shared `_busy()` check — so the invariant lives in the code path, not only in a
    greyed-out control.
  - **A second library draw stranded the first draw's sequence on screen.**
    Repopulating the table in place leaves the selection intact, so re-selecting
    row 0 emitted nothing; the member viewer is now repainted explicitly.
  - **Untrusted service text could rewrite the honesty banner.** ASSP's own error
    text was interpolated unescaped into a RichText label — markup that could hide
    the very "network-derived / UNCALIBRATED / advisory" labels marking it. All
    externally-derived text is now HTML-escaped.
  - **Closing mid-run destroyed a running `QThread`** (pre-existing, but this change
    triples the number of flows that can be in flight). `closeEvent` now cancels
    what is cancelable and gives each live thread a bounded chance to finish.
- **Public expression-backend registry** (`bt4.biomodels.expression.available_backends`
  / `resolve_backend`, re-exported as `bt4.api.available_expression_backends` /
  `resolve_expression_backend`) — the mirror of the splice resolver, so a frontend
  selects an expression head by name through the stable API instead of importing
  `biomodels` across a layer (§3, §10.9). `available_backends()` never raises and
  lists `"ribonn"` only when it can genuinely run; resolution is lazy (no torch
  import, no weight load) and confers **no** calibration.
- **Opt-in, out-of-loop ASSP splice cross-check** (`bt4.api.splice_crosscheck` /
  `bt4.pipeline.run_splice_crosscheck`, `bt4.biomodels.splice.AsspSplicePredictor`)
  — a **network** validator that runs the online ASSP service (Alternative Splice
  Site Predictor, Wang & Marín 2006) over an already-delivered sequence behind the
  existing `SplicePredictor` contract, closing the last non-human-gated gap in the
  splice subsystem (CLAUDE.md §6, §10.15). BT3's fatal splice bug was scraping this
  exact service **in the optimizer's inner loop as its only splice path**; BT4
  inverts every property of that mistake, structurally:
  - **Opt-in and out-of-the-inner-loop.** Requested explicitly by name and gated
    behind the `bt4[assp]` extra (httpx, lazily imported); it runs only as a final
    audit / validation pass on the delivered sequence, never per optimizer move, and
    is **never** returned by `splice.default()` or `available_splice_backends()`.
  - **Never blocking.** Rate-limited with exponential backoff and cached by
    sequence hash; if the service is unreachable or returns a garbled body the raw
    predictor raises an `AsspError`, but `run_splice_crosscheck` catches it and
    reports "unavailable" — a cross-check outage can never fail an optimization. The
    same graceful path covers a wrapped CNN's missing deps.
  - **Network-derived and non-reproducible.** `network_derived` is `True` and
    `calibrated` is `False`; ASSP numbers are excluded from the
    reproducible-from-manifest guarantee and reported as a separate advisory section
    (the CLI prints them to **stderr**, never into the stdout FASTA/JSON artifact or
    a `Result` manifest).
  - **Wired through the CLI** — `bt4 validate --splice-backend assp` and `bt4
    optimize --check-splice assp` (both flags also accept `pwm` / `pangolin` /
    `spliceai` for an offline or installed-CNN cross-check).
  - **CI never makes a live call.** The adapter is driven from committed **offline
    fixtures** (`tests/fixtures/assp/`, `FixtureAsspTransport`, selected via
    `$BT4_ASSP_FIXTURE_DIR`). Honest caveat: the live wire format is *unverified
    against the service* (unreachable during development), so the fixtures are
    *synthetic ASSP-format reports*, not real captures — the same "no bundled panel
    ships" posture as the wrapped CNNs.

### Fixed
- **`MinMaxTerm` is now scale-invariant — `minmax_weight` finally means the same
  thing on every organism.** Its `delta` was a raw frequency *difference*
  (`f(codon) - f_avg(aa)`), so its magnitude tracked the codon table's units: mean
  `|delta|` was ~4.5 on the per-thousand hand-curated tables but ~52,000 on a
  raw-count table — a **~11,700x disparity**. The same `minmax_weight` therefore
  meant four orders of magnitude more on one organism than another, and on a
  raw-count table the term silently swamped CAI, GC and every other frontier axis:
  precisely the magic-scalar failure §10.5 exists to prevent. The term now
  normalizes frequencies to a **within-family fraction** first, which is all it
  ever needed (exactly as CAI's `w = f/f_max` needs only ratios); mean `|delta|` is
  ~0.07 on every organism. Within-family preference order is provably unchanged, so
  a `minmax`-only solve picks the same codons — what changes is that the knob is
  now comparable across organisms. **This was already live before the new tables
  shipped:** `bt4 build-table` emits raw counts, so anyone optimizing with
  `minmax_weight` against their own table was affected. Regression-tested for
  scale-invariance, cross-organism comparability, and order preservation.
- **The table builder could count the wrong species.** Source archives were cached
  by bare filename, and Ensembl reuses the *same* filename across releases, so a
  stale cache entry would be counted silently. Each source archive's expected
  SHA-256 is now pinned in the build spec and checked on every run (cache hits
  included); verified by planting one species' archive under another's filename,
  which now aborts. `--verify` also diffs the provenance **sidecars**, not just the
  TSVs, so a sidecar naming the wrong source can no longer pass clean. And
  `write_table` validates `extra` before writing anything, so a rejected call no
  longer leaves a TSV on disk without its sidecar.
- **RiboNN adapter: correct ensemble aggregation and honest empty-UTR guard.** The
  first end-to-end runs against real RiboNN weights surfaced two integration bugs.
  (1) RiboNN returns the ensemble as **one row per cross-validation model**, so a
  single input yields several rows sharing a `tx_id`; the adapter's `set_index`
  realignment then made `float(ordered[tx_id])` operate on a Series and raised
  `TypeError: cannot convert the series to <class 'float'>`. Realignment now groups
  by `tx_id` and averages (the ensemble mean, also averaging over cell types) via a
  new tested helper `_reduce_te_by_tx_id` — a no-op when rows are already unique.
  (2) Scoring with the default **empty** `utr5`/`utr3` crashed deep inside RiboNN's
  data loader (pandas reads an all-empty UTR column as `NaN` and its `.str`
  preprocessing fails); the adapter now refuses up front with a clear message, since
  the UTRs carry most of RiboNN's signal and an empty-UTR score is not meaningful.

### Changed
- **Python 3.10 is now supported** (was 3.11+). `requires-python` is lowered to
  `>=3.10`, the 3.10 classifier is added, ruff/mypy target 3.10, and CI's quality
  matrix now runs 3.10 alongside 3.11–3.13. The pure core uses no 3.11-only
  features, so this is a compatibility widening with no behavior change. It notably
  lets the wrapped **RiboNN** expression backend be installed into the same
  environment as its own dependency stack, whose pinned `torch==1.13.1` ships only
  CPython ≤3.10 wheels.

### Added
- **BT4 Studio "Candidates & splice audit" tab** — step 5 (final) of the
  expression/splice design flow, surfacing `api.candidates` → `api.splice_audit`
  in the desktop app. A background `CandidatesWorker` (mirroring the known-good
  `OptimizeWorker` `QThread` lifecycle) runs both on a worker thread and hands the
  window the candidate set + splice audit in one signal. The tab renders the
  ranked, honestly-labeled candidate table (delivered pick starred; per-member
  source / CAI / GC / expression+units / calibration / hard-violation / **distinct**
  splice-site counts) with two advisory banners: an *uncalibrated* expression head
  is shown as **discovery order, not a ranking** (solver's pick starred, scores
  annotating only; a calibrated head switches to ranked-by-expression), and the
  splice banner leads with **UNCALIBRATED (advisory)** whenever `all_calibrated` is
  `False`, reporting cross-backend agreement and stating the flags localize sites
  heuristically and edit nothing. Every metric is recomputed per candidate from its
  own DNA (invariant #2); an opt-in toggle routes the installed SpliceAI/Pangolin
  CNNs into the audit. The results area is now a `QTabWidget` (Design | Candidates &
  splice audit); the Design tab is unchanged. No Cancel control on this tab (the
  assemble→audit flow is not point-cancelable), and the cross-flow Optimize/Rank
  gating clears the worker-thread reference so neither button can deadlock.
- **Localize-and-flag splice audit** (`bt4.api.splice_audit` /
  `bt4.biomodels.splice.audit_splice`) — step 4 of the expression/splice design
  flow (`docs/DESIGN_expression_splice_flow.md` Stage C). An **out-of-loop,
  advisory** audit that runs the available `SplicePredictor` backends over a step-3
  candidate set to **localize** residual cryptic splice sites (one flag per
  contiguous above-threshold run, at its peak — non-maximal suppression) and attach
  the whole-panel **backend agreement** (pooled rank + sign) as the authoritative
  cross-backend confidence signal — built from the Delta-splicing values the audit
  already computed (a new shared `agreement_from_deltas` helper), so each backend
  scores every sequence **once**, never twice (§7). **It never edits** the sequences — a targeted synonymous auto-edit at flagged loci is a
  deliberately deferred, calibrated-gated future step. Honesty (CLAUDE.md §6/§10.6):
  every shipped backend is `calibrated=False` today, so `all_calibrated` is `False`
  and every `SpliceFlag` carries its **emitting backend's** `calibrated` flag; the
  site `threshold` is a **heuristic display knob** (not a validated cutoff) and the
  PWM baseline's per-position `score` is an uncalibrated **arbitrary-units**
  pseudo-score. Per-flag `added_risk_vs_reference` is **positive = worse** and
  strictly *intra-backend*, kept distinct from the panel-level `delta_splicing`
  (larger = better). Cross-backend `also_flagged_by` is a **raw positional
  co-occurrence** (±`match_window` nt, sized to the backends' anchor offsets),
  explicitly **not** a kind-level agreement (Pangolin reports one combined
  `P(splice)` and so can never disagree on kind — its flags are labelled `"splice"`,
  never donor-specific). New `biomodels/splice/audit.py` (raw-sequence core, imports
  only `domain` + the splice backends) + `pipeline/splice_audit.py` (the
  `CandidateSet` adapter + `available_splice_backends()`, which adds the wrapped
  SpliceAI/Pangolin CNNs when installed). Deterministic (#7). API-level surface (the
  BT4 Studio annotation UI is step 5).
- **Candidate-set assembly + expression rerank** (`bt4.api.candidates` /
  `assemble_and_rank_candidates`) — step 3 of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). Assembles the finalist set an
  expression head ranks: the **Pareto frontier** plus, when a GLOBAL rule is active
  *and* the delivered exact-DP seed actually violates it, a small **deterministic
  library of repeat-refined variants** (distinct refinement seeds over the delivered
  seed). The set is de-duplicated and scored by an `ExpressionPredictor` — in **one
  batched call** when the backend implements the new `BatchExpressionPredictor`
  contract (`score_many`, e.g. RiboNN), else per sequence — and delivered under the
  same **calibrated-gating** honesty rule as `rerank_by_expression`: an uncalibrated
  head (the default placeholder, and the shipped RiboNN adapter) only *annotates* —
  the set stays in **discovery order** (`order_basis="discovery"`) with the
  solver-delivered sequence `chosen` — while a calibrated head reorders by predicted
  expression (`order_basis="expression_rank"`, total order `(score desc, index asc)`)
  and re-picks the top (CLAUDE.md §10.5/§10.6). Hardened for correctness/honesty: the
  **delivered (`chosen`) sequence is invariant to `n`** (uncalibrated, the
  solver-delivered sequence is pinned first in discovery order; calibrated, the
  head's top pick is the top of the top-n keep — the cap is applied *after* scoring
  so a calibrated reranker never loses its best candidate);
  every member is a full `Result` (round-trips, metrics recomputed, certificate,
  residual GLOBAL violations disclosed); variants are labelled `repeat_refined` (the
  *process*, not a guaranteed fix); and de-dup/cap counts, the batch-path flag, and
  the predictor identity (folded into the manifest, invariant #9) are all reported.
  New `BatchExpressionPredictor` Protocol in `bt4.biomodels.expression`; `_refine`
  gains an optional `seed` (default unchanged). API-level surface (UI wiring is
  step 5). No calibration claim — ranking is a reporting no-op until a head is
  calibrated.
- **Strong splice-consensus motif constraint** (`bt4.constraints.SpliceSiteMotifConstraint`,
  `avoid_splice_sites`) — step 2 of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). A new **LOCAL, exact-in-the-trellis**
  hard constraint that forbids the *strong* splice-consensus **donor** (`GTRAGT`,
  the intronic +1..+6 core) and **acceptor** (`YYYYYYNYAGG`, a polypyrimidine tract
  + `NYAG|G`) motifs on the mRNA **sense strand only** (splicing is strand-specific,
  so — unlike restriction/repeat motifs — there is **no** reverse-complement
  banning). It is an honest **structural heuristic**, not a splice model: it reduces
  only the most *obvious* cryptic-splice risk and makes no calibrated claim; the
  wrapped SpliceAI/Pangolin CNNs do the real audit out of loop (CLAUDE.md §6,
  §10.6). It **never** bans the ubiquitous bare `GT`/`AG` (governing rule 3). The
  default patterns (Shapiro & Senapathy 1987; Zhang 1998; human/mammalian
  major-spliceosome only) are deliberately specific (~1/2048 donor, ~1/8192
  acceptor) so the hard veto rarely over-constrains a design, and are configurable
  via `donor_motifs`/`acceptor_motifs`. `ok_suffix⇔validate` and `context_len`
  sufficiency (5 donor / 10 acceptor) are property-tested (invariant #3). Wired
  through `OptimizeConfig`, the `bt4` CLI (`--avoid-splice-sites`), the `service`
  schema, and BT4 Studio (a checkbox with an explanatory tooltip); off by default.
- **Batched RiboNN scoring** (`RiboNNExpressionModel.score_many` /
  `.delta_logte_many`) — the first step of the expression/splice design flow
  (`docs/DESIGN_expression_splice_flow.md`). RiboNN's cost is dominated by fixed
  *per-invocation* overhead (weight hashing + model load + its DataLoader worker
  spawn), so scoring a whole candidate set one sequence at a time paid that cost
  N times. The new public batch methods route the entire set through the existing
  batched `_predict_te` path (one temporary TSV, one `predict` invocation — RiboNN's
  `top_k`-model ensemble runs inside that single call), so scoring a Pareto frontier
  costs roughly the wall-clock of scoring a single sequence.
  `delta_logte_many` additionally scores the shared **reference once** (appended to
  the batch), not once per design. Both preserve per-input validation (valid DNA,
  length-3N ending in a stop codon, non-empty `utr5`/`utr3`) and the `tx_id`
  realignment; results come back **in input order**. `score_sequence` /
  `delta_logte` now delegate to the batch methods (single source of truth). A
  `num_workers=0` DataLoader path was investigated and **deliberately left out**:
  RiboNN's `predict_using_nested_cross_validation_models` exposes no worker-count
  parameter, so requesting 0 workers would mean patching RiboNN internals (against
  the "wrap, never reimplement" contract), and batching already amortizes the
  one-time worker spawn across the set. `calibrated` stays **`False`** — no
  calibration claim. Tested without torch / pandas / the RiboNN checkout (batch
  ordering, ensemble averaging per `tx_id`, reference-scored-once, and the
  empty-UTR / bad-CDS guards still firing).
- **Wrapped RiboNN expression backend** (`bt4.biomodels.expression.RiboNNExpressionModel`)
  — the Phase-4 learned expression head behind the `ExpressionPredictor` contract
  (CLAUDE.md §6/§9). It runs the published **RiboNN** translation-efficiency CNN
  (Zheng, Persyn, Wang et al., *Nat Biotechnol* 2025; Sanofi / Cenik Lab)
  inference-only as an out-of-loop frontier reranker. **License:** RiboNN's code
  and weights are each **Sanofi non-commercial** (academic/non-commercial only) —
  compatible with BT4's open-source non-commercial scope and, like SpliceAI's
  CC BY-NC weights, **never bundled**: the adapter drives the user's own RiboNN
  clone (lazily importing the repo's `src`, pointed at via `$BT4_RIBONN_DIR`) and
  their Zenodo weights. Every weight it loads is verified against a bundled
  180-entry SHA-256 manifest (`data/ribonn_sha256.json`, 90 human + 90 mouse —
  public content hashes only) **before** `torch.load`. The score is in RiboNN's
  native **CLR-residual TE** units (never exponentiated); `delta_logte(designed,
  reference)` gives the UTR-fixed, CDS-attributable Δ (negative = a CDS change
  predicted to *reduce* expression), analogous to Pangolin's `delta_splicing`.
  Ships **`calibrated=False`** (`default()` still returns `NullExpressionModel`):
  faithful reproduction is not calibration for BT4's CDS-variant regime, so
  promotion requires a passing `verify_expression_gate` on a regime-matched panel
  (human-only, data-gated). New `bt4[expression-ribonn]` extra (torch + pandas),
  lazily imported so `import bt4` stays light.
- **Model-agnostic expression acceptance-gate harness**
  (`bt4.biomodels.expression.gate`) — the honest gate a learned expression head
  must pass to earn `calibrated=True` (CLAUDE.md §6/§8/§10.6, Phase 4). For a
  log-TE regression head it reports **Spearman** (primary), **Pearson**, **R²**,
  and **split-conformal coverage** at a target level (default 90%), evaluated on a
  **group-disjoint split** (homology cluster / chromosome) so no group leaks
  across calibration and test — the distribution-shift-aware check that a head
  validated only on natural-gene TE has *not* earned calibration for BT4's
  CDS-variant regime. `passed` requires both the Spearman threshold **and**
  conformal coverage near target (point accuracy *and* honest uncertainty). The
  gate never flips anything: thresholds are inputs set at gate time, and the
  neutral `NullExpressionModel` provably cannot pass (its zero-variance scores
  give Spearman 0). New `ExpressionEvalCase` / `ExpressionGateReport` and a
  `run_expression_gate(predictor, samples)` wrapper. Fully dependency-free and
  tested without torch or any real model, mirroring how the splice
  fidelity/attestation machinery shipped before a calibrated backend.
- **Shared dependency-free statistics** (`bt4.biomodels._stats`) — `pearson`,
  `spearman` (moved from `splice.agreement`, which now re-exports them), plus
  `r2_score`, `conformal_quantile` (finite-sample split-conformal), and
  `empirical_coverage`. Single well-tested home for the estimators the splice
  agreement report and the expression gate both use.
- **License-clean splice fidelity-attestation layer**
  (`bt4.biomodels.splice.attestation`) — the honest promotion path for the wrapped
  Pangolin / SpliceAI backends (CLAUDE.md §6, §10). A `FidelityAttestation` records
  **only** a passing integration-fidelity gate's derived scalars (`passed`,
  `max_abs_deviation`, `n_cases`, `tolerance`) plus the public pinned weight
  SHA-256s and the tool version — **never** a `FidelityCase` raw per-position score
  (those are the license-encumbered model outputs). The shape is enforced
  structurally (`_ALLOWED_FIELDS` + an honesty test asserting no raw-score field is
  serializable), and `from_dict` refuses any unexpected key. `attest_backend`
  refuses to record a failing or too-loose gate; `verified_predictor(predictor,
  attestation)` is the single seam that flips a backend to `calibrated=True`, and
  only when the attestation passed, clears the `MAX_ATTESTATION_TOLERANCE` floor,
  and its weight SHAs exactly match the adapter's `PINNED_WEIGHT_SHA256` (a
  refusal, never a silent downgrade). A deterministic, timestamp-free
  `content_hash` makes an attestation a provenance-manifest stamp. This layers the
  committed-record / private-execution / user-opt-in / baseline-fallback options;
  no attestation ships, so `default()` still returns the honest PWM baseline. Both
  Pangolin (GPL) and SpliceAI (CC BY-NC) are eligible to certify under BT4's
  open-source, non-commercial scope.

## [0.4.0] - 2026-08-04

First tagged release since 0.3.1, capturing the Phase 1 performance and Phase 3
refinement/splice wave: the full Rust trellis port, richer refinement moves, the
wrapped SpliceAI splice backend, and the last Phase 2 budget item.

### Added
- **Full Rust trellis port** (`bt4_native.trellis_solve`) — the exact-DP inner
  loop of `bt4.optimize.exact_dp.solve_exact` now runs in Rust (Phase 1, CLAUDE.md
  §7), following the existing native-primitive pattern: a PyO3 `#[pyfunction]` with
  a byte-identical pure-Python twin (`bt4._accel._py_trellis_solve`) and a
  Hypothesis equivalence test pinning the two. The DP is callback-driven, so Rust
  never calls back into Python: a **regime gate** restricts the native path to
  position-independent objectives (no `POSITIONAL` term — `CpbTerm` was made
  context-based so PAIRWISE terms stay position-independent), Python **precomputes**
  the reachable-context transition graph and the pre-summed per-transition deltas
  (fixing the float summation order, so the lexicographic tie-break is bit-for-bit
  identical), and the layer DP runs in Rust; it **falls back to the pure-Python
  DP** whenever the regime does not hold, the extension is absent, or a
  context-count cap is exceeded. A single solve is not accelerated (the Python
  precompute costs ~a whole pure DP), so `run_optimize` stays on the pure path; the
  win is the **Pareto frontier**, which builds the transition graph once and reuses
  it across every scalarization point (only the cheap deltas recomputed) with the
  DP in Rust — a measured ~2.7–5.5x `run_frontier` speedup with **byte-identical**
  DNA, objective scalars, and certificates.
- **Block/segment moves + parallel tempering in the SA refinement engine**
  (`bt4.optimize.anneal_refine`, Phase 3 — CLAUDE.md §7, §9). The engine gained
  four opt-in knobs: `block_size` / `block_prob` (coordinated multi-position
  synonymous swaps) and `replicas` / `temps` / `swap_every` (a parallel-tempering
  replica ladder with standard replica-exchange Metropolis swaps). These widen the
  refinement's *reach* so it can cross a barrier that only clears when several
  codons move **together** — a dispersed max-repeat or out-of-frame uORF the
  single-codon chain could leave in place — **without weakening invariant #5**:
  block candidates pass the same local (union-of-windows `ok_suffix`) and global
  (whole-sequence recount) feasibility gates, every replica gates against its own
  current hard-violation count, every visited configuration keeps a global count
  `<=` the seed's, and the delivered result is ranked lower-global-count-first then
  higher-score. All four default off, and with them off the engine reproduces the
  prior single-chain trajectory **byte-for-byte** (invariant #7). Block moves
  always full-`score` re-score (never `delta_score`), since summing per-position
  deltas is only valid for additive disjoint-context terms. The honest
  **feasibility floor** is preserved: a repeat pinned to synonymously-immovable
  bases (Met `ATG` / Trp `TGG`) is unremovable by any move and is still reported as
  a residual, never claimed clean. New Hypothesis tests pin the never-raise-global
  guarantee under block+tempering, determinism/round-trip with replicas and blocks,
  the default-knobs no-op, and the immovable-repeat feasibility floor.
- **Wrapped SpliceAI splice backend** (`bt4.biomodels.splice.SpliceAiSplicePredictor`)
  — the second *wrapped published* splice CNN behind the `SplicePredictor`
  contract, the cross-check to Pangolin (Phase 3, CLAUDE.md §6). It runs the
  published **SpliceAI** model (Jaganathan et al. 2019) inference-only, and its
  3-way per-position softmax (null/acceptor/donor) maps *cleanly* onto
  `SpliceResult.acceptor` and `.donor` (both populated, unlike Pangolin's single
  combined track). **License (verified): SpliceAI code is PolyForm Strict 1.0.0
  and its weights are CC BY-NC 4.0 (noncommercial) — even more restrictive than
  Pangolin's GPL, so no-bundle is mandatory**; the adapter lazily imports the
  user's own installed `spliceai` package + weights, SHA-256 hash-pinning them
  (verified before load). Ships **`calibrated=False`** (`verify_spliceai_fidelity`
  is the gate; no reference panel bundled), so `default()` still returns the PWM
  baseline. With both CNNs installed, the agreement harness now compares two real,
  independently-trained splice models (no harness code change needed — it already
  compares at the pooled-Δ level). New `bt4[splice-spliceai]` extra
  (TensorFlow), lazily imported so `import bt4` stays light.
- **Wrapped Pangolin splice backend** (`bt4.biomodels.splice.PangolinSplicePredictor`)
  — the first *wrapped published* splice model behind the existing
  `SplicePredictor` contract (Phase 3, CLAUDE.md §6). It runs the already-validated
  **Pangolin** CNN (Zeng & Li 2022) as an inference-only backend, feeding its
  per-nucleotide `P(splice)` into the shipped Δsplicing / top-k-log-odds framing.
  **License-clean:** Pangolin is **GPL-3.0** (the earlier roadmap's "MIT" was
  wrong), so — exactly as BT4 wraps GPL ViennaRNA — the adapter **lazily imports
  the user's own installed `pangolin` package and weights and bundles neither**;
  BT4 stays MIT. Weights are **SHA-256 hash-pinned** (the published v1.0.2 digests)
  and verified *before* they are unpickled, keeping runs
  reproducible-from-manifest. The adapter reproduces upstream Pangolin's scores
  **bit-for-bit** yet ships **`calibrated=False`** (no reference panel is bundled;
  `verify_pangolin_fidelity` is the promotion gate), so `default()` keeps returning
  the honest PWM baseline. Heavy deps behind the new `bt4[splice-pangolin]` extra,
  lazily imported so `import bt4` stays light.
- **Two-backend splice agreement harness** — `bt4.biomodels.splice.backend_agreement`
  reports each available backend's Δsplicing ranking, pairwise **Spearman rank
  agreement**, and sign agreement across candidates (the first-class uncertainty
  signal of CLAUDE.md §6/§8); it reports, it does not judge. Exposed as the
  standalone runner `scripts/compare_splice_backends.py` (`--fasta`, `--json`),
  which degrades to the baseline alone — and says so — when neither CNN backend
  (Pangolin nor SpliceAI) is installed.
- **CpG / UpA whole-sequence count budget** (`dinuc_budget` + `dinuc_min` /
  `dinuc_max`; CLI `--cpg-min/--cpg-max` and `--upa-min/--upa-max`) — the last
  Phase 2 item. A dinucleotide count does not decompose per-codon (a 2-mer
  straddles the codon boundary), so the amount-bucketed budget DP
  (`bt4.optimize.lagrangian`) now takes a **context-aware** per-codon amount
  (`bt4.objectives.dinucleotide.dinucleotide_amount`) attributing each occurrence
  to the codon holding its END base, with a new `budget_context` folded into the
  trellis state so a straddling count stays exact. Enforced by the same **exact
  bucketed DP** as the GC budget, with a `proven_optimal` certificate and every
  local constraint still honored. Mutually exclusive with the GC budget, and (like
  it) not combinable with `refine` / `max_repeat_length` / `avoid_uorf`. Wired
  through `OptimizeConfig`, the CLI, and the `service` request schema.
- **Library / degenerate-design mode (opens Phase 5).** `api.library(protein,
  config, n, *, seed, temperature)` and `bt4 library PROTEIN --n N` sample a
  *library* of coding sequences by drawing from each residue's synonymous-codon
  distribution (organism usage frequencies raised to `1/temperature`), keeping
  only codons that pass every LOCAL constraint. This is an honest **stochastic
  sampler, not an optimizer**: every member round-trips and carries metrics
  recomputed from its own DNA, the library is fully deterministic from its seed,
  and each result carries the new **`OptimalityStatus.SAMPLED`** certificate,
  which makes no optimality or expression claim. GLOBAL constraints
  (`max_repeat_length`, `avoid_uorf`) are not enforced during sampling but are
  validated and any residual violation reported honestly per member. New modules
  `bt4.optimize.sample` (deterministic constrained sampler, `domain`-only) and
  `bt4.pipeline.library` (`LibraryResult` + `run_library`).
- **Two more `bt4_native` hot-loop primitives** (`max_gc_run`, `longest_repeat`),
  each with a byte-for-byte pure-Python fallback in `bt4._accel` and a Hypothesis
  equivalence property test that pins the Rust and Python paths together (and, for
  `longest_repeat`, cross-checks `longest_repeat(seq) > m` iff
  `MaxRepeatConstraint(m).validate(seq)` flags a hard violation). This is honest
  incremental native acceleration — **not** a full trellis inner-loop port, which
  still remains (CLAUDE.md §7, §9 Phase 1).

### Changed
- **`GcRunConstraint.ok_suffix` now calls the (optionally Rust-accelerated)
  `bt4._accel.max_gc_run`** on its bounded trailing window, with no change to
  observable behavior (the pure-Python fallback is the same scan as before). The
  `longest_repeat` primitive is added and cross-checked against
  `MaxRepeatConstraint`, but is **deliberately not** placed on the per-SA-move
  `MaxRepeatConstraint.validate` hot path: the whole-sequence longest-repeat is
  O(n²), which is *slower* than the constraint's existing O(n·k) k-mer scan when
  the native extension is absent — so wiring it there would regress the common
  pure-Python path (CLAUDE.md §7, "everything incremental"). Every existing
  `ok_suffix ⇔ validate` and constraint test passes unchanged.

## [0.3.1] - 2026-08-01

BT4 Studio first-run polish: the desktop app now guides a non-technical user
through mistakes with plain-language messages instead of raw Python errors, and
never leaves a stale result behind a failed run.

### Added
- **Cancel button + live progress** for BT4 Studio. The frontier sweep now
  reports per-point progress (`solving frontier point 3 of 9`) and can be stopped
  mid-run; cancelling returns the partial frontier computed so far. `api.frontier`
  / `run_frontier` gained optional `on_progress` and `should_cancel` hooks.
- A one-time **warning before optimizing a very long protein** (it may take a
  while, and the run is cancelable).
- `bt4.api` now re-exports `InfeasibleError`, `validate_protein`, `AMINO_ACIDS`,
  and `available_tai_organisms` so frontends can validate input and translate
  failures without reaching past the API layer.

### Changed
- **Plain-language input handling in BT4 Studio.** Pasting a FASTA record strips
  its header automatically; an empty box, a trailing `*` stop, or non-amino-acid
  characters get a clear, specific message (not a Python `repr`); restriction-
  enzyme names are matched case-insensitively and unknown ones list the valid
  catalog; an infeasible constraint set explains which knobs to relax instead of
  saying "no feasible codon". The **tAI** checkbox is now labelled correctly and
  enabled only for organisms that ship a tRNA table.

### Fixed
- **A failed run no longer leaves a stale, exportable result on screen** — the
  results panel (and the delivered result behind Export) is cleared on failure,
  so Export can't silently write the previous sequence.
- `scripts/sensitivity.py` detected tAI availability via the pre-0.3.0 organism-
  list quirk and silently returned `None` for every organism after that quirk was
  fixed; it now uses `api.available_tai_organisms()`.

## [0.3.0] - 2026-08-01

First release with a **downloadable, double-clickable BT4 Studio app** for
Windows / macOS / Linux, plus a wave of Phase 2/3 objectives, constraints, and
solver backends.

### Added
- **5' translation-ramp objective** (`RampTerm`) -- a heuristic that prefers
  slower codons in the first N codons (`ramp_weight` / `ramp_codons`).
- **CpG / dinucleotide objective** (`DinucleotideTerm`) to deplete (stealth) or
  elevate (immunostimulatory) CpG content (`cpg_weight` / `cpg_mode`).
- **Codon-pair bias** (`CpbTerm` + `build_codon_pair_table`): a pairwise objective
  built from a reference CDS set, solved exactly in the trellis via a new
  `objective_context` on the DP (the state now carries the previous codon).
- **OR-Tools CP-SAT backend** (`bt4.optimize.cpsat.solve_cpsat`, `bt4[ilp]`
  extra): solves the additive objective under a global **GC budget** (`gc_min` /
  `gc_max`) with a proven-optimal / gap-bounded certificate. New `ilp` CI job.
- CLI flags for all of the above (`--ramp-weight`, `--cpg-weight`, `--cpg-mode`,
  `--gc-min`, `--gc-max`) and a CpG control in BT4 Studio.

### Changed
- **Idiot-proof, double-clickable app packaging.** The PyInstaller spec now emits
  a *single* file per desktop OS instead of a one-folder zip: a one-file
  `BT4-Studio-Windows.exe`, a one-file `BT4-Studio-Linux-x86_64`, and (on macOS) a
  `.app` that CI wraps in a drag-to-Applications `BT4-Studio-macOS.dmg`. Verified
  end-to-end on Linux: the one-file build launches BT4 Studio and runs its event
  loop. The README's install section is rewritten for non-technical users
  (download-one-file table + how to click past the unsigned-app OS warnings), with
  the from-source/CLI install moved to a "for developers" section.
- **Release pipeline is now re-drivable and self-healing.** `release.yml` accepts
  a `workflow_dispatch` `ref` input to rebuild an existing tag's source and
  idempotently (re)attach the per-OS app + wheel/sdist to its release — the
  honest, non-destructive way to repair a release that has no assets. The publish
  step now also fails loudly instead of publishing an empty, asset-less release.
  See [`packaging/README.md`](packaging/README.md#repairing-a-release).
- **CI now launches the packaged app.** A `bt4-studio --self-test` hook builds the
  main window (loading the bundled data + Qt/pyqtgraph) and exits without the
  event loop; the release workflow runs it against the freshly built bundle on
  each OS, so a bundle that builds but crashes on first launch fails CI instead of
  shipping. The macOS `.app` also now carries its real version in `Info.plist`,
  the codon/tRNA data dir is a regular package (reliable frozen-bundle resource
  loading), and the Windows asset rename/upload no longer depends on a fragile
  cross-shell absolute path. A full non-technical [`docs/INSTALL.md`](docs/INSTALL.md)
  guide was added.

### Fixed
- The only tagged release (`v0.2.0`) had **no downloadable app**: its publish step
  ran the pre-idempotency workflow and `gh release create` failed on "release
  already exists" (the tag/release were made in the UI first), so the built
  bundles never attached. The pipeline is now idempotent and re-drivable, and the
  docs no longer point users at an empty Releases page.
- **`available_organisms()` listed bogus organisms.** It matched every `*.tsv`,
  so the tAI tRNA tables leaked in as `homo_sapiens.trna`, `mus_musculus.trna`,
  and `saccharomyces_cerevisiae.trna` — visible in the app's organism dropdown and
  `bt4 organisms`, and unloadable as codon tables. The tRNA tables are now
  excluded (they remain available via `available_tai_organisms()`).

## [0.2.0] - 2026-07-31

Richer biology and surfaces on top of the exact-DP core.

### Added
- **Restriction-site constraint** (`bt4.constraints.RestrictionSiteConstraint`,
  `available_enzymes`): an IUPAC-aware matcher and a catalog of common enzymes
  (EcoRI, BamHI, NotI, ...), always avoiding each site's reverse complement.
  Wired into `OptimizeConfig.restriction_enzymes`, the CLI (`--enzyme`,
  `bt4 enzymes`), and BT4 Studio.
- **More organisms**: representative *E. coli* K-12 and *S. cerevisiae*
  codon-usage tables (auto-discovered; clearly labeled representative).
- **`bt4 build-table`** and `bt4.io` FASTA parsing: recompute an authentic codon
  table from a user-supplied CDS FASTA (Laplace-smoothed so the result always
  loads), with a content-hashed provenance sidecar.
- **`bt4.service`**: an optional FastAPI HTTP API (`/optimize`, `/frontier`,
  `/validate`, `/organisms`, `/health`) that calls only `bt4.api`.
- **Benchmark harness** (`scripts/benchmark.py`) and a golden/regression test
  suite pinning current optimizer output.

### Fixed
- BT4 Studio frontier plot now shows raw CAI / GC-fraction axis values instead of
  a rescaled "x0.001" SI-prefix label.

## [0.1.0] - 2026-07-31

First tagged release: an honest exact-DP codon optimizer with a CLI and the BT4
Studio desktop app.

### Added
- **Exact codon-trellis DP solver** (`bt4.optimize`) over the true per-constraint
  context, with an explicit `beam` speed knob and a machine-readable
  `OptimalityCertificate` (`proven_optimal` / `beam_truncated`).
- **Objective terms** (`bt4.objectives`): `CaiTerm` (log relative-adaptiveness)
  and `GcProximityTerm`, both additive with `delta == score` property tests.
- **Constraints** (`bt4.constraints`): `HomopolymerConstraint` and
  `ForbiddenMotifConstraint` (with automatic reverse complements), with
  `ok_suffix ⇔ validate` agreement property tests.
- **Pipeline + stable API** (`bt4.pipeline`, `bt4.api`): `optimize()`,
  `frontier()` (a CAI/GC Pareto frontier), and `validate()`, with metrics
  recomputed from the delivered DNA and a content-hashed provenance manifest.
- **`ObjectiveTerm` / `Constraint` protocols** and the `Scope` enum in the pure
  `domain` layer (the shared vocabulary the optimizer speaks).
- **`bt4` CLI**: `optimize`, `validate`, `organisms`, and `--version`.
- **BT4 Studio** (`bt4.app`): a native PySide6 desktop app calling `bt4.api` on a
  background thread — constraint controls, an honest optimality-certificate
  badge, a recomputed-metrics table, an interactive CAI/GC frontier plot, a
  sequence viewer, and FASTA/JSON export. Offline; nothing leaves the machine.
- **IO** (`bt4.io`): FASTA and versioned, deterministic JSON export.
- **Packaging & distribution**: a `packaging` extra, a PyInstaller spec
  (`packaging/bt4-studio.spec`) that builds a standalone BT4 Studio bundle, and a
  `Release` workflow that publishes the sdist + wheel and per-OS app bundles on a
  version tag.
- **Community health**: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`,
  and GitHub issue / pull-request templates; a landing-page `README.md` with a
  screenshot of BT4 Studio.

### Fixed
- **Wheel/sdist builds** were broken by a `pyproject.toml` `force-include` that
  double-added the codon data files (`pip install .` failed with "a second file
  is being added to the wheel archive at the same path"); replaced with
  `artifacts` so the data and `py.typed` marker ship exactly once.
- **Two import-linter layering violations** (`optimize → constraints`,
  `objectives → biomodels`) that surfaced once `bt4.app` existed — resolved by
  lifting the protocols into `domain` and decoupling `CaiTerm` from the codon
  table, keeping every pure layer importing only `domain`.

### Notes
- Richer objectives (tAI, codon-pair, 5′ ramp), ILP / relaxation backends, and
  the validated splice / folding / expression models are on the roadmap and are
  **not** yet shipped — see [`CLAUDE.md`](./CLAUDE.md) §9.
