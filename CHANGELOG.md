# Changelog

All notable changes to BT4 are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it cuts
its first tagged release.

## [Unreleased]

## [0.5.0] - 2026-08-21

The largest release so far, and the one that starts moving BT4 from *a coding
sequence in isolation* toward *a coding sequence in a construct*: a
`ConstructContext` carrying the 5'UTR and vector backbone, `SeededConstraint` so
every LOCAL rule is junction-correct, one shared junction folding window, and
`api.audit_construct` auditing the assembled construct including restriction-site
uniqueness. Around it, annotated **GenBank I/O** puts residual violations on the
map a user actually opens; the organism set grows to **twelve**, every table
recounted from release-pinned public CDS; and CAI finally **declares its
reference set**, with the highly-expressed set (PaxDb top-300) the default
wherever one exists — 8 of the 12 organisms, and the other four ship none rather
than a guessed one — so `w = 1` marks the codon translation *prefers*, not merely
the one that is most common.

The calibration apparatus also grew up. Both wrapped splice CNNs — **Pangolin**
and **SpliceAI** — passed their integration-fidelity gates at `max_abs_deviation`
exactly 0.0 over the same 18-case panel, and ship with committed, license-clean
attestations honored behind an explicit opt-in. The statistical-calibration gates
that a fidelity pass does *not* provide landed too: the splice gate (Part B) with
its strict panel format, position-convention verification and four permanent
baselines, and the expression gate with `within_group` scoring, conformal
coverage, a cluster bootstrap and five baselines a head must beat. The wrapped
**RiboNN** translation-efficiency head landed with batched scoring and an honest
scope, and BT4 Studio surfaced it — plus the opt-in ASSP cross-check, a
**Candidates & splice audit** tab, **Library (sampled)** mode, menus, shortcuts
and runtime theming.

Running through the release is a series of measurements that were allowed to
contradict this project's own claims, and did. BT4's splice risk pooling was
**structurally mute in BT4's own regime** — one uncalibrated background constant
flooring the opt-in path (6 of 93 designed sequences above it) while saturating
the shipped one (93 of 93) — so a floored zero is now distinguishable from a
measured one, with a background-free ranking statistic beside it and the
threshold deliberately left where it was. The adapters' `N`-padding was shown
**not to be neutral**, making an `N`-padded splice number a lower bound rather
than an estimate. The models were shown **not to be blind** on designed CDS (a
planted donor lifts the local peak ~11x at exactly the predicted anchor) while
the detection floor was shown to be **high** enough that intermediate-strength
sites clear nothing. And the candidate ranking was shown **reliable** —
Eρ² = 0.90–0.96, so "the ordering is ensemble noise" is excluded without needing a
single label — while the *delivered pick* was shown unstable, its argmax moving with
fold or tissue in 2 of 3 proteins. Reliability is not validity: nothing here measured
whether the ranking is **right**, and a reproducible ordering can be reproducibly
wrong. That question needs labels BT4 does not have.

**No model became statistically calibrated in this release.** RiboNN and every
splice backend report `calibrated=False` by default, and `default()` still returns
the PWM baseline. The one exception is exactly scoped and must not be read as more:
under the explicit `BT4_SPLICE_USE_ATTESTED=1` opt-in a wrapped CNN reports
`calibrated=True`, and that asserts **integration fidelity** — BT4's adapter
reproduces the published model bit-for-bit — not that the model's scores are
calibrated probabilities for designed coding sequence, which is a separate and
still-unmet gate. The splice operating point remains a convention rather than
evidence, and the frontier reranker still refuses to re-pick the delivered sequence
from an uncalibrated score. What shipped is the
apparatus that could *earn* a calibration claim, and the surfaces that keep an
uncalibrated number from being read as one.

> The expression-promotion seam below and its two follow-ups (#131, #132, #133) landed
> on `main` *after* the version bump and originally sat under `[Unreleased]`. Because
> the `v0.5.0` tag had not been pushed yet, they ship in this release — so they are
> documented here rather than left describing an unreleased state they are not in.
> `[Unreleased]` is therefore empty at the tag: everything on `main` at that commit is
> accounted for below. Tag **this** commit rather than whatever `main` reaches later,
> so the section stays exactly true of what it labels.

### Added
- **The expression head can finally be promoted -- and BT4 Studio can use it.**
  `verified_predictor` has existed since the attestation apparatus landed, but **nothing
  in `src/` called it**: a maintainer could run the acceptance gate, earn a claim and
  commit a record, and every user would still get an uncalibrated head, a candidate set
  in discovery order, and the solver's pick delivered. The gate, the record format and
  the calibrated-gating on the consuming side were all built; the seam between them was
  missing. It now exists, opt-in and scope-bound:
  - **`bt4.biomodels.expression.attestations`** (mirroring the splice side's
    `promote_if_attested`): `BT4_EXPRESSION_USE_ATTESTED=1`, or an explicit
    `use_attested=True` on `api.resolve_expression_backend`, honours a resolvable
    attestation. **Nothing auto-promotes**, `default()` still returns the neutral
    placeholder, and **no attestation is bundled** -- none has been earned, and shipping
    one that had not would be the fabricated artifact CLAUDE.md §10.6 forbids.
  - **A maintainer's own attestation is first-class.** An expression attestation is
    earned against a *measured panel*, frequently unpublished, so
    `$BT4_EXPRESSION_ATTESTATION` points at a local record and needs no commit. A
    mis-pointed path **refuses** rather than falling back, because a typo would otherwise
    look exactly like "no attestation".
  - **BT4 Studio surfaces it.** A per-run toggle on the Candidates tab (passed
    explicitly, never by mutating process env), disabled with a tooltip naming what is
    missing when nothing resolves; the attestation's **scope on the page** (species /
    cell types / readout / `top_k` / UTR contexts / panel hash), not just in a tooltip;
    the head **pinned to that scope** while it is honoured; and the banner flipping from
    *discovery order, NOT a ranking* to a ranking **with its scope named**. A head the
    attestation does not cover stops the run with the mismatch named, rather than being
    silently downgraded to uncalibrated.
  - **The manifest records which claim authorized a calibrated pick** (invariant #9): a
    promoted head carries the attestation's content hash, and `pipeline/candidates.py` /
    `pipeline/rerank.py` fold it into the stamp, so two runs steered by different
    attestations cannot share one.

### Fixed
- **The expression attestation's scope was decorative, and is now the run's.** Three
  measured holes, all closed before any real gate is ever run -- the gate is a one-shot
  pre-registered procedure, so the record has to be complete and unfakeable when it is
  made, not afterwards:
  - **`attest_expression`'s `species` and `cell_types` were caller-declared free text,
    never checked against how the gate actually scored.** Run the gate averaging all 78
    cell types, then declare `cell_types=("HEK293T",)`, and the record -- and every later
    check -- accepted the lie. The scope is now **derived** from the run
    (`GateComparison.scope`); anything a caller declares is a **cross-check that refuses
    on mismatch**, never an override. Where the panel's own `species` / `cell_type` /
    `readout` columns declare the same fact it is checked against the panel bytes too,
    and `verified_against_panel` records exactly which fields got that second check --
    so a verified scope is distinguishable from a merely-declared one.
  - **The record never named the model that was scored.** `attest_expression` took
    `backend` as free text and wrote it straight through, so a gate run against the
    neutral placeholder -- or any other head -- could be filed as a RiboNN result, and
    `verified_predictor` would then promote a real RiboNN head against it, because it
    only ever compared the *label the caller typed*. It is now derived from the `name` of
    the predictor the gate actually constructed, and a run whose scores were **handed in**
    rather than computed (`scoring_source != "gate"`) cannot be attested at all: nothing
    ties supplied numbers to the model the record names, and no after-the-fact check
    recovers it. The maintainer path that needed `head_scores` -- gate once, attest from a
    second run -- is served by `--attest`, which attests from the same comparison.
  - **`verified_predictor` bound neither `top_k` nor the UTR context**, so an
    attestation earned at `top_k=5` under one transcript context would promote a head
    configured with a different ensemble size or different UTRs entirely. Both are now
    bound. `batch_size` / `num_workers` deliberately are **not**: RiboNN pads to a fixed
    width and does not shuffle when predicting, so neither can change a score, and
    binding them would be false precision.
  - **The gate's JSON record omitted `cell_types`, `top_k` and the UTRs**, so a finished
    run was not reconstructable from its own output. `scripts/run_expression_gate.py` and
    `bt4 expression-gate` now emit and print the full `GateScope`, including whether the
    scores were computed by the gate or supplied by the caller (`scoring_source`) -- the
    step at which the link between the named backend and the numbers stops being
    mechanical.
  - Two of these are now caught **before** the scoring pass: `run_panel_gate` refuses a
    head configured to average every cell type against a panel that declares one, or a
    species the panel contradicts. Leaving `--cell-type` off used to run cleanly to a
    wrong verdict with no error and no warning.
- **A corrupt attestation file took BT4 Studio down at startup.**
  `load_expression_attestation` let `json.JSONDecodeError` escape, and every caller that
  handles "this attestation is unusable" catches `ExpressionAttestationError` -- so a
  hand-edited record raised out of window construction instead of greying out one
  checkbox. Unparseable JSON is now wrapped, and refused rather than read as absent.
- **`--attest` truncated the run's own record when it refused.** The report is
  documented as `> gate_result.json` and was printed *after* the attestation step, so a
  refusal returned early and left a zero-byte file -- losing the result over a failure of
  the optional step that follows it. The report is now printed first, unconditionally.
- **A panel mixing two assays could be filed under one readout.** The check was
  membership, so naming either one passed while rows measured by the other were still
  scored. A multi-readout panel is now refused: gate one assay at a time.
- **An explicit `use_attested=True` that could not be fulfilled returned an uncalibrated
  head silently.** It now raises for an attestable head when nothing resolves -- a
  per-call argument is a request about that call, so answering it with an uncalibrated
  head is the failure this layer exists to prevent. A standing `$BT4_EXPRESSION_USE_ATTESTED`
  stays a no-op, because a forgotten export is a preference, not a request.
- **`scripts/ribonn_sensitivity.py` honoured a standing promotion opt-in**, unlike the
  gate it feeds. These checks decide whether a panel is worth acquiring; a prior
  attestation must not colour them.
- **A `schema_version` 1 record was refused with a bare field list.** It now says which
  fields schema 2 added and why they cannot be filled in afterwards.
- **Studio downgraded silently if the attestation stopped resolving mid-session.** The
  toggle was read at window-open; if the record was deleted or corrupted before Run, the
  head was built uncalibrated with no warning -- answering "give me a calibrated ranking"
  with an uncalibrated one, which is the failure the layer exists to prevent. It now
  refuses the run and names what happened.
- **`promote_if_attested` gated on a class name where `verified_predictor` uses
  `isinstance`**, so a *subclassed* head would be skipped here and accepted there --
  a silent downgrade in the one place that must never produce one.
- **`content_hash` moved when a record's JSON lists were reordered**, which a content
  hash must not do; `from_dict` now sorts the same tuples `attest_expression` writes.
- **The docs claimed the environment variable alone made `api.candidates` rank.** It does
  not, and should not: the opt-in governs *promotion*, not *selection*, so a script that
  never asked for RiboNN cannot have its behaviour changed by an exported variable. The
  head must be handed to `api.candidates(..., predictor=head)`.
- **`bt4 expression-gate` ignored `--top-k`** (it had no such flag and always used 5,
  while the script did forward it). The flag now exists and is forwarded; since `top_k`
  is part of the scope an attestation binds, a promotion refuses a differently-sized
  ensemble.

### Fixed
- **The packaged BT4 Studio app did not start at all, and only launching it said so.**
  Cutting this release meant building the bundle and opening it, which is how a
  release-blocking defect surfaced that every from-source gate had passed over:
  `packaging/bt4-studio.spec` collected `**/*.tsv` and `**/*.provenance.json`, and two
  kinds of data file added since have neither shape. The frozen app died before its
  first window with `FileNotFoundError: .../bt4/biomodels/expression/data/ribonn_sha256.json`
  -- RiboNN's public weight-hash pin is read at *import* time, and `bt4.api` imports
  `bt4.biomodels.expression`, so nothing about the failure was optional or lazy. The
  committed splice fidelity attestations
  (`biomodels/splice/data/{pangolin,spliceai}.attestation.json`) were missing from the
  bundle for the same reason, and would have failed **quieter and worse**: an absent
  attestation file legitimately reads as "no attestation ships for this backend", so a
  packaged BT4 would have disclaimed two gates that had in fact passed.
  - The spec now collects `**/*.json`. A pattern that matches the *kind* of file cannot
    drift behind the tree the way an enumeration of today's files did.
  - `tests/test_bundle_spec.py` pins it: every non-Python file under `src/bt4` must be
    matched by the spec's own `includes` patterns, resolved exactly the way PyInstaller
    resolves them (`Path.glob`), with a companion test asserting the check is not
    vacuous. This defect was invisible to the whole suite because the suite runs from
    a source tree where the files are simply *there*.
- **`bt4-studio --self-test` proved the app opened, not that it worked.** It built
  `StudioWindow` and exited, so a bundle that opens and then cannot *design* would have
  passed the release workflow's only per-OS check on the real artifact and shipped,
  surfacing on the user's first click. It now designs the same protein twice inside the
  frozen bundle -- once plainly, once under a restriction rule -- and **requires the
  rule to have changed the answer**.
  - That last requirement is the point, and the first version of this check did not
    have it: it asked only whether the delivered DNA lacked an EcoRI site, and the
    unguarded optimum for its protein never contained one, so it passed identically
    whether or not the restriction rule was applied at all. A check that cannot fail is
    not a check. The pair now used (a 74-residue albumin leader, SacII) is live -- the
    unguarded solve carries `CCGCGG`, the guarded one does not -- and the self-test
    re-establishes that on every run, failing loudly with "gone vacuous" rather than
    quietly proving nothing if a future codon table changes the unguarded answer.
  - The same review corrected the justification originally written here: the REBASE
    catalog is **not** loaded lazily. `bt4.constraints.restriction` reads it at *import*
    time, so window construction already covered it, and citing it as "a data file
    nothing loads until a run" was wrong. What the engine check actually adds is proof
    that the frozen app can *solve* and that a rule the user configures reaches the
    delivered sequence.
  - It stays silent on success: only `cli` prints (CLAUDE.md section 3).
- **A version bump falsely invalidated the committed splice attestations.**
  `test_attestation_records_the_producing_version` asserted
  `attestation.bt4_version == bt4.__version__`, which is strictly stronger than the
  provenance it documents: the fidelity gates ran on a maintainer machine holding the
  licensed weights, at the version current then, and cutting a new BT4 does not re-run
  them. So the assertion left only two ways to release -- re-run a gate that needs
  weights no CI can hold, or rewrite the one field whose job is to record *which build
  measured this*. The check now pins what it says: the recorded version parses and is
  not from the future. Promotion was never affected and is unchanged -- `verified_predictor`
  gates on the backend, `passed`, the tolerance floor and an exact weight-SHA match, never
  on the version -- so a release neither re-certifies nor silently invalidates an
  attestation. **What this does not do is tie a promoted `calibrated=True` to the adapter
  code that earned it** -- an attestation pins the *weights* by SHA-256, not the adapter,
  so a future change to the wrapper would still promote under an attestation captured
  against the old one. The version equality did not really cover that either (it forced a
  re-run only when someone happened to bump the version, which is not when the adapter
  changes), so nothing was lost; but the gap is real and is recorded here rather than
  papered over by a check that looked like it addressed it. Two *synthetic*-attestation
  fixtures that hardcoded `"0.4.0"`
  (`tests/test_splice_attestation.py`, `tests/test_calibration_flag_leaks.py`) now read
  `bt4.__version__` instead, which `docs/DESIGN_splice_cnn_calibration.md` had already
  asked for "so it does not rot at the next bump" -- this was that bump.
- **The released Linux app was not standalone, and the fix belongs in the build, not the
  docs.** Opening the packaged app on a real X display (offscreen proves nothing about
  the X path) aborted it with *"Could not load the Qt platform plugin xcb"* and nothing
  else. The cause is that **PyInstaller embeds the shared libraries present on the build
  machine**, and both workflows installed only the five that make a build succeed and an
  offscreen `--self-test` pass — so every `libxcb-*` the real `xcb` plugin needs was
  simply absent from the shipped binary, and the app worked only on machines that
  happened to have them.
  - Both Linux build jobs now install the full set, so the download carries its own
    copies. Measured in both directions, which is what pins the cause rather than the
    correlation: a bundle built **with** those packages present embeds 20 `libxcb-*` /
    `libxkbcommon-*` libraries and still launches and designs after every one of them is
    *removed* from the machine running it; a bundle built **without** them embeds none
    and dies at startup on a machine that has a working display. What a full build needs
    from the system drops to base graphics — libc, `libEGL`/`libGL`, `libdrm`, `libxcb1`
    — which any desktop has.
  - `docs/INSTALL.md`'s Linux troubleshooting now says that, and leads with the cause it
    cannot fix by installing anything: no display at all (plain `ssh`, a container, bare
    WSL). The twenty-package workaround is kept, scoped to **releases before v0.5.0**,
    whose binaries genuinely do expect the system to supply them.
  - An earlier draft of this entry claimed the list was "every X11 library the app links
    against, verified by launching it on a bare system". Neither was true: the list was a
    subset, and no bare-system check had been run. It is withdrawn rather than softened.
- **A GC budget crashed the packaged app with `ModuleNotFoundError: No module named
  'ortools'`.** A GC-count budget routes to CP-SAT only when *no* local rule and no
  pairwise term is active (otherwise the pure-Python Lagrangian backend takes it, which
  is why this never showed in normal use). CP-SAT lives behind the optional `bt4[ilp]`
  extra and the BT4 Studio bundle deliberately does not carry it — so a user of the
  shipped app who set GC bounds with every local rule switched off got a raw Python
  traceback out of a control the app offers them, inside a frozen bundle where they can
  install nothing. `_solve_with_gc_budget` now falls back to the Lagrangian backend when
  OR-Tools cannot be imported; it solves that case exactly, so the feature works instead
  of failing. The fallback fires only on a genuinely absent OR-Tools, and the certificate
  reported is whichever backend actually ran.
- **The packaged app was only ever built and opened *after* a tag was pushed.** The
  per-OS build-and-launch check lived in `release.yml` alone — i.e. it could only tell
  you the bundle was broken once the release was already being cut. `ci.yml` gains a
  `bundle` job that builds the app on Linux, macOS and Windows and opens it on every PR,
  so this whole class fails at review time. `tests/test_bundle_spec.py` covers the
  cheap half (data-file coverage) in the fast job; only freezing and launching covers
  the rest.

- **BT4 Studio silently dropped a preset field it could not apply.** `_on_preset_chosen`
  guards against exactly this -- a preset that sets something with no visible control is
  named in the status bar -- but `refine` was registered with a `lambda v: None` setter,
  which kept it *out* of the unmapped list. So choosing the **IVT mRNA** preset, whose
  published spec asks for refinement, produced a design without it and said nothing. The
  no-op entry is gone (a field with no control now lands in the list the guard reports),
  and the wording changed from "Not shown in this panel" to "**Not applied** - this panel
  has no control for", because the run is built from the controls and an unmapped field
  is dropped entirely, not merely hidden.
- **The packaged app could not say which version it was.** The About box described what
  BT4 Studio does but never named a release, so a user reporting behaviour gave no way to
  tell which build they were on -- and 0.5.0 changes what the engine delivers by default
  (the CAI reference set). It now shows `BT4 Studio <version>`.
- **`--self-test` now also requires the bundled splice attestations to load.** The
  packaging defect's quiet half was a missing attestation reading as "never gated"; a
  failure whose signature is *silence* is only caught by asking for the file by name.
- **Nothing checked that a release tag matches the version inside the package.** A tag is
  a label a human types: `git tag v0.5.1` on a tree that still says `0.5.0` would publish
  assets whose filenames, release page, `bt4 --version` and run manifests all disagree,
  and it cannot be fixed afterwards without deleting a published release. `release.yml`
  now refuses the mismatch before building anything.
- **Four honesty defects the release audit found in claims BT4 makes about itself.**
  - `CLAUDE.md` section 1 still opened by saying BT4 optimizes a coding sequence "in
    isolation: no field carries the 5'UTR, the vector backbone, or any sequence outside
    the CDS" — contradicted by `ConstructContext`, which ships, and by this file's own
    preamble. Rewritten to say what is true now *and* to keep the limit that did not
    change: context is optional, the bare-CDS path still `N`-pads the splice CNNs, and
    "expression-relevant objectives" is still the strongest honest claim.
  - Both wrapped CNN adapters' docstrings asserted "no reference panel ships, so the
    shipped adapter is **always** `calibrated is False`" — false since the attestations
    landed: under `BT4_SPLICE_USE_ATTESTED=1` a shipped adapter reports `True`. Corrected
    to name the opt-in and to keep the distinction it turns on (integration fidelity, not
    statistical calibration). The reference *panel* still does not ship, and that part
    was right.
  - `docs/INSTALL.md` promised users "BT4 Studio is 100% local and offline. Nothing you
    paste in is ever uploaded anywhere" — contradicted by the app's own ASSP button, and
    by the About box, which says so. Now states the exception, that it asks first, and
    that nothing else in the app ever transmits.
  - `docs/NEXT_SESSION.md` queue item 10 was still titled "SpliceAI still human-gated"
    while its own body and the status board record that SpliceAI passed.

### Fixed
- **The plain-English RiboNN guide, corrected by an adversarial pass over its own claims.**
  A four-lens review (executable accuracy / overclaim audit / instructional quality /
  independent statistical re-derivation) found real defects in the version that shipped,
  now fixed:
  - **`tar -xf weights.zip -C models` fails on Linux.** It was lifted from the runbook's
    Windows `cmd.exe` block into a `bash` block; GNU `tar` cannot read a zip archive, so
    the reader hit `This does not look like a tar archive` at the download step that gates
    all of Parts 2-6. Now `unzip -q weights.zip -d models`, with
    `python -m zipfile -e weights.zip models` as the dependency-free fallback. Both
    verified; so was the failure.
  - **It stated the exact misreading of RiboNN's attribution that `CLAUDE.md` §9 warns
    against** -- "most of RiboNN's signal comes from the UTRs", dropping the load-bearing
    *per-nucleotide* qualifier. Per nucleotide the split is 67/31/2; **length-integrated it
    is 22/73/5, so the CDS is the majority.** Both numbers now stated.
  - **The panel-size Monte Carlo was under-powered.** 92 rows is **51.1% +- 2.5%**
    (1,500 trials), not the 44% first reported from 200 -- a coin flip, not majority-fail.
    Table corrected throughout, ~500 rows added as what actually buys ~94%, and the
    "agrees with the research doc's ~102 rows" claim withdrawn: that figure counts only
    test-fold noise and ignores the calibration fold's own uncertainty in the conformal
    width, so it is a lower bound rather than a matching computation.
  - **The weights are ~200 MB, not ~3 GB** -- an invented figure that appeared in no
    source.
  - **The Ranaghan panel was misdescribed** as "31 different ways by real commercial
    codon-optimization tools". It is 1 native CDS + **three anonymized algorithms run ten
    times each**, and the panel's own `.LICENSE.md` explicitly forbids mapping them to
    named commercial tools.
  - **Two beginner hard-stops:** four `.fa` UTR files were referenced as if they already
    existed, with no guidance on obtaining them (RiboNN refuses to run without real UTRs);
    and the decision point gating an afternoon against a research grant gave no numbers.
    Both fixed -- a "getting your sequences" step with an Ensembl walkthrough, and explicit
    pre-committed thresholds (`within_over_between` >= 0.2, `median_abs_gc3_spearman`
    <= 0.7) labelled as BT4's own pre-commitment, not a standard.
  - Plus: the licence is purpose **and** affiliation (not affiliation instead of purpose);
    the `max_shift` check now honours `$BT4_RIBONN_WEIGHTS` instead of hard-coding a path
    it had just told the reader they could change; a one-protein panel now carries a
    pre-flight stop (the gate exits 2 with `needs at least two distinct groups`, so the
    outcome row describing it was unreachable); the attestation is **content-hashed, not
    "signed"**; the r2 = 0.17-0.19 figure regained its "0.49-0.50 after fine-tuning"
    clause (a reverse overclaim -- it made a null result look more terminal than the
    evidence supports); the splice side was wrongly described as having "exactly the same
    gap" when its `verified_predictor` *is* called behind an opt-in; three off-by-one step
    cross-references and one dead anchor; and glossary entries for r2, GC3, CAI, tAI and
    `top_k`, each of which had been used as load-bearing content while undefined.

### Added
- **A plain-English walkthrough for the RiboNN calibration job**
  (`docs/GUIDE_ribonn_calibration.md`). `DESIGN_ribonn_calibration.md` is an expert
  runbook; this is the same procedure for someone who can use a terminal but has never
  met conformal prediction — 21 numbered steps, a glossary, a 15-row trap table, and the
  **free, weights-free steps pulled to the front** (the 195-test apparatus and a
  full-pipeline dry run against the `null` backend, neither of which the runbook
  mentions). Every command was checked against the source and the runnable ones were
  executed; the test counts (90 / 195) and the report field names are measured, not
  transcribed.

### Fixed
- **The runbook claimed a fidelity result BT4 does not have.** Its job table marked job 1
  *"✅ done — the adapter reproduces upstream bit-for-bit"*. There is **no** RiboNN
  fidelity gate, capture script or attestation anywhere in the tree — "bit-for-bit" is the
  *splice* side's vocabulary, and `expression/attestation.py` explicitly contrasts itself
  with it. What is established is that real end-to-end runs produce numbers, which is why
  the runbook's own Stage 1.2 (fold semantics) exists at all. Corrected in place; §10.6
  applies to BT4's own docs as much as to its models.

### Measured
- **The runbook's own panel-size floor fails a good head more often than it passes one.**
  Stage 2 says "≥ 4 proteins and ≥ ~90 rows". The gate splits by whole protein at 50%, so
  90 rows leaves ~45 to judge coverage on — against a 90% ± 5% band. Simulated against the
  real `verify_expression_gate` with a head at median within-protein ρ = 0.96: **44%** pass
  at 92 rows, 82% at 200, 79% at 200/10 proteins, **99%** at 900. The rank half passed
  **100% of the time at every size** — coverage was the entire failure mode. So a negative
  at the documented floor is not evidence about RiboNN, and Stage 5 maps it to "RiboNN does
  not do BT4's job". Size for **~200 rows**, which is also what the research doc's own
  ±0.05 arithmetic (n ≈ 102, post-split) has always implied.
- **Six further runbook defects, all verified in source** and recorded in the guide's
  Appendix B: the `bt4 expression-gate` shortcut offered at Stage 4 has **no `--json`** (nor
  `--baselines`/`--top-k`/`--batch-size`), so it cannot produce the `gate_result.json` the
  same stage tells you to keep; `harness_ok` is a not-equal test at a **1e-9** floor, not
  the "differ substantially" the prose promises; `$BT4_RIBONN_WEIGHTS` exists and is
  documented in **no** markdown file in the repo; the Stage 0.3 `max_shift` snippet prints
  nothing and exits 0 on a wrong path, indistinguishable from "all clear"; RiboNN's
  `torch.load` carries no `map_location`, so the CPU-only environment the runbook offers may
  simply not load (recorded in `NEXT_SESSION.md`, absent from the runbook); and the CI claim
  that `NullExpressionModel` "cannot pass in either mode" is checked in pooled mode only.

### Fixed
- **Overstated the `N`-padding effect by roughly 2×, via a set mismatch.** `CLAUDE.md` §6,
  `NEXT_SESSION.md` and the review doc all read *"median peak inside the CDS 0.276 →
  0.462"*. That paired arm A's median over the **9-sequence** main set with **0.4622**,
  which is *one protein's* value from the **3-sequence replication** — different sets, and
  the pairing ran in the flattering direction. The same-set figure is **0.2757 → 0.3691,
  +0.093**. The finding survives; its magnitude did not.

### Measured
- **These models do not need 5 kb of flank, and the construct-size objection does not
  disqualify them.** Receptive field is not requirement: detection of a planted textbook
  donor **saturates at ~100–250 nt** of real flank (0.5447 at zero → 0.6224 at 250 nt →
  0.6561 at 5,000 nt; twenty times more context buys ~5%). A 1.5 kb CDS in a 4.7 kb AAV
  payload sits well past saturation. Independently corroborated by published context
  ablations — Jaganathan 2019 (top-k 0.57/0.90/0.93/0.95 at 40/200/1,000/5,000 nt per side)
  and OpenSpliceAI's retrained series, where 80→400 nt is the large gain and 400 nt→10 kb
  is worth only a few percent.
- **…and the earlier flank effect was an extreme-value artifact.** On the same unmodified
  sequences, the model's response at a *given position* is flat to four decimals across
  the whole flank range (0.0536 at every length from 0 to 5,000), while the maximum over
  the *whole CDS* climbs 0.2757 → 0.4622. BT4's `pooled_risk` is top-3 over the whole CDS,
  so **BT4's aggregation is what is flank-sensitive, not the model's detection** — the same
  shape as the pooling-hinge and saturating-baseline findings.

### Measured
- **The splice candidate ranking is reliable; the delivered pick is not.** Before asking
  whether a splice model ranks synonymous candidates *correctly* — which needs labels nobody
  has — there is a prior question needing none: is the ranking stable, or ensemble noise? A
  ranking that changes with the training seed cannot be right even in principle. Pangolin's
  12 members (3 CV folds × 4 tissues) were retained separately and analysed as a **two-facet
  generalizability study**; a naive split-half was designed, critiqued, and **discarded**,
  because folds are re-training replicates (noise) while tissues are different biological
  targets (signal), and a random split averages the two into an uninterpretable number.
  Result: **Eρ² = 0.959 / 0.901 / 0.942** (tissue-general universe) across the three
  proteins, with true candidate variance exceeding every error term by 5–10×; fold-vs-fold
  Spearman +0.861 to +0.970, tissue-vs-tissue median +0.827 to +0.891, and σ²(candidate ×
  tissue) about **2×** σ²(candidate × fold) — heterogeneity behaving like heterogeneity.
  **One failure mode is excluded with no assay required.**
- **…but a reliable ordering does not give a stable winner.** The top candidates are
  near-ties, and in **2 of 3 proteins the argmax changes** with the fold or tissue used;
  Beclin1's worst case delivers a sequence the full 12-member ensemble ranks **7th of 30**.
  If splice Δ were routed into candidate *selection*, the delivered sequence would depend on
  Pangolin's configuration. This also sharpens the earlier low agreement with the PWM
  baseline: Pangolin's own ranking is stable, so the disagreement is substantive rather than
  Pangolin being noisy, and one of the two is wrong. **Reliability is not validity** — a
  ranking can be perfectly reproducible and perfectly wrong; every member shares an
  architecture and most of its training data, so a shared blind spot is invisible to all of
  this. A **floor census** justifies the statistic switch: the 0.5 hinge destroys 38–67% of
  all (candidate, member) cells and floors 5–7 of 30 candidates on **every** member, so
  nothing can be ranked on the hinged risk.

### Measured
- **The splice CNNs are not blind on designed coding sequence — measured with a graded
  implantation ladder.** The open question after the `N`-padding result was whether
  Pangolin is *blind* inside a designed CDS or *correctly silent* because there is no
  strong site to find; those are indistinguishable on a label-free panel. Planting a
  9-mer donor consensus by **substitution** (length and reading frame preserved) into two
  third-party designed hosts at three positions each: host baseline **0.0524** →
  full consensus **0.5700**, about **11×**, with the peak landing on exactly the base
  `CNN_ANCHOR_OFFSETS` predicts in 5 of 6 plants. Three controls establish the response is
  the splice *signal* rather than a reaction to an edit — composition-matched scramble
  **0.0525**, a weaker motif **0.0547**, and decisively a `GT`→`CT` ablation that keeps 7
  of the 9 bases and destroys only the invariant dinucleotide, **0.0543**. Graded
  throughout: 0.570 → 0.357 → ~0.053.
- **…and the detection floor is high enough to matter.** A *weakened but real* donor
  scores **0.357** and clears nothing, so BT4's 0.5 cutoff sits **above the
  intermediate-strength sites cryptic splicing actually uses**; on the `N`-padded shipped
  path even a full consensus is missed in 1 of 6 plants (4/6 vs 5/6). This **retires the
  "train the models" framing**: they demonstrably detect a strong site in this exact
  regime, so silence on clean designed CDS is not an inability to see, and what remains is
  an operating point — derived on labelled data, never trained. It does **not** establish
  correct silence: detecting a site BT4 planted is not evidence about sites nobody put
  there, and no label-free panel can supply that.

### Measured
- **The splice adapters' `N`-padding is not neutral — it systematically deflates scores
  inside the CDS.** Both wrapped CNNs pad with 5,000 literal `N` (upstream's own
  convention), and nothing had separated "the model sees little on a designed CDS" from
  "the input is off-distribution", since a broken input produces a low structureless band
  too. Four arms over the designed-CDS panel, differing only in the flank, scores sliced
  back to coding positions: `N`-padding median peak **0.2757**, random uniform ACGT
  **0.1731**, real human chr1 **0.3691**, composition-matched shuffle **0.4944**. Under
  real flanks several designed CDSs reach or cross 0.5 that never did padded — so the
  "everything floors below 0.5" picture is **partly an input artifact**, not purely a
  misplaced threshold. Replicated over three *different* real regions and three
  *independent* shuffles because the surprising half rested on one draw: real regions
  agree to **three decimal places** on every protein (real context is a stable
  background), while shuffles scatter and run higher in **9 of 9** comparisons
  (distribution shift, not restored function) — and random ACGT scoring *below* `N` rules
  out "any real bases beat `N`". **Licenses no threshold change**: there are no labels
  here, a higher score is not a more correct score, and the honest consequence is that a
  number from the `N`-padded path is a lower bound while passing a real `ConstructContext`
  changes the answer rather than refining it.

### Fixed
- **Corrected a false claim this project had just written into its own constitution.**
  Nine places across six files asserted that on the designed-CDS panel *"no position on
  any of 93 sequences reached 0.5"* — `CLAUDE.md` §6, `docs/NEXT_SESSION.md` (×2), the
  `CHANGELOG`, `docs/REVIEW_splice_calibration.md`, and three docstrings in
  `biomodels/splice/base.py`. It was **generalized from a four-sequence spot check to all
  93**, and the review document already contradicted it two sections later: PDE3A's Δ
  spread of `1.0885` is arithmetically impossible under a hinge at 0.5 unless something
  cleared 0.5. Re-measured against the hash-verified weights, the true count is **6 of
  93** (all six designs of one protein; no native clears it in any group). The finding
  the claim was offered as evidence for — that the hinge discards the CNN's signal —
  survives; the universal quantifier did not.
- **And the same measurement surfaced the inverse defect in the shipped default, which
  nothing had recorded.** The PWM baseline `bt4.biomodels.splice.default()` returns
  clears 0.5 on **93 of 93** designed sequences, native included, at peaks of
  **0.981–1.000** — so its hinge never binds, its risk and response coincide to four
  decimals, and it **flags a splice site on 100% of designed coding sequences**. A
  detector that fires on everything is uninformative in exactly the way one that fires on
  nothing is, and this is the path most users are on: it drives Studio's "distinct splice
  sites" column. So `DEFAULT_SITE_PROBABILITY = 0.5` is simultaneously **too high for the
  opt-in CNNs and too low for the shipped baseline** — one constant standing in for two
  different score scales, which is the concrete reason an operating point has to be
  derived per backend on data rather than nudged.

### Fixed
- **`--use-attested-splice` was a dead flag on both `optimize` and `validate`.**
  `_enable_attested_splice` was defined and **never called**, so the flag parsed,
  appeared in `--help`, and did nothing: a user who asked for the attested backend
  still got the uncalibrated one, and was told so only by the very tag they had just
  tried to change. Only `BT4_SPLICE_USE_ATTESTED=1` in the environment worked. A
  control that no-ops is worse than an absent one, because it reports an opt-in that
  did not happen. Wired in `main()` rather than per-subcommand — a per-command call
  site is exactly how it came to be dead on both. Regression-tested through the
  observable effect (verified to fail before the fix), not by asserting a call.
- **The cross-check printed the bare word `calibrated`.** Wiring the flag made that tag
  reachable from the CLI for the first time, and it named the stronger of the two claims
  BT4 keeps apart: the splice `calibrated` flag is set by a **fidelity** attestation —
  the adapter reproduces the published model bit-for-bit — and asserts nothing about
  whether a score of 0.5 means a 50% chance of splicing. It now reads
  `fidelity-attested (reproduces upstream; NOT statistically calibrated)`. The adjacent
  `(top-3 log-odds; uncalibrated)` line was already correct and is unchanged: that
  describes the pooling scale, which is uncalibrated either way.

### Added
- **Three industrial expression hosts: CHO, *B. subtilis* and *K. phaffii*** — taking
  BT4 from nine selectable organisms to **twelve**. Each ships a genome-wide codon
  table recounted from a release-pinned Ensembl CDS set (`Cricetulus_griseus_chok1gshd`
  CHOK1GS_HDv1, Ensembl 116; *B. subtilis* subsp. subtilis str. 168 ASM904v1, Ensembl
  Bacteria 63; *K. phaffii* GS115 GCA_000027005.1, Ensembl Fungi 63 — a division the
  builder did not previously reach) **and** a GtRNAdb tRNA table, because a codon table
  without tRNA data would make tAI silently unavailable exactly where a user asked for
  it — a shipped invariant, and the one that decided the scope of this change. All six
  tables reproduce their committed bytes under `--verify`.

  Checked against external ground truth rather than self-consistency: CHO's GC3 lands
  within 0.6 points of mouse and rat (it is a rodent, counted from a separate assembly);
  *B. subtilis* separates from *E. coli* by ~11 GC3 points through the identical Ensembl
  Bacteria path, reproducing their known ~7-point genome-GC difference; both AT-rich
  hosts take AAA/GAA for Lys/Glu and CHO takes AAG/GAG/CTG, as the literature reports.
  *B. subtilis*'s tRNA set is **86 genes, exactly the published count** for strain 168
  (Kunst et al. 1997) — corroboration, not just a pinned number; the *K. phaffii* count
  is the source's own and its note says so.

  Two limits are recorded in the sidecars rather than smoothed over. **CHO is the one
  organism whose two inputs are not assembly-matched**: GtRNAdb serves Chinese hamster
  tRNAs on the more fragmented CriGri_1.0 while the codon table is counted on
  CHOK1GS_HDv1, so its tRNA copy numbers are the weaker input. And **none of the three
  ships a highly-expressed reference set**, for three *different* measured reasons —
  PaxDb v6.1 has no dataset at all for *K. phaffii* (taxon 644223); for CHO it has only
  a single study (`PXD014877_Mueller_Nature_2020`), not the whole-organism *integrated*
  set every shipped table uses, and one built from a single study would carry the same
  `highly_expressed` label while meaning materially less; and for *B. subtilis* the
  integrated set exists and joins at **4,042/4,052 = 99.8%** via a declared `^BSU` →
  `BSU_` locus-tag rewrite, which is admissible (derivable from the two pinned files, not
  a third-party mapping) but needs the builder to support a per-spec identifier rewrite.
  That is now the next queued task rather than an unexplained absence.

  **The *B. subtilis* start-codon gap was measured for this organism, not inherited.**
  The shared validity filter requires an `ATG` start and drops **954 of 4,237 CDS
  (22.5%)** — TTG 553, GTG 387, ATT 8, CTG 5, ATC 1 — more than double the 9.6% it costs
  *E. coli*, because *B. subtilis* genuinely uses alternative starts. Counting the
  dropped genes back in (skipping the initiator, which is not a codon *choice*) moves
  **no** amino acid's most-used codon and shifts relative adaptiveness by at most 0.023,
  so the filter costs precision rather than correctness — the same verdict as *E. coli*,
  now established for the organism it actually applies to. Pinned by a test that asserts
  the drop tally from the sidecar and was verified to fail when the number is wrong.

### Fixed
- **BT4's splice risk pooling discarded its entire signal in BT4's own regime, silently.**
  `pool_log_odds` sums `max(0, logit(p) - logit(background))` with `background =
  DEFAULT_SITE_PROBABILITY = 0.5`, so only positions *above* 0.5 contribute. Measured
  against the hash-verified Pangolin weights on the designed-CDS panel, **no position on
  any of the 93 sequences reached 0.5** — peak scores ran 0.128–0.445 and differed more
  than twofold between a native CDS and its synonymous redesigns, and every one of them
  pooled to a risk of exactly `0.0`. *(Corrected below: that count was generalized from a
  four-sequence spot check and is wrong — the true figure is 6 of 93.)* So `delta_splicing` was identically zero for every
  candidate, the cross-backend rank agreements computed from those deltas were Spearman
  correlations of constants, and none of it was visible: a pooled risk of `0.0` meant
  either "no risk" or "nothing cleared an uncalibrated cutoff", with nothing telling the
  two apart. The constant is documented as *"a display / localization knob, not a
  calibrated cutoff"* but was wired in as a hard gate inside risk pooling, where instead
  of shifting a display it zeroed the output. **The background was deliberately not
  lowered** — that is the same uncalibrated knob pointed somewhere more flattering, and
  deriving a real operating point is the statistical-calibration gate's job, on data.
  Instead: new `pool_top_k_logit`, the same top-k log-odds with the hinge *and the
  background* removed (no operating point at all, monotone in the scores everywhere, so
  it separates what the risk has flattened — and **not a risk**: it goes negative and has
  no calibrated zero); new `PooledRisk` / `pooled_risk_detail`, carrying
  `n_above_background`, `max_score` and `below_background` so a zero is attributable; and
  every consumer that reports a risk now reports which zero it is — `bt4 designed-probe`,
  `bt4 validate --splice-backend`, the Studio ASSP banner, `BackendCandidateAudit`,
  `SpliceCrossCheck`, and `AgreementReport.degenerate`. Re-measured, Pangolin's
  background-free response spread across designs of one protein is 3.9–5.9 log-odds: it
  responds to synonymous change, and BT4 was throwing that away. `pool_log_odds` and
  `pooled_risk` are byte-identical — no shipped number moved.
- **The splice audit and cross-check pooled against a different cutoff than they
  localized at.** Both passed the caller's `threshold` to site localization and let
  pooling keep the default 0.5, so `--threshold 0.2` reported sites at 0.35 while the
  pooled risk meant to summarize them stayed exactly zero. They now pool against the
  threshold they localize at; at the default the two are the same number, so the shipped
  path is unchanged.
- **`bt4 designed-probe` rounded away the distinction its own conclusion rests on.** At
  `%.4f` an exactly-zero Δ spread and one of 1e-7 both print `0.0000`, and the two
  readings are opposite: exactly zero means the backend gave the native and every design
  the same pooled risk, while 1e-7 means it separates them and the signal is merely tiny.
  On the first real run Pangolin displayed `0.0000` for two of three proteins, so the
  formatting was deciding the finding. Zero now prints as `0` and sub-1e-4 values in
  scientific notation. *(The explicit line that change added — "this backend cannot rank
  the candidates at all" — was itself wrong, and is corrected by the pooling fix above:
  the zero was BT4's hinge, not the backend's silence.)*
- **`bt4 designed-probe` printed nothing while it ran** — the same silence-looks-like-a-hang
  defect fixed for `splice-gate` one change earlier, reintroduced in the new command:
  `probe_designed_cds` accepted a `progress` callback and the CLI did not pass one. Now
  wired, with `--quiet` to suppress it.
  - Reported per **protein group** rather than per sequence, and the docstring says why:
    the probe hands a whole group to `backend_agreement`, which scores its members
    internally, so the group is the finest boundary this layer can honestly report.
    Claiming per-sequence progress it cannot observe would be worse than coarse progress.
- Two guarantees the probe already had are now pinned: it reports the backends' **own**
  names (`consensus-pwm-baseline`, not the `pwm` alias — the bug found by running it),
  and `DesignedCdsProbe` has **no** `passed` / `promotable` field, structurally, because
  the panel it reads has no labels.

### Added
- **The designed synonymous CDS panel — BT4's own regime, and the one nothing had
  measured.** Every splice measurement so far is *recall on natural sites in natural
  genes*; BT4 emits a synonymous re-encoding of one protein's CDS, where the question is
  *specificity*.
  - **It carries no splice labels, and the reader refuses to accept one.** Designed
    coding sequence has no splice ground truth — nothing is assayed, none of it is
    annotated, and a motif is not a site. A `label` column is rejected **by name**, with
    the reason, because anyone adding one believed the panel has truth it does not.
  - **`bt4.api.read_designed_cds_panel`** (`biomodels/splice/designed_panel.py`) verifies
    the defining property rather than trusting it: within a group, every member must
    translate to the **same protein**, and the refusal names the offending member. A
    panel whose members are not synonymous is not a weaker synonymous panel — it is a
    different experiment, and every number from it would be about the wrong thing.
  - **`scripts/make_designed_cds_panel.py`** builds it from the already-committed
    `ranaghan2021_tab4.fasta` (CC BY 4.0): **93 sequences, 3 proteins, each the native
    human CDS plus 30 designs** from three anonymized commercial optimizers over ten
    repeat runs. Third-party by default — generating the designs with BT4 would make the
    answer partly a fact about BT4, which is the tool that would consume it. Verified:
    3 groups, one protein each, 93 distinct DNAs. `--include-bt4` adds BT4's own design
    when that comparison is wanted explicitly.
  - **`bt4 designed-probe`** reports what is measurable without labels: the **Δ spread
    across synonymous variants** — the number that decides whether these models are
    usable in BT4's loop at all, since synonymous positions are the only thing BT4
    changes — plus cross-backend rank agreement and sign agreement. It states that it is
    **not a gate** and has no threshold, and that its agreement figure is **not** the
    site panel's Jaccard and must not be set beside it.

### Added
- **`docs/REVIEW_splice_calibration.md` records the SpliceAI run and the first
  two-backend agreement measurement**, so all three panel20 results are in the repo
  rather than a chat log.
  - **Like-for-like on the shared task** (`--combined-track on`): Pangolin **0.983** vs
    SpliceAI **0.965** skill, 0.940 vs 0.907 top-k. The **top-k gap (0.033) reproduces
    the published one (0.040)** even though both absolute levels are inflated by an easy
    panel; the AP gap is compressed because both sit near ceiling there.
  - **Separating cost SpliceAI nothing** — its combined figures are the exact mean of its
    separated pair to three decimals, so its kind discrimination is effectively perfect.
    The comparability caveat added in the previous change is right in principle and its
    measured magnitude here is **nil**; an earlier draft implied otherwise.
  - **Agreement: Jaccard 0.855**, and of 333 sites both find 300, only Pangolin 15, only
    SpliceAI 6, **neither 12**. Running both is not redundant — 21 sites (6.3%) are found
    by exactly one model — but the 12 both miss are a **correlated blind spot**, which is
    the standing limit on agreement as an uncertainty signal: it bounds independent
    error, not shared error.
  - Records an open question rather than resolving it: **7 positions both models call are
    not annotated sites**, and on a MANE-Select-only panel those are as plausibly genuine
    non-MANE isoform sites as shared false positives. The two readings have opposite
    implications and this panel cannot separate them.

### Fixed
- **Two consecutive `splice-gate` runs invited a comparison their numbers do not
  support, and nothing said so.** Pangolin emits one combined `P(splice)` track and is
  scored as a single `splice` stratum; SpliceAI emits separate donor and acceptor tracks
  and is scored as two. On the same GENCODE panel that reads as Pangolin **0.983** against
  SpliceAI **0.969 / 0.961** — and the natural conclusion is wrong, because a
  kind-separated stratum treats **the other kind's sites as negatives**: the backend must
  locate the site *and* get its kind right. Separated is the harder task, so setting it
  beside a combined figure understates it.
  - A separated run now carries a note saying exactly that, mirroring the one a combined
    run already had, and naming the flag that fixes it.
  - **`--combined-track {auto,on,off}`** exposes the stratification the pipeline already
    supported but the CLI could not reach, so a two-track backend can be scored on the
    same task a one-track backend solves and the two numbers become like-for-like.

### Added
- **`bt4 splice-agreement` — do two backends call the same *positions*?** §6 names
  cross-backend agreement a first-class uncertainty signal, and `backend_agreement`
  provided it only for the design flow (ranking candidate CDSs by Δsplicing). A
  site-prediction panel has no candidates to rank, so the question was unanswerable
  there — which is where it matters most now that both CNNs are attested.
  - **Comparing two gate reports is not a substitute, and a test proves it.** Two
    backends can each recover every annotated site while pointing at entirely
    different bases; constructed as `test_two_backends_can_score_identically_and_agree_on_nothing`,
    where both score perfectly and Jaccard is 0.
  - Reports **Jaccard of called positions**, **Spearman over the union of calls**
    (not over every base — at 1-in-2,600 prevalence a whole-panel correlation is
    dominated by the ~99.96% where both are ~0 and reads near 1.0 regardless), and the
    2×2 of annotated sites recovered by **both / only one / neither**.
  - Calls are the top-*k* **positive** scores per window, `k` = that window's site
    count. Both details were forced by failing tests: pooling across windows measures
    window length rather than agreement, and taking top-*k* unconditionally pads with
    zero-scored positions chosen by the index tie-break — manufacturing agreement
    between two backends that found nothing in common, in the flattering direction.
  - No anchor shift is applied *between* the backends (both CNNs anchor identically);
    the shift places the annotated sites in their shared frame. A combined-track
    backend collapses exactly as the gate collapses it.

### Added
- **SpliceAI passed its integration-fidelity gate — `max_abs_deviation` exactly 0.0 over
  18 cases — so BOTH wrapped CNNs now reproduce their published models bit-for-bit and
  both attestations ship.** Part A is complete.
  - Run against the CC BY-NC weights, all five of which matched
    `PINNED_WEIGHT_SHA256`, on the **same panel as Pangolin's**
    (`content_hash f3589fd1e10ffb73e…`) — so the two are directly comparable rather than
    approximately so, which is what a real cross-backend agreement figure needs.
  - **The pass is not the vacuous kind.** Peak site probability across the panel spans
    **0.029–0.925** (mean 0.442), so a transposed one-hot, a swapped acceptor/donor
    channel, or a wrong ensemble would have shown up. The runner reports the spread and
    says so itself rather than leaving it to be asserted.
  - BT4's SpliceAI adapter had **never been executed against real weights** — every
    prior test drove a fake predictor — and it was correct on the first run.
  - Promotion is unchanged: still opt-in via `BT4_SPLICE_USE_ATTESTED=1`, `default()`
    still returns the PWM baseline, and a real-flank score still clears `calibrated`.
    An attestation is a claim about the **wrapper**, not about the scores; statistical
    calibration for designed coding sequence remains a separate, unmet gate.

### Fixed
- **`--no-deps` alone does not make the SpliceAI capture work, and the docs said it
  did.** `spliceai/utils.py` imports `pandas`, `numpy`, `pyfaidx` and `keras`, all of
  which `--no-deps` skips. The failure is late and misleading: `available()` still
  reports `True` (it needs only Keras and the weight files) while the capture dies
  importing `one_hot_encode`. Both docs now give the second `pip install` line, verified
  with **pysam still absent**.
- **The capture script's refusal named a cause it had not established.** It led with
  "pin `numpy<2` — `spliceai/utils.py` still calls `np.fromstring`" and then printed an
  underlying error that contradicted it (the real cause was a missing `pandas`). It now
  leads with what actually failed, names the three dependencies `--no-deps` skips, and
  offers the numpy story as one possibility rather than the diagnosis. The refusal
  itself was correct throughout — it declined to substitute its own encoder, which is
  the property that keeps the gate capable of failing.

### Changed
- **The splice-CNN install route is documented from a real Windows setup**, in a new
  "Splice CNN environment gotchas" section of `NEXT_SESSION.md` mirroring the RiboNN
  one, plus corrections to the runbook.
  - **`pip install spliceai` fails on Windows and it does not matter.** It depends on
    **pysam**, which has no Windows wheels and cannot build. `--no-deps` is the answer:
    pysam serves SpliceAI's own VCF command line, and BT4 never uses it — the adapter
    resolves weights with `importlib.util.find_spec` (which does not execute the module)
    and loads the `.h5` files with Keras directly. Verified with no pysam present:
    `available()` is `True` and upstream's `one_hot_encode`, which the capture script
    requires, imports fine.
  - **`TF_ENABLE_ONEDNN_OPTS` must be held constant across capture and gate.**
    TensorFlow warns that oneDNN can change numerical results by reordering
    computation; capturing with one setting and gating with the other would produce a
    deviation caused by TensorFlow rather than by BT4's adapter — the one thing a
    fidelity gate must not confuse.
  - **The two model stacks share one environment**, contrary to the plan's assumption:
    TF 2.15 forces `numpy<2`, and `torch 2.13.0+cpu` runs fine on numpy 1.26.4. Also
    recorded: Windows has no `conda`, so `python -m venv`, and the activated prompt is
    the check that a `pip install` lands where intended.
  - Two now-stale runbook claims corrected: the `tensorflow>=2.6` pin note (done), and
    the assertion that `_import_keras` "already falls back to `tf_keras`" — it did not
    effectively, which is what the same-day fix addressed.

### Fixed
- **The `tf_keras` fallback for SpliceAI's weights was unreachable, and the extra's
  TensorFlow pin steered installs straight into the case it existed for.** SpliceAI's
  weights are 2019 **Keras 2** `.h5` graphs and Keras 3 cannot load them; from
  TensorFlow 2.16 `tensorflow.keras` **is** Keras 3.
  - `_import_keras` ordered its candidates by module *availability* — try
    `tensorflow.keras`, fall back on `ImportError` — so under TF ≥ 2.16 the first
    import succeeded and returned the wrong Keras. The shim was only ever reached when
    TensorFlow was absent entirely, which is the one situation it cannot help. Order is
    now decided by the installed Keras **generation** (`_ambient_keras_is_v3`), read
    from the packages rather than the TensorFlow version so an explicit `keras<3` or
    `TF_USE_LEGACY_KERAS` is respected.
  - `bt4[splice-spliceai]` declared `tensorflow>=2.6` with **no upper bound**, so a
    fresh install pulled a TensorFlow that installs cleanly and then fails at load with
    an opaque deserialization error about a file that is not corrupt. Now
    `>=2.6,<2.16`, with the reason recorded inline and the `tf_keras` route named for
    anyone who wants a newer TensorFlow.
  - Three tests pin the resolution across environments (Keras 3 + shim, Keras 2, and
    Keras 3 with no shim — which must degrade rather than refuse to import, so the
    honest failure happens at the weights). The first was verified to fail unfixed.

### Added
- **`docs/REVIEW_splice_calibration.md` now records the site-prediction half too**, so
  both halves of Part B are measured rather than planned. Run against BT4's own wrapped
  Pangolin and the hash-verified GPL weights, on a GENCODE v44 / GRCh38 panel of 20
  held-out MANE windows (861,096 positions, 333 sites, motif consistency 100%).
  - **The per-kind anchors are confirmed on real data** — donors peak at −1 for 100% of
    sites, acceptors at +1 for 99%. Those offsets were derived in #102 from upstream
    source rather than measured, and this was the largest correctness risk in this half.
  - Pangolin scores **skill 0.983 / top-k 0.940** against the `pwm` baseline's 0.096, and
    with `--min-pr-auc-skill 0.75` declared beforehand the run reports **`PROMOTABLE on
    this panel: True`** — a first for any BT4 splice backend. Stable across a 2.3× panel
    size change (skill moved −0.005).
  - Recorded as **not a promotion**: the figures sit **+0.133 / +0.150 above** Zeng & Li's
    published 0.85 / 0.79, consistently, which says the panel is easier than the
    genome-wide benchmark rather than that the model is better. Pangolin's combined track
    also means there is **no exonic/intronic split** here, so the penalty that matters most
    to BT4 is not checkable on this panel shape. `calibrated` is unchanged and `default()`
    still returns the PWM baseline.

- **`bt4 splice-gate` now reports which window it is scoring.** A CNN-backed run over a
  real GENCODE panel takes tens of minutes — the wrapped models read ~10 kb of context
  per position — and the command printed *nothing* until it finished, which is
  indistinguishable from a hang. Two real runs were interrupted or queried on exactly
  that ambiguity before this was added.
  - `score_splice_panel` / `run_splice_panel_gate` take an optional `progress`
    callback, `(index, total, window_id, length)`, invoked **before** each window so the
    one being waited on is the one named. `None` is the default, keeping the API
    print-free (§3: only `cli` prints).
  - It carries the window's **length** because that, not the count, is what the
    remaining wait is proportional to: windows are whole gene spans and vary by more
    than an order of magnitude, so "12/20" alone predicts nothing.
  - The CLI writes it to **stderr**, so the report on stdout stays pipeable, and
    `--quiet` turns it off.

### Fixed
- **An ECE ceiling was accepted as a pre-registered bar, and it is a bar nothing can
  fail.** The first real GENCODE site-prediction run exposed it: at that panel's
  prevalence (129 sites in 372,634 positions) the `constant` baseline — which predicts
  the base rate — scored ECE **0.000000**, the identical figure a *perfect* classifier
  gets, while the backend under test scored 0.050. The no-information control was better
  calibrated than the model.
  - `thresholds_declared` now requires a **discrimination** threshold (`min_pr_auc` or
    `min_pr_auc_skill`). `max_ece` alone leaves a run un-promotable and says why. This is
    the same hole as dropping the baseline a backend would lose to: a condition that
    cannot fail is not a pre-registration.
  - A run whose baselines **match or beat the head's ECE** now says so, naming them,
    wherever the ECE column could be read as evidence — not only when it was declared as
    a threshold. ECE rewards predicting the base rate, so at splice prevalence it
    describes the score distribution and never the quality; the skill column carries the
    verdict.
  - The `--max-ece` CLI help now states this rather than presenting the flag as a bar.
  - A pre-existing test asserted the disproved belief — that any one of the three
    thresholds counts as pre-registration. It is corrected, with the measurement that
    superseded it named in its docstring.

### Fixed
- **`read_splice_panel` could not read a single panel the GENCODE builder produces.**
  Python's `csv` caps one field at 131,072 characters; a splice window is a gene span
  plus 5,000 nt of flank each side, which routinely exceeds that. So BT4's own writer
  and its own reader were mutually incompatible for **every real gene**, and the failure
  surfaced as a bare `_csv.Error` from inside the stdlib that named nothing about the
  format. Every fixture in the panel tests is a few hundred bases, which is exactly why
  it survived to be found by running it on chr1.
  - New `biomodels/_csv.py`: `relaxed_field_size()`, a context manager that raises the
    cap to `2**31 - 1` for one parse and **restores the caller's value afterward**,
    including on an exception — the limit is process-global module state, so a library
    must not raise it permanently. `2**31 - 1` rather than `sys.maxsize` because
    `csv.field_size_limit` takes a C `long`, which is 32-bit on Windows even in a 64-bit
    build.
  - Applied to the expression panel reader too, for a **different and smaller** reason,
    stated rather than blurred: a valid row there can never approach the cap
    (`MAX_CDS_UTR3_LEN` holds CDS+3'UTR to 11,937 nt), so this only decides *which*
    error an over-long row gets — the panel's own message naming RiboNN's limit, instead
    of a stdlib error citing a limit the format does not have.
  - Two regression tests, each verified to fail against the unfixed reader, and each
    using a field genuinely over the cap — the only size that reaches the bug.

### Fixed
- **`scripts/make_gencode_splice_panel.py` ignored `--limit` after the first window.**
  `build_windows` rebound its own `limit` parameter to the window's sequence length
  partway through the loop, so from the second iteration the guard read
  `len(windows) >= 35000` and could never fire. The flag exists so a user can cost a
  trial run before committing to a real one; instead `--limit 5` built every MANE
  transcript on five chromosomes — a multi-hour build and a panel two orders of magnitude
  too large to score with a CNN. Found by running it. The local is renamed, and a
  regression test asserts the built count equals the requested one.

### Added
- **`docs/REVIEW_splice_calibration.md` — the first measured results Part B has ever
  had.** The variant-effect half was run on a maintainer machine against
  `splicebench2023`'s own pre-computed scores (no weights involved, so this measures
  BT4's gate against a published benchmark rather than BT4's adapters against a model).
  Panels are recorded by content hash.
  - **Removing training-chromosome genes costs a third of intronic skill** (0.724 →
    0.480) while exonic *rises* slightly (0.365 → 0.419). The exonic/intronic
    skill gap collapses from 0.359 to **0.061** on chr3 alone. The write-up names the
    confound rather than claiming leakage: holding out chr1/3/5/7/9 *replaces* the gene
    set instead of filtering a fixed one, and gene and chromosome are collinear, so no
    split of this benchmark can separate difficulty from training overlap.
  - **Held-out exonic ECE is 0.345.** These delta scores cannot be read as probabilities
    of disruption, which matters because BT4's own operating point treats scores as if
    they were. Direct evidence against importing published cutoffs.
  - Records what it does **not** establish: nothing about BT4's adapters (these are the
    benchmark's numbers), and nothing about BT4's regime (natural variants in natural
    genes, not designed synonymous CDS). No flag changed; `default()` still returns the
    PWM baseline.

### Fixed
- **`bt4 variant-gate` compared strata on raw average precision, and told the reader a
  panel was broken when it was not.** Its guidance said the exonic figure "should sit
  well BELOW the intronic one. If it does not, suspect the panel build before the
  model." A real held-out run on `splicebench2023`'s chr3 genes produced exactly that
  inversion — exonic AP **0.665** above intronic **0.598** — and the panel was fine. The
  cause was prevalence: **42.3% exonic against 22.7% intronic**, so average precision's
  floor differed by nearly double and raw AP was never comparable between the strata.
  On `pr_auc_skill`, which exists for precisely this, the expected ordering held
  (0.419 < 0.480), as it did on ROC-AUC.
  - The report now prints **`prev`** beside every AP, because average precision without
    its floor is not interpretable — a point the gate's own docstrings make and the CLI
    was not surfacing.
  - The ordering is read on **skill**, never on raw AP.
  - An AP inversion that skill contradicts is **detected and explained** as the
    prevalence artifact it is, rather than left to send a reader hunting a panel bug.
    A genuine skill inversion still warns, because that one really is worth suspecting
    the panel over.

### Changed
- **`bt4 variant-gate` no longer claims the exonic/intronic *gap* is comparable to the
  published one.** The previous wording said "only the GAP is meant to be compared
  directly". Measured against the real benchmark, it is not: a single tool sits above a
  median-over-tools on both strata and by more on the harder one, which **compresses**
  the gap — 0.266 observed against 0.354 published. Only the **ordering** is comparable,
  and that is the finding that matters.
  - The panel-composition hypothesis was tested and **refuted**: matching the published
    six-dataset composition with `--include-mlh1` moved both figures *further* from the
    published pair (exonic 0.535 → 0.575, intronic 0.777 → 0.841), because MLH1's
    variants are easier than average. The output now says so rather than implying the
    flag would close the gap.
  - Notes that an exact reproduction *is* possible — the archive carries all eight tools'
    scores — and that it needs the other six mapped in.

### Changed
- **`bt4 variant-gate` now states the published anchor's panel composition.** It printed
  Smith & Kitzman's 0.419 / 0.773 beside the run's own figures without saying that the
  published number is a *median across tools* pooled over **all six** of their datasets
  — while the converter excludes MLH1 by default, so a default run is five datasets and
  3,616 variants against an anchor computed on 3,912. Two numbers printed side by side
  imply a comparability that was not there. The output now names both reasons a single
  tool will not match, points at `--include-mlh1` for a like-for-like composition, and
  says plainly that only the **gap** is meant to be compared directly.

### Added

- **`bt4 variant-gate`** — runs the splice gate on a measured variant panel. The variant
  half was API-only, so following the runbook meant writing Python; it is now a command.
  It lists the panel's score columns when none is chosen (masked and unmasked answer
  different questions), warns when genes sit on chromosomes both models trained on, and
  prints the published **0.419 exonic / 0.773 intronic** anchor beside the run's own
  numbers — because reproducing that gap is what the benchmark is for, and a run that
  does not should be read as a panel-build problem before a model one.

- **The variant half of the splice gate is reachable from data.**
  `SpliceVariantCase` shipped with no way to construct one, which left the exonic /
  intronic task — the one that matters most for BT4, since it designs coding sequence —
  unusable.
  - **`bt4.api.read_variant_panel`** (`biomodels/splice/variant_panel.py`) — a
    tab-separated format where every column beyond `variant_id`/`group`/`region`/`label`
    is read as a **named score column**, so a benchmark's own pre-computed predictions
    come through as data rather than needing a schema change. `panel.cases(column)`
    selects one.
  - A score column with **gaps is refused**, not silently scored on its covered subset:
    a tool that did not cover every variant is a real situation, and quietly dropping the
    rows it missed while reporting the panel's name is the dishonest response to it.
  - **Held-out status is checkable, and usually fails.** Each row may declare its gene's
    chromosome; `training_overlap` / `held_out` report what the run can support, and a
    gene with no declared chromosome is *unknown* rather than clean.
  - **`scripts/make_splicebench_variant_panel.py`** converts `kitzmanlab/splicebench2023`
    (MIT, 3,616 variants, already scored by eight tools). It maps `sdv_fc2` → label and
    `exon` → region, keeps masked and unmasked scores as **separate** columns so that
    choice stays at gate time, stamps each gene's chromosome, and excludes the separate
    296-variant MLH1 clinical set unless asked — 3,616 and 3,912 are both true about
    different things.
  - This makes the cheapest real check available: gate the benchmark's own scores, with
    no model installed, and confirm BT4 reproduces the published 0.419/0.773 split. If it
    does not, the defect is in BT4.

- **`scripts/make_gencode_splice_panel.py` — the site-prediction panel, built rather than
  hand-assembled.** Turns a pinned GENCODE v44 + GRCh38 into the format
  `bt4.api.read_splice_panel` reads, with the position convention correct by construction.
  - **The arithmetic was executed, not reasoned about.** Over 1,206 annotated sites from
    64 chr1 MANE Select transcripts it comes out **99.42% canonical GT/AG on both
    strands**, while the two plausible wrong conventions score 0.08% and 44.2%. The
    residual ~0.6% is the real minor spliceosome, which is why the reader's floor is 90%.
  - **Two traps that silently relabel true positives as negatives**, neither catchable by
    a motif check: a ±5,000 nt window contains a **median of 8** annotated sites (only
    2.8% contain the centre site alone), and **27%** of gene-body windows contain
    opposite-strand sites. The script collects every overlapping MANE transcript's sites,
    and skips antisense-overlapping windows by default (`--keep-antisense` opts in and
    records the count per window).
  - MANE Select filtering is what keeps a panel stable: from GENCODE v44 to v50 the
    protein-coding transcript count on the held-out chromosomes grows **4.1×** while MANE
    Select grows **1.3%**.
  - Tested against a synthetic genome for both strands, with BT4's own reader
    independently verifying the output at `motif_consistency: 1.0` — so the composition of
    GTF convention, strand handling, reverse-complementing and index mapping is evidence
    rather than assertion.

- **SpliceAI's integration-fidelity tooling, completing Part A's scripts.** Only the
  licensed weights step is left for that backend.
  - **`scripts/capture_spliceai_panel.py`** — the sibling of the Pangolin capture, with
    the same runtime independence guard and the same statically-enforced rule that it
    never imports `bt4`.
  - **It imports upstream's own `one_hot_encode` rather than re-deriving it.** Pangolin's
    CLI encodes inline, which forced that capture to reimplement it; SpliceAI ships the
    encoder as a reusable function, so importing it is strictly stronger evidence — a
    transposed layout, a wrong base order, or a mishandled `N` in BT4's own
    `_one_hot_rows` now surfaces as a gate **failure** instead of being reproduced
    identically on both sides.
  - **No fallback encoder, deliberately**, and a test asserts the script defines none.
    The obvious fix for the NumPy 2 `np.fromstring` breakage — adding a local encoder —
    would silently destroy the independence while leaving the gate passing. The script
    refuses and names the `numpy<2` pin instead.
  - **`run_splice_fidelity_gate.py` now dispatches on the capture payload's own
    `backend`**, so a Pangolin capture can never be checked against the SpliceAI adapter:
    the numbers would still be numbers, they would just describe a different model. A
    payload with no `backend` field is treated as Pangolin, so existing captures keep
    working. Its panel-strength warning reads both capture shapes — one combined
    `P(splice)` track for Pangolin's binary head, two for SpliceAI's 3-way softmax — so a
    two-track capture is no longer reported as a flat, zero panel.
- **An annotated splice-panel format, and the baselines a splice backend must beat.**
  The runnable half of Part B of `docs/DESIGN_splice_cnn_calibration.md`: everything
  needed to gate a wrapped CNN on real data except the data.
  - **`bt4.api.read_splice_panel`** (`biomodels/splice/panel.py`) — a tab-separated
    window format (`window_id`/`group`/`sequence`/`donors`/`acceptors`), a strict
    reader, and an order-independent content hash so a gate result is bound to exact
    bytes (invariant #7). `group` is the chromosome, and overlap with the models'
    training set (chr 2, 4, 6, 8, 10–22, X, Y) is reported.
  - **The position convention is verified, not trusted.** A splice panel has one
    catastrophic failure mode and it is silent: annotate one base off and every score
    is misaligned, the model looks incompetent, and nothing says why. BT4 pins the
    anchor its own PWM baseline already uses — a donor is the `G` of the
    intron-opening `GT`, an acceptor the `G` of the intron-closing `AG` — and, since
    ~99% of human introns are canonical, **refuses** a panel below 90% consistency
    while naming the shift that would have worked. A panel built to the
    exonic-boundary convention (what the GENCODE recipe produces) is rejected with
    *"move each donor +1 and each acceptor -1"* rather than silently scored.
  - **`bt4.api.splice_panel_gate` / `bt4 splice-gate`** (`pipeline/splice_gate.py`) —
    runs the acceptance gate on a backend and **the same gate on four permanent
    baselines**: `permutation` (the null), `gt_ag` (the canonical dinucleotide rule
    ~99% of introns follow), `pwm` (BT4's shipped `ConsensusPwmSplicePredictor`, the
    free incumbent) and `constant` (the per-stratum base rate — perfectly calibrated,
    completely useless, so its excellent ECE is visible rather than a trap). The
    comparison is **per stratum**, so beating the motif on donors cannot excuse losing
    to it on acceptors.
  - **BT4's own default cannot be evidence for itself.** Run the PWM backend as the
    head and it ties the `pwm` baseline exactly, so `beats_every_baseline` is `False`
    — the structural counterpart of the expression gate's "the null model provably
    cannot pass".
  - **Both alignment traps are reported, never assumed.** `anchor_offset` is an
    explicit input and an `AlignmentDiagnostic` shows where the backend's score
    actually peaked around each declared site: a perfect oracle shifted two bases reads
    as PR-AUC 0.006, and the diagnostic says the anchors disagree instead of leaving it
    to look like a bad model. And a **combined-track** backend (Pangolin emits one
    `P(splice)` and leaves `acceptor` all-zero) collapses to a single `"splice"`
    stratum rather than being scored as hopeless at acceptors — an artifact of the
    wrapper, not a finding about the model.
  - Promotion needs three conditions at once, reported separately: the gate's
    thresholds, beating every baseline in every stratum, and the panel being **held
    out**. A panel drawn from training chromosomes can never be `promotable`.
  - **Nothing here flips a flag**, and it answers a different question from the
    fidelity attestation: that one proves BT4's wrapper reproduces the published
    model, this one asks whether the numbers mean what they claim. ASSP is deliberately
    not an option — it is network-derived and outside the reproducible-from-manifest
    guarantee, so it cannot support a gate result.
- **Pangolin passed its integration-fidelity gate — the first wrapped model in BT4
  to pass one.** On a maintainer machine holding the GPL weights: 18 cases,
  tolerance `1e-3`, **max absolute deviation exactly 0.0**. BT4's adapter
  reproduces upstream Pangolin's per-position scores bit-for-bit across 6,549
  positions. The `FidelityAttestation` is committed at
  `src/bt4/biomodels/splice/data/pangolin.attestation.json` — eight license-clean
  scalars plus the public weight SHA-256s, never a raw model score.
  - **The pass is not vacuous.** All three donor probes peaked at position 302 with
    the consensus planted at 300–309, across three different random flank sets, so
    the model was tracking sequence rather than the `N`-padding boundary. Designed
    CDSs averaged 0.085 peak P(splice); consensus probes 0.61–0.71.
  - **Promotion is opt-in** (`BT4_SPLICE_USE_ATTESTED=1`, or
    `--use-attested-splice` on `optimize` / `validate`). `default()` still returns
    the honest PWM baseline, and a real-flank score still reports uncalibrated: the
    gate was captured on the `N`-padded path, and that scoping survives promotion.
  - **It certifies the wrapper, not the biology.** These models score median prAUC
    **0.419 on exonic** variants vs 0.773 intronic (Smith & Kitzman 2023), and BT4
    designs coding sequence. A statistical-calibration gate remains unmet.
  - SpliceAI is unchanged at `calibrated=False`; its gate has not been run.
- **A reproducible three-script workflow for the gate**, whose separation is the
  point: `make_splice_panel.py` (may use `bt4`) → `capture_pangolin_panel.py`
  (**never** imports `bt4`, enforced by a runtime guard and a static AST test) →
  `run_splice_fidelity_gate.py`. The runner reports panel *strength* and warns when
  a pass rests on a panel too flat to discriminate.
- `scripts/check_splice_weights.py` — verify installed weights against BT4's pins
  before investing in a capture, plus a public `weights_dir()` on both adapters.


- **`docs/DESIGN_ribonn_calibration.md`** — the step-by-step runbook for `NEXT_SESSION.md`
  item 11, mirroring the splice side's `DESIGN_splice_cnn_calibration.md`. It carries the
  operational half that the research doc deliberately left out: the Windows (`cmd.exe`) /
  WSL install, the `max_shift` determinism check, the free Stage-1 sensitivity checks (run
  these *before* buying any data — several can end the project decisively), the ranked
  panel hunt with the panel-TSV format, the pre-registration template, the exact
  `run_expression_gate` / `bt4 expression-gate` commands, and the Stage-5 decision tree
  with the `attest_expression` → `verified_predictor` promotion snippet. Every command and
  API reference is checked against the code that shipped in the previous change; the
  "wiring caveat" (nothing in `src/` calls `verified_predictor` yet) and the scoped nature
  of any pass are stated plainly. `NEXT_SESSION.md` item 11 and the research doc now
  cross-link it.
- **`docs/RESEARCH_ribonn_calibration.md`** — the evidence behind `NEXT_SESSION.md` item
  11, assembled from primary sources (the paper's full text, the upstream repository, the
  dataset records) with an adversarial verification pass over every dataset and
  statistical claim. It records what is now verified rather than assumed:
  - **The paper never tests synonymous CDS variants of one protein under a fixed UTR.**
    "synonymous" appears only in its introduction and reference titles; its codon
    analysis is *insertional* and changes the protein. The one CDS-attributable number is
    **r² = 0.11**, and zero-shot on designed reporters it scores **0.17–0.19** against
    **0.62** on natural genes. That is the regime gap, stated from the source.
  - **CORRECTED — the homology-leakage argument does not survive.** The paper reports a
    homology control ("removing highly homologous test mRNAs led to highly similar r²").
    It is qualitative and weak, but the charge as previously framed is refuted; the
    argument that stands is simply that the synonymous-variant regime was never
    evaluated.
  - **CORRECTED — conformal coverage stays *valid* under arbitrary unit mismatch.** The
    guarantee holds for any score function; unit mismatch destroys *sharpness*, not
    validity. So "different units, therefore coverage is meaningless" would have been
    wrong — the honest statement is that unrecalibrated intervals are valid and useless,
    which is why `width_over_iqr` exists.
  - **CORRECTED — there is no community-standard Spearman cutoff.** Any bare 0.4 or 0.5
    would be invented, so `min_spearman` is documented everywhere as a pre-commitment.
  - **CORRECTED — reporting a calibration slope of 1 after fitting the link is
    circular**, which is why `link_slope_spread` is reported instead.
  - **No public dataset fully qualifies.** A ranked table records each candidate with its
    licence, size, readout and verdict, so the panel hunt starts from evidence rather
    than optimism — including that the two Ranaghan panels already in this repo carry
    **no measurements** (that paper measured one sequence, in *E. coli*), making them a
    sensitivity resource and never a validation panel.

- **The construct-context gap is now evidence-backed, not just architectural.**
  - **N-padding the splice CNNs is an artifact generator, not a neutral default.**
    OpenSpliceAI documents that SpliceAI predicts donors and acceptors at N-padded
    boundaries "with an extremely high signal that disappears when the sequence is
    padded with the actual genomic sequence" -- so BT4 hallucinates sites at exactly
    the two positions a designer cares about. Its re-benchmark shows accuracy rises
    steeply only through ~400 nt of context, so **+/-400 nt of real flank beats
    +/-5 kb of N**.
  - **Every documented cassette-splicing failure has the acceptor in the vector**:
    Cheng 2022 (13/17 genes at their own retained exon-exon junction, 17/17 at a V5
    tag, acceptors from the mPGK-PuroR linker / SV40 / Neo-KanR), Kowarz 2022
    (codon optimization *created* donors; acceptor was adenoviral pIX past the
    poly(A)), De Ravin 2022 (a cHS4-insulator acceptor drove clonal expansion in a
    clinical trial; fixed by a 2-bp AG->TG change).
  - **oORF pairing is now the highest value-per-cost item in the roadmap.** An
    out-of-frame AUG in the user's 5'UTR whose in-frame stop lands *inside* the CDS
    represses significantly more than a non-overlapping uORF (Johnstone 2016,
    P = 1.23e-3) -- and the stop position is a function of BT4's synonymous choices,
    so BT4 can move it. Deterministic, no ML, no calibration gate, and no shipping
    optimizer does it.
  - **`RampTerm` implements a falsified mechanism.** It rewards *lower* codon
    adaptiveness over the first ~35 codons; Goodman, Church & Kosuri (*Science* 342,
    >14,000 reporters) found "reduced RNA structure and not codon rarity itself is
    responsible." The 5' effect is real; the lever is wrong.
  - **The highly-expressed reference premise is falsified in *E. coli***: Welch 2009
    found the favourable codons were those read by tRNAs most charged under
    starvation, *explicitly not* the codons abundant in highly expressed genes.
    Keeping highly-expressed as the default stays defensible, but no tooltip may
    imply it makes CAI predictive.
  - **Radrizzani 2024 (*Nat Rev Genet*) is good news for BT4's direction**: human
    synonymous-site selection is not translational selection but avoidance of
    spurious transcripts, mis-splicing and cryptic splice sites, plus a
    high-GC / low-CpG "self" signature -- which promotes things BT4 already ships to
    first class while demoting CAI.
  - **Prior art recorded so BT4 does not overclaim**: TIsigner already takes
    `-u/--UTR` with mammalian accessibility windows, and DNA Chisel (MIT) already
    evaluates constraints over an entire plasmid record. Context-as-constraint-scope
    is not novel; context-as-optimization-substrate is.
- **Vendor manufacturability rules BT4 lacks, exactly as published**: windowed GC
  **range** (max - min across 50 bp windows <= 50 points, Twist DOC-001081 REV4), a
  **Tm-based repeat trigger** (any repeat with Tm >= 60 C regardless of length --
  length-only rules miss GC-rich short repeats), and **per-base homopolymer limits**
  (IDT: A/T >= 10, G/C >= 6). Report a profile, not pass/fail, and do not hard-code
  a threshold BT4 cannot cite.
- **An expected-effect-size table** the app and CLI should render beside every
  result, with each range cited and the vendor-authorship caveat attached.
- **A "what NOT to do" section** in the survey -- fourteen plausible-sounding moves
  the evidence refuses, including not building a single global MFE objective, not
  asserting ensemble-free-energy design is NP-hard (it is explicitly open), not
  flipping a model to `calibrated=True` because its *input* improved, and not
  transmitting a user's vector backbone anywhere.
- **A measured product review, and a re-pointed queue**
  (`docs/REVIEW_2026-08_expression_and_context.md`). BT4 was audited against what
  it is *supposed* to be — a codon optimizer for **protein expression**, aware of
  the **5′UTR** and the **vector backbone**, honouring user-specified
  manufacturability rules, in a **simple honest UI**. Every number in the review
  was produced by running the code in this tree, not read from a docstring, and
  §10 carries the reproduction commands. **No source file changed; this is
  documentation and a queue re-point.** Headline findings:
  - **Verdict.** BT4 is an unusually well-engineered, unusually honest **CAI
    optimizer with manufacturability constraints** — not yet an expression
    optimizer, because **the optimizer only ever sees the CDS**. Folding is scored
    on `CDS[0:48]` with no leader; the splice CNNs are scored on the CDS floating
    in 5,000 literal `N` bases; RiboNN can take UTRs but sits where it is
    structurally forbidden from influencing delivery. `OptimizeConfig` has 40
    fields and not one is sequence outside the CDS.
  - **Four defects, measured.** (1) `run_frontier` — the path **BT4 Studio's
    Optimize button runs** — only *reports* GLOBAL rules, so setting *Max repeat
    length* returns a green `proven_optimal` badge on a sequence with a 58-nt
    repeat and 96 violations. (2) `run_validate` / `POST /validate` silently drop
    every GLOBAL constraint, reporting zero violations on a sequence with a 31-nt
    exact repeat; no test covers it. (3) `folding_dg` is computed whole-sequence
    but labelled `5' dG` (optimized `-39.0`, reported `-138.0` on a 156-nt CDS) —
    reported ≠ computed, invariant #2. (4) `avoid_internal_start` is infeasible on
    82–100% of 400–700 aa proteins (Met is single-codon; `MAAMG` reproduces it),
    and `InfeasibleError` names every *active* constraint rather than the culprit;
    the `relax()` promised by §4.2 does not exist and `OptimalityStatus.RELAXED` is
    defined but never used.
  - **The default operating point is out of family.** BT4's own
    `scripts/compare_tools.py` puts nine shipping tools at CAI 0.63–0.83 / GC
    42–54% and BT4 alone at **CAI 1.000 / GC 62.4%**. Relatedly, **no GC-content
    constraint exists at all**: the soft term is separable and saturates at weight
    2 without reaching its target, and the hard count budget cannot control
    clustering (74% worst 50-nt window at 50% total GC). A windowed-GC rule is
    listed in CLAUDE.md §6 and was never built, though `bt4 tracks` already
    computes the window.
  - **The engine is strong; nothing leads the user to it.** CpG 47 → 0 costs 0.077
    CAI and keeps a proven-optimal certificate, and a nine-flag AAV/LVV profile
    produces a clean 700-aa design (CAI 0.882, GC 53.0%, 50-nt window 42–64%, CpG
    5, zero hard violations) in 20.7 s — a combination no user will ever find.
    Hence named application presets.
- **Constitutional wording corrected to match the code** (§10.6 applied to BT4's
  own framing): §1 now claims "optimized for **expression-relevant objectives**"
  rather than "for real expression outcome", and states the CDS-in-isolation scope
  explicitly; the `io/` architecture diagram and the BT4 Studio export bullet no
  longer list **GenBank**, which is a design target that was never built (the
  README was already accurate on both counts).

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
    design target, so their codon usage is not evidence about the nuclear
    translation a design will meet. Both numbers are stamped —
    `genes_excluded_organelle_encoded` (what this filter removed) and
    `organelle_records_in_cds_source` (how many were in the annotation at all) —
    because most organelle CDS are dropped as invalid under the standard code
    before ever reaching the filter.
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
  - **Surfaces:** `--reference-set` on `bt4 optimize`, `library`, `validate` and
    `tracks` (there is no `bt4 frontier` subcommand; the frontier is reached via
    `api.frontier`); `bt4 organisms` now prints each organism's default and
    available sets; a **Reference set** picker in BT4 Studio that repopulates from the engine
    per organism (and disables itself, with a reason, when only one exists);
    `reference_set` on the service's optimize request and in `/organisms`;
    `OptimizeConfig.reference_set`; and `api.available_reference_sets` /
    `api.default_reference_set`.
  - **Honest scope, unchanged.** A highly-expressed reference set makes CAI a
    better-founded proxy, not a validated expression predictor — Welch et al.
    (PLoS ONE 2009) found an *E. coli* variant built by maximizing exactly this
    quantity expressed at a fraction of alternatives. CAI stays one axis of the
    objective vector (§10.7).

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

### Changed

- **The runbook's B1 now carries verified acquisition recipes for both panels** — pinned
  URLs with md5s, which GTF variant and why, exact column names, and the
  `for_zenodo` → `data` rename `splicebench2023`'s notebooks require.
  - **Recorded: 53% of `splicebench2023` is not held out.** BRCA1 (chr17), FAS (chr10) and
    WT1 (chr11) are on chromosomes both models trained on; only the chr3 genes are held
    out. That includes BRCA1, otherwise the closest public thing to BT4's regime.
  - **Recorded: SpliceAI's held-out split is second-hand.** Pangolin's is in its own
    paper; SpliceAI's Cell paper is paywalled, so chr1/3/5/7/9 comes from OpenSpliceAI
    (eLife 2025), a reimplementation that rebuilt its data pipeline.

- **Corrected a false claim BT4 made about its own metric.** `pr_auc_skill` was
  documented as "the only PR figure comparable across panels of differing prevalence".
  It is not. Rescaling average precision's floor to 0 does not remove its prevalence
  dependence: measured at **fixed model quality**, skill falls 0.589 → 0.152 as negatives
  go from 1k to 30k while ROC-AUC holds at 0.91. The metric is unchanged and still worth
  reporting — its floor and ceiling really are fixed — but the comparability claim was
  wrong, and a declared `negative_construction` is the only thing that makes two panels'
  numbers mean the same thing. A test now pins the prevalence-dependence so the claim
  cannot drift back.
- `SpliceSiteCase.group` no longer claims cases "are never split across folds"; this gate
  has no folds. It now says what the group actually does.

- `docs/NEXT_SESSION.md` item 11 now records that the calibration **machinery** has
  landed and only the data step remains; `CLAUDE.md`'s Phase 4 paragraph says the same,
  and both restate that RiboNN is still `calibrated=False` — the apparatus is not a claim.
- `THIRD_PARTY_DATA.md` records the attestation shape: derived scalars and public content
  hashes only, with RiboNN's non-commercial per-sequence outputs structurally excluded.
- **`ExpressionAttestation` -- the single seam that can flip an expression head to
  `calibrated=True`, replacing a bare `dataclasses.replace`.** Mirrors
  `bt4.biomodels.splice.attestation`, with the differences a *usefulness* claim needs
  that a *fidelity* claim does not.
  - **Licence-clean by construction.** RiboNN's weights are Sanofi non-commercial, so its
    raw per-sequence outputs must never enter MIT-licensed BT4. An attestation carries
    only derived scalars plus public content hashes (the 90 pinned weight SHA-256s for
    its species, and the panel's hash). `_ALLOWED_FIELDS` pins the shape, a module-level
    assert fails the *import* if the dataclass drifts from it, `from_dict` refuses any
    unexpected key, and a test asserts no field name could hold an array of scores.
  - **Four floors, so a run configured to pass cannot self-certify**: the run must have
    been **within-group** (a pooled run credits between-protein skill and cannot certify
    BT4's regime, however good it looks), it must have **beaten every baseline**, its
    interval must be informative (`width_over_iqr` < 1.0), and its `min_spearman` /
    coverage tolerance must clear `MIN_ATTESTATION_SPEARMAN` / 
    `MAX_ATTESTATION_COVERAGE_TOLERANCE`. All four are re-checked at promotion, so
    hand-editing the JSON afterwards buys nothing.
  - **Scope is part of the claim.** The record carries species, cell-type selection and
    readout, and `verified_predictor` refuses a predictor whose species or cell types
    differ — an attestation earned on HEK293T does not certify a head averaging all 78
    tissues, because those are different quantities. It also refuses unless the weight
    hashes match the adapter's own pins, binding the claim to the same bytes the adapter
    hash-verifies before loading.
  - `content_hash()` is timestamp-free (invariant #7) and scope-sensitive, so it is a
    stable provenance stamp and two different scopes cannot share one.
  - The default head is unchanged: `default()` still returns the uncalibrated
    placeholder, and a bare `RiboNNExpressionModel()` is still `calibrated=False`.
- **The gate is now a supported surface, not a snippet: `bt4.pipeline.expression_gate`,
  `api.expression_gate`, `bt4 expression-gate PANEL.tsv`, and
  `scripts/run_expression_gate.py`.** The orchestration lives in `pipeline/` (per §3, so
  the CLI reaches it through `api` without crossing a layer) and the script and the CLI
  both render the same `GateComparison` -- they cannot drift about what a result means.
  - **Every run is scored against five permanent baselines**: `permutation` (the null --
    the head's own predictions against a deterministic shuffle), `cai`, `gc3`, `length`,
    and `constant`. The reason they are not optional: a within-protein Spearman of 0.3
    is worthless if plain CAI scores 0.35, because BT4 already computes CAI *inside* the
    optimizer loop and for free. `constant` is there because split conformal is valid for
    any score function, so a predictor with no information achieves exactly valid
    coverage -- its "pass" belongs in the same table rather than being a trap the reader
    has to remember.
  - **`promotable` requires three things at once** and reports each separately, so a
    failure says which one failed: the gate's own thresholds, the head's bootstrap CI
    lower bound above *every* baseline's estimate, and an interval narrower than the
    label IQR.
  - **One backend invocation per UTR context, not per row.** A predictor carries its
    fixed UTRs on the model, so a panel spanning transcripts genuinely needs one
    predictor each -- and no more; each bucket is scored in a single batched call, and
    panel order is restored afterwards so a measurement can never be paired with another
    sequence's score.
  - Pooled mode still works (it is a useful contrast) but **warns on stderr** that it
    credits between-protein skill and is not the regime BT4 deploys in.
  - The report carries the panel's `content_hash`, the settings, and an explicit honesty
    note: nothing here flips a flag, `promotable` means "the pre-registered conditions
    held on this panel", and `min_spearman` is a pre-commitment rather than a community
    standard -- no such standard exists.
- **A measured CDS-variant panel format and a strict reader**
  (`bt4.biomodels.expression.panel`). The gate consumes in-memory triples, which is
  right for the gate and no help to a maintainer turning a published supplementary table
  into something runnable and provable months later. Tab-separated, `group` /
  `variant_id` / `cds` / `measured` / `utr5` / `utr3` required, with optional `readout` /
  `cell_type` / `species` carried through so a number is never separated from the
  question it answers.
  - **It refuses rather than copes, and that is the point.** RiboNN *silently drops* any
    row whose 5′UTR exceeds 1381 nt or whose CDS+3′UTR exceeds 11937 nt -- the caps are
    applied inside its data module, which filters the frame. A quietly-shortened panel is
    the worst possible gate input: the gate would report an honest `n_test` while
    answering a question about a dataset nobody chose. Such a row is now a hard error
    naming the row, as are a non-3N or stop-less CDS, a non-ACGT or empty UTR, a
    non-finite `measured`, a duplicate `variant_id`, an unknown species, and — so a
    mislabelled column cannot sit unused while the gate runs on nothing — an
    **unrecognised column**.
  - **`content_hash()`** is order-independent and timestamp-free (invariant #7), so a
    panel's identity can be pre-registered *before* a gate run and compared afterwards.
    Re-ordering or re-quoting a file does not change it; changing any value does.
  - **`describe()`** surfaces the sizing facts the gate's own arithmetic makes
    load-bearing: rows, groups, and **groups with 2+ members**, since a 90% conformal
    interval needs ≥ 9 calibration rows for a finite half-width, within-group scoring
    needs groups with 2+ members, and a grouped split needs 2+ groups. An unfit panel is
    visible before a single model runs.
  - **`contexts()`** buckets rows by their `(utr5, utr3)` pair. A predictor carries its
    UTR context on the model, so a panel spanning transcripts with different UTRs cannot
    be scored in one invocation; this makes that split explicit rather than accidental.
- **The expression gate can now judge a CDS-variant panel honestly: `within_group`,
  `recalibrate`, a cluster-bootstrap CI, and a vacuity check.** As written the gate
  could hand out a **false pass**, and each addition closes one way that happens.
  Defaults are unchanged, so every previous call behaves exactly as before.
  - **`within_group=True` -- the strict bar.** Pooled scoring computes one Spearman over
    the whole test fold; when the groups are proteins that fold mixes proteins with
    wildly different baselines, so a head that knows nothing about codons but recognises
    "this is a highly-expressed gene" scores well. That is precisely what training
    across natural genes teaches and precisely what BT4 cannot use. Within-group mode
    centres predictions and measurements inside each group and aggregates a per-group
    Spearman **unweighted across groups** (ProteinGym's aggregation), so a protein with
    30 variants cannot outvote one with 4. A regression test pins the defect: a
    gene-identity-only head **passes pooled and fails within-group**.
  - **`recalibrate=True` -- the fitted link.** A head whose output is in arbitrary units
    (RiboNN reports a CLR compositional residual) cannot be compared to an assay's units
    by subtraction. The affine link `measured ≈ slope × predicted + intercept` is fitted
    on the **calibration fold only** -- fitting it on the fold that is then conformalized
    would break the independence split conformal requires. `link_slope_spread` (the link
    refitted inside each calibration group) is reported *instead of* a calibration slope,
    which is 1.0 by construction once fitted and would be a circular pass.
  - **The rank metric stays on the head's *raw* predictions, deliberately.** Rank
    correlation needs no link, and it must describe what BT4 would really do: BT4 ranks
    candidates by the raw score and never applies a fitted link at design time. Scoring
    linked predictions would let a head that ranks **backwards** be rescued by a negative
    fitted slope and reported as passing, while a deployed BT4 handed the user the worst
    candidate. Pearson, R², the conformal residuals and the interval width *do* use the
    link, because they live on the measurement scale.
  - **`width_over_iqr` -- the vacuity check.** Split conformal is valid for *any* score
    function, so a **constant predictor achieves exactly valid coverage** with a useless
    interval. Median interval width over the label IQR is the number that exposes it, and
    a test pins that such a predictor is caught on both the rank and the width axis.
  - **A cluster-bootstrap CI on the primary metric**, resampling **whole groups** --
    variants of one protein are a dependent cluster, and resampling individual cases
    would treat 30 variants as 30 independent observations and produce a CI far too
    narrow. Seeded and deterministic (invariant #7); reports `nan` with
    `bootstrap_resamples=0` rather than a CI computed from too little.
  - **`coverage_conditional_on_group_anchor`.** In within-group mode the target is a
    variant's offset from its protein's own baseline, so the interval is only achievable
    at design time when a member of that protein has already been measured to anchor it
    -- BT4's `delta_logte` framing, but a **narrower claim** than an unconditional
    interval, so it is stamped as one rather than quietly conflated.
  - Also reported: `per_group_spearman`, `n_groups_test`, and `n_groups_ranked` -- the
    last being the effective sample size for a cross-group claim, which is the number of
    *proteins*, never the number of rows.
  - **`run_expression_gate` now uses `score_many`** where the backend offers it. Scoring
    a panel row-by-row through RiboNN would multiply the wall clock by the row count and
    re-hash 90 weight files each time.
  - New shared estimators in `bt4.biomodels._stats`: `linear_fit` (least squares, with
    the honest intercept-only answer when the predictor has no variance) and `iqr`.
- **`scripts/ribonn_sensitivity.py` -- the zero-data checks that decide whether a
  RiboNN calibration panel is worth acquiring at all.** Four checks, no measured data
  required, driven entirely from a RiboNN checkout plus sequences already in this
  repository. Every report is stamped `calibrated: False` and can promote nothing.
  - **`utr-control`** scores one CDS under two different UTR pairs. It is the positive
    control that makes a *null* result elsewhere trustworthy: if swapping both UTRs
    leaves the score untouched, the sequences are not reaching the model and no other
    number means anything. It refuses a backend with no UTR context rather than
    crashing on `dataclasses.replace`.
  - **`cds-spread`** is the decisive check and covers the GC-confound check in the
    same invocation. Holding UTRs fixed, it reports the spread of scores *within* each
    protein's synonymous variants against the spread *between* proteins, plus the rank
    correlation of the within-protein response against CAI, tAI, GC, GC3 and length --
    every feature recomputed here by BT4's own functions. `within_over_between` near
    zero means the backend reads gene identity rather than codon choice, which is what
    RiboNN was trained to do and precisely what BT4 cannot use. It runs on the in-tree
    `ranaghan2021_tab4.fasta` (93 records, three human proteins x 31 real
    codon-optimizer outputs, CC BY 4.0) -- which carries **no measurements**, so it is
    a sensitivity resource and never a validation panel.
  - **`direction`** builds a max-CAI and a min-CAI design per protein and runs an
    exact two-sided binomial sign test (`math.comb`, no scipy). **Ties are counted and
    excluded**, as the sign test requires: scoring them as failures would report a
    blind backend as "0/N prefer the optimized design", which reads as a strong
    preference for the deoptimized one.
  - **`ladder`** walks a real BT4 Pareto frontier for one protein and reports the
    Spearman of score against CAI along it -- a coherence check, since a jagged
    response is unusable for ranking even when it is nonzero.
  - House style throughout: build a dict, render it as a table or `--json`; UTRs are
    **required** rather than defaulted (a bundled UTR would be a hidden modelling
    choice) and are identified in the report by content hash rather than printed; the
    reference set travels with every CAI number; batched backends are driven through
    `score_many`. Tests pin the verdicts that matter -- a blind backend is reported as
    blind, a gene-identity backend is *not* credited with synonymous skill, and a
    GC3-only backend is exposed by the confound correlation -- and an end-to-end CLI
    rehearsal runs on the `null` placeholder, which needs no weights and is the
    reference for what "blind" looks like.
- **RiboNN gains cell-type selection and a fold-resolved read, closing two scope
  errors that would have corrupted a calibration gate before it ever ran.** Both are
  diagnostics and neither touches `calibrated`, which stays `False`.
  - **`cell_types`** picks which of RiboNN's per-cell-type outputs to average (78
    human / 68 mouse). Empty stays the default and averages all of them, which is the
    right summary for a generic design -- but comparing the mean of 78 tissues against
    a measurement from *one* cell line is a scope error, not a rounding error, so a
    HEK293T panel is now scored with `cell_types=("HEK293T",)`. An unmatched name
    **raises** and lists what is available rather than quietly averaging the wrong
    set. The units label names the selection, so "mean over all human cell types" and
    "mean over HEK293T" can never share a label in a report or manifest.
  - **`predict_folds()`** returns `RiboNNFoldPrediction(index, fold, te)` per
    (input, fold) instead of the fold mean. RiboNN emits one row per input per outer
    fold, each already that fold's `top_k`-model mean; averaging all ten is correct
    for the novel designed sequences BT4 produces -- no fold saw them -- and **wrong**
    for a natural transcript, where nine folds trained on its own label, making the
    averaged number optimistic and incomparable to RiboNN's published held-out
    accuracy. Keeping the fold identity is what makes the zero-cost
    adapter-validation check possible: score RiboNN's own published labels, keep the
    holdout fold, and confirm the held-out r² lands near the published value while
    the other nine sit visibly higher. If they are indistinguishable, the fold
    semantics are wrong and every downstream number is uninterpretable.
  - `_run_predict` now returns the raw table, and the fold-averaged and fold-resolved
    views are two consumers of that **one** invocation and code path, so they cannot
    drift; a test pins that the averaged view is exactly the mean of the resolved one.

- **The RiboNN adapter now forwards `batch_size` and `num_workers`, and a claim that
  it could not has been corrected in three places.** `ribonn.py`'s own comment,
  `CLAUDE.md` and `docs/DESIGN_expression_splice_flow.md` all stated that
  `predict_using_nested_cross_validation_models` "exposes no worker-count parameter"
  and that a `num_workers=0` path was therefore deliberately left out. Verified
  against upstream `src/predict.py`, that is **wrong**: the signature is
  `(input_path, species, run_df, top_k_models_to_use=5, batch_size=1024,
  num_workers=4)`. `RiboNNExpressionModel` gains both as validated fields
  (`batch_size=64`, `num_workers=0`), threaded through `resolve_backend` and into both
  `predict` call sites.
  - **Neither knob can change a score,** which is why lowering the defaults is safe
    rather than a silent behaviour change: RiboNN pads every transcript to a *fixed*
    width (`max_utr5_len + max_cds_utr3_len` = 13318, set in `predict.py`), not to a
    batch's longest member, and `RiboNNDataModule.make_dataloader` sets
    `shuffle=False` for every non-training stage (`reorder = stage == "train"`), so
    batch composition affects neither the one-hot encoding nor row order. They are
    memory and throughput only.
  - **`num_workers=0` is a correctness requirement, not tuning.** The adapter scores
    from a mutated `sys.path` and a temporary working directory
    (`_run_predict_with_models_layout`), neither of which a *spawned* worker inherits
    -- so `num_workers>0` hangs or fails wherever the multiprocessing start method is
    spawn (Windows, macOS). RiboNN also rebuilds the predict dataloader once per
    ensemble member (`top_k` x folds, up to 50 times), paying the spawn cost each
    time.
  - **`batch_size=1024` OOMs an ordinary CPU box** -- 1024 fixed-width
    `(channels, 13318)` float32 tensors at once, before worker prefetch.
  - `calibrated` is untouched and remains `False` for every configuration; a knob is
    not a gate (CLAUDE.md §10.6). New tests pin the passthrough, the validation, the
    defaults, and that the helper restores the working directory even when the
    upstream call raises.
- **`docs/NEXT_SESSION.md`'s RiboNN environment notes gain three upstream findings**:
  the licence is an *affiliation* grant ("any person from academic research or
  non-profit organizations"), not merely a non-commercial one; `max_shift` in the
  shipped `runs.csv` MLflow params is a determinism hazard because
  `_stochastic_shift` is not gated on `self.training` and uses an unseeded
  `torch.randint`; and native Windows is viable but unsupported upstream (no `make`
  needed, weights folder must be named `models`, `torch.load` is called with no
  `map_location`). `setuptools<81` corrected to `<82`, and the weights-extraction
  target clarified (`-C models`, since the zip root holds `human/` and `mouse/`).
- **Two load-bearing claims in the docs did not survive a fact-checking sweep, and
  are corrected everywhere they appeared.** A six-lens literature review (joint
  codon+structure design, 5'UTR-aware expression models, vector/AAV/LVV sequence
  hazards, synthesis-vendor thresholds, the tool landscape, expression
  determinants), each lens followed by an adversarial verification pass, re-checked
  every claim in `docs/RESEARCH_codon_optimization_SOTA.md`. **No source behaviour
  changed** -- the only file under `src/` touched is a module docstring.
  - **"The CDS is a minority of the expression signal (~31%)" was a misreading of
    RiboNN**, and it was repeated in seven places (`CLAUDE.md`, `README.md`,
    `COMPARISON.md`, `DESIGN_expression_splice_flow.md`, `NEXT_SESSION.md`, the
    survey, and `biomodels/expression/ribonn.py`). RiboNN reports **two**
    attributions: per-nucleotide information density 67 / 31 / 2 and
    **length-integrated total attribution 22 / 73 / 5** (5'UTR / CDS / 3'UTR,
    human; mouse 23 / 73 / 4). Per base the 5'UTR is denser; **integrated over
    length the CDS carries ~73% of the attributed translation-efficiency signal**.
    Quoting only the first pair argued against BT4's own existence. The honest
    ceiling is a different fact and is now stated as such: mRNA *abundance* rather
    than translation rate is the majority channel for protein abundance (Li, Bickel
    & Biggin 2014, >=56%), integration site alone spans ~1,000x (Akhtar 2013), and
    Kozak context alone spans ~100x (Shukla 2026). The sharper reason RiboNN is not
    calibrated for BT4 is now stated too: it has **never been shown to discriminate
    synonymous CDS variants of the same protein under a fixed UTR**.
  - **LinearDesign's "up to 128x" is not an expression result and has not
    replicated.** The 57-128x is anti-spike IgG in n = 6 mice against a *vendor
    codon optimizer's* output; HEK293 protein was **2.9x**, and the durable claim is
    stability. Against it: an independent 2026 mammalian bake-off in which
    "strategies prioritizing RNA stability consistently reduced expression" and
    LinearDesign gave the lowest yields; monosome-dominated polysome profiles from
    MFE-optimized mRNA (Leppek 2022); "folding free energy shows only a weak
    correlation with in-cell lifetime and protein expression" (Jin 2025); and the
    observation that MFE minimization is substantially a **GC-maximization proxy**,
    so folding dG and GC are **not** independent Pareto axes. Joint codon+structure
    optimization is therefore **demoted** from the survey's #1 recommendation and
    from Tier 5 of the review; the part with uncontested evidence (position-
    dependent structure, which needs a real cap distance) moves into the context
    work instead.

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

- **Python 3.10 is now supported** (was 3.11+). `requires-python` is lowered to
  `>=3.10`, the 3.10 classifier is added, ruff/mypy target 3.10, and CI's quality
  matrix now runs 3.10 alongside 3.11–3.13. The pure core uses no 3.11-only
  features, so this is a compatibility widening with no behavior change. It notably
  lets the wrapped **RiboNN** expression backend be installed into the same
  environment as its own dependency stack, whose pinned `torch==1.13.1` ships only
  CPython ≤3.10 wheels.

### Fixed

- **A second adversarial review, over the ~1,200 lines written since the first, found ten
  more defects — including two the per-kind anchor rewrite introduced.**
  - **Declaring an anchor crippled every sequence-derived control.** Case building moved
    each site's *label* into the backend's frame, but `gt_ag` and `pwm` were still read at
    the case's own index — the backend's anchor, not the panel's. Measured: with the head
    identically perfect, the `pwm` control fell from **0.853 skill to 0.0001** and `gt_ag`
    to exactly 0, making `beats_every_baseline` near-automatic for precisely the real
    backends that need an offset. All four review lenses found it independently. The
    controls now invert the shift, unioning both kinds for a combined stratum.
  - **The `"splice"` offset key was inert but reported as applied.** `{"splice": -5}`
    scored identically to `{"splice": 0}` while the diagnostic claimed each value was in
    force and recommended a correction derived from that fiction — so the advertised
    round-trip destroyed a correctly aligned run. Offsets are now keyed by *site* kind
    throughout and the key is refused.
  - **A shifted-out-of-window site silently left the positive class**, lowering prevalence
    and *raising* every metric — rewarding a backend that structurally cannot find an
    annotated site. It is clamped back to the nearest scoreable index, preserving the
    forced miss the old formulation produced.
  - **`aligned` could be carried by a single peak.** `n_flat < n_sites` let one site
    outvote any number of silent ones; a kind where 7 of 8 produced nothing reported as
    aligned at 12% agreement. It now needs a majority.
  - **The GENCODE builder's motif tolerance discarded small genes.**
    `len(bad) > len(same) * 0.1` reads as 10% slack, but a k-intron gene contributes only
    2k sites, so the allowance was 0.2k — below **1** for any gene with four introns or
    fewer. One legitimate GC-AG intron discarded every site of the gene, the opposite of
    absorbing the minor spliceosome it was meant for.
  - **A window running off a contig end was written silently.** Python slices truncate,
    and minus-strand indices are measured from `w_end`, so every one shifted. The window
    is now shrunk to what the assembly has and the indices recomputed; a truncation that
    cuts into the transcript is skipped and counted.
  - Neighbouring sites with no room for their own dinucleotide at a window edge are
    dropped rather than written into a panel BT4's own reader would refuse.
  - **A ragged row crashed the splicebench converter.** `csv.DictReader` fills a short
    row's missing fields with `None`, so `_boolean(None)` raised `AttributeError` instead
    of counting the row unparseable — and a 972-column published supplement can easily
    have one. Null score spellings (`NA`, `NaN`, `None`, `-`, `.`) are now all matched.
  - A whitespace-only key defeated the `variant_id` fallback (truthiness was tested before
    stripping), emitting an empty id that made the **whole** panel unreadable.
  - A genuine score of `0.0` was replaced by a falsy fallback, and `describe()` published
    a measured-looking prevalence of `0.0` for an empty region.
- **`anchor_offset` was a single scalar, and no scalar is correct for a real backend.**
  Research into SpliceAI's and Pangolin's source established that **both anchor on the
  exonic boundary base**, one base *before* BT4's donor position and one base *after* its
  acceptor position — opposite directions. Measured with a perfect exonically-anchored
  backend: at `-1` donors score AP 1.000 and acceptors 0.006; at `+1` the reverse.
  - **And the alignment diagnostic endorsed the broken setting.** Under either value half
    the sites aligned and half landed two bases off, and the modal tie-break resolved
    that `{0: N, ±2: N}` split to `0` — printing **"anchors agree" at 50% alignment**
    while a perfect model read as hopeless on half the panel. The one check meant to
    catch misalignment was confirming the wrong value.
  - `anchor_offset` now accepts a **per-kind mapping**, `CNN_ANCHOR_OFFSETS` records the
    verified values, and the diagnostic is **per site kind** so the false-agreement state
    is unrepresentable. Case building shifts each site's label by *its own kind's* offset
    before any union, which is the only formulation that can handle Pangolin's combined
    track (its union drops the kind).
  - `bt4 splice-gate` gains `--cnn-anchors`, `--donor-offset` and `--acceptor-offset`.
    A scalar remains valid: a kind-separated panel is a legitimate way to run a real
    backend, and `read_splice_panel` already accepts donors-only or acceptors-only
    windows.

- **A second round of adversarial review found five more ways the splice gate's verdict
  could be won without winning.** All measured, all in the Part B gate.
  - **Dropping the `pwm` control let BT4's own default certify itself.**
    `SPLICE_BASELINES` says the controls are "kept permanently: a control that
    disappears once it is inconvenient was never a control" — and nothing enforced it.
    Running the PWM backend with `baselines=("permutation","gt_ag","constant")` reported
    `promotable=True`, one keyword away from the module's headline property. Three
    review lenses found this independently. A subset run is still allowed; it just
    cannot recommend anything.
  - **An uncontested stratum counted as a beaten one.** With `baselines=()` no comparison
    ran, `beats_every_baseline` kept its initial `True`, and `best_baseline` reported
    `-inf` as though it were a control's skill.
  - **A GENCODE/Ensembl-named training panel reported itself held out.** The check
    matched only `chr`-prefixed names, so a panel drawn *entirely* from the models'
    training chromosomes but named `2`, `4`, `X` reported `held_out=True` — the naming a
    builder following the runbook's own GENCODE recipe is most likely to produce. Names
    now normalize across conventions, and an **unrecognisable** group (a RefSeq
    accession, a scaffold) is reported as *unknown* rather than silently clean.
  - **Float dust defeated combined-track detection.** `_is_combined` tested exact
    `== 0.0`, so an acceptor channel carrying `1e-12` was scored with a donor/acceptor
    split — the exact artifact the collapse prevents, and silently.
  - **Sites too close to a window edge were scored as silent forced misses.** A donor at
    position 0 carries a real `GT`, so the motif check passes and the panel is accepted —
    but no backend has flanking sequence there, so the PWM returns exactly `0.0`. That is
    a `label=1` case the model structurally cannot get right, depressing every metric
    through no fault of its own. `SplicePanel.edge_sites()` now reports them and the
    runner notes the count, because nothing is wrong with the panel's *labels* and the
    cure is a wider window rather than a dropped site.
  - **A silent backend was credited with perfect alignment.** The alignment probe's
    tie-break resolved a flat window to offset 0, so a backend emitting no signal read as
    "anchors agree" — a positive claim from an absence of evidence, in the one diagnostic
    whose job is to separate "misaligned" from "scoring near zero".

- **One calibrated window certified a whole run.** On the caller-supplied-scores path
  `backend_calibrated` was an `any()` over per-window results, so a single
  `SpliceResult` carrying `calibrated=True` reported the entire run as calibrated. It is
  now `all()`, guarded by `bool(scored)` because `all(())` is `True`. A calibration flag
  is the last place to take the generous reading (§10.6).
- **`promotable` was reachable without anyone setting a bar.** The gate ships permissive
  defaults on purpose — a threshold this module blessed would be one a weak backend
  could be pointed at — but that left `passed` trivially true for a bare call, making
  "the pre-registered conditions held" vacuous. `promotable` now also requires
  `thresholds_declared`, and a bare call says so in its notes. The defaults are
  unchanged; what changed is that they can no longer be mistaken for a verdict.
- **The `pwm` control went blind to every acceptor site when the strata collapsed.**
  `_tracks` collapses donor and acceptor into one `"splice"` stratum for a
  combined-track backend, which is right for the *head* — Pangolin's acceptor track is
  identically zero, so its donor track carries everything. But the **baselines** are
  scored through the same function, and BT4's PWM baseline is a real two-track
  predictor. Taking its donor track alone made half the positives unreachable for the
  control: measured, **0.344 skill instead of 0.654** on a mixed panel, so a
  Pangolin-class head was easier to beat than it should have been — in exactly the mode
  Pangolin runs in, and against the very claim the module is built on. The collapse now
  **unions** the two tracks, which is a no-op for a genuinely combined backend and
  restores the control otherwise. Found by adversarial review.
- **The alignment diagnostic reported a residual offset as if it were absolute.** The
  probe runs on the already-aligned track, so its modal offset is the shift *remaining*
  after `anchor_offset` is applied. The note printed that number as the value to
  declare, so a user who had already passed `+1` against a true `+3` was told to pass
  `+2` — further from the truth, in a message they had every reason to trust. It now
  reports `recommended_offset` (applied + residual) and names the offset currently in
  force, and the agreeing message states which anchor it agreed at.

- **Three confirmed-open findings from the completed verification workflow.** The
  exhaustive multi-lens review (seven lenses × three refuters per finding) was
  re-run to completion after its first pass hit a usage limit; its honest tally
  was **6 confirmed, 36 refuted, 0 unverified** (the earlier run's "38 refuted"
  had silently counted never-run verifiers as refuted — the aggregation now
  distinguishes *unverified* from *refuted*). A completeness critic then
  re-checked the six confirmed against what actually shipped. Four were already
  fixed; three were genuinely open and are closed here:
  - **The `load_provenance` fail-open guard was untested — and the test named for
    it was vacuous.** A test titled "the provenance guard must fail closed" only
    checked `write_table`/`load_table_from_file` and never called `load_provenance`
    with a keyless sidecar, so reverting the guard to its fail-open form left the
    whole suite green. It is renamed to what it actually checks, and real tests now
    exercise the guard at the loader — each verified to fail under the exact
    one-line revert.
  - **`load_table` mislabelled a custom table's reference set.** The label was
    derived from the *filename* only in `load_table`, and cross-checked against the
    sidecar only in `load_provenance` — so a user's `build-table` output (honestly
    stamped `custom`) dropped in at `<organism>.tsv` loaded as `genome_wide`, and
    `bt4 tracks` (which reads `load_table`, never `load_provenance`) would print
    that false label while `bt4 optimize` correctly refused. The sidecar check is
    now a shared helper both loaders call, so `tracks` is as strict as `optimize`:
    the two paths can no longer disagree about a sequence's declared reference set.
  - **The candidate-assembly reference-set threading was untested.** The library
    and tracks revert-detectors were effective, but the candidates one was not: its
    only test used a config with no GLOBAL rule, so zero repeat-refined variants
    were produced and the reference-set-carrying table was never consumed — a
    revert left CI green while a `genome_wide` run with `max_repeat_length` active
    would have scored its repeat-refined variants against the highly-expressed
    table (invariant #2). A new revert-detector forces the variants and asserts
    invariant #2 for every member.

  Every new test was verified to **fail under the exact mutation it guards** —
  the lesson from the vacuous test above.

- **Second review round: seven more lenses, seven more defects.** A follow-up
  seven-lens pass (data, plumbing, claims, test quality via source mutation,
  honesty invariants, behavior-change blast radius, packaging/CI) ran over the
  fixes above. It hit a usage limit partway, so its own "refuted" tally conflates
  *disproved* with *never checked* — every finding below was therefore re-verified
  by hand before being acted on.
  - **`load_table("homo_sapiens.highly_expressed")` returned the right counts
    under the wrong label.** The genome-wide suffix is the empty string, so a
    dotted key resolved straight onto the highly-expressed *file* and came back
    stamped `reference_set="genome_wide"`. An organism key is now required to be a
    bare stem, exactly as `available_organisms()` lists it — nothing downstream
    could have detected that mislabel.
  - **The non-optimize half of the reference-set axis had no test at all.**
    Deleting `reference_set=` from `api.library`, `api.candidates`, `api.tracks`
    and all three comparison scripts left the suite green — so the invariant-#2
    fix above could have been reverted without CI noticing. New
    `tests/test_reference_set_reaches_every_surface.py` is shaped as a revert
    detector: one failing test per dropped call site.
  - **BT4 Studio's reference set stuck after visiting *A. thaliana*.** Selecting
    the one organism with a single reference set forced the combo to
    `genome_wide`; switching back kept it, because a forced value was
    indistinguishable from a user preference. Two clicks reproduced the same
    "silently hands out a codon-commonness index" failure the organism default
    had. The app now remembers the user's explicit pick separately.
  - **The delivered result was labeled from the live controls.** Runs are
    asynchronous, so changing the organism mid-run relabeled the result with
    tables it was not built from. The pair is frozen at run start now.
  - **The "peptide FASTA must come from the CDS release" rule was prose only.**
    The two URLs are pinned in different files and both digest checks pass after a
    divergence, while versioned gene IDs from two releases silently fail to join.
    It is a structural check now.
  - **The machine-readable comparison boards carried unlabeled CAI.** Only the
    human-readable branch named its tables. Both scripts now emit
    `{organism, codon_reference_set, rows}` — a shape change to `--json` — and
    `compare_reproducibility.py` gained the `--reference-set` flag it was missing.
  - **`THIRD_PARTY_DATA.md` did not record PaxDb** (CC BY 4.0 — attribution is a
    license obligation), and its codon-table row still described the tables as
    "Kazusa-style representative published frequencies", stale since they were
    recounted from Ensembl. Both rows are now accurate and split by reference set.
  - **Three of my own CHANGELOG claims were wrong.** The "abundant enough to reach
    a top-300 list" rationale survived in the Added section while the Fixed section
    said it had been removed; the `%MinMax 56.3` figure was not reproducible (the
    real value is **57.81**, and the entry now names the exact config that shows
    it); and the Welch "worst-performing" correction was listed as done when an
    exception had aborted that edit before it landed.
- **Adversarial-review fixes, before the change shipped.** Three independent
  reviewers (data/provenance, `reference_set` plumbing, and a numeric audit of
  every prose claim) found nine real defects; all are fixed here, and the
  shipped tables are byte-identical before and after, which is what shows none of
  them silently altered a number.
  - **The ambiguity guard did not gate the join.** An identifier resolving to two
    genes was flagged, but only consulted *after* both lookups missed — so a key
    ambiguous in the exact map yet singular in the version-stripped one was joined
    anyway and stamped as clean. Ambiguity is now judged across both maps (pooled
    per key, which is what actually catches the WormBase-style `ZK1010.1.1` /
    `ZK1010.1` collision the docstring described) and checked *before* either
    lookup. Latent for the eight shipped tables, provably: their TSVs are
    unchanged by the fix.
  - **`filters.cds_counted` was the genome-wide gene count**, shipped beside
    `cds_count: 300` under a key that means "the number counted into this table"
    in the genome-wide sidecars — one key name carrying two quantities. The block
    is now `cds_source_filters` with `genes_with_a_representative_cds`, naming
    what it actually describes.
  - **The join tally's key names implied arithmetic that does not close.**
    `rows_matched_*` and `rows_unmatched_*` now partition the abundance rows
    exactly, `rows_matched_whose_gene_has_no_counted_cds` is named as the subset
    it is, and `genes_*` are marked as a different unit.
  - **`excluded_organelle_encoded: 0` was misreadable** as "no organelle-encoded
    genes here". Most organelle CDS never reach that filter — under the standard
    code they read as having internal stops and are dropped as invalid first — so
    the source's own tally is now stamped alongside it
    (`organelle_records_in_cds_source`, 13 for human: exactly the human
    mitochondrial protein-coding gene count).
  - **The comparison scripts recomputed against the wrong table.**
    `compare_tools.py`, `compare_reproducibility.py` and `benchmark.py` all
    dropped `cfg.reference_set` when re-scoring, so every row on the published
    board was measured with a different table than the one that produced BT4's
    sequence (invariant #2). Reproducible case: `benchmark.py` on protein
    `MKTAYIAKQRQISFVKSHFSRQ` under
    `OptimizeConfig(organism="saccharomyces_cerevisiae", reference_set="genome_wide")`
    reported a %MinMax mean of **57.81** where the correct value is **100.0**.
    Only callers passing a non-default `reference_set` could hit it — which is
    every caller from the moment this same change added the `--reference-set`
    flag, so it is a defect introduced and fixed in one release, not a latent
    one that ever shipped.
  - **`load_provenance` failed open** on a sidecar with no `reference_set` key,
    validating it as whatever it was loaded as. The field is now required, and
    `write_table` gained it as a reserved parameter so it can no longer be
    smuggled through the free-form `extra` dict.
  - **BT4 Studio launched on the one organism with no highly-expressed table.**
    The organism combo took the alphabetically first entry, *A. thaliana*, so a
    freshly opened Studio silently handed out a codon-commonness index. It now
    starts on the engine's default organism.
  - **A stale tooltip contradicted the control it sat on**, and renders read the
    live combos instead of what the delivered run used — so a theme switch after
    changing the organism would recompute an old result's tracks under a new
    table. Renders now read the delivered run's own labels.
  - **Prose claims that outran the data:** "both are GC-richer" (the Ser change
    `AGC`→`TCC` is GC-neutral; the whole GC rise is the Arg change), "the genome
    at large does not" (false for three of the thirteen asserted optimal codons),
    "abundant enough to reach a top-300 list" (the exclusion is justified by the
    genetic code, not by an unmeasured ranking claim), "BT4 never fetches anything
    at runtime" (the opt-in ASSP cross-check does), and Welch et al.'s result
    stated as "worst-performing" rather than what the paper supports.
- **`bt4 --help` and `bt4 tracks --help` crashed** with
  `ValueError: unsupported format character 'M'`. argparse %-formats help strings
  *while rendering help*, and the `tracks` parser wrote a literal percent sign as
  `%MinMax` instead of `%%MinMax` — in both the subcommand's own `--organism`
  help and its one-line summary, which the top-level listing reprints. So the
  first command a new user types has been broken, and every other test passed
  because nothing in the suite invoked `--help`. Fixed, and every command's
  `--help` is now a parametrized regression test.

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
