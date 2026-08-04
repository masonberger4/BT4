# Changelog

All notable changes to BT4 are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) once it cuts
its first tagged release.

## [Unreleased]

### Added
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
  which degrades to the baseline alone — and says so — when Pangolin is absent.
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
