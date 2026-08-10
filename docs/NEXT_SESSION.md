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

**Phases:** 0–2 complete · Phase 3 groundwork landed · Phase 4 not started
(scaffold only) · Phase 5 opened. All merged and green on `main`.

**Status vocabulary:** `DONE` · `GROUNDWORK` (contract + baseline shipped,
calibration pending) · `BLOCKED-data` (needs a matched-regime panel) ·
`BLOCKED-human` (needs licensed weights / a maintainer machine) · `NOT-STARTED`.

| Component | State | Calibrated? | Primary file(s) |
|---|---|---|---|
| Exact-DP codon trellis + certificate | DONE | n/a | `optimize/exact_dp.py` |
| Rust trellis port (`trellis_solve`, regime-gated) | DONE | n/a | `rust/bt4_core`, `bt4_native` |
| Objectives: CAI, tAI, GC, ramp, CpG, %MinMax, codon-pair | DONE | n/a | `objectives/` |
| tAI (real GtRNAdb, 8 organisms) | DONE | n/a | `biomodels/codon/tai.py` |
| Constraints: homopolymer, GC-run, max-repeat, tandem/inverted, forbidden+presets, restriction, Kozak-ATG, uORF, splice-motif | DONE | n/a | `constraints/` |
| Budget backends: CP-SAT, Lagrangian, dinucleotide-count | DONE | n/a | `optimize/{cpsat,lagrangian}.py` |
| SA refinement + block moves + parallel tempering | DONE | n/a | `optimize/anneal_refine.py` |
| Folding (ViennaRNA + labeled baseline) | GROUNDWORK | ViennaRNA=yes, baseline=no | `biomodels/folding/` |
| Splice PWM baseline | GROUNDWORK | no (baseline) | `biomodels/splice/` |
| Splice CNNs: Pangolin (GPL) / SpliceAI (CC BY-NC) | GROUNDWORK | **no** (fidelity gate pending) | `biomodels/splice/{pangolin,spliceai}.py` |
| Splice audit (localize-and-flag) + backend agreement | DONE | advisory (`all_calibrated=False`) | `biomodels/splice/audit.py` |
| Splice fidelity-attestation layer | DONE (unused until a gate runs) | n/a | `biomodels/splice/attestation.py` |
| **ASSP cross-check (opt-in, out-of-loop network validator)** | DONE | `network_derived`, not calibrated | `biomodels/splice/assp.py`, `pipeline/splice_crosscheck.py` |
| Expression: `ExpressionPredictor` + `NullExpressionModel` + rerank hook | GROUNDWORK | placeholder=no | `biomodels/expression/`, `pipeline/rerank.py` |
| Expression: wrapped RiboNN (Sanofi non-commercial) | GROUNDWORK | **no** (acceptance gate pending) | `biomodels/expression/ribonn.py` |
| Candidate-set assembly + expression rerank | DONE | calibrated-gated | `pipeline/candidates.py`, `bt4.api.candidates` |
| Library / degenerate-design (SAMPLED) mode | DONE | n/a (sampler, not optimizer) | `optimize/sample.py`, `pipeline/library.py` |
| Surfaces: `bt4.api`, `bt4` CLI, FastAPI service, provenance | DONE | n/a | `api/`, `cli/`, `service/`, `provenance/` |
| **BT4 Studio** incl. Candidates & splice-audit tab | DONE | n/a | `app/studio.py`, `app/worker.py` |
| Packaged installers (PyInstaller/Briefcase) | NOT-STARTED | n/a | `packaging/` |

The **expression/splice design flow**
([`DESIGN_expression_splice_flow.md`](DESIGN_expression_splice_flow.md)) is
**complete through step 5** — batched RiboNN scoring, the splice-consensus motif
constraint, candidate assembly + rerank, the localize-and-flag splice audit, and
the BT4 Studio UI that surfaces them. Only step 6 (auto-edit / auto-select)
remains, and it is calibration-gated (see the queue). The opt-in **ASSP**
cross-check (the last named Phase-3 "still ahead" item) has also landed as a
network validator.

---

## Next-task queue

Ordered. Each item is tagged by precondition. **Pick the first `self-contained`
item unless you have a reason not to.**

1. **[START HERE · self-contained] Phase-5 organism breadth** — add more organisms
   with authoritative provenance via `bt4 build-table` on public CDS sets
   (content-hash the sources; never fabricate a table), and grow the
   restriction-enzyme / REBASE catalog. Fully autonomous, always welcome. Tables
   feed CAI (and tAI where GtRNAdb data exists); wire each new organism through
   `available_organisms()` and add its provenance sidecar.
2. **[self-contained] BT4 Studio polish** — light/dark theming, accessibility,
   responsive layout (Phase 4 polish, §6.6); a library-mode control (`api.library`)
   is a natural add.
3. **[self-contained] External-validation report** — compare BT4 output
   codon/GC/CpG distributions against real highly-expressed gene panels (§8), using
   public data and BT4's own recompute functions.
4. **[self-contained → then human] Packaged installers** — PyInstaller/Briefcase
   for macOS/Windows/Linux. Advance up to the point where signing / tag-pushing /
   release-cutting is needed; those steps are human-only (HTTP 403 in the sandbox).
5. **[BLOCKED-human] Promote the splice CNNs to `calibrated=True`** — capture
   reference panels and run `verify_pangolin_fidelity` / `verify_spliceai_fidelity`,
   then commit a `FidelityAttestation`. Needs the licensed weights on a maintainer
   machine. Never assign `calibrated=True` by hand — it is earned on data
   (§10.5/§10.6).
6. **[BLOCKED-data · human] Promote RiboNN to `calibrated=True`** — assemble a
   license-clean, regime-matched **CDS-variant** TE panel and run
   `verify_expression_gate` (Spearman + split-conformal coverage on a group-disjoint
   split). Reproducing RiboNN faithfully is **not** calibration for BT4's
   CDS-variant regime (its ablation puts only ~31% of per-nt signal in the CDS). Do
   not relabel a hand-weighted composite as "calibrated".
7. **[BLOCKED until #5/#6] Design-flow step 6** — targeted synonymous splice
   **auto-edit** and RiboNN **auto-select**, each unlocked only once its backend
   passes its gate.

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
`numpy<2` (torch 1.13.1 ABI) and `setuptools<81` (its older `pytorch_lightning`
calls `pkg_resources`), and the Zenodo `weights.zip` extracted to a directory
literally named `models/` under `$BT4_RIBONN_DIR` (so the hard-coded
`models/<species>/<run_id>/state_dict.pth` path resolves without a symlink).
Scoring requires **non-empty** `utr5`/`utr3` (empty → refused; the UTRs carry most
of RiboNN's signal). Weights are non-commercial — never bundled or CI-run.

---

*History (what each session shipped) lives in
[`../CHANGELOG.md`](../CHANGELOG.md) — this file carries only current state and the
forward queue.*
