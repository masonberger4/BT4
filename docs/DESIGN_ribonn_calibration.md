# Calibrating the RiboNN expression head

> **Status: machinery landed (PR #85), data step NOT STARTED.** This is the
> procedure for `NEXT_SESSION.md` queue item **#11**. The evidence behind every
> claim here lives in
> [`RESEARCH_ribonn_calibration.md`](RESEARCH_ribonn_calibration.md); this file is
> the runbook. Live status stays in [`NEXT_SESSION.md`](NEXT_SESSION.md), not here.
>
> **New to this job? Read [`GUIDE_ribonn_calibration.md`](GUIDE_ribonn_calibration.md)
> instead** — the same procedure in plain language, step by step, with the free checks
> pulled to the front. Every command in it was verified against the code (several were
> executed), and its **Appendix B** lists the places where this file is wrong or
> under-specified — including a `bt4 expression-gate` invocation that cannot produce the
> record it asks you to keep, and a panel-size floor at which a *good* head fails more
> often than it passes.

## Context

BT4 wraps **RiboNN** (Zheng et al., *Nat Biotechnol* 2026; Sanofi / Cenik Lab)
behind its `ExpressionPredictor` contract. The adapter is written, tested, and
drives *your own* installed copy of the model, hash-verifying all 90 per-species
weight files before loading. But it ships **`calibrated=False`**, and per CLAUDE.md
§10.6 that flag is *earned on data, never assigned*. Until it flips:

- `bt4.biomodels.expression.default()` keeps returning the neutral placeholder;
- the Candidates tab and `api.candidates` show RiboNN's scores as **discovery
  order, not a ranking**, with the solver's pick delivered;
- `rerank_by_expression` annotates the frontier but never re-picks the delivered
  point;
- design-flow **step 6** (RiboNN auto-select) stays blocked.

The apparatus to earn that flip is all on `main` after PR #85. What remains is a
**human-run, data-gated** procedure against the licensed weights. This document is
that procedure.

## The one thing to understand before starting

The word "calibration" is doing three different jobs, and separating them is most
of the work. RiboNN already does job 1. This runbook is about job 2, and job 3 only
if job 2 succeeds.

| # | Question | Plain meaning | Status |
|---|---|---|---|
| 1 | **Runnable?** | Does BT4's wrapper drive RiboNN and get numbers out of it? | ✅ real end-to-end runs happen against the licensed weights. **Not** a fidelity claim — unlike the splice side there is no RiboNN fidelity gate, capture script or attestation, so "reproduces upstream bit-for-bit" is *not* established (Stage 1.2 exists precisely because the fold semantics are unproven) |
| 2 | **Discriminating?** | Within one protein, UTRs fixed, does its ordering of synonymous variants match reality? | ❓ **the real question** |
| 3 | **Calibrated (strict)?** | When it emits a number, is its error bar honestly right ~90% of the time? | ❓ needs #2 first |

Two things people mistake for calibration, and are not:

- **Reproducing the published model is not calibration for BT4.** It proves the
  wrapper is correct. The research doc puts it plainly: *"Feeding RiboNN a real UTR
  makes it runnable. Neither is a gate."*
- **Doing well across different genes is not calibration for BT4.** That is the job
  RiboNN was built and validated for (r² ≈ 0.62 on natural genes), and it is a
  *different job* from choosing between synonymous variants of **one** protein under
  a fixed UTR — which is exactly BT4's regime, and which the RiboNN paper never
  tests.

**Doing job 1 and calling it calibration is how you ship a confident wrong answer.**

### Why job 3 needs an extra step

RiboNN outputs **CLR-residual TE** — a centred, effectively unitless quantity. Any
panel you find measures something else (log₂ protein output, ribosome load, an
mRNA-normalised luciferase ratio). Those are different rulers.

- **Ranking (job 2) survives a ruler change.** Spearman only cares about order, so it
  is the honest primary metric.
- **Error bars (job 3) do not.** So on the *calibration* half of the data you fit the
  straight line `measured ≈ a·predicted + b`, then judge error bars on the *test*
  half through that line. Fitting that line **is** the literal act of calibrating
  RiboNN to your assay — and it must be fit only on calibration proteins, and the `a`
  and `b` reported. The gate does this when you pass `--recalibrate`.

## Findings that shape the job

These are the ones that change what you actually do; the full evidence and the
corrections it forced on BT4's own docs are in
[`RESEARCH_ribonn_calibration.md`](RESEARCH_ribonn_calibration.md).

1. **The honest prior is that this gate fails.** RiboNN scores r² = 0.62 on natural
   genes, **0.17–0.19 zero-shot on designed reporters**, and **0.11** on the only
   CDS-attributable test its own paper reports. Design the work so a *null result is
   a cheap, recorded, publishable outcome* — several free checks below can end the
   project decisively before you spend anything.
2. **No public dataset fully qualifies.** There is no mammalian panel of synonymous
   CDS variants of one protein under a fixed UTR, with per-variant sequences *and*
   measurements downloadable at usable n. The two Ranaghan panels already in this
   repo carry **no measurements** (that paper measured one sequence, in *E. coli*),
   so they are a sensitivity resource only — never a validation panel.
3. **The licence is an *affiliation* grant, not merely non-commercial.** RiboNN's
   `LICENSE` and weights licence grant use "to any person from academic research or
   non-profit organizations." Non-commercial *intent* may not qualify an unaffiliated
   maintainer. Resolve this in writing **before** downloading the weights (Stage 0.1).
4. **`max_shift` is a determinism hazard.** RiboNN's `_stochastic_shift` is not gated
   on `self.training` and uses an unseeded `torch.randint`. A nonzero `max_shift` in
   the shipped `runs.csv` makes *inference* non-deterministic and breaks invariant #7.
   Check it before trusting any number (Stage 0.3).
5. **Windows works but is not upstream-supported.** You never need `make`; BT4 imports
   `src.predict` in-process. But the weights folder must be named `models`, and
   `num_workers` must be `0` (now the adapter default) or spawned workers hang.

---

# Stage 0 — Prerequisites

## 0.1 — Resolve the licence *before downloading anything*  ⚠️ blocking

Email `patent.gos@sanofi.com`, state your affiliation and intended use, and get the
answer in writing before you download the weights. If you are not affiliated with an
academic or non-profit institution, you may have no licence at all, and everything
downstream depends on this.

## 0.2 — Install (Windows `cmd.exe`)

RiboNN is a conda/mamba repo. You never run its `Makefile` — its `install` target is
only `mamba env create -f environment.yml -y`, and BT4 puts your checkout on
`sys.path` and imports `src.predict` directly. RiboNN documents no Windows support;
its one tested environment is Ubuntu 20.04 + one NVIDIA GPU. The native path below
works because BT4 bypasses the shell tooling.

```bat
git clone https://github.com/Sanofi-Public/RiboNN.git C:\RiboNN
cd /d C:\RiboNN

REM --- pick ONE env line ---
REM (a) You have an NVIDIA GPU: use the upstream pins as-is.
mamba env create -f environment.yml -y
REM (b) No NVIDIA GPU, or (a) fails to solve on win-64: CPU-only, same pins minus CUDA.
mamba create -n RiboNN -c conda-forge -y python=3.10.13 pytorch=1.13.1 ^
    pytorch-lightning=1.8.5 torchmetrics=1.3.1 lightning-utilities=0.10.1 ^
    mlflow=2.18.0 numpy=1.22.4 pandas=2.2.3 scikit-learn=1.0.2

conda activate RiboNN
pip install "setuptools<82"          REM upstream issue #10 / PR #11

REM --- weights: Zenodo record 17258709. Extract INTO a folder named models\ ---
mkdir models
curl.exe -L -o weights.zip "https://zenodo.org/records/17258709/files/weights.zip?download=1"
tar -xf weights.zip -C models
dir models\human                     REM expect run-id folders, each with state_dict.pth

REM --- BT4 goes INTO the same env (it imports src.predict in-process) ---
set BT4_RIBONN_DIR=C:\RiboNN
cd /d C:\path\to\BT4
pip install -e ".[expression-ribonn,dev]"
python -c "from bt4 import api; print(api.available_expression_backends())"
REM expect ('null', 'ribonn')
```

`setx BT4_RIBONN_DIR "C:\RiboNN"` makes the variable stick in *new* prompts. If
`ribonn` is missing from that tuple, `available()` failed — check `C:\RiboNN\src\`
and `C:\RiboNN\models\human\` both exist.

**Windows landmines (all verified against the upstream sources):**

1. **The weights folder must be named `models`.** RiboNN hard-codes the relative path
   `models/<species>/<run_id>/state_dict.pth`. When it is named anything else, BT4
   falls back to `os.symlink`, which on Windows needs Developer Mode or an elevated
   prompt. The `tar -xf ... -C models` above sidesteps that. The zip's root holds
   `human/` and `mouse/`, so it must be extracted *into* a folder called `models`.
2. **`num_workers` must be `0`** — it is now the adapter default, but pass it
   explicitly on the CLI. Windows *spawns* DataLoader workers, which re-import the
   module after BT4 has mutated `sys.path` and `chdir`'d into a temp directory; a
   spawned worker inherits neither and hangs or fails.
3. **`batch_size` defaults to 64**, not RiboNN's 1024, which OOMs a CPU box. Neither
   knob changes a score (RiboNN pads to a fixed width and does not shuffle when
   predicting).
4. **No heredocs.** `python - <<'PY'` does not work in `cmd.exe`; write each inline
   snippet below to a `.py` file first.
5. **Quoting** uses `"..."`, not `'...'`; **line continuation** is `^`, not `\`.

**If any of this fights you, WSL2 is smoother** (`wsl --install`, then the upstream
Linux instructions verbatim). Not required, but RiboNN is only tested on Linux.

## 0.3 — Check `max_shift` (determinism)  ⚠️ blocking

Save as `check_max_shift.py` and run `python check_max_shift.py`:

```python
import glob
import pandas as pd

for f in glob.glob(r"C:\RiboNN\models\*\runs.csv"):
    df = pd.read_csv(f)
    cols = [c for c in df.columns if "max_shift" in c]
    print(f, {c: sorted(df[c].unique()) for c in cols} or "NO max_shift COLUMN")
```

Expect all zeros. **If any value is nonzero, stop** — the adapter needs a seeding fix
before any number it produces is reproducible (invariant #7).

---

# Stage 1 — The free checks (no measured data needed)

**Do these before hunting for a dataset.** They cost an afternoon, need only your
checkout plus sequences already in this repo, and several can end the project with a
decisive, honest answer. The driver is `scripts/ribonn_sensitivity.py` — every report
is stamped `calibrated=False` and can promote nothing. Note `--utr5`/`--utr3` are
**required**: a bundled UTR would be a hidden modelling choice, so you supply real
ones and they are recorded in the report by content hash.

## 1.1 — Positive control: prove the harness works

Score one CDS under two *different real* UTR pairs; the scores **must** differ
substantially — the paper's own attribution puts most of the per-nucleotide signal in
the UTRs. **Do not skip this:** without it, a "no effect" result below could be a
wiring bug misread as biology.

```bat
python scripts\ribonn_sensitivity.py --check utr-control ^
    --utr5 utr5.fa --utr3 utr3.fa --utr5-alt alt5.fa --utr3-alt alt3.fa ^
    --num-workers 0
```

`harness_ok: true` means the UTR context reaches the model. `false` means nothing
else in this stage is interpretable.

## 1.2 — Adapter validation: prove the fold semantics (free)

RiboNN returns one row per input **per outer fold**; averaging all ten is right for a
novel design and wrong for a natural transcript (nine of ten folds trained on its
label). `RiboNNExpressionModel.predict_folds()` exposes the fold identity so you can
keep only the holdout fold. Score RiboNN's own published labels
(`data_with_human_TE_cellline_all_NA_plain.csv`, CenikLab `TE_classic_ML`, GPL-3.0),
join each transcript to the prediction row whose `fold` matches its own, and confirm
the held-out r² lands near the published **0.62** while the other nine sit visibly
higher. If they are indistinguishable, the fold semantics are wrong and everything
downstream is uninterpretable. This proves wiring, not calibration. (This is a small
Python script of your own using `predict_folds`; there is no `--check` for it.)

## 1.3 — The decisive check: does it respond to synonymous change at all?

Use the in-tree panel `scripts/data/ranaghan2021_tab4.fasta` (93 records, CC BY 4.0,
three human proteins × 31 real codon-optimizer outputs — genuine synonymous variants,
no measurements needed). Hold one UTR pair fixed and score all 93 in **one** batched
call; the report gives, per protein, the spread of scores across its variants next to
the spread *between* proteins, plus the GC/CAI/tAI/GC3/length confound correlations.

```bat
python scripts\ribonn_sensitivity.py --check cds-spread ^
    --fasta scripts\data\ranaghan2021_tab4.fasta ^
    --utr5 utr5.fa --utr3 utr3.fa --num-workers 0 --json > stage1_spread.json
```

**How to read `responds_to_synonymous_change`:**

- **False** (within-protein spread at the noise floor) → RiboNN is effectively blind
  to synonymous CDS change under a fixed UTR. It cannot rank BT4's candidates, and no
  panel will change that. **Stop.** Item 11 closes honestly with measured evidence —
  a real, valuable, publishable outcome.
- **True** → continue, but read `median_abs_gc3_spearman` first (Stage 1.4).

## 1.4 — Is it just GC?

`cds-spread` reports the within-protein rank correlation of RiboNN's score against GC,
GC3, CAI, tAI and length in the same run. If GC3 alone explains nearly all of the
within-protein variation, RiboNN is a **GC detector** in this regime — and BT4 already
has GC as an objective term for free, so any later apparent skill must be reported as
such.

## 1.5 — Direction, dose-response, known biology

- **`--check direction`** builds a max-CAI and a min-CAI variant per protein (from the
  in-tree panel's proteins) and runs an exact two-sided sign test on which RiboNN
  prefers. Crude — optimized/deoptimized pairs also differ in GC — so a clear
  preference proves nothing alone; a coin flip is the informative failure.
- **`--check ladder`** walks a real BT4 Pareto frontier for one protein and reports the
  Spearman of score vs CAI along it — a coherence check; a jagged response is unusable
  for ranking even when nonzero.

**Stage 1 gate:** proceed to Stage 2 only if 1.3 shows real within-protein spread and
1.4 shows it is not purely GC.

---

# Stage 2 — Assemble the best panel available

No public dataset fully qualifies, so this is a "best partial panel, reported
honestly" exercise. Work in this order; stop when you have **≥ 4 proteins and ≥ ~90
rows**.

| # | Dataset | System | Readout | n / proteins | Verdict |
|---|---|---|---|---|---|
| 2.1 | **Mauger 2019** Dataset S1 (`pnas.1908052116.sd01.xlsx`) | HeLa + mouse | luminescence / ELISA | ~86 / **4** | **Highest value, ~10 min to check.** Open S1: if it holds per-variant measurements (S2 already has the sequences), this is the best mammalian, fixed-UTR panel anywhere. |
| 2.2 | **PERSIST-seq** (`morrislab/mrl-hl-lbkwk`, HuggingFace) | HEK293T | **mean ribosome load** | 203 / 1 per arm | Readout is the *right* one (never raw protein output). One protein per arm → only a `group = design-family` split, a weaker claim than CLAUDE.md's. CC BY 4.0. |
| 2.3 | **iCodon** (GEO `GSE207584`) | Zebrafish | mRNA decay | ~1,395 / **100** | Only dataset with enough *groups*, but zebrafish + decay → scope-limited secondary result ("validated on regime X, not BT4's Y"). RiboNN has no zebrafish model. CC BY 4.0. |
| — | Mordstein 2020 / CodonBERT MLOS | — | — | — | **Skip unless needed:** SRA raw reads only (days of reprocessing), or UTRs never disclosed (unusable for a full-length model). |

Then write the panel as a TSV in the format below (`api.read_panel` validates it and
**refuses** any row RiboNN would silently drop), add a `.LICENSE.md` sidecar next to
it (the `scripts/data/*.LICENSE.md` convention), add a row to `THIRD_PARTY_DATA.md`,
and record its `content_hash` (`api.read_panel(path).content_hash()`).

**Panel format** — one row per measured sequence:

```
group	variant_id	cds	measured	readout	cell_type	utr5	utr3	species
NLUC	persist_001	ATGGTC...TAA	1.42	mean_ribosome_load	HEK293T	ACATTT...	GCTCGC...	human
NLUC	persist_002	ATGGTG...TGA	1.19	mean_ribosome_load	HEK293T	ACATTT...	GCTCGC...	human
```

- `group` = **the protein** — both the leakage-control unit and the centring unit.
- `measured` = higher means more expression; log-transform ratios first, and prefer
  `log(protein/mRNA)` or ribosome load over raw protein output (protein output
  conflates translation with mRNA stability).
- `readout` / `cell_type` = recorded so the attestation's scope is explicit.
- `utr5` / `utr3` = the experiment's real UTRs, non-empty; `cds` = ACGT, length-3N,
  ending TAA/TGA/TAG.

**Sizing, from the gate's own arithmetic:** ≥ 9 rows in the conformal split (below it
the half-width is `+inf`), ≥ ~100 for a ±0.05 coverage claim, ≥ ~89 for 80% power at
ρ = 0.3, ≥ 3 variants per protein, ≥ 4 proteins. The split assigns whole groups sorted
**alphabetically** — no RNG — so group naming determines the fold.

---

# Stage 3 — Pre-register, then don't touch it

Commit this file **before** running the gate, and do not edit it afterward. A
few-hundred-row panel provably degrades under repeated threshold-tuning: if you run,
see a fail, nudge a threshold and re-run, you have silently turned a validation into a
search.

```
# docs/ribonn_gate_preregistration.md
panel_sha256:         <content_hash of the panel TSV>
species:              human
cell_type:            <the single column, e.g. HEK293T>
readout:              mean_ribosome_load
group_key:            protein
within_group:         true
recalibrate:          true
target_coverage:      0.90
coverage_tolerance:   0.05
min_spearman:         0.30      # a PRE-COMMITMENT, not a community standard
                                # (verified: no such standard exists)
decision_rule:        RiboNN's cluster-bootstrap 95% CI lower bound on within-protein
                      Spearman must exceed EVERY baseline's point estimate
baselines:            permutation, CAI-only, GC3-only, length-only, constant
width_rule:           median interval width / label IQR must be < 1.0
runs_permitted:       1
```

---

# Stage 4 — Run the gate

The gate scores the panel in one batched call per UTR context, runs the same gate on
every baseline, and reports a `promotable` verdict that requires the thresholds **and**
beating every baseline **and** an informative interval — each reported separately.

```bat
python scripts\run_expression_gate.py ^
    --panel panel.tsv --backend ribonn --species human --cell-type HEK293T ^
    --within-group --recalibrate ^
    --min-spearman 0.30 --target-coverage 0.90 --coverage-tolerance 0.05 ^
    --baselines permutation,cai,gc3,length,constant ^
    --num-workers 0 --json > gate_result.json
```

Or, without leaving the CLI: `bt4 expression-gate panel.tsv --within-group
--recalibrate --cell-type HEK293T --num-workers 0`. **Omit `--within-group` and it
warns loudly** that pooled mode credits between-protein skill, which is not BT4's
regime. Keep `gate_result.json` and `stage1_spread.json` with the pre-registration
file — together they are the record.

---

# Stage 5 — The decision

| Result | Meaning | Action |
|---|---|---|
| `promotable: true` — within-protein bootstrap CI lower bound beats **every** baseline, coverage within tolerance, width/IQR < 1 | RiboNN genuinely ranks synonymous variants in this scope | Commit an `ExpressionAttestation` (below); update CLAUDE.md §6/§9, `NEXT_SESSION.md` item 11, CHANGELOG |
| Ranking passes but does not beat **CAI-only** | No value over what BT4 already computes for free, in-loop | Stays `False`; documented with the numbers |
| Ranking fails | RiboNN does not do BT4's job | Stays `False` **for this regime**; evidence into `RESEARCH_ribonn_calibration.md`; item 11 closes as **answered**, not pending. Fine-tuning is a separate, much larger project |
| Ranking passes, coverage fails | Orders variants usefully; error bars lie | Report both. The gate requires both halves, so the default is "stays `False`". A rank-only promotion would be a deliberate contract change, not a slip-in |
| Only a single-protein panel is ever available | `group = design-family` is the only runnable split | A narrower claim than CLAUDE.md's; report it as such, do not promote |

## Recording a pass (only if `promotable`)

`ExpressionAttestation` is the single seam that flips `calibrated=True`. It is
licence-clean by construction (derived scalars + public content hashes only — never a
raw per-sequence score), carries four floors so a self-serving threshold cannot
self-certify, and carries its **scope** (species, cell types, readout, panel hash).

```python
# promote.py  — run once, on the machine that ran the gate.
import json
from bt4 import api, __version__
from bt4.biomodels.expression import attest_expression

panel = api.read_panel("panel.tsv")
comparison = api.expression_gate(
    panel, "ribonn",
    settings=api.GateSettings(within_group=True, recalibrate=True,
                              min_spearman=0.30, coverage_tolerance=0.05),
    species="human", cell_types=("HEK293T",),
)
assert comparison.promotable, "gate did not pass; nothing to attest"
att = attest_expression(
    comparison, species="human", cell_types=("HEK293T",),
    readout="mean_ribosome_load", bt4_version=__version__,
)
open("ribonn_attestation.json", "w").write(json.dumps(att.to_dict(), indent=2))
```

**Simpler, and preferred:** pass `--attest ribonn_attestation.json --readout
mean_ribosome_load` to `scripts/run_expression_gate.py` itself. That records the
attestation from the *same* comparison the verdict was read off — no second scoring pass,
and no chance of configuring the second invocation differently from the first.

A user then opts in for their own session:

```bash
export BT4_EXPRESSION_ATTESTATION=/path/to/ribonn_attestation.json
export BT4_EXPRESSION_USE_ATTESTED=1
```

or explicitly, per call:

```python
from bt4 import api
model = api.resolve_expression_backend(
    "ribonn", species="human", utr5="...", utr3="...",
    cell_types=("HEK293T",), use_attested=True,
)
assert model.calibrated  # refuses on species / cell-type / top_k / UTR / weight mismatch
```

## The wiring (landed 2026-08) and the scope it binds

**Superseded caveat, kept for the record:** this section used to read *"`verified_predictor`
is the seam, but nothing in `src/` calls it yet"* — so even a committed, passing
attestation changed nothing for users. That is no longer true.
`bt4.biomodels.expression.attestations.promote_if_attested` is called from
`resolve_backend`, gated on `BT4_EXPRESSION_USE_ATTESTED` (or an explicit
`use_attested=`), and BT4 Studio carries a per-run toggle on the Candidates tab that
displays the attestation's scope and pins the head to it. A promoted head reorders the
candidate set and re-picks the delivered sequence, and its attestation's content hash
enters the run manifest.

`$BT4_EXPRESSION_ATTESTATION` reads a maintainer's own record, because an expression
attestation is earned against a measured panel that is often unpublished — and the record
carries the panel hash plus a hash per UTR context, from which a *short* UTR is
brute-forceable. **Nothing is bundled**: no expression head has passed its gate, and
shipping a record for one that has not would be CLAUDE.md §10.6's fabricated placeholder.

The scope is strict by design, and is now taken from the run rather than declared:
`attest_expression` derives species / cell types / readout from `GateComparison.scope`
and refuses a declaration that disagrees, while `verified_predictor` binds `top_k` and
the UTR context in addition to species, cell types and the weight hashes. An attestation
earned on HEK293T does **not** certify a head averaging all 78 cell types, and the honest
claim is always scoped — "calibrated for ranking synonymous variants of a known protein,
cell type X, readout Y, UTR context Z". The cross-protein / novel-CDS regime cannot be
honestly gated at all: it needs ~100 held-out proteins that do not exist.

---

# Verification

**Already green in CI (no weights, and must stay so):** the gate's honesty tests
(`test_expression_gate.py`) — a gene-identity-only head passes pooled and *fails*
within-group; a rank-backwards head is not rescued by a negative fitted slope; a
constant predictor is caught on both the rank and the width axis; `NullExpressionModel`
cannot pass in either mode. The panel reader (`test_expression_panel.py`) refuses an
over-length row rather than dropping it. The attestation layer
(`test_expression_attestation.py`) refuses a failing gate, a pooled run, a below-floor
threshold, and a scope or weight-hash mismatch, and the dataclass shape is asserted
licence-clean at import. The sensitivity script (`test_ribonn_sensitivity.py`) reports
a blind backend as blind and a GC3-only backend as a GC detector.

**On the machine with the weights (human-run):**

1. `python -c "from bt4 import api; print(api.available_expression_backends())"` →
   `('null', 'ribonn')`.
2. Stage 0.3 prints all-zero `max_shift`.
3. Stage 1.1 `harness_ok: true`; Stage 1.2 held-out r² ≈ 0.62.
4. Determinism (invariant #7): the same panel scored twice is byte-identical.
5. If `promotable`, `--attest` writes the attestation and the opt-in flips the flag; a
   deliberately wrong `cell_types`, `top_k` or UTR context on the model makes it
   *refuse* rather than hand back an uncalibrated head.

# What not to do

- **Do not run the gate more than once** against the same panel with tweaked
  thresholds. Pre-register, run once.
- **Do not promote on a pooled result**, however high the Spearman — pooled credits
  between-protein skill, which is not BT4's regime, and `attest_expression` refuses it.
- **Do not compare RiboNN's CLR-residual TE to raw protein output.** Use ribosome load
  or `log(protein/mRNA)`; raw protein output re-introduces the mRNA-abundance term TE
  divides out.
- **Do not relabel a hand-weighted CAI+GC composite as "calibrated"** to force a pass
  (the §10.5 magic-scalar trap).
- **Do not flip `calibrated=True` by hand.** `dataclasses.replace(..., fidelity_verified=True)`
  bypasses every floor and scope check; `verified_predictor` is the only sanctioned path.
- **Do not commit an attestation for a panel you cannot publish** — the record binds the
  panel's hash and each UTR context's hash. Use `$BT4_EXPRESSION_ATTESTATION` instead.
- **Do not bundle an attestation that has not been earned** to demonstrate the calibrated
  path. Test doubles belong in tests; a committed artifact is a claim.
