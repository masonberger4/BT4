# Calibrating the splice CNNs (SpliceAI + Pangolin)

> **Status: NOT STARTED.** This is the procedure for `NEXT_SESSION.md` queue item
> **#10**. The evidence behind every upstream claim here lives in
> [`RESEARCH_splice_cnn_calibration.md`](RESEARCH_splice_cnn_calibration.md);
> this file is the runbook. Live status stays in
> [`NEXT_SESSION.md`](NEXT_SESSION.md), not here.

## Context

BT4 wraps two published splice models — **SpliceAI** (Illumina) and **Pangolin**
(Zeng & Li) — behind its `SplicePredictor` contract. Both adapters are written,
tested, and drive *your own* installed copy of each model, hash-verifying the
weights before loading. But both ship **`calibrated=False`**, and per CLAUDE.md
§10.6 that flag is *earned on data, never assigned*. Until it flips:

- `bt4.biomodels.splice.default()` keeps returning the honest PWM baseline;
- every splice audit is banner-led **"UNCALIBRATED (advisory)"**;
- design-flow **step 6** (targeted synonymous splice auto-edit) stays blocked;
- the Studio Candidates tab shows flags that advise but never edit.

## The one thing to understand before starting

"Calibration" means **two different things**, and BT4's code only implements the
first. This plan does both, in order.

| | **Part A — integration fidelity** | **Part B — statistical calibration** |
|---|---|---|
| Question | *Does BT4's wrapper reproduce the published model's own numbers?* | *Do those numbers mean what they claim on real splice sites?* |
| Analogy | Checking your scale reads the same as the factory's scale | Checking the scale is actually accurate in grams |
| Machinery in BT4 | **Exists** — `verify_*_fidelity`, `FidelityAttestation`, `verified_predictor` | **Does not exist** — no PR-AUC / MCC / ECE / Brier anywhere |
| Flips `calibrated=True`? | Yes — this is literally the gate | No, not by BT4's current definition |
| What you need | The licensed weights, on your machine | An MIT-licensed public panel (see B1) |
| Rough effort | 1–2 days | 2–4 days |

Part A unblocks the codebase. Part B is the honest science that says whether the
unblocked feature deserves trust. **Doing A without B is how you ship a confident
wrong answer** — and finding 5 below is exactly why that matters here.

## Five findings that shape the job

Verified against primary sources and against BT4's own code. Details and citations
in [`RESEARCH_splice_cnn_calibration.md`](RESEARCH_splice_cnn_calibration.md).

1. **All 17 pinned weight hashes are already correct.** The upstream bytes were
   downloaded and re-hashed: 5/5 SpliceAI `.h5` and 12/12 Pangolin `.3.v2` match
   `PINNED_WEIGHT_SHA256` exactly. Step A3 should pass on the first try — which
   de-risks the whole of Part A.

2. **The attestation layer is wired to nothing.** `FidelityAttestation`,
   `attest_backend`, `load_attestation` and `verified_predictor` all exist and are
   well tested — but outside `tests/test_splice_attestation.py` **nothing in BT4
   ever calls them.** The three places that build a CNN backend —
   `pipeline/splice_audit.py:40`, `pipeline/splice_crosscheck.py:149`, and
   `scripts/compare_splice_backends.py:72` — all construct the predictors bare. So
   even after the gate runs and an attestation is committed, **nothing changes for
   anyone** until the wiring lands (step A7). That is the largest code task here.

3. **Pangolin ships two model sets, and the wrong one is in the example script.**
   `pangolin/models/` holds **64** files: 40 base `final.{1..5}.{0..7}.3` and 24
   fine-tuned `final.{1..3}.{0..7}.3.v2`.

   | | weights | folds | channels |
   |---|---|---|---|
   | `pangolin/pangolin.py` (the CLI — what BT4 mirrors) | `final.{1,2,3}.{0,2,4,6}.**3.v2**` | 3 | `[1,4,7,10]` |
   | `scripts/custom_usage.py` (the example) | `final.{1..5}.{i}.**3**` | 5 | `INDEX_MAP` |

   Same repo, different numbers. BT4 correctly tracks the CLI (the `.v2` models are
   the paper's fine-tuned ones, used for its variant-effect figures). **Capturing
   reference expectations by copying `custom_usage.py` makes the gate fail and look
   like a BT4 bug when it is not.** The capture script in A4 mirrors the CLI.

4. **Illumina's SpliceAI repo was archived 2026-04-20 and is read-only.** No
   upstream fix will ever land — including for the Keras-3 breakage. Also:
   **`pip install pangolin` installs an unrelated package** (a probabilistic
   programming language); Pangolin is GitHub-only. Both belong in BT4's docs.

5. **These models are weakest exactly where BT4 operates.** Smith & Kitzman,
   *Genome Biol* 24:294 (2023): across 8 predictors on 3,616 variants, median
   **prAUC 0.773 intronic vs 0.419 exonic** — improvements needed *"especially
   within exons"*, with lower concordance for *"missense or **synonymous**"*
   variants. BT4 designs coding sequence: its entire regime is the exonic half.
   This does not block calibration, but it **must be stated wherever BT4 reports
   splice risk**, and it is the strongest argument for not stopping at Part A.

   *Related, and also structural:* BT4 pads each sequence with 5,000 literal `N`
   per side (`pangolin.py:484`, `spliceai.py:387`) — upstream's own documented
   convention. But `pipeline/splice_audit.py::_FlankedPredictor` forwards
   `calibrated` **unchanged**, so a promoted backend would silently claim
   calibration on the real-flank path too, which no BT4 doc allows. Fix in A7.4.

   *The clean way out, which also simplifies capture:* `validate_dna` rejects `N`,
   so every panel sequence must be pure ACGT. Make each case a **≥10,001 nt real
   genomic window with the site in the middle** — exactly what Pangolin's docs
   prescribe ("5000 bases before the site, base at the site, 5000 bases after").
   The N-padding then only touches the window's outer edges, never the site.

---

# Part A — Integration fidelity (flips `calibrated=True`)

### A0 — Scratch area, outside the repo

```bash
mkdir -p ~/bt4-splice-calib && cd ~/bt4-splice-calib
```

Everything licensed lives here, **never inside the BT4 checkout**. `.gitignore`
already blocks `*.h5`, `*.pth`, `external/`, but distance is safer than a rule.
**Use two separate environments** — Pangolin needs PyTorch, SpliceAI needs an old
TensorFlow, and they do not coexist. The two tracks share nothing but this plan.

### A1 — Install the two upstream models

**Pangolin** — GPL-3.0 (the LICENSE file is verbatim GNU GPLv3), v1.0.2, actively
maintained.

```bash
conda create -n pangolin python=3.10 -y && conda activate pangolin
pip install torch                       # upstream pins NO version; CPU wheel is fine
pip install gffutils biopython pandas pyfastx
git clone https://github.com/tkzeng/Pangolin.git ~/bt4-splice-calib/Pangolin
pip install ~/bt4-splice-calib/Pangolin
pip install -e '.[splice-pangolin]'     # from the BT4 checkout
```

> **Do not `pip install pangolin`** — that is an unrelated probabilistic
> programming package. GitHub install only. `conda install -c conda-forge pyvcf` is
> needed only for the variant CLI, not for anything in this plan. Weights ship
> inside the repo at `Pangolin/pangolin/models/` (64 files, 2,877,321 bytes each).

**SpliceAI** — code PolyForm Strict 1.0.0, weights CC BY-NC 4.0 (from the LICENSE
file; note `setup.py` still says "GPLv3", which is stale and wrong).

```bash
conda create -n spliceai python=3.10 -y && conda activate spliceai
pip install "tensorflow==2.15.*"        # last TF whose bundled Keras is 2.x
pip install "numpy<2" "pandas<2.2" "setuptools<81"
pip install spliceai==1.3.1
pip install -e '.[splice-spliceai]'     # from the BT4 checkout
```

Weights ship inside the pip package at `spliceai/models/spliceai{1..5}.h5`
(3,131,720 bytes each). Every pin is version-critical because `spliceai` 1.3.1
declares **no upper bounds and no `requires_python`**:

| Pin | Why |
|---|---|
| `tensorflow==2.15.*` | TF ≥ 2.16 defaults to **Keras 3**, which cannot load these 2019 Keras-2 `.h5` graphs. Alternative: modern TF + `pip install tf_keras` + `export TF_USE_LEGACY_KERAS=1` — BT4's `_import_keras` already falls back to `tf_keras`, but TF 2.15 is the verified-safe path. |
| `numpy<2` | `spliceai/utils.py` uses the long-deprecated `np.fromstring` |
| `setuptools<81` | the package imports `pkg_resources`, removed in setuptools ≥ 81 |

`WARNING:absl:No training configuration found in the save file` on load is
**benign** for inference.

> **Note for BT4's own metadata:** `bt4[splice-spliceai]` currently declares
> `tensorflow>=2.6` with **no upper bound**, so a fresh install pulls a TF that
> cannot load the weights. Tighten that pin as part of this work.

> **Because Illumina's repo is archived (finding 4)**, the maintained forks
> `bw2/SpliceAI` (v1.3.4) and `bw2/Pangolin` (v1.0.5) — which back
> `spliceailookup.broadinstitute.org` — are the pragmatic long-term install path,
> and the Broad's Docker images (`docker.io/weisburd/spliceai-38`, `pangolin-38`)
> skip the dependency pain entirely. **For the fidelity gate, use the upstream
> weights BT4 actually pins.**

### A2 — Point BT4 at the weights

Resolution order: explicit `model_dir=` → environment variable → the installed
package's own `models/` directory.

```bash
export BT4_PANGOLIN_MODEL_DIR=~/bt4-splice-calib/Pangolin/pangolin/models
export BT4_SPLICEAI_MODEL_DIR=$(python -c "import spliceai,os;print(os.path.join(os.path.dirname(spliceai.__file__),'models'))")
```

Sanity check (each should print `True`):

```bash
python -c "from bt4.biomodels.splice import PangolinSplicePredictor as P; print(P().available())"
python -c "from bt4.biomodels.splice import SpliceAiSplicePredictor as S; print(S().available())"
```

`available()` is a cheap existence check — it does **not** hash. That is next.

### A3 — Verify the weights match BT4's pins

Per finding 1 this should pass cleanly, but run it — it is the cheapest possible
check and everything downstream depends on it.

```bash
python - <<'PY'
import hashlib, os
from pathlib import Path
from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256 as PANGOLIN_PINS
from bt4.biomodels.splice.spliceai import PINNED_WEIGHT_SHA256 as SPLICEAI_PINS

def check(label, pins, env):
    d = Path(os.path.expanduser(os.environ[env]))
    print(f"\n=== {label} ({d}) ===")
    ok = True
    for name, want in sorted(pins.items()):
        p = d / name
        if not p.is_file():
            print(f"  MISSING   {name}"); ok = False; continue
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                h.update(blk)
        got = h.hexdigest()
        if got == want:
            print(f"  ok        {name}")
        else:
            ok = False
            print(f"  MISMATCH  {name}\n      pinned {want}\n      actual {got}")
    print(f"  -> {label}: {'ALL PINS MATCH' if ok else 'DOES NOT MATCH - STOP'}")

check("Pangolin", PANGOLIN_PINS, "BT4_PANGOLIN_MODEL_DIR")
check("SpliceAI", SPLICEAI_PINS, "BT4_SPLICEAI_MODEL_DIR")
PY
```

If anything mismatches, do **not** edit the pins to silence it. Either the download
came from a mirror/fork (re-download from the official source), or upstream
re-released — in which case update `PINNED_WEIGHT_SHA256` as a **separate,
clearly-described commit** naming the release, then redo A3. A subset of files can
never satisfy `verified_predictor`, which compares the full sorted tuple with `!=`.

### A4 — Capture the reference panel

The scientific core of Part A, and the one place it is easy to cheat by accident.

**The rule: expected numbers must come from upstream's own code, not from BT4's
adapter.** Capturing with BT4 and then checking BT4 makes the gate pass trivially.
**Neither script below imports `bt4`.**

**First, a free sanity check (SpliceAI only).** Upstream ships `examples/input.vcf`
(10 records) and `examples/output.vcf` with exact delta scores — e.g.
`19:38958362 C>T` → `T|RYR1|0.00|0.00|0.91|0.08|-28|-46|-2|-31`. Reproduce it
before capturing anything:

```bash
spliceai -I examples/input.vcf -O out.vcf -R hg19.fa -A grch37
diff <(grep -v '^##' out.vcf) <(grep -v '^##' examples/output.vcf)
```

(hg19 FASTA: `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz`
— note `.soe.`, not `.cse.`, which has a certificate mismatch.) This exercises the
*variant* path, not `score_sequence`, so it does not replace the capture below —
but a clean diff proves the install is sound. **Pangolin ships no equivalent
oracle** (`examples/` holds inputs only), so its panel must be captured.

**Panel composition** — 10–20 cases. This gate is about wrapper faithfulness, not
statistical power:

| Kind | Count | Why |
|---|---|---|
| Windows centred on well-known real donor/acceptor sites | 3–5 | the signal the model exists to detect |
| Windows from deep inside a large intron, no annotated site | 3–5 | true negatives |
| Real CDSs of the kind BT4 designs | 2–3 | BT4's actual regime |
| BT4-designed synonymous variants of those CDSs | 2–3 | the *exact* regime — designed, not natural |
| Edge cases: ~100 nt, and GC-extreme | 2 | guards the padding/crop arithmetic |

Every sequence **pure ACGT, no `N`** (`validate_dna` rejects `N`).

**Pangolin capture** — mirrors `pangolin/pangolin.py::compute_score` (the CLI),
*not* `custom_usage.py` (finding 3):

```python
# capture_pangolin.py  — run in the pangolin env. Does NOT import bt4.
# Weights are located from $BT4_PANGOLIN_MODEL_DIR, NOT pkg_resources: upstream's
# own scripts use `resource_filename`, which setuptools >= 81 removed along with
# pkg_resources. torch pulls in a modern setuptools, so that import fails on a
# fresh install.
import json, os
import numpy as np, torch
from pangolin.model import Pangolin, L, W, AR

MODEL_DIR = os.environ["BT4_PANGOLIN_MODEL_DIR"]

IN_MAP = np.asarray([[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
CHANNEL = {0: 1, 2: 4, 4: 7, 6: 10}        # the CLI's [1,4,7,10] P(splice) heads

def one_hot(seq):                          # '+' strand only, matching BT4
    s = (seq.upper().replace('A','1').replace('C','2')
                    .replace('G','3').replace('T','4').replace('N','0'))
    return IN_MAP[np.asarray(list(map(int, s))).astype('int8')]

models = {}
for i in (0, 2, 4, 6):                     # heart, liver, brain, testis
    fold_models = []
    for j in (1, 2, 3):                    # the CLI's three production folds
        m = Pangolin(L, W, AR)
        w = torch.load(os.path.join(MODEL_DIR, f"final.{j}.{i}.3.v2"),
                       map_location=torch.device("cpu"))
        m.load_state_dict(w); m.eval()
        fold_models.append(m)
    models[i] = fold_models

def site_scores(seq):
    padded = "N"*5000 + seq.upper() + "N"*5000
    x = torch.from_numpy(np.expand_dims(one_hot(padded).T, axis=0)).float()
    per_tissue = []
    with torch.no_grad():
        for i in (0, 2, 4, 6):
            folds = [m(x)[0, CHANNEL[i], :].numpy() for m in models[i]]
            per_tissue.append(np.mean(folds, axis=0))
    return np.mean(per_tissue, axis=0).tolist()      # tissue-agnostic average

json.dump([{"sequence": s, "expected": site_scores(s)} for s in SEQUENCES],
          open("pangolin_panel.json", "w"))
```

Note the `.T` — Pangolin is **channel-major `[4][L]`**. Its P(splice) head is a
*binary* softmax, so donor and acceptor are **not** separated; that is why BT4 puts
the combined track in `SpliceResult.donor` and leaves `acceptor` all-zero.

**SpliceAI capture — `scripts/capture_spliceai_panel.py` (landed).** Same three-step
workflow as Pangolin, in the *other* virtualenv (TensorFlow 2.15 and PyTorch do not
coexist):

```
python scripts\make_splice_panel.py --out panel_sequences.json
python scripts\capture_spliceai_panel.py --panel panel_sequences.json --out expected_spliceai.json
```

`make_splice_panel.py` is shared — the *same* panel feeds both backends, so the two
gates are run on identical sequences.

Two things differ from the Pangolin capture, both deliberate:

- **It imports upstream's own `one_hot_encode`** rather than re-deriving it. Pangolin's
  CLI encodes inline, which forced the Pangolin capture to reimplement it; SpliceAI
  ships the encoder as a reusable function, so importing it is strictly stronger
  evidence — a transposed layout, a wrong base order, or a mishandled `N` in BT4's own
  `_one_hot_rows` shows up as a gate **failure** instead of being reproduced identically
  on both sides. There is deliberately **no fallback encoder**: if `spliceai.utils` will
  not import, the script refuses and names the cause (NumPy 2 removed `np.fromstring`,
  which `spliceai/utils.py` still calls — pin `numpy<2`). A "helpful" local fallback
  would silently destroy the independence, so a test asserts the script defines no
  encoder of its own.
- **It records two tracks, not one.** SpliceAI's 3-way softmax (null / acceptor / donor)
  genuinely separates them, so `expected_acceptor` and `expected_donor` must *both*
  match. Pangolin's binary head emits one combined `P(splice)`, which is why its capture
  carries a single track and BT4 puts it in `SpliceResult.donor` with `acceptor`
  all-zero.

SpliceAI is **position-major `[L][4]`** (Keras channels-last) — the opposite of
Pangolin's channel-major `[4][L]`. Save both captures **in the scratch directory only**:
those arrays *are* the licensed model outputs (GPL-derived for Pangolin, CC BY-NC for
SpliceAI) and must never enter the repo.

### A5 — Run the fidelity gate

```
python scripts\run_splice_fidelity_gate.py --panel panel_sequences.json --captured expected_pangolin.json
python scripts\run_splice_fidelity_gate.py --panel panel_sequences.json --captured expected_spliceai.json
```

The runner reads **which backend from the capture payload itself**, so a Pangolin
capture can never be checked against the SpliceAI adapter — the numbers would still be
numbers, they would just be describing a different model. Passing `--panel` binds the
capture to the panel it came from by content hash, so a stale capture is caught rather
than silently compared against regenerated sequences.

It also reports the panel's **peak-score spread** and warns when that spread is too
narrow: a gate that passes on a panel where the model scores ~0 everywhere is nearly
vacuous, because a wrong channel, a wrong fold set, or a transposed one-hot would match
within tolerance too. A pass on a flat panel is reported as passing *and* as weak.

`passed=True` requires `max_abs_deviation <= 1e-3`, which is also
`MAX_ATTESTATION_TOLERANCE`. The tolerance cannot be relaxed to force a pass —
`attest_backend` refuses and says why.

| `max_abs_deviation` | Likely cause |
|---|---|
| ~1e-7 | float32 summation-order noise — arguably fine, but understand it before accepting |
| ~1e-2 – 1e-1 | wrong ensemble members or channels (Pangolin heart 0→1, liver 2→4, brain 4→7, testis 6→10) |
| Large, scores look shuffled | one-hot layout (Pangolin `[4][L]` vs SpliceAI `[L][4]`), or SpliceAI channels swapped |
| Wildly off on every case | captured with `custom_usage.py`'s base `.3` weights (finding 3) |
| Matches mid-sequence, diverges at the ends | padding/crop misalignment |

A genuine failure is a **bug in BT4's adapter**, to be fixed in the adapter.

### A6 — Build and commit the attestation

```bash
python - <<'PY'
import json, bt4
from bt4.biomodels.splice import attest_backend
from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256 as PINS
att = attest_backend("pangolin", report, PINS, bt4_version=bt4.__version__)   # report from A5
print(json.dumps(att.to_dict(), indent=2, sort_keys=True))
print("content_hash:", att.content_hash())
PY
```

The file may contain exactly eight fields — `backend`, `passed`,
`max_abs_deviation`, `n_cases`, `tolerance`, `weight_sha256`, `bt4_version`,
`schema_version` — enforced by `_ALLOWED_FIELDS`, an import-time assertion, and a
test. The captured per-position arrays from A4 **are the licensed outputs** and
must never appear.

Commit to a new package data directory:

```
src/bt4/biomodels/splice/data/pangolin.attestation.json
src/bt4/biomodels/splice/data/spliceai.attestation.json
```

and add to `pyproject.toml` under `[tool.hatch.build.targets.wheel].artifacts`:

```toml
  "src/bt4/biomodels/splice/data/*.json",
```

(the equivalent entry already exists for `expression/data/*.json`; without it the
files reach the wheel but not the sdist).

### A7 — Wire the attestation in *(the part that makes it real)*

Design: **committed attestation + explicit opt-in.** The attestation proves
fidelity; the *user* decides to run a licensed model. Keeping promotion opt-in also
means **no existing test breaks** — `default()` still returns the PWM baseline and
the `all_calibrated is False` assertions stay green.

1. **New module `src/bt4/biomodels/splice/attestations.py`**, mirroring
   `expression/ribonn.py`'s `load_pinned_sha256`:
   ```python
   @lru_cache(maxsize=None)
   def bundled_attestation(backend: str) -> FidelityAttestation | None: ...
   def promote_if_attested(predictor, *, enabled: bool): ...
   ```
   `promote_if_attested` returns the predictor unchanged when `enabled` is False or
   no attestation is bundled; otherwise it calls `verified_predictor`, which raises
   `AttestationError` on mismatch rather than silently downgrading.

2. **One opt-in switch**, three spellings: env var `BT4_SPLICE_USE_ATTESTED=1`;
   CLI flag `--use-attested-splice` on `optimize` / `validate`; a Studio checkbox
   beside the existing "run installed splice CNNs" toggle.

3. **Apply it at the three construction sites** from finding 2:
   `pipeline/splice_audit.py::available_splice_backends`,
   `pipeline/splice_crosscheck.py::resolve_splice_backend`,
   `scripts/compare_splice_backends.py::available_backends`.

4. **[DONE] The flank-regime leak is fixed at both seams.** An audit of the
   flag's propagation found this was two places, not one, and that the root had a
   second caller the plan had missed:
   - `score_in_context` (`biomodels/splice/base.py`) now returns
     `calibrated=False` whenever flanks are applied. This is the root — and it also
     feeds `pipeline/tracks.py`, which set `TracksResult.splice_calibrated` and the
     track's unit label from it, so a flanked track would have been labelled
     calibrated site probabilities.
   - `_FlankedPredictor.calibrated` (`pipeline/splice_audit.py`) needs its own
     guard regardless, because `audit_splice` reads `predictor.calibrated`
     directly — for each `SpliceFlag`, each `BackendCandidateAudit`, and the
     report-level `all_calibrated` that drives Studio's banner — never the
     `SpliceResult`'s flag.

5. **[DONE] The unbound-configuration gap is closed, at the promotion seam.**
   `verified_predictor` now refuses a Pangolin predictor whose `tissues` are not
   the full default set. The reasoning is sharper than "bind the config": an
   honored attestation is *required* to claim the adapter's full 12-file pinned map
   (a subset fails the equality check), yet a gate run at one tissue loads only 3
   of those files and reads 1 of 4 output channels. Enforcing at
   `verified_predictor` rather than in `promote_if_attested` means it holds for
   every caller, including a user promoting by hand.

   *Not an issue:* `top_k` only affects pooling of already-computed per-position
   scores, and SpliceAI always loads all five ensemble members, so neither needs
   binding.

6. **Stamp it into provenance.** Pass the attestation's `content_hash()` through
   `build_manifest(extra=...)` wherever a calibrated splice backend influenced a
   result, mirroring `candidates.py:379`'s `predictor` / `predictor_calibrated`.
   Invariant #9: a calibrated audit must be reproducible from its stamp.

### A8 — Tests, docs, changelog

**Tests** (all CI-runnable, none needing weights):

- the committed attestation JSON loads and round-trips;
- its `weight_sha256` equals the adapter's `PINNED_WEIGHT_SHA256` exactly;
- its `content_hash()` equals a pinned constant (catches silent edits);
- `promote_if_attested(..., enabled=False)` is a no-op;
- `_FlankedPredictor.calibrated` is `False` when flanks are present;
- `default()` **still** returns `ConsensusPwmSplicePredictor` — the test that proves
  promotion stayed opt-in.

`tests/test_splice_attestation.py` hardcodes `_VERSION = "0.4.0"`; make it read
`bt4.__version__` so it does not rot at the next bump.

**Docs** — CLAUDE.md §6 and [`NEXT_SESSION.md`](NEXT_SESSION.md) in the *same*
change (§10.11): item #10 comes off `[BLOCKED-human]`, the CNN status-board row
changes, item #12 unblocks. Add a "Splice CNN environment gotchas" section to
`NEXT_SESSION.md` mirroring the existing RiboNN one, covering the `.3`-vs-`.3.v2`
trap, the Keras-3/TF-2.15 pin, `pip install pangolin` being the wrong package, and
the two-env split. Document `BT4_PANGOLIN_MODEL_DIR` / `BT4_SPLICEAI_MODEL_DIR` —
which today appear in **no** markdown file. Also record that **Illumina/SpliceAI
was archived 2026-04-20** and name `bw2/SpliceAI` + the Broad Docker images as the
live install path. Update `README.md`'s wrapped-models table and add a
`CHANGELOG.md` entry.

---

# Part B — Statistical calibration (the honest science)

Part A proves the wrapper is faithful. It says **nothing** about whether a score of
0.5 means "50% chance this is a splice site" — and BT4 treats 0.5 as a localization
threshold (`audit.py:74`) explicitly documented as *"a display / localization knob,
not a calibrated cutoff."* The moment `calibrated=True` reaches users, that knob
starts reading like a decision boundary.

**There is published evidence it is the wrong number.** Walker et al., *AJHG* 2023
(ClinGen SVI Splicing Subgroup), calibrated SpliceAI as a likelihood ratio on 2,736
non-canonical-splice-site variants with in-vitro splicing assays across 8 genes:

| SpliceAI Δ | Likelihood ratio | Interpretation |
|---|---|---|
| ≤ 0.1 | 0.17 [0.14–0.21] | moderate evidence **against** spliceogenicity |
| 0.1 – 0.2 | 1.00 [0.71–1.39] | **uninformative** |
| ≥ 0.2 | 15.99 [13.23–19.32] | moderate evidence **for** |

They conclude **"SpliceAI score ≥ 0.5 may be calibrated too high"** and recommend
0.2 / 0.1 operating points. Separately, OpenSpliceAI (Chao et al., *eLife* 2025) is
the only published reliability/ECE analysis: SpliceAI-architecture models are
**slightly overconfident**, and class-wise **temperature scaling** on a held-out
split fixes it. Its `calibrate` subcommand is a ready-made recipe to mirror. *(No
published work applying isotonic regression to SpliceAI deltas was found.)*

### B1 — Get a panel *(start with the MIT one — it saves days)*

**Start here: `github.com/kitzmanlab/splicebench2023`** (Smith & Kitzman, *Genome
Biol* 2023). 3,616 variants across 5 genes (POU1F1 ex2, MST1R/RON ex11, FAS ex6,
WT1 ex9, BRCA1 11 exons), **already scored by 8 tools including SpliceAI and
Pangolin**, plus 500,000 random background exonic SNVs. **Explicitly MIT-licensed.**
This gives labels *and* a reference scoring to cross-check the install against,
before building anything. Its BRCA1 saturation-genome-editing subset is the closest
public thing to BT4's synonymous-CDS regime.

Then add, as budget allows: **MFASS** (27,733 ExAC variants, 83% outside canonical
sites; no LICENSE file — check before redistributing), **MaPSy** (4,964 exonic
disease variants), **SpliceVarDB** (50,715 assessed variants), or build labels
directly from a pinned **GENCODE** release. Full table in
[`RESEARCH_splice_cnn_calibration.md`](RESEARCH_splice_cnn_calibration.md) §5.

**Match the published test split.** Both models trained on human chr **2, 4, 6, 8,
10–22, X, Y** — so **chr 1, 3, 5, 7, 9 are held out**. Building from training
chromosomes produces flattering nonsense. *(Caveat worth recording: OpenSpliceAI
found 0.71% of MANE transcripts on the test chromosomes are paralogous to training
sequences — the split is not perfectly clean.)*

If building labels from GENCODE directly: for a `+`-strand transcript each internal
exon's end+1 is a **donor** and each exon's start−1 an **acceptor** (flip for `−`);
drop the first acceptor and last donor. Negatives are the other positions in the
same gene bodies — exactly the top-k denominator both papers use. SpliceAI's own
bundled `spliceai/annotations/grch38.txt` is already a flat `EXON_START`/`EXON_END`
table if you want to skip GTF parsing.

#### The file format to build (and the off-by-one it refuses)

**Landed:** `bt4.api.read_splice_panel` reads a small tab-separated format, so B1 has
a defined target rather than an ad-hoc script. Three required columns —
`window_id`, `group`, `sequence` — and four optional: `donors`, `acceptors`
(comma-separated 0-based positions), `strand`, `note`.

```
window_id	group	sequence	donors	acceptors
BRCA1_ex11	chr1	ACGT…	1204,3391	880,2755
deep_intron_1	chr3	ACGT…
```

* **One row is a window**, not a site: one window yields `len(sequence)` scored
  positions, of which only the annotated ones are positive. A window with no sites at
  all is a legitimate pure-negative control.
* **`group` is the chromosome** — the leakage-control unit. Overlap with the training
  set (chr 2, 4, 6, 8, 10–22, X, Y) is reported, and a panel that overlaps can never be
  `promotable`.
* **Minus-strand windows are reverse-complemented *before* they reach the file.** The
  sequence is stored in the orientation the sites are annotated on, so every position
  indexes the string directly; `strand` is provenance only.
* **No `N` anywhere.** An `N` inside a scored window is an unscoreable position
  masquerading as a real negative. (The models' own N-padding is applied by BT4's
  adapter at the window's outer edges, which is a different thing.)

**The position convention is the whole ballgame, so it is verified, not trusted.**
A splice panel has exactly one catastrophic failure mode and it is silent: annotate
one base off and every score is misaligned, the model looks incompetent, and nothing
in the numbers says why. BT4 pins the anchor its own PWM baseline already uses — the
one convention verifiable from this repository rather than assumed about someone
else's:

| Site | Position is | Check |
|---|---|---|
| **donor** | the `G` of the intron-opening `GT` (first **intronic** base) | `sequence[i:i+2] == "GT"` |
| **acceptor** | the `G` of the intron-closing `AG` (last **intronic** base) | `sequence[i-1:i+1] == "AG"` |

Note this is **not** the exonic-boundary convention the GENCODE recipe above produces.
Since ~99% of human introns are canonical `GT-AG`, a correct panel matches almost
everywhere and a mis-anchored one matches almost nowhere — so the reader **refuses**
below 90% and names the fix:

```
only 0.0% of annotated sites carry their canonical dinucleotide at the declared
position, under the 90% floor. […] shifting donors by +1 would reach 100.0% and
acceptors by -1 would reach 100.0%; that pattern is the exonic-boundary convention:
your donors look like the LAST EXONIC base and your acceptors like the FIRST EXONIC
base. BT4 anchors on the intronic dinucleotide -- move each donor +1 and each
acceptor -1
```

Lower `--min-motif-consistency` **only** for a deliberately non-canonical panel (a U12
`AT-AC` set), and say so in `--annotation`. Never to quiet a failure — that is the
off-by-one this check exists to catch.

### B2 — The classification metrics *(landed)*

`biomodels/_stats.py` had only regression estimators (`pearson`, `spearman`,
`r2_score`, `conformal_quantile`, `empirical_coverage`). CLAUDE.md §6 names the
metrics a splice model must report — **PR-AUC / MCC / ECE / Brier, never bare
accuracy** — and none existed. Now shipped, keeping the module's constraints (**pure standard
library, no numpy** — BT4's core has zero dependencies):

```python
def pr_auc(labels, scores) -> float
def roc_auc(labels, scores) -> float
def mcc(labels, predictions) -> float
def brier_score(labels, probs) -> float
def expected_calibration_error(labels, probs, *, n_bins=10) -> float
def reliability_bins(labels, probs, *, n_bins=10) -> list[tuple[float, float, int]]
def top_k_accuracy(labels, scores) -> float
```

`top_k_accuracy` is the metric both papers headline — take the top *k* scored
positions where *k* is the number of true sites, report the fraction correct — so it
is directly comparable to published numbers. Follow the module's existing
convention: return an honest `0.0` on degenerate input rather than raising.

### B3 — The splice acceptance gate *(landed)*

`src/bt4/biomodels/splice/gate.py` and `src/bt4/pipeline/splice_gate.py` ship, and
building them against the evidence changed the sketch this section originally carried.
Four differences, each because the obvious choice would have certified the wrong thing:

**Two case types, never mixed.** `SpliceSiteCase(predicted, label, kind, group)` for
site prediction and `SpliceVariantCase(predicted, label, region, group)` for variant
effect. The exonic/intronic split that matters most for BT4 is a *variant-effect*
result; applying it to site prediction is near-degenerate, because annotated sites sit
at exon/intron boundaries by construction. So `region` lives only on the variant case,
and pooling the two types raises.

**Spearman is excluded.** On a binary label it is an exact affine function of ROC-AUC,
so it adds nothing — and at splice prevalence it is unusable as a bar: a *perfect*
classifier scores **0.055 at 0.1% prevalence**, far under the expression gate's
`min_spearman=0.30`. Transplanting that threshold would have failed a flawless model.

**The verdict is per stratum, and `overall` has no pass authority.** A blended figure
lets a backend certify on intronic strength while failing exactly where BT4 operates.

**`negative_construction` is a required argument.** Average precision's floor is the
prevalence, which is a construction choice, so a threshold without a pinned denominator
is passable by sampling fewer negatives. The report also carries
`pr_auc_skill = (AP − p) / (1 − p)`, which is 0 at no-skill and 1 at perfect for every
prevalence — deliberately not the lift ratio `AP / p`, whose ceiling drifts with
prevalence too.

Also reported: ROC-AUC (for prevalence-stability, never as a pass axis), top-k accuracy
(for comparability with the anchors in B4), MCC at BT4's own operating point, and
Brier + ECE + reliability bins — plus **Brier skill**, because ECE alone cannot
separate an informative model from a vacuous one.

#### The baselines a backend must beat

`pipeline/splice_gate.py` runs the same gate on four permanent controls
(`SPLICE_BASELINES`) and reports them in one table. On this task the dumb predictors
are unusually strong, so the comparison is the point:

| Baseline | What it is | Why it is permanent |
|---|---|---|
| `permutation` | the backend's own scores, shuffled within each stratum | the null: what this panel yields from no relationship at all |
| `gt_ag` | 1.0 where the canonical dinucleotide sits | ~99% of human introns follow it; a CNN must beat "look for GT" |
| `pwm` | BT4's shipped `ConsensusPwmSplicePredictor` | the free incumbent — the direct analogue of `cai` in the expression gate. A backend that cannot beat what `default()` already returns has not earned a PyTorch dependency, a hash-pinned weight set, or a non-commercial licence term |
| `constant` | the per-stratum base rate | perfectly calibrated, completely useless — its excellent ECE is visible in the same table instead of being a trap the reader must remember |

The comparison is **per stratum**, so beating the motif on donors cannot excuse losing
to it on acceptors. And the structural counterpart of the expression gate's "the null
model provably cannot pass": **run the PWM backend as the head and it ties the `pwm`
baseline exactly**, so `beats_every_baseline` is `False` and BT4's own default can
never be evidence for itself.

#### Two alignment traps the runner reports rather than assumes

* **Anchor offset.** A backend anchors its per-position score somewhere; one base of
  disagreement turns a good model into a hopeless one, silently. `anchor_offset` is an
  explicit input, and the report's `AlignmentDiagnostic` shows where the backend's
  score actually *peaked* around each true site. Measured: a perfect oracle shifted two
  bases reads as PR-AUC 0.006, and the diagnostic says
  `the backend's score peaks +2 from the declared position […] the anchors DISAGREE`.
* **Combined tracks.** Pangolin emits **one** `P(splice)` track and leaves `acceptor`
  all-zero. Scoring that with a donor/acceptor split would report it as perfectly
  hopeless at acceptors — an artifact of the wrapper, not a finding about the model. A
  combined track collapses to a single `"splice"` stratum, detected from the output and
  recorded in the report.

Promotion needs three conditions at once, reported separately so a failure says which:
the gate's own thresholds, beating every baseline in every stratum, and the panel being
**held out** — a panel overlapping chr 2/4/6/8/10–22/X/Y can never be `promotable`.

### B4 — Run it and compare against published anchors

```
bt4 splice-gate panel.tsv ^
  --negative-construction "all other positions in the same gene bodies" ^
  --annotation "GENCODE v44 / GRCh38" ^
  --backend pangolin
```

(`^` continues a line in Windows Command Prompt.) It prints the panel's provenance and
motif consistency, the alignment note, then the head and every baseline side by side,
then the three promotion conditions. Compare the head's row against:

| Benchmark | Pangolin | SpliceAI |
|---|---|---|
| Splice-site prediction, top-1 | **79%** | **75%** |
| Splice-site prediction, top-0.5 | **94%** | **87%** |
| Splice-site AUPRC | **0.85** | **0.77** |
| MFASS variant effect, AUPRC | **0.56** | **0.47** |
| MaPSy, Pearson r | **0.61** | **0.50** |

*(Zeng & Li 2022, all on the same held-out chromosomes. SpliceAI's own paper
headlines **top-k accuracy 0.95** on its GENCODE test set; the widely-quoted PR-AUC
0.98 could not be confirmed against the primary source, so treat it as unverified.)*

Numbers landing far below these mean the panel extraction is wrong, not the models.
Two outcomes, both fine, neither optional to report:

- **Well-calibrated** → record the reliability curve, keep 0.5, say so with evidence.
- **Systematically over-confident** (the expected result, per OpenSpliceAI) → do
  **not** silently rescale inside the adapter; that would break Part A's bit-for-bit
  fidelity claim. Report the miscalibration and set BT4's *threshold* from the data,
  replacing `DEFAULT_SITE_THRESHOLD = 0.5`. Note that SpliceAI's README cutoffs
  (0.2 / 0.5 / 0.8) are for **delta** scores, a different quantity from the raw
  per-position probability BT4 thresholds — do not import them without deriving them.
  Smith & Kitzman also found optimal thresholds **vary widely by exon and variant
  class**, and recommend normalizing by genome-wide call rate rather than nominal
  score; their 500k background SNV set makes that directly computable.

Deliverable: `docs/REVIEW_splice_calibration.md` — panel provenance, per-backend
metrics split exonic/intronic, the reliability curve, the chosen threshold and why,
and an explicit statement of what the numbers do **not** cover.

### B5 — The limits to state plainly

Two honest scope limits, both to appear wherever BT4 reports splice risk:

1. **The exonic penalty.** Median prAUC **0.419 exonic vs 0.773 intronic** (Smith &
   Kitzman 2023) across 8 predictors including these two. BT4 designs coding
   sequence — its entire regime is the weaker half. A calibrated flag must not be
   allowed to imply otherwise.
2. **The regime gap.** The panel is natural genomic sequence and natural variants.
   BT4 designs **synonymous variants of a coding sequence in a vector**. A model
   validated on natural splice sites has not been shown to discriminate *designed
   synonymous variants of the same protein* — precisely the gap CLAUDE.md already
   documents for RiboNN, for exactly the same reason.

Also worth recording: gene-model annotation choice altered SpliceAI's predictions
for **> 10% of variants** in some genes (Smith & Kitzman), projected to affect ~1 in
5 human genes — so which annotation the panel used belongs in the provenance.

---

# Verification

### After Part A — on a machine with the weights

```bash
# 1. Pins match (the A3 script) — all 12 Pangolin + all 5 SpliceAI
# 2. SpliceAI examples/output.vcf diff is clean
# 3. Gate passes: passed=True, max_abs_deviation <= 1e-3, both backends

# 4. Promotion actually flips the flag
BT4_SPLICE_USE_ATTESTED=1 python -c "
from bt4.api import available_splice_backends
print([(b.name, b.calibrated) for b in available_splice_backends()])"
# expect: consensus-pwm-baseline False, pangolin[...] True, spliceai True

# 5. And the opt-in is genuinely opt-in
python -c "
from bt4.api import available_splice_backends
print([(b.name, b.calibrated) for b in available_splice_backends()])"
# expect: every entry False

# 6. End-to-end through the real CLI
bt4 optimize MAAAGGKLQ --check-splice pangolin --use-attested-splice
```

### After Part A — in CI (no weights, must stay green)

```bash
python -m ruff check src tests scripts
python -m mypy
lint-imports
QT_QPA_PLATFORM=offscreen python -m pytest tests/ -p no:cacheprovider
```

The tests that prove nothing regressed:
`tests/test_splice.py::test_default_returns_working_labeled_model`,
`tests/test_splice_attestation.py::test_default_returns_uncalibrated_baseline`,
`tests/test_splice_pangolin.py::test_adapter_is_uncalibrated_by_default`, and the
`all_calibrated is False` assertions in `tests/test_splice_audit.py`,
`tests/test_construct_context.py`, `tests/test_app_smoke.py`. **All must still
pass.** If any turns red, promotion leaked into the default path — the failure mode
this design exists to prevent, not a test to update.

### After Part B

```bash
python -m pytest tests/test_splice_gate.py tests/test_stats.py -q
python scripts/run_splice_gate.py --panel <panel.json> --backend pangolin
```

### Determinism (invariant #7)

`FidelityAttestation` carries no timestamp by design, so re-running A6 on the same
weights must produce a **byte-identical** file and the same `content_hash()`.

# What not to do

- **Do not hand-set `fidelity_verified=True`.** `verified_predictor` is the only
  seam, and it exists so the flag cannot be assigned (§10.6).
- **Do not loosen the tolerance to force a pass.** `attest_backend` refuses anything
  looser than `1e-3`: *"not a bit-for-bit fidelity claim."*
- **Do not capture the panel with BT4's own adapter** — the gate would pass
  trivially and prove nothing.
- **Do not capture Pangolin with `custom_usage.py`** — wrong weights, wrong fold
  count, wrong channel map (finding 3).
- **Do not `pip install pangolin`** — wrong package entirely.
- **Do not commit the captured panel, the weights, or any raw per-position score** —
  the entire reason `FidelityAttestation` stores eight scalars and nothing else.
- **Do not edit `PINNED_WEIGHT_SHA256` to silence a mismatch.**
- **Do not let a real-flank score inherit the attestation** (A7.4).
- **Do not report one pooled PR-AUC.** Exonic and intronic separately, always.
- **Do not treat better input as calibration.** Real flanks fix an input defect;
  that is not a gate.
  [`RESEARCH_codon_optimization_SOTA.md`](RESEARCH_codon_optimization_SOTA.md):
  *"Do not flip any model to `calibrated=True` because context improved."*
- **Do not enable design-flow step 6 (splice auto-edit) in this change.** It
  unblocks once a backend is calibrated; shipping both at once means one change both
  grants trust and spends it.
