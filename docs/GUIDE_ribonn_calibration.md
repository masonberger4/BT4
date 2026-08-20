# Calibrating RiboNN — a plain-English walkthrough

Welcome. This guide walks you through the whole RiboNN calibration job from start to
finish, in ordinary language, one step at a time. It is the friendly companion to the
technical runbook in [`DESIGN_ribonn_calibration.md`](DESIGN_ribonn_calibration.md)
and the evidence file [`RESEARCH_ribonn_calibration.md`](RESEARCH_ribonn_calibration.md)
— same job, fewer assumptions about what you already know. Where this guide and the
runbook differ on a command, **this guide is the one that was checked against the
code**; the differences are listed in [Appendix B](#appendix-b--where-this-guide-corrects-the-runbook).

A few reassurances before we start:

- **A lot of this is free.** Steps 1–4 cost you nothing — no downloads, no licence, no
  money. They are worth doing first, and one of them might tell you to stop.
- **"It didn't work" is a real, good result here.** You are running an experiment to
  find out whether a model is any good at a job it was never built for. The honest
  expectation is that it isn't. Finding that out cheaply and writing it down *is* the
  deliverable. You have not failed if the answer is no.
- **Nothing you do here changes BT4 for users by accident.** There is no switch you
  can trip. Turning the model on for real needs a deliberate, separate step, and even
  that isn't wired up yet (see [Step 20](#step-20--the-part-that-still-wont-work)).

---

## What you are actually doing

BT4 can call a published AI model called **RiboNN**, which predicts how well a piece of
mRNA gets translated into protein. BT4 currently labels every number RiboNN produces as
**uncalibrated** — meaning "we show you this, but we don't trust it enough to let it
choose anything for you."

Your job is to run a test that decides whether that label can honestly be removed.

### The word "calibration" means three different things here

This trips up almost everyone. Keep them separate:

| | The question | In plain words | Where you stand |
|---|---|---|---|
| **1** | Does it run? | When BT4 calls RiboNN, does RiboNN actually produce numbers? | Yes — verified on real hardware |
| **2** | **Does it discriminate?** | For **one protein**, with everything else held still, does RiboNN correctly say which version is better? | **This is the real question. Nobody knows.** |
| **3** | Are its error bars honest? | When it says "3.2, give or take 0.4", is it right about the 0.4? | Only worth asking if #2 passes |

**Question 2 is the whole exercise.** Getting question 1 to work feels like progress and
is not. A model that runs perfectly and ranks randomly is worse than no model, because
it looks like an answer.

### Why the answer will probably be "no"

RiboNN was trained to compare *different genes* to each other. BT4 needs something else
entirely: comparing *different spellings of the same gene*. Every protein can be written
many ways in DNA — the protein comes out identical, only the spelling changes. BT4's job
is to pick the best spelling. RiboNN was never tested on that.

The published numbers say what to expect:

- On the job it was built for (comparing natural genes): **r² = 0.62** — good.
- On designed, made-up sequences: **r² = 0.17–0.19** — weak.
- On the one test in its own paper that isolates the coding region: **r² = 0.11** — near
  useless.

So plan for a "no", and make getting there cheap. Several steps below can end the
project in an afternoon.

---

# Part 1 — Free. Do this first.

**Nothing in this part needs the licence, the model, or any money.** It is the cheapest
way to learn how the machinery works and to catch your own mistakes before they cost
anything.

## Step 1 — Send the licence email

⚠️ **Do this before you download anything.** It is the one step that can invalidate
everything else.

RiboNN's licence does **not** simply say "free for non-commercial use". It grants use
"to any person from academic research or non-profit organizations" — that is a rule about
**who you are**, not about what you intend to do with it. If you are an independent
researcher with no university or non-profit affiliation, you may not be covered at all,
however non-commercial your motives.

Email **`patent.gos@sanofi.com`**. State your affiliation and what you intend to do. Get
the answer **in writing** before downloading the weights.

While you wait, do Steps 2–4 — they need none of this.

## Step 2 — Get BT4 itself running

You need a copy of BT4 and a terminal. From the folder where BT4 lives:

```bash
pip install -e ".[dev]"
```

> **Note:** [`INSTALL.md`](INSTALL.md) is *not* the guide for this job — that one is for
> people downloading the BT4 Studio desktop app, and never mentions Python. This job is
> a command-line job.

## Step 3 — Prove the machinery works (2 minutes, free)

BT4 ships the entire test apparatus, and **all of it runs without RiboNN**. Run it:

```bash
pytest tests/test_expression_gate.py tests/test_expression_panel.py \
       tests/test_expression_attestation.py tests/test_ribonn_sensitivity.py -q
```

You should see **90 tests pass in about a second**. For the full expression surface
(195 tests, ~8 seconds):

```bash
pytest tests/test_expression_gate.py tests/test_expression_panel.py \
       tests/test_expression_attestation.py tests/test_ribonn_sensitivity.py \
       tests/test_expression.py tests/test_expression_ribonn.py \
       tests/test_run_expression_gate.py tests/test_calibration_flag_leaks.py \
       tests/test_candidates.py -q
```

These tests are not busywork — they are the honesty guarantees you are about to rely on.
Among other things they prove that a model which has learned nothing except *"which
protein am I looking at"* sails through the easy version of the test and **fails** the
strict version. That single fact is why Step 17 insists on one particular flag.

> **A warning about green tests:** every RiboNN test uses a stand-in, not the real model.
> A green test suite tells you BT4's side is sound. It tells you **nothing** about whether
> a real RiboNN run works on your machine.

## Step 4 — Take the whole test for a dry run (free)

This is the single most useful free step, and the technical runbook doesn't mention it.

BT4 includes a **fake model** (called `null`) that always returns zero. You can push a
practice data file through the *entire* real scoring machinery using it. Nothing about
the result is meaningful — but the *plumbing* is completely real, so you learn the file
format, see the report layout, and find your typos now instead of after you've spent
money.

First make a practice data file:

```bash
python -c "
lines = ['group\tvariant_id\tcds\tmeasured\tutr5\tutr3']
for g in range(4):
    for v in range(25):
        lines.append(f'P{g}\tg{g}v{v}\tATG{\"AAA\" * (v + 1)}TAA\t{10.0 * g + v}\tGCCACC\tGCTAAT')
open('practice.tsv','w').write('\n'.join(lines) + '\n')"
```

Then run the real test against the fake model:

```bash
python scripts/run_expression_gate.py --panel practice.tsv --backend null \
       --within-group --recalibrate --json > dryrun.json
```

You will get a full report. The verdict will be a resounding "no" — correct, since the
fake model predicts nothing. What you're checking is that **it ran at all** and that you
can read the output.

Spend a few minutes reading `dryrun.json`. The part that decides everything is `verdict`:

```json
"verdict": {
  "gate_passed": false,            // did it clear the two numeric thresholds?
  "beats_every_baseline": false,   // did it beat all five dumb competitors?
  "interval_is_informative": false,// are the error bars narrower than the data itself?
  "promotable": false              // all three of the above, at once
}
```

**`promotable` is the answer.** It is true only when all three are true.

---

# Part 2 — Setting up the real model

From here on you need the licence from Step 1 and about 3 GB of downloads.

## Step 5 — Install RiboNN

RiboNN is a research codebase from Sanofi. You do **not** need to run its `Makefile`;
BT4 loads it directly.

```bash
git clone https://github.com/Sanofi-Public/RiboNN.git ~/RiboNN
cd ~/RiboNN
mamba env create -f environment.yml -y     # if you have an NVIDIA GPU
conda activate RiboNN
pip install "setuptools<82"                # works around an upstream bug
```

**No NVIDIA GPU?** Use this instead of the `mamba env create` line:

```bash
mamba create -n RiboNN -c conda-forge -y python=3.10.13 pytorch=1.13.1 \
    pytorch-lightning=1.8.5 torchmetrics=1.3.1 lightning-utilities=0.10.1 \
    mlflow=2.18.0 numpy=1.22.4 pandas=2.2.3 scikit-learn=1.0.2
```

> **Heads-up for CPU-only machines:** RiboNN loads its saved model files without telling
> PyTorch "there is no graphics card here". If the released files were saved from a GPU,
> loading them on a CPU-only box can fail outright. That is RiboNN's behaviour, not a BT4
> bug — but it means a CPU-only setup may simply not work, and you want to discover that
> at Step 8, not at Step 17.

## Step 6 — Download the model weights

The weights are ~3 GB on Zenodo. **The folder they live in must be named exactly
`models`** — RiboNN looks for that literal name.

```bash
cd ~/RiboNN
mkdir models
curl -L -o weights.zip "https://zenodo.org/records/17258709/files/weights.zip?download=1"
tar -xf weights.zip -C models
ls models/human            # expect a list of folders, each holding state_dict.pth
```

The zip contains `human/` and `mouse/` at its top level, which is why it is extracted
**into** `models` rather than alongside it.

> **If you must keep the weights elsewhere** (a different drive, say), point BT4 at them
> with the `BT4_RIBONN_WEIGHTS` environment variable. The folder still has to be *named*
> `models`, but it can live anywhere. This escape hatch is real but appears in no other
> document.

## Step 7 — Install BT4 into the *same* environment

⚠️ **Order matters and cannot be worked around.** BT4 doesn't launch RiboNN as a separate
program — it loads RiboNN's code directly into itself. If BT4 and RiboNN live in two
different Python environments, they can never see each other.

With the `RiboNN` environment still active:

```bash
export BT4_RIBONN_DIR=~/RiboNN          # Windows: set BT4_RIBONN_DIR=C:\RiboNN
cd /path/to/BT4
pip install -e ".[expression-ribonn,dev]"
```

> **Never add the `[ml]` extra here.** It requires newer PyTorch and NumPy than RiboNN
> can tolerate and will quietly break the environment you just built.

Now check:

```bash
python -c "from bt4 import api; print(api.available_expression_backends())"
```

You want `('null', 'ribonn')`. If you only see `('null',)`, BT4 can't find your files —
check that `~/RiboNN/src/` and `~/RiboNN/models/human/` both exist.

> **This check is weaker than it looks.** It confirms that two folders exist and that
> PyTorch is installed. It does **not** open a single model file. Empty folders with the
> right names will pass it. Step 8 is the real test.

## Step 8 — Prove it *actually* works

This is the first command that genuinely proves your setup. It checks all 90 model files
against known fingerprints and then runs a real prediction.

You need three real sequences: a 5′UTR, a 3′UTR (the untranslated bits at each end of an
mRNA), and a coding sequence. **The coding sequence must have a length divisible by 3 and
end in `TAA`, `TGA`, or `TAG`.**

```bash
python -c "
from bt4.biomodels.expression import RiboNNExpressionModel
m = RiboNNExpressionModel(species='human', utr5='PUT_REAL_UTR5_HERE',
                          utr3='PUT_REAL_UTR3_HERE', cell_types=('HEK293T',))
r = m.score_sequence('ATG...TAA')
print(r.score, '|', r.units, '| calibrated =', r.calibrated)
"
```

Success looks like a number, a description of its units, and `calibrated = False`.

**`calibrated = False` is correct and expected.** It stays False for this entire guide
until Step 19, and possibly forever.

Common refusals, all of them loud and clearly worded:

| What you'll see | What it means |
|---|---|
| `RiboNN requires non-empty 5' and 3' UTR context` | You left a UTR blank. They're mandatory. |
| `RiboNN requires each CDS to be length-3N ending in a stop codon` | Your sequence's length isn't divisible by 3, or lacks a stop codon. |
| `RiboNN weight file missing:` | The download is incomplete. |
| `sha256 ... != pinned ...; refusing to load` | A file is corrupted or the wrong version. BT4 refuses before loading it. |
| `RiboNN clone not found; set $BT4_RIBONN_DIR` | The environment variable isn't set in *this* terminal. |

### Also check it gives the same answer twice

```bash
python -c "
from bt4.biomodels.expression import RiboNNExpressionModel
m = RiboNNExpressionModel(species='human', utr5='...', utr3='...')
a = [x.score for x in m.score_many(['ATG...TAA'])]
b = [x.score for x in m.score_many(['ATG...TAA'])]
print(a, b, 'OK' if a == b else 'NOT REPRODUCIBLE — see Step 9')
"
```

## Step 9 — The randomness check

⚠️ **Blocking.** RiboNN has a setting called `max_shift` that nudges sequences by a random
amount. It was meant for training only, but it isn't switched off when making predictions,
and it isn't seeded — so if it's on, **the same sequence gives different answers each
time**, and every number below becomes meaningless.

```bash
python -c "
import glob, sys, pandas as pd
files = glob.glob('$HOME/RiboNN/models/*/runs.csv')
if not files: sys.exit('FAIL: no runs.csv found — the weights are not where you think')
for f in files:
    df = pd.read_csv(f)
    cols = [c for c in df.columns if 'max_shift' in c]
    print(f, {c: sorted(df[c].unique()) for c in cols} or 'NO max_shift COLUMN')
"
```

**Expect all zeros.** If anything is non-zero, stop — the adapter needs a fix first.

> The version of this check in the technical runbook has a flaw worth knowing about: if
> the file path is wrong it finds nothing, prints nothing, and exits successfully —
> which looks exactly like "all clear". The version above fails loudly instead.

---

# Part 3 — The cheap experiments that can end the project

These use only sequences already in BT4. **No purchased data.** They can save you months.

The tool is `scripts/ribonn_sensitivity.py`. It can't promote anything; every report it
writes is stamped uncalibrated.

**You must supply your own UTRs** with `--utr5` and `--utr3` (either the sequence itself
or a path to a FASTA file). This is deliberate — BT4 won't bundle UTRs, because choosing
them for you would be a hidden decision affecting every result. **Keep those files**;
only a short hash of them is recorded, and two runs with different UTRs cannot be
compared.

## Step 10 — Positive control: is the wiring even live?

Before trusting any "no effect" result, prove the machinery can detect an effect at all.
Score one sequence under two genuinely different real UTR pairs. Because most of RiboNN's
signal comes from the UTRs, the scores should move.

```bash
python scripts/ribonn_sensitivity.py --check utr-control \
    --utr5 utr5.fa --utr3 utr3.fa --utr5-alt alt5.fa --utr3-alt alt3.fa \
    --cell-type HEK293T --json > stage1_control.json
```

Look for **`harness_ok: true`**.

> **What that flag really checks:** only that the two scores are not *identical* — the
> threshold is one part in a billion. It is a wiring check, not a biology check. A
> genuine result should move the score by a visible amount; if `harness_ok` is true but
> the two numbers differ in the eighth decimal place, treat that as a failure, because
> the UTRs are barely reaching the model.

**If this fails, stop.** Nothing below is interpretable.

## Step 11 — The decisive experiment

Now the real question: **does RiboNN notice at all when you rewrite a gene's spelling?**

BT4 ships 93 real sequences — three human proteins, each written 31 different ways by
real commercial codon-optimization tools. No measurements needed; you're only asking
whether RiboNN's scores *move*.

```bash
python scripts/ribonn_sensitivity.py --check cds-spread \
    --fasta scripts/data/ranaghan2021_tab4.fasta \
    --utr5 utr5.fa --utr3 utr3.fa --cell-type HEK293T --json > stage1_spread.json
```

Read the results **in this order**:

1. **`within_over_between`** — how much scores vary *between spellings of one protein*
   compared with *between different proteins*. This is the number that matters. A tiny
   ratio means RiboNN is mostly just recognising which protein it's looking at.
2. **`median_abs_gc3_spearman`** — see Step 12.
3. **`responds_to_synonymous_change`** — read this **last**, and don't over-trust it.

> **About `responds_to_synonymous_change`:** it is true whenever the variation exceeds
> one part in a billion. It does *not* mean "responds usefully" — it means "isn't
> literally frozen". A `true` here is permission to keep reading, not evidence of skill.

**If the spread is at the noise floor: stop, and write it up.** You have measured that
RiboNN is blind to the only thing BT4 can change. No dataset on earth fixes that, and
this is a genuine, publishable finding that closes the question honestly.

## Step 12 — Is it just measuring GC content?

The same command already answered this. `cds-spread` reports how closely RiboNN's scores
track five simple properties: GC content, GC3, CAI, tAI, and length.

If **`median_abs_gc3_spearman`** is high, RiboNN is largely acting as a GC detector in
this regime. That matters enormously, because **BT4 already optimizes GC directly, for
free, instantly**. A model whose only skill is GC adds nothing at all.

## Step 13 — Two more cheap sanity checks

```bash
# Does it prefer a highly-optimized spelling over a deliberately bad one?
python scripts/ribonn_sensitivity.py --check direction \
    --utr5 utr5.fa --utr3 utr3.fa --cell-type HEK293T --json

# Does its score change smoothly as sequences change gradually?
python scripts/ribonn_sensitivity.py --check ladder \
    --utr5 utr5.fa --utr3 utr3.fa --cell-type HEK293T --json
```

Read these asymmetrically: a clear preference proves little (optimized and deoptimized
sequences differ in many ways at once), but **a coin-flip result is genuinely bad news**,
and a jagged, erratic response means the scores can't be used for ranking even if they
aren't zero.

### 🚦 Decision point

Continue to Part 4 **only if** Step 11 shows real spread and Step 12 shows it isn't
purely GC. Otherwise stop, write down the numbers, and close the question. **That is a
successful outcome**, reached for the price of an afternoon instead of a research grant.

---

# Part 4 — Getting real measurements

To judge whether RiboNN's ordering is *correct*, you need real experimental measurements
of different spellings of the same protein. This is the expensive part.

## Step 14 — Find the best available dataset

**No public dataset fully fits.** This is a "best partial answer, honestly labelled"
exercise. In priority order:

| Dataset | Why | Catch |
|---|---|---|
| **Mauger 2019** (PNAS, Dataset S1) | Best fit if it works: human cells, 4 proteins, fixed UTRs | **Unknown whether it lists per-variant measurements. ~10 minutes to check — do this first.** |
| **PERSIST-seq** (HuggingFace `morrislab/mrl-hl-lbkwk`) | Exactly the right kind of measurement (ribosome load); 203 rows; CC BY 4.0 | Only one protein per arm, so the claim you can make is narrower |
| **iCodon** (GEO `GSE207584`) | 100 proteins — by far the most | Zebrafish, and it measures mRNA decay rather than translation. RiboNN has no zebrafish model. |

Skip Mordstein 2020 (raw sequencing data only — days of reprocessing) and CodonBERT MLOS
(its UTRs were never published, and RiboNN needs them).

⚠️ **The sequences already in BT4 have no measurements.** The Ranaghan panel you used in
Step 11 is sequence-only — that paper measured exactly one sequence, in bacteria. It is
perfect for Part 3 and useless here.

> **How much would making your own cost?** Roughly **$15–19k in DNA synthesis alone** for
> ~300 sequences, before any lab work; 6–12 months. Out of scope, but worth knowing when
> weighing whether to continue.

### Choosing the right kind of measurement

Not all "expression" measurements can be compared to RiboNN's output. Use **ribosome
load** or **log(protein ÷ mRNA)**. Do **not** use raw protein output — it mixes together
how well the mRNA is translated with how long it survives, and RiboNN's output has
already divided that second part out. Comparing them is comparing different quantities.

## Step 15 — Write the data file

A tab-separated file, one row per measured sequence. **Six columns are required:**

```
group	variant_id	cds	measured	utr5	utr3
```

Three more are optional but strongly recommended: `readout`, `cell_type`, `species`.

```
group	variant_id	cds	measured	readout	cell_type	utr5	utr3	species
NLUC	persist_001	ATGGTC...TAA	1.42	mean_ribosome_load	HEK293T	ACATTT...	GCTCGC...	human
NLUC	persist_002	ATGGTG...TGA	1.19	mean_ribosome_load	HEK293T	ACATTT...	GCTCGC...	human
```

What each column means:

- **`group`** — **the protein.** This is the most important column. It's how the test
  keeps different spellings of one protein from leaking between halves of the experiment,
  and how it knows to compare like with like.
- **`variant_id`** — a unique name for this row. Duplicates are refused.
- **`cds`** — the coding sequence: `A`/`C`/`G`/`T` only, length divisible by 3, ending in
  a stop codon.
- **`measured`** — the experimental result. **Bigger must mean more expression.** Take
  logs of ratios first.
- **`utr5` / `utr3`** — the experiment's real UTRs, never blank.

The reader is deliberately strict, and **refuses bad rows rather than skipping them** —
so your row count always means what you think it means. It will reject: a misspelled
column name, a missing column, a blank required field, a non-numeric measurement, a
sequence with any letter other than ACGT, a length not divisible by 3, a missing stop
codon, a duplicate `variant_id`, or a species other than human/mouse.

**Two size limits nothing else warns you about:** the 5′UTR must be ≤ **1381** letters,
and the coding sequence plus 3′UTR together ≤ **11937**. Over either, the row is refused.

Now validate it — free, no model needed:

```bash
python -c "
from bt4 import api
p = api.read_panel('panel.tsv')
print(p.content_hash())
print(p.describe())"
```

Save that `content_hash` — it's a fingerprint proving which exact data you used.

Also add a `.LICENSE.md` file next to your data (copy the shape of the existing
`scripts/data/*.LICENSE.md` files: `## Source` with the full citation and DOI,
`## License` with the licence name, URL and an explicit "Changes made:", and
`## Contents`), and add a row to `THIRD_PARTY_DATA.md`.

## Step 16 — How much data do you actually need? (read this before buying anything)

The technical runbook says "at least 4 proteins and about 90 rows". **That floor is too
low, and the failure it causes looks exactly like a real negative result.**

Here's why. The test splits your data in half by protein. Half calibrates, half tests. So
90 rows leaves only ~45 rows to measure error-bar accuracy on — and it demands that
accuracy land within 90% ± 5%. On 45 rows, ordinary statistical luck misses that band
about half the time **even when the model is excellent**.

Measured directly against BT4's real gate, using a simulated model with near-perfect
within-protein ranking (median ρ = 0.96):

| Data size | Chance a genuinely good model passes | What actually failed |
|---|---|---|
| 4 proteins, **92 rows** | **44%** | error-bar band — never the ranking |
| 4 proteins, **200 rows** | **82%** | " |
| 6 proteins, **198 rows** | **82%** | " |
| 10 proteins, **200 rows** | **79%** | " |
| 20 proteins, **900 rows** | **99%** | " |

The ranking part passed **100% of the time** at every size. Coverage was the entire story.

**Aim for ~200 rows, not 90.** This agrees with BT4's own research doc, which computes
that a ±0.05 coverage claim needs about 102 rows — and after the 50% split, 200 rows is
what leaves you 100. At 90 rows you are more likely than not to record "RiboNN can't do
this" about a model that can.

Other floors worth knowing: at least **9 rows** in the calibration half (below that the
error bars are reported as infinite), at least **2 proteins**, and at least one protein
with 2+ variants in the test half.

> **A quirk to plan around:** the split isn't random — it sorts your protein names
> **alphabetically** and puts the first half in calibration. So what you *name* your
> proteins determines which ones get tested. Name them before you look at any results.

---

# Part 5 — Locking it in and running the test

## Step 17 — Write down your thresholds *before* running

This step is pure self-discipline, and it is the difference between a validation and a
fishing expedition. If you run the test, see a fail, nudge a threshold, and re-run, you
have quietly converted an experiment into a search for a number you like. On a few
hundred rows that is enough to manufacture a false pass.

**Nothing in the code enforces this. It is entirely on you.**

Create `docs/ribonn_gate_preregistration.md` and **commit it before Step 18**:

```
panel_sha256:         <the content_hash from Step 15>
species:              human
cell_type:            HEK293T
readout:              mean_ribosome_load
group_key:            protein
within_group:         true
recalibrate:          true
target_coverage:      0.90
coverage_tolerance:   0.05
min_spearman:         0.30      # our own pre-commitment, NOT an industry standard
decision_rule:        the model's cautious (lower-bound) score must beat every baseline's
                      best score
baselines:            permutation, cai, gc3, length, constant
width_rule:           median error-bar width ÷ spread of the data must be < 1.0
runs_permitted:       1
```

Because the gate's own record is thin, **also write down here**: which cell types you
used, your `top_k` setting, and which UTR files. None of those appear in the output file,
and you will not remember them later.

## Step 18 — Run the test. Once.

```bash
python scripts/run_expression_gate.py \
    --panel panel.tsv --backend ribonn --species human --cell-type HEK293T \
    --within-group --recalibrate \
    --min-spearman 0.30 --target-coverage 0.90 --coverage-tolerance 0.05 \
    --baselines permutation,cai,gc3,length,constant \
    --num-workers 0 --json > gate_result.json
```

### The three flags you must not forget

**`--within-group`** — ⚠️ **off by default, and forgetting it wastes the entire run.**
Without it, the test measures whether RiboNN can tell *different proteins* apart. It can;
that's what it was trained for; and it is not what BT4 needs. With it, the test asks the
real question: within a single protein, can it rank the spellings? Forgetting it produces
a flattering, meaningless number — and the promotion step will refuse it afterwards
anyway, so you'd have to run again. The only warning is one line on stderr, easy to miss
if you're redirecting output to a file.

**`--recalibrate`** — also off by default. RiboNN's output is on its own arbitrary scale;
your measurements are in real units. This fits the straight line connecting them. Without
it, the error-bar half of the test is comparing two different rulers.

**`--cell-type`** — ⚠️ **the silent one.** Leave it out and RiboNN averages **all 78 human
cell types**. Against a panel measured in one cell line, that is a different quantity
entirely — and there is no error, no warning, and no cross-check against your file's own
`cell_type` column. It runs cleanly to a wrong verdict.

> **Use the script, not the `bt4 expression-gate` shortcut**, if you want the record file.
> The shortcut command exists and works, but it has **no `--json` flag**, so it cannot
> produce `gate_result.json` — and it also can't set `--top-k` or `--baselines`.

Keep `gate_result.json`, `stage1_spread.json` and your pre-registration file together.
Those three files *are* the record.

## Step 19 — Read the verdict

Open `gate_result.json` and look at `verdict.promotable`. It is true only when three
separate things all hold:

1. **`gate_passed`** — the ranking cleared 0.30 **and** the error bars landed within 90% ± 5%.
2. **`beats_every_baseline`** — RiboNN beat all five stand-ins. **The hard one is CAI**, a
   simple arithmetic formula BT4 already computes instantly, for free, *inside* the
   optimizer. A model that can't beat it adds nothing. Note the bar is deliberately
   demanding: RiboNN's *cautious lower estimate* must exceed each baseline's *best*
   estimate.
3. **`interval_is_informative`** — the error bars are narrower than the spread of the data
   itself. This catches a genuinely sneaky failure: a model that always predicts the same
   number gets its error bars technically "right" while saying nothing at all.

### What each outcome means

| Result | What it means | What to do |
|---|---|---|
| **`promotable: true`** | RiboNN really does rank spellings of a protein | Go to Step 20 |
| Ranking works, but doesn't beat **CAI** | No gain over what BT4 already does for free | Stays uncalibrated. Publish the numbers. |
| **Ranking fails** | RiboNN doesn't do BT4's job | Stays uncalibrated **for this regime**. Write it up — the question is now *answered*, not pending. |
| Ranking works, error bars don't | It orders things usefully but misstates its confidence | Report both halves. Default is still uncalibrated. |
| Only one protein available | Only a much narrower claim is possible | Report the narrower claim; do not promote |

**Note the asymmetry:** a pass proves something specific and scoped. A fail at 90 rows
might just be Step 16's coin flip. Check your row count before concluding "no".

---

# Part 6 — If it passed

## Step 20 — Record it properly

The only sanctioned way to flip the switch is an **attestation** — a small signed record
saying "this model passed this test, on this data, in this scope."

```python
# promote.py — run once, on the machine that ran the test.
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

⚠️ **Two things to be careful about here.**

First, this script **runs the whole test a second time**, which costs another full scoring
pass and can silently disagree with Step 18 if you type different settings. Make the
`species` and `cell_types` match Step 18 **exactly**.

Second — and this is the sharpest edge in the whole procedure — **the `species` and
`cell_types` you type here are not checked against how the test actually ran.** They are
free text, copied straight into the permanent record. Run the test across all 78 cell
types, then type `cell_types=("HEK293T",)` here, and you produce a committed record
whose stated scope is simply false, which every later check will accept. Nothing catches
this but you.

The attestation layer will, however, refuse outright to record: a failing test, a run
done without `--within-group`, a model that didn't beat every baseline, error bars as
wide as the data, or thresholds set below the built-in floors.

Then update [`CLAUDE.md`](../CLAUDE.md) §6/§9, item 11 in
[`NEXT_SESSION.md`](NEXT_SESSION.md), and [`CHANGELOG.md`](../CHANGELOG.md).

## Step 21 — The part that still won't work

**Be careful what you promise anyone.** Even after a passing test and a committed
attestation, **nothing changes for BT4's users.**

The function that flips the switch is `verified_predictor`, and **no part of BT4 calls it**
outside its own tests. A user must invoke it by hand:

```python
from bt4.biomodels.expression import (
    RiboNNExpressionModel, load_expression_attestation, verified_predictor,
)
model = RiboNNExpressionModel(species="human", utr5="...", utr3="...",
                              cell_types=("HEK293T",))
model = verified_predictor(model, load_expression_attestation("ribonn_attestation.json"))
assert model.calibrated
```

Wiring that into BT4's normal paths is **separate work — budget it separately.** (The
splice side has exactly the same gap, recorded as item 10.)

And keep the claim honest and scoped. An attestation earned on HEK293T does **not**
certify a model averaging all 78 cell types. The honest sentence is always something like:
*"calibrated for ranking synonymous variants of a known protein, in cell type X, measured
by readout Y."* The broader claim — ranking sequences for a protein nobody has measured —
**cannot be honestly tested at all**, because it would need ~100 held-out proteins that
don't exist.

---

# Appendix A — The traps, in one place

| # | Trap | Why it hurts |
|---|---|---|
| 1 | Two Python environments | BT4 loads RiboNN's code *into itself*. Separate environments can never see each other. |
| 2 | Adding the `[ml]` extra | Pulls in newer PyTorch/NumPy than RiboNN tolerates. |
| 3 | Weights folder not named `models` | RiboNN hard-codes that name; the fallback needs admin rights on Windows. |
| 4 | Trusting `('null', 'ribonn')` | Only proves two folders exist. Empty folders pass. |
| 5 | Forgetting `--within-group` | Measures the wrong thing entirely and wastes the run. Warns on stderr only. |
| 6 | Forgetting `--recalibrate` | Compares two different scales; error bars become garbage. |
| 7 | Forgetting `--cell-type` | Silently averages 78 cell types. No error, no warning, wrong answer. |
| 8 | Mismatched FASTA headers in Step 11 | Grouping is "text before the first `|`". Get it wrong and every sequence becomes its own group; the crash message never mentions headers. |
| 9 | Reading `responds_to_synonymous_change` as biology | It's a one-part-in-a-billion check. |
| 10 | A 90-row dataset | A good model fails ~56% of the time. See Step 16. |
| 11 | Re-running with adjusted thresholds | Turns validation into a search. Nothing enforces the "run once" rule. |
| 12 | Typing a scope in `promote.py` that doesn't match the run | Produces a permanent record that lies, which every later check accepts. |
| 13 | Expecting a green test suite to mean a real run works | Every RiboNN test uses a stand-in. |
| 14 | Comparing against raw protein output | Re-introduces the exact quantity RiboNN's output divides out. |
| 15 | Expecting a pass to change anything for users | Nothing calls `verified_predictor`. Separate work. |

---

# Appendix B — Where this guide corrects the runbook

Every command in this guide was checked against the source, and several were executed.
Where [`DESIGN_ribonn_calibration.md`](DESIGN_ribonn_calibration.md) differs:

1. **"The adapter reproduces upstream bit-for-bit."** The runbook's summary table marks
   this ✅ done. **No RiboNN fidelity gate, capture script or attestation exists in the
   repo** — "bit-for-bit" is the *splice* side's vocabulary, and the expression
   attestation module explicitly contrasts itself with it. What has been established is
   that real end-to-end runs happen and produce numbers. That is why the runbook's own
   Step 1.2 (fold semantics) exists at all. This guide says "it runs", not "it matches".
2. **`bt4 expression-gate` cannot produce `gate_result.json`.** The runbook offers it as
   an equivalent to the script at Stage 4, but it has **no `--json` flag** (nor
   `--baselines`, `--top-k`, `--batch-size`). Use the script for the record.
3. **`harness_ok` does not mean "the scores differ substantially".** It is a
   not-equal test at a 1e-9 floor — about six orders of magnitude away from the runbook's
   wording.
4. **The ~90-row panel floor is too small** — measured above; it fails a good model more
   often than it passes one. The runbook's own research doc computes ~102 rows for the
   coverage claim it asks for.
5. **`$BT4_RIBONN_WEIGHTS` exists** and lets the weights live outside the checkout. It is
   documented in no markdown file in the repo, only in source and two error messages.
6. **The `max_shift` check passes silently on a wrong path** — an empty match prints
   nothing and exits 0, indistinguishable from "all clear". Fixed in Step 9.
7. **CPU-only machines may not work at all.** RiboNN calls `torch.load` without
   `map_location`. Recorded in `NEXT_SESSION.md` but absent from the runbook, which offers
   a CPU-only environment as a supported option.
8. **The claim that CI proves `NullExpressionModel` "cannot pass in either mode"** is not
   backed by a test — it is checked in pooled mode only, and the fixture's all-singleton
   groups would make the within-group version raise rather than fail.
9. **Sizing framing.** The runbook presents "≥3 variants per protein, ≥4 proteins" as
   following from the gate's arithmetic. Only the ≥9-row conformal floor does; the gate's
   real minimums are 2 groups, 2 test cases, and one test group with 2+ members.
10. **A missing quote attribution.** "Feeding RiboNN a real UTR makes it runnable. Neither
    is a gate." is credited to `RESEARCH_ribonn_calibration.md`; it is in
    `RESEARCH_codon_optimization_SOTA.md`.

---

# Glossary

| Term | Plain meaning |
|---|---|
| **CDS** | The coding sequence — the stretch of DNA that spells out the protein. |
| **UTR** | Untranslated region. The bits at each end of an mRNA that aren't the protein recipe but strongly affect how well it's read. |
| **Synonymous variant** | A different spelling of the same protein. Identical product, different DNA. |
| **TE (translation efficiency)** | How hard the cell's machinery works on one mRNA. RiboNN's output. |
| **Calibrated** | In BT4, a flag meaning "this model earned the right to influence decisions, by passing a test on real data." |
| **The gate** | The pass/fail test itself. |
| **Panel** | Your data file of measured sequences. |
| **Baseline** | A deliberately dumb competitor the model must beat to prove it's adding something. |
| **Spearman (ρ)** | A score from −1 to +1 for whether two rankings agree. Only cares about order. |
| **Coverage** | Of everything the model said "I'm 90% sure", how often was it right? |
| **Conformal interval** | The error bar. Built from how wrong the model was on data it hadn't seen. |
| **Within-group** | Comparing only spellings of the *same* protein — BT4's actual job, and the strict version of the test. |
| **Pooled** | Comparing everything at once, including across different proteins — the easy version, and the wrong question. |
| **Attestation** | The signed record saying a model passed, and in exactly what scope. |
