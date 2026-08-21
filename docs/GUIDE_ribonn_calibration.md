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
  can trip. Turning the model on for real needs a deliberate, separate opt-in — an
  environment variable or a checkbox — and it only works inside the exact scope your
  test measured (see [Step 21](#step-21--using-it-this-part-now-works)).
- **Don't recognise a word?** There is a one-line-per-term [Glossary](#glossary) at the
  very end. CDS, UTR, Spearman, GC3, r², coverage, conformal interval and attestation are
  all in it.

**The shape of the job:** **Part 1** free checks · **Part 2** install the model ·
**Part 3** cheap experiments that can end the project · **Part 4** get real measurements ·
**Part 5** run the test, once · **Part 6** if it passed. Parts 1 and 3 are where most
projects stop, and that is by design.

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
- On designed, made-up sequences: **r² = 0.17–0.19** — weak. (Re-training it on such
  data recovers 0.49–0.50 — but BT4 uses the model frozen, exactly as published, so that
  number is not available to you here.)
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

RiboNN's licence does **not** simply say "free for non-commercial use". It restricts
**both** what you may do with it (academic-research, non-commercial purposes only) **and
who you may be** — it grants use "to any person from academic research or non-profit
organizations". Most non-commercial licences only do the first; this one does both. So if
you are an independent researcher with no university or non-profit affiliation, you may
not be covered at all, however non-commercial your motives.

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
strict version. That single fact is why Step 18 insists on one particular flag.

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

From here on you need the licence from Step 1. The weights themselves are about
**200 MB**; the conda environment that PyTorch lives in is a few GB on top of that.

## Step 5 — Install RiboNN

RiboNN is a research codebase from Sanofi. You do **not** need to run its `Makefile`;
BT4 loads it directly.

> **You need conda, and you probably already have it.** Miniconda, Anaconda and Miniforge
> all work. Check with `conda --version`; if that errors, install
> [Miniforge](https://github.com/conda-forge/miniforge) or Miniconda, close your terminal,
> and open a new one.
>
> **Which one you have changes two small things.** Miniforge sets conda-forge as its
> default *and only* channel and ships `mamba` already, so a Miniforge user can skip the
> next paragraph. On Miniconda/Anaconda the `defaults` channel is still consulted when
> creating an environment from a file — RiboNN's `environment.yml` lists `conda-forge` but
> no `nodefaults` — so if you would rather not involve Anaconda's channel at all, add a
> `- nodefaults` line to its `channels:` list before running `conda env create`. Either
> way the packages you get come from conda-forge.
>
> **`mamba` is optional.** It is a faster drop-in for `conda`, and RiboNN's own docs
> assume it. If `mamba --version` errors, either use `conda` wherever this guide says
> `mamba` — the commands are otherwise identical — or install it once with
> `conda install -n base -c conda-forge mamba -y`. You are not missing speed by using
> `conda`: since **conda 23.10** the libmamba solver is conda's default ("With this
> 23.10.0 release we are changing the default solver of conda to conda-libmamba-solver!"),
> so a current conda already solves the way mamba does. Only on an older conda is it worth
> running `conda install -n base -c conda-forge conda-libmamba-solver -y` and
> `conda config --set solver libmamba`.

> ⚠️ **Give RiboNN its own environment — do not reuse one you built for Pangolin.**
> RiboNN's `environment.yml` pins **torch 1.13.1** and **numpy 1.22.4** exactly, while
> BT4's `splice-pangolin` extra declares **torch >= 2.2**, so installing that extra here
> would upgrade torch straight off RiboNN's pin. That is what the `-n RiboNN` below is
> for, not tidiness. **This is the whole scope of the claim** — `torch>=2.2` is BT4's own
> floor (shared with `[ml]`); upstream Pangolin declares no torch requirement at all, and
> nothing checks a torch version at runtime. Whether Pangolin would actually *run* on
> 1.13.1 is untested, so treat separate environments as the safe arrangement rather than a
> measured impossibility.
>
> **SpliceAI is a different question and is not answered here.** Its extra is
> `tensorflow>=2.6,<2.16` with no torch at all, so the reason above does not reach it. The
> plausible collision is numpy — TensorFlow 2.15 requires `numpy>=1.23.5`, above RiboNN's
> 1.22.4 pin — but nobody has tried it, and `NEXT_SESSION.md` records a *measured* finding
> that Pangolin and SpliceAI share one environment more easily than their metadata
> suggested. Do not assume; if you need both, measure. You will install BT4 into this new environment as well as
> the old one (Step 7); an editable install points both at the same source folder, so
> nothing is duplicated except the dependency set. The consequence to know: inside the
> `RiboNN` environment the wrapped splice CNNs are simply not available and BT4 falls back
> to its PWM baseline, so a single run cannot use Pangolin *and* RiboNN. Run them
> separately and compare.

```bash
git clone https://github.com/Sanofi-Public/RiboNN.git ~/RiboNN
cd ~/RiboNN
mamba env create -f environment.yml -y     # if you have an NVIDIA GPU
conda activate RiboNN
pip install "setuptools<82"                # works around an upstream bug
```

> ⚠️ **Windows: clone into your home folder, not `C:\`.** Use
> `cd %USERPROFILE%` then `git clone https://github.com/Sanofi-Public/RiboNN.git RiboNN`,
> giving `C:\Users\<you>\RiboNN`. Step 9's check defaults to `~/RiboNN/models`, so a
> clone at `C:\RiboNN` makes it look in the wrong place — and that check is the one whose
> whole job is to fail loudly rather than find nothing. Use **Anaconda Prompt** (or
> Miniforge Prompt) from the Start menu, not PowerShell: `conda activate` is wired up
> there and may not be elsewhere. The `\` line continuations below are bash — in
> `cmd.exe`, put the command on one line.

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
> at Step 8, not at Step 18.

## Step 6 — Download the model weights

The weights are about **200 MB** on Zenodo. Name the folder they go into `models` —
that is the name RiboNN looks for, and using it avoids a fallback that needs
administrator rights on Windows.

```bash
cd ~/RiboNN
mkdir models
curl -L -o weights.zip "https://zenodo.org/records/17258709/files/weights.zip?download=1"
unzip -q weights.zip -d models
ls models/human            # expect a list of folders, each holding state_dict.pth
```

The zip contains `human/` and `mouse/` at its top level, which is why it is extracted
**into** `models` rather than alongside it.

> ⚠️ **Use `unzip`, not `tar`.** On Linux, `tar -xf weights.zip` fails — GNU `tar` cannot
> read zip archives, and you get `This does not look like a tar archive` after the
> download. (It happens to work on macOS and Windows, which ship a different `tar`.) No
> `unzip`? Install it (`sudo apt install unzip`), or use
> `python -m zipfile -e weights.zip models`, which needs nothing extra.

> **Windows form of the same block**, in `cmd.exe`. `tar` and `curl` both ship with
> Windows 10 and later, so the Linux `tar` caveat above does not apply to you:
>
> ```
> cd %USERPROFILE%\RiboNN
> mkdir models
> curl -L -o weights.zip "https://zenodo.org/records/17258709/files/weights.zip?download=1"
> tar -xf weights.zip -C models
> dir models\human
> ```
>
> If `tar` complains, `python -m zipfile -e weights.zip models` works anywhere Python does.

> **If you must keep the weights elsewhere** (a different drive, say), point BT4 at them
> with the `BT4_RIBONN_WEIGHTS` environment variable — an escape hatch that appears in no
> other document. On Linux and macOS the folder can then have any name; BT4 bridges the
> difference itself. Only on Windows does the name need to be literally `models`, because
> the bridge it uses there needs Developer Mode or an elevated prompt.

## Step 7 — Install BT4 into the *same* environment

⚠️ **Order matters and cannot be worked around.** BT4 doesn't launch RiboNN as a separate
program — it loads RiboNN's code directly into itself. If BT4 and RiboNN live in two
different Python environments, they can never see each other.

With the `RiboNN` environment still active:

```bash
export BT4_RIBONN_DIR=~/RiboNN   # Windows: set BT4_RIBONN_DIR=%USERPROFILE%\RiboNN
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

> **Every step from here on runs inside this environment.** If you open a new terminal,
> run `conda activate RiboNN` and set `BT4_RIBONN_DIR` again before doing anything else —
> otherwise you are silently running against a different Python. On Windows the variable
> is set with `set BT4_RIBONN_DIR=%USERPROFILE%\RiboNN` (no spaces around the `=`), and
> `set` lasts only for that window; `setx BT4_RIBONN_DIR "%USERPROFILE%\RiboNN"` makes it
> permanent but does **not** affect the window you type it in, so open a new one.

> **If a tool or agent is running these commands for you**, note that `conda activate`
> changes shell state, and anything that starts a fresh shell per command will silently
> lose it — running Step 8 against the wrong Python and reporting a confusing error. The
> stateless form works everywhere and is worth preferring in that setting:
>
> ```bash
> BT4_RIBONN_DIR=~/RiboNN conda run -n RiboNN --no-capture-output python -c "..."
> ```
>
> **Both** pieces of state have to ride along, which is easy to get half-right: `conda
> run` fixes the interpreter, but it inherits the parent process environment, so a shell
> that lost the `conda activate` lost the `export` too — and the command then runs the
> right Python against an unset `BT4_RIBONN_DIR`, failing with "RiboNN clone not found"
> for a reason that looks nothing like the cause. On Windows, set the variable
> permanently once (`setx BT4_RIBONN_DIR "%USERPROFILE%\RiboNN"`, then open a new
> terminal) so every shell inherits it. `--no-capture-output` is what lets you see
> progress while it runs rather than only at the end.

### First: getting the sequences you'll need

Every remaining step needs **real UTRs** — the untranslated stretches at each end of an
mRNA. RiboNN refuses to run without them, and BT4 deliberately doesn't supply any, because
picking them for you would be a hidden choice affecting every number you produce.

You need **two pairs**: a main pair, and a genuinely different second pair for the control
in Step 10.

**Where to get them.** Pick a well-studied human transcript — human beta-globin (*HBB*) is
the classic reporter choice, and something structurally unlike it (say *ACTB*) makes a good
second pair. On [Ensembl](https://www.ensembl.org): search the gene → pick the main
transcript → **Sequence** → the 5′UTR and 3′UTR are marked in the exon table.

**Put each in its own file** in FASTA format — a header line starting with `>`, then the
letters:

```
>utr5
ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC
```

Save four files: `utr5.fa`, `utr3.fa`, and `alt5.fa`, `alt3.fa` for the second pair.
Wherever a command below says `--utr5 utr5.fa`, you can also paste the sequence directly
instead of a filename.

⚠️ **Keep these files.** Only a 12-character hash of them is recorded in any report, so
two runs made with different UTRs are not comparable and you will have no way to tell them
apart afterwards.

**You also need one coding sequence** for Step 8 — any real CDS will do, including one
from `scripts/data/ranaghan2021_tab4.fasta` in the BT4 repo. It must be `A`/`C`/`G`/`T`
only, a length divisible by 3, and end in `TAA`, `TGA`, or `TAG`.

## Step 8 — Prove it *actually* works

This is the first command that genuinely proves your setup. It checks all 90 model files
against known fingerprints and then runs a real prediction.

```bash
python -c "
from bt4.biomodels.expression import RiboNNExpressionModel
m = RiboNNExpressionModel(species='human', utr5='PUT_REAL_UTR5_HERE',
                          utr3='PUT_REAL_UTR3_HERE', cell_types=('HEK293T',))
r = m.score_sequence('PUT_REAL_CDS_HERE')
print(r.score, '|', r.units, '| calibrated =', r.calibrated)
"
```

Replace all three `PUT_REAL_..._HERE` strings before running — anything that isn't real
`A`/`C`/`G`/`T` is refused.

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
m = RiboNNExpressionModel(species='human', utr5='PUT_REAL_UTR5_HERE',
                          utr3='PUT_REAL_UTR3_HERE')
a = [x.score for x in m.score_many(['PUT_REAL_CDS_HERE'])]
b = [x.score for x in m.score_many(['PUT_REAL_CDS_HERE'])]
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
import glob, os, sys, pandas as pd
root = os.environ.get('BT4_RIBONN_WEIGHTS') or os.path.expanduser('~/RiboNN/models')
files = glob.glob(os.path.join(root, '*/runs.csv'))
if not files: sys.exit(f'FAIL: no runs.csv found under {root} — check where Step 6 put the weights')
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
Score one sequence under two genuinely different real UTR pairs — the scores should move.

Why they should: **per nucleotide**, the UTRs carry most of RiboNN's information density
(the paper's split is 67/31/2 across 5′UTR/CDS/3′UTR). That qualifier is load-bearing and
often dropped: *integrated over their length*, the split is 22/73/5, so the CDS is
actually the **majority** of the total attributed signal. Both numbers are real. Swapping
UTRs is simply the biggest change you can make in one step, which is what makes it a good
control — not evidence that the CDS doesn't matter.

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

BT4 ships 93 real sequences: three human proteins, each with its **natural** coding
sequence plus 30 optimizer outputs — three codon-optimization algorithms run ten times
each (31 records per protein). No measurements needed; you're only asking whether RiboNN's
scores *move*.

> Two honest caveats about this panel. The ten runs of one algorithm measure that tool's
> run-to-run variability, so you are seeing **three** design processes, not 31. And the
> algorithms are **anonymized** in the source paper — the repo's licence file explicitly
> forbids mapping them to named commercial tools, because that mapping was never
> published.

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
    --utr5 utr5.fa --utr3 utr3.fa --cell-type HEK293T

# Does its score change smoothly as sequences change gradually?
python scripts/ribonn_sensitivity.py --check ladder \
    --utr5 utr5.fa --utr3 utr3.fa --cell-type HEK293T
```

Leave `--json` off for these two — the plain output is easier to read and ends with a
`note:` line telling you how to interpret it.

- **`direction`** prints how many non-tied proteins it preferred the max-CAI design in.
  About half is the coin flip.
- **`ladder`** prints `spearman(score, CAI)`. Near zero, or a score that jumps around
  while CAI climbs steadily, means the scores can't be used for ranking.

Read these asymmetrically: a clear preference proves little (optimized and deoptimized
sequences differ in many ways at once), but **a coin-flip result is genuinely bad news**,
and a jagged, erratic response means the scores can't be used for ranking even if they
aren't zero.

### 🚦 Decision point

Continue to Part 4 **only if both** of these hold:

- **`within_over_between` ≥ 0.2** (Step 11) — RiboNN's response to rewriting one protein is
  at least a fifth of its response to changing protein entirely.
- **`median_abs_gc3_spearman` ≤ 0.7** (Step 12) — its within-protein response is not almost
  entirely explained by GC3.

**Those two numbers are our own pre-commitment, not an industry standard** — no standard
exists. Write them down *before* you look at the output, for the same reason Step 17 makes
you pre-register the gate's thresholds.

If `within_over_between` comes back as `null`, the *between*-protein spread was itself at
the floor — that means the harness isn't working, so re-check Step 10 rather than reading
it as a result.

Otherwise stop, write down the numbers, and close the question. **That is a successful
outcome**, reached for the price of an afternoon instead of a research grant.

---

# Part 4 — Getting real measurements

To judge whether RiboNN's ordering is *correct*, you need real experimental measurements
of different spellings of the same protein. This is the expensive part.

## Step 14 — Find the best available dataset

> **Read [Step 16](#step-16--how-much-data-do-you-actually-need-read-this-before-buying-anything)
> first.** It sets the size floor you should judge these candidates against — and it is
> roughly double what the technical runbook says. Two of the three options below are
> chosen partly on size, so knowing the floor first changes the decision.

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
> ~300 sequences — and **mid-five to low-six figures all-in** once transcription,
> transfection and sequencing are counted; 6–12 months. Out of scope, but worth knowing
> when weighing whether to continue.

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
within-protein ranking (median ρ = 0.96), 1,500 trials per row:

| Data size | Chance a genuinely good model passes | What actually failed |
|---|---|---|
| 4 proteins, **92 rows** | **51%** — a coin flip | error-bar band — never the ranking |
| 4 proteins, **200 rows** | **76%** | " |
| 6 proteins, **198 rows** | **77%** | " |
| 10 proteins, **200 rows** | **75%** | " |
| 10 proteins, **500 rows** | **94%** | " |
| 20 proteins, **900 rows** | **99%** | " |

The ranking part passed **100% of the time at every size.** Coverage was the entire story.

**Treat ~200 rows as the pragmatic floor, and ~500 as what actually buys you confidence.**
At 90 rows a genuinely good model is a coin flip away from being recorded as a failure —
and Step 19's table maps that failure to "RiboNN does not do BT4's job."

> **Why not the research doc's ~102 rows?** That figure counts only the *test* fold's
> sampling noise. Split conformal has a second, roughly equal source of uncertainty — the
> calibration fold's own estimate of the error-bar width — so ~102 is a lower bound, not a
> guarantee. The measured table above is the honest version.

**Two caveats on the table.** It assumes each protein contributes roughly the same number
of variants; with lopsided proteins the two folds diverge and the pass rate drops below
these numbers. And it assumes the model's error is a similar size across proteins — if it
isn't (the report tells you, as `link_slope_spread`), then adding *proteins* helps more
than adding variants. Maximise both where you can.

Other floors worth knowing: at least **9 rows** in the calibration half (below that the
error bars come back as infinite), at least **2 proteins**, and at least one protein with
2+ variants in the test half.

> **A quirk to plan around, and it is not cosmetic:** the split isn't random. It sorts your
> protein names **alphabetically** and puts the first half into calibration. So what you
> *name* your proteins decides which ones get tested — and the alphabetically-first half
> must hold at least 9 rows, ideally about half the panel. Fix your names before you look
> at any results.

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
                      point estimate
baselines:            permutation, cai, gc3, length, constant
width_rule:           median error-bar width ÷ spread of the data must be < 1.0
runs_permitted:       1
```

The gate now records the rest of the scope itself — `gate_result.json` carries a `scope`
block with the cell types, `top_k`, the species and a content hash per UTR context, and
`bt4 expression-gate` prints the same line. **Where the two disagree, `scope` is the
record**: it is what the run was configured with, and it is what an attestation is
derived from, so a note here that contradicts it is refused rather than believed.

One thing the record genuinely cannot hold: it stores a *hash* of each UTR context, not
a filename. **Write down which UTR files you used**, so that hash can be resolved back to
something you can find again.

## Step 18 — Run the test. Once.

> ⚠️ **Stop here if your panel has only one protein.** The gate refuses to run — you get
> `error: expression gate needs at least two distinct groups for a leakage-free split,
> got 1` and no output file at all. This is not a bug: with one protein there is no way to
> split calibration from test without the two halves leaking into each other. If that is
> the only data you have (PERSIST-seq is shaped this way), the honest move is to group by
> *design family* instead of by protein — a real but distinctly narrower claim, and one
> you should describe as such rather than promote.

```bash
python scripts/run_expression_gate.py \
    --panel panel.tsv --backend ribonn --species human --cell-type HEK293T \
    --within-group --recalibrate \
    --min-spearman 0.30 --target-coverage 0.90 --coverage-tolerance 0.05 \
    --baselines permutation,cai,gc3,length,constant \
    --num-workers 0 --json > gate_result.json
```

(Step 20 adds one flag to this same command, so you may prefer to jump ahead and run the
final form once rather than scoring twice.)

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

**`--cell-type`** — ⚠️ **it used to be the silent one.** Leave it out and RiboNN averages
**all 78 human cell types**; against a panel measured in one cell line that is a
different quantity entirely. **If your panel fills in its `cell_type` column, this is now
caught before anything is scored** — the gate refuses and tells you which `--cell-type`
to pass. If your panel leaves the column blank there is still nothing to check against
and it will run cleanly to a wrong verdict, so fill it in. (Same for `species`.)

> **Use the script, not the `bt4 expression-gate` shortcut**, if you want the record file.
> The shortcut command exists and works, but it has **no `--json` flag**, so it cannot
> produce `gate_result.json` — and it also can't set `--baselines` or `--attest`.
> (It *can* now set `--top-k`.)

Keep `gate_result.json`, `stage1_spread.json` and your pre-registration file together.
Those three files *are* the record — four, once Step 20 adds the attestation.

## Step 19 — Read the verdict

Open `gate_result.json` and look at `verdict.promotable`. It is true only when three
separate things all hold:

1. **`gate_passed`** — the ranking cleared 0.30 **and** the error bars landed within 90% ± 5%.
2. **`beats_every_baseline`** — RiboNN beat all five stand-ins. **The hard one is CAI**, a
   simple arithmetic formula BT4 already computes instantly, for free, *inside* the
   optimizer. A model that can't beat it adds nothing. Note the bar is deliberately
   demanding: RiboNN's *cautious lower bound* must exceed the best baseline's plain
   mid-range score (its point estimate) — not the other way round.
3. **`interval_is_informative`** — the error bars are narrower than the spread of the data
   itself. This catches a genuinely sneaky failure: a model that always predicts the same
   number gets its error bars technically "right" while saying nothing at all.

### What each outcome means

| Result | What it means | What to do |
|---|---|---|
| **`promotable: true`** | RiboNN really does rank spellings of a protein | Go to Step 20 |
| Ranking works, but doesn't beat **CAI** | No gain over what BT4 already does for free | Stays uncalibrated. Publish the numbers. |
| **Ranking fails** | RiboNN doesn't do BT4's job | Stays uncalibrated **for this regime**. Write it up — the question is *answered* for the frozen model, not pending. (Re-training RiboNN on such a panel is a separate, far larger project, not a re-run of this one.) |
| Ranking works, error bars don't | It orders things usefully but misstates its confidence | Report both halves. Default is still uncalibrated. |

**Note the asymmetry:** a pass proves something specific and scoped. A fail at 90 rows
might just be Step 16's coin flip. Check your row count before concluding "no".

---

# Part 6 — If it passed

## Step 20 — Record it properly

The only sanctioned way to flip the switch is an **attestation** — a small,
content-hashed record saying "this model passed this test, on this data (bound by the
panel's SHA-256), in this scope." It is not signed and proves nothing about *who* ran it;
it binds the claim to exact bytes, not to an identity.

**Add one flag to the Step 18 command.** That is the whole procedure now:

```bash
python scripts/run_expression_gate.py \
    --panel panel.tsv --backend ribonn --species human --cell-type HEK293T \
    --within-group --recalibrate \
    --min-spearman 0.30 --target-coverage 0.90 --coverage-tolerance 0.05 \
    --baselines permutation,cai,gc3,length,constant \
    --num-workers 0 \
    --readout mean_ribosome_load \
    --attest ribonn_attestation.json \
    --json > gate_result.json
```

`--attest` writes the record **from the comparison the verdict was read off**. That
matters for two reasons, and both used to be real hazards:

- It does not re-score. The old two-script procedure ran the whole test a second time —
  another full RiboNN pass — purely to rebuild an object it had already had.
- It cannot be configured differently the second time, because there is no second time.

If the run is not promotable, `--attest` refuses with exit code 3 and writes nothing.
That is a result, not a failure — go to the outcome table in Step 19.

> ⚠️ **This is not the sharpest edge in the procedure any more.** Earlier versions of
> this guide warned that the `species` and `cell_types` you typed into `promote.py` were
> free text, never checked against how the test actually ran — so a run averaging all 78
> cell types could be filed as a HEK293T result, and every later check would accept it.
> **That hole is closed.** The scope is now taken from the run itself. If you pass
> `species=` or `cell_types=` and they disagree with what was scored, you get a refusal,
> not a record. And where your panel's own `species` / `cell_type` / `readout` columns
> declare the same fact, that is checked against your file too — the attestation's
> `verified_against_panel` field lists exactly which parts of the scope got that second
> check, so a reader can tell a verified scope from one taken on your word.

The attestation layer still refuses outright to record: a failing test, a run done
without `--within-group`, a model that didn't beat every baseline, error bars as wide as
the data, or thresholds set below the built-in floors.

**If you prefer to do it in Python** (e.g. you already have the comparison in a session):

```python
# promote.py — only if you are not using --attest.
import json
from bt4 import api, __version__

panel = api.read_panel("panel.tsv")
comparison = api.expression_gate(
    panel, "ribonn",
    settings=api.GateSettings(within_group=True, recalibrate=True,
                              min_spearman=0.30, coverage_tolerance=0.05),
    species="human", cell_types=("HEK293T",),
)
assert comparison.promotable, "gate did not pass; nothing to attest"
att = api.attest_expression(
    comparison, readout="mean_ribosome_load", bt4_version=__version__,
)
open("ribonn_attestation.json", "w").write(json.dumps(att.to_dict(), indent=2))
```

Note there is no `species=`/`cell_types=` on `attest_expression` here: they come from
the `expression_gate` call above. You *may* pass them, and then they are checked — but
you cannot use them to say something the run did not do.

Keep `gate_result.json`, `ribonn_attestation.json`, `stage1_spread.json` and your
pre-registration file together. Those four files *are* the record.

Then update [`CLAUDE.md`](../CLAUDE.md) §6/§9, item 11 in
[`NEXT_SESSION.md`](NEXT_SESSION.md), and [`CHANGELOG.md`](../CHANGELOG.md).

## Step 21 — Using it (this part now works)

**Nothing auto-promotes, and nothing changes until you ask.** That is deliberate. But
asking is now one environment variable or one checkbox, where it used to be a Python
script you had to write yourself.

```bash
export BT4_EXPRESSION_ATTESTATION=/path/to/ribonn_attestation.json
export BT4_EXPRESSION_USE_ATTESTED=1
```

With those set, `api.resolve_expression_backend("ribonn", ...)` returns a **calibrated**
head. Hand that head to the design flow and the set becomes a real ranking:

```python
from bt4 import api
head = api.resolve_expression_backend(
    "ribonn", species="human", utr5="...", utr3="...", cell_types=("HEK293T",),
)
result = api.candidates("MKAYVQTL...", predictor=head)   # ranked, not discovery order
```

⚠️ **The environment variable governs promotion, not selection.** `api.candidates` with
no `predictor=` still uses the neutral placeholder however many variables you export —
exporting them does not silently change what a script that never asked for RiboNN does.
The run manifest records the attestation's content hash, so a design produced with a
promoted head is distinguishable from one produced without it.

In **BT4 Studio**: Candidates tab → *Expression head* → tick **honor expression
attestation**. The box is greyed out with an explanatory tooltip when no attestation
resolves, so it is never a control that silently does nothing. When it is on, the panel
shows the scope you are trusting — species, cell types, readout, `top_k`, how many UTR
contexts, and the panel hash — and the head is **pinned** to that species and cell-type
selection, so the form cannot display one scope while the run uses another.

You can still do it explicitly in code:

```python
from bt4 import api
model = api.resolve_expression_backend(
    "ribonn", species="human", utr5="...", utr3="...",
    cell_types=("HEK293T",), use_attested=True,
)
assert model.calibrated
```

**A configuration the attestation does not cover is refused, not downgraded.** Species,
cell-type selection, `top_k`, the UTR context, and the adapter's pinned weight hashes are
all bound: change any of them and you get an error naming the mismatch. You do *not* get
an uncalibrated head handed back quietly, which would be the worst outcome — a user who
asked for a calibrated ranking and cannot tell they did not get one. (`batch_size` and
`num_workers` are deliberately *not* bound: RiboNN pads to a fixed width and does not
shuffle when predicting, so neither can change a score.)

**The UTR context is part of the claim.** Your gate measured ranking *inside* the
transcript context your panel used. A head configured for a different 5′/3′ UTR is
outside what was measured, so it is refused. If you want to design in a different
context, that is a different question and needs its own panel.

**Should you commit the attestation to the repo?** Only if your panel is public. The
record carries the panel's SHA-256 and a hash of each UTR context — and a *short* UTR is
recoverable from its hash by brute force, so committing one publishes that context. The
`$BT4_EXPRESSION_ATTESTATION` path exists precisely so a maintainer with unpublished data
can use their own result without publishing anything.

**And keep the claim scoped.** An attestation earned on HEK293T does **not** certify a
model averaging all 78 cell types — which is why BT4 now refuses that rather than
trusting your label. The honest sentence is always something like: *"calibrated for
ranking synonymous variants of a known protein, in cell type X, measured by readout Y,
in UTR context Z."* The broader claim — ranking sequences for a protein nobody has
measured — **cannot be honestly tested at all**, because it would need ~100 held-out
proteins that don't exist.

---

# Appendix A — The traps, in one place

| # | Trap | Why it hurts |
|---|---|---|
| 1 | Two Python environments | BT4 loads RiboNN's code *into itself*. Separate environments can never see each other. |
| 1b | Reusing your Pangolin environment | Installing `[splice-pangolin]` there would upgrade torch off RiboNN's pin. Build RiboNN its own and install BT4 into both. Versions and the exact scope of the claim: Step 5. |
| 2 | Adding the `[ml]` extra | Pulls in newer PyTorch/NumPy than RiboNN tolerates. |
| 3 | Weights folder not named `models` | RiboNN hard-codes that name; the fallback needs admin rights on Windows. (On Linux/macOS `$BT4_RIBONN_WEIGHTS` lets any name work.) |
| 3b | `tar -xf weights.zip` on Linux | GNU `tar` cannot read zip. Use `unzip -q weights.zip -d models`. (Windows `tar` reads zip fine.) |
| 3c | Cloning to `C:\RiboNN` on Windows | Step 9 defaults to `~/RiboNN/models`, so it looks somewhere else and reports it cannot find `runs.csv`. Clone into `%USERPROFILE%` instead. |
| 4 | Trusting `('null', 'ribonn')` | Only proves two folders exist. Empty folders pass. |
| 5 | Forgetting `--within-group` | Measures the wrong thing entirely and wastes the run. Warns on stderr only. |
| 6 | Forgetting `--recalibrate` | Compares two different scales; error bars become garbage. |
| 7 | Forgetting `--cell-type` | **Now caught, if your panel says so.** A panel with a `cell_type` column makes the gate refuse before it scores. With no such column there is still nothing to check against, and it silently averages 78 cell types — so fill the column in. |
| 8 | Mismatched FASTA headers in Step 11 | Grouping is "text before the first `|`". Get it wrong and every sequence becomes its own group; the crash message never mentions headers. |
| 9 | Reading `responds_to_synonymous_change` as biology | It's a one-part-in-a-billion check. |
| 10 | A 90-row dataset | A good model fails about half the time. See Step 16. |
| 11 | Re-running with adjusted thresholds | Turns validation into a search. Nothing enforces the "run once" rule. |
| 12 | Typing a scope that doesn't match the run | **Fixed.** The scope is taken from the run; a disagreeing declaration is a refusal, not a record. Use `--attest` and there is nothing to type. |
| 13 | Expecting a green test suite to mean a real run works | Every RiboNN test uses a stand-in. |
| 14 | Comparing against raw protein output | Re-introduces the exact quantity RiboNN's output divides out. |
| 15 | Expecting a pass to change things for users *automatically* | It doesn't, by design — but the opt-in exists now (`BT4_EXPRESSION_USE_ATTESTED`, or the Studio checkbox), and a head outside the attested scope is refused rather than quietly downgraded. |
| 15b | An agent or script running the commands for you | `conda activate` **and** `BT4_RIBONN_DIR` are both shell state; a fresh shell per command loses both, and `conda run` alone only replaces the first. See Step 7. |
| 16 | Committing an attestation for a private panel | It carries the panel hash and a hash per UTR context; a *short* UTR is brute-forceable from its hash. Use `$BT4_EXPRESSION_ATTESTATION` locally instead. |

---

# Appendix B — Where this guide corrects the runbook

Every command in this guide was checked against the source, and several were executed.
Where [`DESIGN_ribonn_calibration.md`](DESIGN_ribonn_calibration.md) differs:

1. **"The adapter reproduces upstream bit-for-bit."** The runbook's summary table *used
   to* mark this ✅ done — **corrected in place** when this guide landed. **No RiboNN
   fidelity gate, capture script or attestation exists in the repo**; "bit-for-bit" is the
   *splice* side's vocabulary, and the expression attestation module explicitly contrasts
   itself with it. What has been established is that real end-to-end runs happen and
   produce numbers — which is why the runbook's own Stage 1.2 (fold semantics) exists at
   all. This guide says "it runs", not "it matches".
2. **`bt4 expression-gate` cannot produce `gate_result.json`.** The runbook offers it as
   an equivalent to the script at Stage 4, but it has **no `--json` flag** (nor
   `--baselines`, `--batch-size`, `--attest`). Use the script for the record. *(It
   silently ignored `--top-k` too — it had no such flag and always ensembled 5, while the
   script forwarded it. The flag now exists, because `top_k` is part of the scope an
   attestation binds.)*
3. **`harness_ok` does not mean "the scores differ substantially".** It is a
   not-equal test at a 1e-9 floor — about six orders of magnitude away from the runbook's
   wording.
4. **The ~90-row panel floor is too small** — measured above: it is a coin flip for a
   genuinely good model, and every failure is the coverage band rather than the ranking.
   The runbook's own research doc computes ~102 rows, but that counts only test-fold
   noise; measured end-to-end, ~200 is the pragmatic floor and ~500 is what buys ~94%.
5. **`$BT4_RIBONN_WEIGHTS` exists** and lets the weights live outside the checkout. Before
   this guide it appeared in no markdown file in the repo — only in source and two runtime
   error messages; it is now recorded here and in `CHANGELOG.md`.
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
11. **`tar -xf weights.zip` does not work on Linux.** The runbook gives it in a Windows
    `cmd.exe` block, where it is correct — Windows and macOS ship a `tar` that reads zip.
    GNU `tar`, the Linux default, does not, and Linux is the only platform RiboNN is
    tested on. Use `unzip`.
12. **The runbook's Stage 1 "free checks" are not free.** All four `--check` modes score
    through the backend, so every one needs the licensed weights. The genuinely free work
    — the test suite and a `--backend null` dry run — appears nowhere in the runbook. This
    guide leads with it (Steps 3–4).
13. **The scope you typed was never checked against the run** — this guide's own earlier
    version called it "the sharpest edge in the whole procedure", and it was. **Fixed in
    the code, not in prose (2026-08):** `attest_expression` now derives species / cell
    types / readout from the gate comparison, refuses a declaration that disagrees, and
    cross-checks against the panel's own columns where it has them. The gate additionally
    refuses a cell-type or species mismatch *before* it scores, so the run-once budget is
    not spent on a wrong-scope answer. Steps 18 and 20 are rewritten accordingly.
14. **Windows guidance was thinner than the Linux path, and one line of it was wrong.**
    Verified against the code while walking a maintainer through a real Windows install:
    Step 7's inline comment suggested `C:\RiboNN`, which puts the weights outside the
    `~/RiboNN/models` that Step 9's check defaults to — so the very check written to fail
    loudly would instead report that it could not find `runs.csv`. Steps 5–7 now carry
    `cmd.exe` forms, the `%USERPROFILE%` clone path, `set` vs `setx`, and the
    `conda run -n RiboNN --no-capture-output` form for anything running commands in a
    fresh shell each time.
15. **Miniforge was presented as a prerequisite.** It is not — Miniconda and Anaconda work
    too, and `mamba` is a speed convenience rather than a requirement. One difference does
    survive, and an earlier draft of this very item wrongly denied it by claiming "every
    command names its channel": `conda env create -f environment.yml` names none, so on
    Miniconda/Anaconda the `defaults` channel is consulted unless `nodefaults` is added to
    RiboNN's `environment.yml`. Step 5 has the accurate version.
16. **The Pangolin/RiboNN environment conflict was recorded nowhere.** A maintainer who
    already had a Pangolin environment would reasonably try to reuse it; Step 5 now says
    not to, why, and what the separation costs (no CNN splice audit in the same run).
    Version numbers and the scope of the claim live in Step 5 **only** — an earlier draft
    restated them here and in `NEXT_SESSION.md`, and the three copies had already drifted
    apart within a single commit, which is §10.11's hazard demonstrated rather than
    described.
17. **"Nothing calls `verified_predictor`"** — true when this guide landed, and the reason
    old Step 21 was titled "the part that still won't work". No longer: the promotion seam
    ships (`BT4_EXPRESSION_USE_ATTESTED`, `$BT4_EXPRESSION_ATTESTATION`, a BT4 Studio
    toggle), and `run_expression_gate.py --attest` writes the record from the same
    comparison the verdict came from, so the second scoring pass is gone. What has **not**
    changed: nothing is bundled, nothing auto-promotes, and a head outside the attested
    scope is refused rather than downgraded.

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
| **r²** | 0 to 1: how much of the real variation the model's predictions account for. ~0.6 is a useful model; ~0.1 is nearly none. |
| **GC3** | The G/C content of just the *third* letter of each codon — the letter you are usually free to change without changing the protein. |
| **CAI** | A simple arithmetic score for how closely a sequence uses an organism's preferred codons. BT4 computes it instantly, in-loop, for free. |
| **tAI** | A second codon-quality score, weighted by how much of each matching tRNA the cell actually has. |
| **`top_k`** | How many of RiboNN's ensemble members to average (default 5). It changes the score, so it is part of the scope: the gate records it, and a promotion refuses a head that ensembles a different number. |
| **Coverage** | Of everything the model said "I'm 90% sure", how often was it right? |
| **Conformal interval** | The error bar. Built from how wrong the model was on data it hadn't seen. |
| **Within-group** | Comparing only spellings of the *same* protein — BT4's actual job, and the strict version of the test. |
| **Pooled** | Comparing everything at once, including across different proteins — the easy version, and the wrong question. |
| **Attestation** | The content-hashed record saying a model passed, and in exactly what scope. Not a signature — it binds to bytes, not to a person. |
