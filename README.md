# BT4 — honest protein → mRNA back-translation

[![CI](https://github.com/masonberger4/BT4/actions/workflows/ci.yml/badge.svg)](https://github.com/masonberger4/BT4/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

**BT4 turns a protein into an optimized coding DNA / mRNA sequence** for a target
organism (default *Homo sapiens*) — choosing synonymous codons to maximize codon
adaptation while honoring hard biological constraints, and it is **honest** about
exactly how optimal the answer is and how it was produced.

It treats codon optimization as what it really is — *constrained combinatorial
optimization over a codon trellis*, solved **exactly** — not greedy per-codon
substitution. Every result recomputes its own metrics from the delivered DNA,
carries a machine-readable **optimality certificate**, and ships a content-hashed
**provenance manifest** so it is reproducible from its stamp.

![BT4 Studio — the desktop app](docs/images/bt4-studio.png)

*BT4 Studio: paste a protein, set your constraints, and get an exactly-optimized
coding sequence with an honest optimality badge and a CAI/GC trade-off frontier.*

---

## What it does today

- **Exact codon-trellis DP** over the *true* per-constraint context — no silent
  window cap, no greedy shortcut — with an explicit `beam` speed knob.
- **A real optimality certificate** on every result: `proven_optimal` when the
  full state space was explored, or honestly `beam_truncated` when you traded it
  for speed. BT4 never claims optimality it didn't earn.
- **Objectives:** CAI, **tAI** (real tRNA copy numbers via the dos Reis wobble
  model), GC-target proximity, a 5′ translation ramp, CpG deplete/elevate, a
  **%MinMax** codon-commonness term, and **codon-pair bias** (built from a
  user-supplied reference CDS — no default table is fabricated) — returned as a
  **multi-objective Pareto frontier** (a simplex sweep over every active axis, not
  just CAI/GC), never a single magic-weighted number.
- **Global GC budget, two honest backends:** an OR-Tools **CP-SAT** backend for
  the pure-additive case (proven-optimal), and a **Lagrangian relaxation** that
  dualizes the budget into the exact DP so — unlike CP-SAT — it keeps local
  constraints and pairwise terms honored, with a gap-bounded certificate.
- **CpG / UpA count budget:** cap or floor the whole-sequence CpG (`--cpg-max` /
  `--cpg-min`) or UpA (`--upa-*`) count — CpG depletion for stealth, elevation for
  immunogenicity. Enforced *exactly* by an amount-bucketed DP whose per-codon count
  is boundary-aware (a CpG can straddle two codons), with a proven-optimal
  certificate and every local constraint still honored.
- **Hard constraints:** maximum homopolymer run, **max GC-run** (the "max GC
  length"), a whole-sequence **max repeat length** (direct/inverted/palindromic
  repeats anywhere, reverse-complement aware — enforced by refinement and reported
  honestly, since it is genuinely non-local), forbidden motifs with named
  **forbidden-sequence presets** (poly-A signal, TATA box, telomere repeat, …),
  **tandem & inverted-repeat** (hairpin) bans, an **internal strong-Kozak ATG**
  guard, an **out-of-frame uORF** suppressor (a structural flag, refinement-
  enforced, not a calibrated expression claim), and a **restriction-enzyme
  catalog** (IUPAC-aware, auto reverse-complement).
- **Multiple organisms:** human, *E. coli*, and *S. cerevisiae* codon tables out
  of the box, plus real **tAI** tables for eight organisms (human, mouse, rat,
  zebrafish, *Drosophila*, *C. elegans*, *Arabidopsis*, and *S. cerevisiae*) from
  GtRNAdb tRNA counts; `bt4 build-table` builds an authentic codon table from your
  own CDS FASTA.
- **Benchmarked against real tools:** `scripts/compare_tools.py` places BT4 next
  to GeneOptimizer / IDT / Twist / GenScript on a cited, CC BY 4.0 panel — every
  metric recomputed from the sequence, and BT4 never claimed "better", just placed.
  `scripts/compare_reproducibility.py` adds a run-to-run **variability** view over
  three proteins × three *anonymized* algorithms × ten repeat runs (kept honestly
  separate from the named-tool board), with deterministic BT4 as a zero-spread
  reference.
- **Honest metrics:** every reported number (CAI, GC, violations) is recomputed
  from the delivered sequence, never trusted from the solver.
- **Reproducible provenance:** a content-hashed manifest (codon-table SHA-256 +
  config + seed + version) — a swapped table produces a different stamp.
- **Three ways to use it:** the **BT4 Studio** desktop app, a `bt4` **CLI**, and a
  stable **Python API** (`bt4.api`). FASTA / JSON export everywhere.
- **100% local and offline.** Nothing leaves your machine.
- **Optional Rust accelerator** (PyO3) for the hot-loop primitives, with an
  identical pure-Python fallback so nothing is required to run.

- **Optional 5′-folding refinement (`--refine`):** a simulated-annealing pass
  (synonymous-only, never adds a hard violation, honest `heuristic` certificate)
  that opens up start-proximal mRNA structure via a `FoldingModel`. With
  **ViennaRNA** installed (`bt4[fold]`) this is real ΔG; otherwise it uses a
  clearly-labeled uncalibrated proxy and says so — BT4 never presents the
  baseline as calibrated thermodynamics.

> **Honest about scope.** BT4's design aims wider still — SpliceAI/Pangolin-class
> splice models and a learned expression head. **Both the Pangolin and SpliceAI
> backends have now landed** as inference-only wrappers: each drives your own
> installed model (BT4 bundles neither code nor weights — Pangolin is GPL-3.0;
> SpliceAI's code is PolyForm Strict and its weights CC BY-NC 4.0, noncommercial)
> and hash-pins the weights it loads, with a two-backend agreement harness that
> turns running both into an uncertainty signal. But they ship
> **`calibrated=False`** until they pass an integration-fidelity gate, so BT4 still
> refuses to present them as validated results; the learned expression head is
> **not shipped yet**. See [`CLAUDE.md`](./CLAUDE.md) §9 for the full plan.

---

## Install BT4 Studio (no coding required)

> New to installing apps? See the **[full step-by-step install & troubleshooting guide](docs/INSTALL.md)**.

BT4 Studio is a normal desktop app: **you download one file and open it.** No
terminal, no Python, no setup. Everything runs on your own machine and nothing is
uploaded anywhere.

**Get the app:** open the
[**Releases page**](https://github.com/masonberger4/BT4/releases), open the latest
release, and under **Assets** download the one file for your computer:

| Your computer | Download this | Then |
|---|---|---|
| **Windows** | `BT4-Studio-Windows.exe` | Double-click it. |
| **Mac** | `BT4-Studio-macOS.dmg` | Double-click it, then drag **BT4 Studio** into **Applications**. Open it from Applications. |
| **Linux** | `BT4-Studio-Linux-x86_64` | Right-click → **Properties → Permissions → Allow executing as program**, then double-click. |

That's it — BT4 Studio opens, you paste a protein, pick an organism, and click
**Optimize**.

> **Mac note:** the `.dmg` is built for **Apple Silicon** Macs (M1/M2/M3/M4 — 2020
> or later; check **Apple menu → About This Mac**). On an older **Intel** Mac, use
> the developer install below for now.

<details>
<summary><b>"Windows protected your PC" / "macOS cannot verify the developer" — how to open it anyway</b></summary>

The app is safe, but it is **not code-signed**, so Windows and macOS show a
one-time warning for apps they don't recognize. This is expected. To open it:

- **Windows:** on the blue "Windows protected your PC" box, click **More info**,
  then **Run anyway**.
- **Mac:** if you see *"'BT4 Studio' cannot be opened because Apple cannot check
  it for malicious software"*, **right-click (or Control-click) the app → Open →
  Open**. You only need to do this the first time. (On recent macOS you can also
  go to **System Settings → Privacy & Security** and click **Open Anyway**.)
- **Linux:** if double-click does nothing, your file manager may not run
  executables directly — open a file manager, right-click the file, and enable
  "Allow executing as program" (wording varies by desktop).

These apps are intentionally unsigned -- code-signing would mean paying for Apple
and Windows certificates just to remove a warning you click through once.

</details>

> **Note:** packaged downloads appear on a release only after its build has run.
> If the latest release shows no `.exe` / `.dmg` / Linux file under **Assets**
> yet, it hasn't been built — a maintainer can produce them with the
> [packaging guide](packaging/README.md#repairing-a-release), or use the developer
> install below.

---

## Install for developers (from source)

BT4 is also a `bt4` CLI and a `bt4.api` Python library. It needs **Python 3.10+**.

```bash
git clone https://github.com/masonberger4/BT4
cd BT4
pip install -e '.[app]'     # core + PySide6 + pyqtgraph (drop [app] for just the CLI/API)
bt4-studio                  # launch the desktop app  (or:  python -m bt4.app)
bt4 --help                  # or use the command line
```

On **Linux**, Qt needs a few system libraries for the GUI to open. If launching
`bt4-studio` fails with `libEGL.so.1: cannot open shared object file`, install
them once (Debian/Ubuntu shown):

```bash
sudo apt-get install -y libegl1 libgl1 libglib2.0-0 libxkbcommon0 libdbus-1-3
```

macOS and Windows need no extra system packages. Prefer an isolated install that
puts both the app and CLI on your `PATH`? Use [`pipx`](https://pipx.pypa.io):

```bash
pipx install "bt4[app] @ git+https://github.com/masonberger4/BT4"
```

---

## Use it

### Desktop app

`bt4-studio` (or `python -m bt4.app`): paste a protein, pick the organism, set a
GC target / max-homopolymer / max GC length / max repeat length / forbidden
motifs (and tick any forbidden-sequence presets), and click **Optimize**. **Hover
any control for a tooltip explaining what it does.** You get the optimality badge,
a recomputed-metrics table, the interactive CAI/GC frontier (the delivered point
starred), the coding sequence, and one-click FASTA/JSON export. The optimization
runs on a background thread, so the window never blocks.

### Command line

```bash
bt4 optimize MAALKHETQW --max-homopolymer 5 --enzyme EcoRI    # summary
bt4 optimize MAALKHETQW --max-gc-run 5 --max-repeat-length 10 # GC-run + repeat caps
bt4 optimize MAALKHETQW --forbid-preset poly_a_signal         # ban a preset's motifs
bt4 optimize MAALKHETQW --avoid-uorf                          # suppress out-of-frame uORFs
bt4 optimize MAALKHETQW --cpb-weight 1 --cpb-cds ref.fasta    # codon-pair bias (your CDS)
bt4 optimize MAALKHETQW --fasta                               # FASTA to stdout
bt4 optimize MAALKHETQW --json                                # JSON + manifest
bt4 library MAALKHETQW --n 20 --temperature 1.0               # sample a library (SAMPLED, not optimized)
bt4 validate ATGGCC...TAA --max-homopolymer 6                 # audit a sequence
bt4 tracks ATGGCC...TAA --nt-window 50                        # per-site GC/CpG/%MinMax tracks
bt4 organisms   # codon tables    bt4 enzymes   # enzymes    bt4 presets   # forbidden presets
bt4 build-table my_cds.fasta --organism my_species --out .    # table from real CDS
```

An optional headless HTTP API is available too (`pip install -e '.[service]'`,
then `uvicorn bt4.service.api:app`) exposing `/optimize`, `/frontier`,
`/validate`, `/organisms`, and `/health`.

### Python API

```python
from bt4 import api

result = api.optimize("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5))
print(result.dna, result.certificate.status.value)   # ...TAA  proven_optimal
print(result.audit["cai"], result.metrics.gc)        # recomputed from the DNA

frontier = api.frontier("MAALKHETQW", steps=11)       # CAI vs GC Pareto frontier
report = api.validate(result.dna, api.OptimizeConfig(max_homopolymer=5))

# Library / degenerate-design mode: sample a library from the codon distribution.
# This is a stochastic sampler, not an optimizer -- every member is SAMPLED, not
# proven optimal, and makes no expression claim.
library = api.library("MAALKHETQW", n=20, seed=1, temperature=1.0)
print(library.distinct, library.results[0].certificate.status.value)  # ...  sampled
```

---

## How it works (60 seconds)

The core is a **codon trellis**: one layer per residue, each holding that amino
acid's synonymous codons. A dynamic program walks the trellis, carrying exactly
the sequence context each constraint declares it needs, and finds the assignment
that maximizes the (weighted) additive objective while every hard constraint's
`ok_suffix` veto holds. Because the state carries the *true* per-constraint
context, the result is **exactly optimal** for that objective — and the certificate
says so. Sweeping the objective weights over the unit simplex — CAI vs GC, and
any ramp/CpG/%MinMax axes you activate — traces the **Pareto frontier**; the
delivered point is the trade-off you chose.

The whole system is a **strict, acyclic layering** (`domain` → pure layers →
`pipeline` → `api` → `cli`/`app`/`service`) enforced in CI by import-linter, and
a set of **honesty invariants** (round-trip, reported == computed,
`ok_suffix ⇔ validate`, `delta == score`, certificate honesty, determinism) that
are **property-tested**, not just documented. The full design rationale lives in
[`CLAUDE.md`](./CLAUDE.md).

---

## Develop

```bash
pip install -e '.[dev]'   # or '.[dev,app]' to include the desktop app

ruff check src tests      # lint
mypy                      # strict type-check
lint-imports              # enforce the layering contract
pytest                    # property + integration + determinism tests
```

Optional Rust accelerator (a pure-Python fallback is used when it isn't built):

```bash
pip install maturin
maturin develop -m rust/bt4_core/Cargo.toml   # builds the `bt4_native` extension
```

New to the codebase? Read [`CLAUDE.md`](./CLAUDE.md) first (the design
constitution), then [`CONTRIBUTING.md`](./CONTRIBUTING.md). Adding a constraint or
objective is a new file plus its honesty property test — never an engine edit.

## Roadmap

Phases 0–2 are complete (and Phase 3 groundwork has landed; Phase 5 is opened):
the multi-objective
frontier, the desktop app, **tAI** (real GtRNAdb tRNA data), 5′-ramp / CpG /
%MinMax / **codon-pair-bias** objectives, tandem & inverted-repeat, internal-ATG,
max-GC-run, max-repeat-length and **out-of-frame uORF** constraints, forbidden-
sequence presets, a cited tool benchmark, and two GC-budget backends — OR-Tools
CP-SAT and an honest **exact budget DP** (which, unlike CP-SAT, keeps local
constraints and pairwise terms under the budget) — have landed. Phase 3 groundwork
is in (the `FoldingModel` and `SplicePredictor` contracts with honest baselines,
the SA refinement engine, plotted per-site tracks), **both wrapped published
splice backends** (Pangolin + SpliceAI, inference-only, hash-pinned,
`calibrated=False` until their fidelity gates) with a two-backend agreement
harness, and the **`ExpressionPredictor` contract is scaffolded** for Phase 4 (a
neutral, honestly-uncalibrated placeholder until a validated head passes its
gate). Recording the splice fidelity gates, the validated expression model, and a
Rust trellis port are next. See [`CLAUDE.md`](./CLAUDE.md) §9.

## How BT4 compares

Wondering how BT4 stacks up against IDT, Twist, GeneArt, ATUM, or DNA Chisel?
See **[`docs/COMPARISON.md`](docs/COMPARISON.md)** for an honest, sourced
positioning — including where BT4 genuinely differs (exact multi-objective
optimization, optimality certificates, reproducible provenance, validated ML with
honest calibration) **and where it does not** (empirical expression grounding,
synthesis-manufacturability), plus the peer-reviewed caveat that the CDS is only a
minority of the expression signal.

## Contributing

Contributions are welcome — see [`CONTRIBUTING.md`](./CONTRIBUTING.md) and the
[`Code of Conduct`](./CODE_OF_CONDUCT.md). Please keep the layering contract and
the honesty invariants intact (CI enforces both).

## License

[MIT](./LICENSE).
