# BT4

**Constrained, multi-objective protein → mRNA back-translation.** BT4 turns a
protein into a coding DNA / mRNA sequence optimized for real expression outcome
in a target organism, subject to biological constraints — and it is *honest*
about how optimal the answer is and how it was derived.

BT4 is a from-scratch successor to BT3. The design, contracts, invariants, and
roadmap live in [`CLAUDE.md`](./CLAUDE.md) — read it first.

> **Status: honest exact core + first multi-objective app.** On top of the
> Phase 0 foundations (pure `domain` layer, content-addressed provenance, the
> optional Rust/PyO3 accelerator, layering contract, CI) BT4 now has a working,
> honest vertical slice: an **exact codon-trellis DP** that optimizes CAI under
> local hard constraints (max-homopolymer, forbidden motifs incl. reverse
> complements) with true per-constraint context and a real **optimality
> certificate**; a **CAI/GC Pareto frontier**; a stable `bt4.api`; a `bt4` CLI;
> and **BT4 Studio**, the PySide6 desktop app (§6.6). Richer objectives (tAI,
> codon-pair, ramp), ILP/relaxation backends, and the validated splice/folding/
> expression models follow per the [`CLAUDE.md`](./CLAUDE.md) roadmap.

## What makes BT4 different (vs BT3)

- **Multi-objective, not a magic scalar.** CAI is one weak prior among many
  (tAI, codon-pair bias, 5′ ramp, folding ΔG, splice Δrisk, CpG budget, learned
  expression). BT4 returns a **Pareto frontier** and tells you where your answer
  sits on it.
- **Auditable optimality.** Exact windowed DP with the *true* per-constraint
  context (no silent cap); ILP/CP-SAT or Lagrangian relaxation with a real gap
  bound when needed. Every result carries an **optimality certificate**.
- **Real, validated ML — or it says so.** SpliceAI/Pangolin-class Δsplicing,
  ViennaRNA folding, an optional learned expression head — gated on real data, or
  loudly labeled a baseline. An opt-in **ASSP** cross-check is available for
  validation (never in the optimization loop).
- **Reproducible from the stamp.** Provenance hashes actual table/model
  *contents* + config + seed + git commit, so nothing lies about where a number
  came from.

## Install (development)

```bash
pip install -e '.[dev]'          # pure-Python core + dev tooling
```

Optional extras: `ilp` (OR-Tools), `fold` (ViennaRNA), `ml` (torch),
`app` (PySide6 desktop app), `service` (FastAPI HTTP API), `assp` (online check).

## Usage

**Python API** (stable, print-free — every frontend goes through it):

```python
from bt4 import api

result = api.optimize("MAALKHETQW", api.OptimizeConfig(max_homopolymer=5))
print(result.dna, result.certificate.status.value)   # e.g. ...TAA proven_optimal
print(result.audit["cai"], result.metrics.gc)         # recomputed from the DNA

frontier = api.frontier("MAALKHETQW", steps=11)        # CAI vs GC Pareto frontier
report = api.validate(result.dna, api.OptimizeConfig(max_homopolymer=5))
```

**CLI:**

```bash
bt4 optimize MAALKHETQW --max-homopolymer 5 --forbid GAATTC   # summary
bt4 optimize MAALKHETQW --fasta                               # FASTA out
bt4 optimize MAALKHETQW --json                                # JSON + manifest
bt4 validate ATGGCC...TAA --max-homopolymer 6                 # audit a sequence
bt4 organisms                                                 # list codon tables
```

**BT4 Studio** (native desktop app): paste a protein, pick constraints and a GC
target, and watch the optimality-certificate badge, metrics, sequence, and the
interactive Pareto frontier. It runs each solve on a background thread and calls
only `bt4.api`; nothing leaves the machine.

```bash
pip install -e '.[app]'
bt4-studio            # or:  python -m bt4.app
```

### Optional Rust accelerator

The hot-loop primitives have a Rust/PyO3 implementation in `rust/bt4_core`.
It's optional — without it, identical pure-Python fallbacks are used.

```bash
pip install maturin
maturin develop -m rust/bt4_core/Cargo.toml   # builds the `bt4_native` extension
```

## Develop

```bash
ruff check src tests      # lint
mypy                      # strict type-check
lint-imports              # enforce the layering contract
pytest                    # tests (property + determinism)
```

## License

MIT — see [`LICENSE`](./LICENSE).
