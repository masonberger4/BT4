# BT4

**Constrained, multi-objective protein → mRNA back-translation.** BT4 turns a
protein into a coding DNA / mRNA sequence optimized for real expression outcome
in a target organism, subject to biological constraints — and it is *honest*
about how optimal the answer is and how it was derived.

BT4 is a from-scratch successor to BT3. The design, contracts, invariants, and
roadmap live in [`CLAUDE.md`](./CLAUDE.md) — read it first.

> **Status: Phase 0 (foundations).** The pure `domain` layer (genetic code,
> sequence validation, multi-objective vector + Pareto frontier, optimality
> certificate, result types), content-addressed provenance, the optional
> Rust/PyO3 accelerator with a pure-Python fallback, packaging, the layering
> contract, and CI are being laid down. The optimizer, constraints, ML models,
> and the BT4 Studio desktop app come in later phases.

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
