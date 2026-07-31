# Contributing to BT4

Thanks for your interest in BT4! BT4 back-translates a protein into a coding
DNA / mRNA sequence optimized for a target organism, subject to biological
constraints — and it aims to be *honest* about how optimal each answer is and
how it was derived.

Before you write any code, please read **[`CLAUDE.md`](./CLAUDE.md)**. It is the
design constitution for this repository: the architecture, the contracts, the
honesty invariants, and the roadmap all live there, and contributions are
expected to respect it.

> A note on scope: BT4 today does exact-DP CAI optimization with GC-target
> proximity, hard constraints (max-homopolymer and forbidden motifs including
> reverse complements), a CAI/GC Pareto frontier, an optimality certificate, and
> content-hashed provenance manifests. Richer objectives (tAI, codon-pair, 5′
> ramp), ILP/relaxation backends, and validated splice/folding/expression models
> are on the roadmap in `CLAUDE.md` — they are **not** shipped yet. Please don't
> describe roadmap items as if they already work.

## Development setup

BT4 targets **Python 3.11+**. Install the package in editable mode with the dev
extra:

```bash
pip install -e '.[dev]'          # pure-Python core + dev tooling
```

The core is pure Python and stdlib-light; heavy dependencies live behind
optional extras (`ilp`, `fold`, `ml`, `app`, `service`, `assp`) and are lazily
imported behind contracts, so `import bt4` stays lightweight.

## Required checks

Every pull request must pass these four checks — CI runs them and merges are
blocked on failure, so please run them locally first:

```bash
ruff check src tests      # lint
mypy                      # strict type-check
lint-imports              # enforce the layering contract (import-linter)
pytest                    # tests (property + determinism)
```

## Architecture and layering

BT4 uses a **strict acyclic layering**, enforced by `import-linter` (the
`lint-imports` check above), not by good intentions:

```
domain  ->  biomodels | objectives | constraints | optimize | io  ->  pipeline  ->  api  ->  cli | app | service
```

- `domain` imports **nothing** from `bt4` (and none of the heavy deps).
- The pure layers (`biomodels`, `objectives`, `constraints`, `optimize`, `io`)
  import only `domain`.
- `pipeline` composes them; `api` composes `pipeline`; `cli`, `app`, and
  `service` import only `api`.
- No private symbol (`_foo`) crosses a layer boundary — registries and exports
  are public.

If your change needs to cross a layer boundary the wrong way, that's a signal to
rethink the design, not to loosen the contract.

## Honesty invariants (please don't break these)

BT4's core value is honesty, and the load-bearing invariants from `CLAUDE.md` §5
are **property-tested** — they are enforced in CI, not just documented:

- **Round-trip:** `translate(result.dna) == protein (+ stop)`.
- **Reported == computed:** every metric in a result is recomputed from the
  result's DNA by the owning model, never trusted from an accumulator.
- **`ok_suffix` ⇔ `validate` agreement:** every hard `Constraint`'s streaming
  `ok_suffix` veto must agree with its whole-sequence `validate`, and its
  declared `context_len` must actually be sufficient for that veto.
- **`delta == score`:** an `ObjectiveTerm`'s per-codon `delta`s must sum to its
  whole-sequence `score` (no non-additive term masquerading as additive).
- **Certificate honesty:** the optimality certificate must not overclaim — a
  `proven_optimal` result must agree with an independent exact solve, and any
  cap/prune/relaxation must be reflected.

## Adding a constraint or objective

Adding biology to BT4 is a **new file + an export + its property test** — never
an engine edit. Concretely:

1. Add your concrete class:
   - a `Constraint` (see `src/bt4/constraints/`), or
   - an `ObjectiveTerm` (see `src/bt4/objectives/`).
   Implement the full contract from `CLAUDE.md` §4 (declare `scope()` and
   `context_len()` honestly; implement `delta`/`score` for terms, and
   `ok_suffix`/`penalty`/`validate` for constraints).
2. **Export it publicly** — add it to the module's `__all__` and to the package
   `__init__.py` `__all__`, so it is a public, discoverable symbol (no private
   cross-layer wiring).
3. **Add the corresponding property test.** A PR that adds a constraint MUST add
   the `ok_suffix ⇔ validate` / sufficient-`context_len` test for it; a PR that
   adds an objective MUST add the `delta == score` test for it. See the existing
   tests under `tests/` (e.g. `test_constraints.py`, `test_objective.py`) for the
   pattern.

The solver, pipeline, and API should not need to change to gain a new
constraint or objective — if they do, please raise it in the PR so we can keep
the contracts doing the work.

## Commit and PR flow

- Branch off `main`, and open your pull request **into `main`**. BT4 is
  single-trunk — please don't strand work on long-lived divergent branches.
- Keep PRs **focused**: one logical change per PR is much easier to review.
- Make sure the four required checks are green locally before you push.
- If your change affects the architecture or contracts, update `CLAUDE.md` (and
  the `README.md` if user-facing) in the same PR — docs are kept in sync, not
  left to drift.
- The repository ships a pull-request template; please fill it in.

## Running the desktop app

BT4 Studio is the native PySide6 desktop app. Install the `app` extra and launch
it:

```bash
pip install -e '.[app]'
bt4-studio            # or:  python -m bt4.app
```

The app calls only `bt4.api` (on a background thread) and runs entirely offline.
Its smoke tests run headless:

```bash
QT_QPA_PLATFORM=offscreen pytest tests/test_app_smoke.py
```

## Building the optional Rust accelerator

The trellis hot loop has an optional Rust/PyO3 implementation in
`rust/bt4_core`. It is **not required** — an identical pure-Python fallback is
used when the extension isn't built — but you can build it to run the accelerated
path:

```bash
pip install maturin
maturin develop -m rust/bt4_core/Cargo.toml   # builds the `bt4_native` extension
```

## Code of Conduct

This project follows a [Code of Conduct](./CODE_OF_CONDUCT.md). By participating,
you are expected to uphold it.

## Reporting security issues

Please report vulnerabilities privately — see [`SECURITY.md`](./SECURITY.md).
Don't open a public issue for a security problem.

---

Thanks again for contributing. If you're unsure whether an idea fits, open an
issue and let's talk it through first.
