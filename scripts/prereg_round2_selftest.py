#!/usr/bin/env python
"""Ground-truth self-test for the round-2 gate. Run this before trusting a real result.

A gate is only worth pre-registering if it can be shown to fail when it should. This
script feeds `prereg_round2_analyze.analyze` three synthetic regimes whose right
answers are known by construction, and asserts the gate gets all three right:

===========  ====================================================  ============
regime       what the synthetic score is                           gate must
===========  ====================================================  ============
A            a pure blend of the free features (3*GC + 0.5*CAI)    **not pass**
B            that blend plus a per-variant CDS term, identical     **pass**
             in every UTR context
C            that blend plus fresh per-context noise               **not pass**
===========  ====================================================  ============

B is the only regime with information that is both non-free and a property of the
coding sequence. A is exactly the hypothesis "RiboNN adds nothing over what BT4
already computes"; C is "RiboNN responds, but the response is a CDS x UTR
interaction that no fixed-UTR panel could calibrate".

**This test earned its place by failing.** Two earlier versions of the gate passed
regime A -- i.e. would have green-lit a five-figure panel for a model that is a
GC/CAI blend:

1. *Stability alone.* Once the fit removes the free features, a pure blend leaves
   only numerical dust -- but that dust is deterministic, so it is **identical** in
   every UTR context and its cross-context correlation is exactly 1.0.
2. *Stability x size, against single-feature baselines.* The fit runs in rank space,
   and a monotone transform of a *blend* is not linear in the ranked features, so a
   blend leaves a small stable residual that a ``score = GC`` baseline does not.
   The baseline has to be "any blend", which is why the analysis compares against a
   family of seeded random blends and takes the worst case.

Neither flaw is visible by reading the definition; both are obvious the moment a
known-negative regime is scored. Re-run this after **any** change to the analysis.

Usage::

    python scripts/prereg_round2_selftest.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from prereg_round2_analyze import analyze  # noqa: E402

N_PROTEINS = 16
N_VARIANTS = 40
CONTEXTS = ("c1", "c2", "c3", "c4")
SEED = 7


def synth(regime: str, rng: np.random.Generator) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for p_i in range(N_PROTEINS):
        protein = f"P{p_i:02d}"
        gc = rng.uniform(0.40, 0.62, N_VARIANTS)
        gc3 = gc + rng.normal(0, 0.01, N_VARIANTS)
        cai = rng.uniform(0.60, 0.95, N_VARIANTS)
        tai = rng.uniform(0.30, 0.70, N_VARIANTS)
        cpg = rng.uniform(0.01, 0.05, N_VARIANTS)
        length_nt = float(rng.integers(400, 2000))
        cds_term = rng.normal(0, 1.0, N_VARIANTS)  # a stable property of each variant
        for c_i, ctx in enumerate(CONTEXTS):
            base = 3.0 * gc + 0.5 * cai + c_i * 0.1  # a free-feature blend
            if regime == "A":
                score = base
            elif regime == "B":
                score = base + 0.6 * cds_term
            elif regime == "C":
                score = base + 0.6 * rng.normal(0, 1.0, N_VARIANTS)
            else:  # pragma: no cover - programmer error
                raise ValueError(regime)
            for v_i in range(N_VARIANTS):
                rows.append(
                    {
                        "protein": protein,
                        "variant_index": v_i,
                        "utr_context": ctx,
                        "score": float(score[v_i]),
                        "gc": float(gc[v_i]),
                        "gc3": float(gc3[v_i]),
                        "cai": float(cai[v_i]),
                        "tai": float(tai[v_i]),
                        "cpg": float(cpg[v_i]),
                        "length_nt": length_nt,
                    }
                )
    return rows


def main() -> int:
    rng = np.random.default_rng(SEED)
    expected = {"A": False, "B": True, "C": False}
    failures: list[str] = []
    for regime, should_pass in expected.items():
        report = analyze(synth(regime, rng))
        gate = report["gate"]
        got = bool(gate["passes"])
        ok = got == should_pass
        print(
            f"regime {regime}: gate {gate['median']:+.4f} "
            f"CI [{gate['ci95'][0]:+.4f}, {gate['ci95'][1]:+.4f}]  "
            f"passes={got}  expected={should_pass}  {'OK' if ok else 'WRONG'}"
        )
        if not ok:
            failures.append(f"regime {regime}: expected passes={should_pass}, got {got}")

    print()
    if failures:
        print("SELF-TEST FAILED -- do not trust a real result from this analysis:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("SELF-TEST PASSED: the gate fails both known-negative regimes and passes the positive.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
