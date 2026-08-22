#!/usr/bin/env python
"""Compute the round-2 gate and magnitudes from an already-scored panel.

Implements the *analysis* half of `docs/PREREG_ribonn_part3_round2.md`, and is frozen
with it: re-analysing these scores with different code is a new pre-registration, not a
correction. It is deliberately separate from the scoring script so the endpoints can be
recomputed without re-spending the CPU, and so the scoring code never sees them.

The chain, in order:

1. **Free-feature fit** (§3.1). Within each protein x UTR context, regress RiboNN's
   score on the six free features *rank-transformed, plus their squares*. The squares
   matter: a rank-linear fit removes only monotone dependence, so a non-monotone
   response to GC would survive as a residual that looks like new information while
   still being GC.
2. **Void check** (§3.1). If the residual is still correlated with GC, the fit failed
   and the run is void -- neither pass nor fail.
3. **The gate** (§3.2). Per protein, combine the residual's **size** and its
   **stability across UTR contexts** --
   ``stable_non_free = residual_fraction x max(0, cross-UTR rho)`` -- and require it to
   beat the worst case over a family of free-feature baseline scorers (pure GC, pure
   CAI, and 24 seeded random blends) pushed through this same pipeline.

   Size and stability are both required, and that is not a matter of taste:
   ``scripts/prereg_round2_selftest.py`` scores three known-answer regimes and shows
   stability alone passes a pure GC/CAI blend (its post-fit residual is deterministic
   dust, hence identical in every context, hence correlated at 1.0), while size alone
   passes per-context noise. The blend family is required for the same reason -- a
   ``score = GC`` baseline is too weak, because the fit is in rank space and a monotone
   transform of a *blend* is not linear in the ranked features.
4. **Magnitudes** (§3.3), reported with CIs for the human spending judgement, never
   gating on their own.

Uncertainty everywhere is a cluster bootstrap resampling **whole proteins**, because
one protein's variants are a dependent cluster and resampling rows would understate the
interval badly.

Every score consumed here is an uncalibrated model output in arbitrary CLR-residual
units. Nothing this script prints can promote a backend.

Usage::

    python scripts/prereg_round2_analyze.py --scores round2_scores.jsonl --json
"""

from __future__ import annotations

import argparse
import json
import statistics
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from bt4.biomodels._stats import spearman

FEATURES = ("gc", "gc3", "cai", "tai", "cpg", "length_nt")
VOID_RESIDUAL_GC_MAX = 0.10
SANITY_WITHIN_OVER_BETWEEN_MIN = 0.20
BOOTSTRAP_DRAWS = 10_000
BOOTSTRAP_SEED = 20260821
N_BLEND_BASELINES = 24


def _ranks(values: np.ndarray) -> np.ndarray:
    """Average-rank transform (ties shared), then centred."""
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(1, len(values) + 1, dtype=float)
    # average ties
    uniq, inverse = np.unique(values, return_inverse=True)
    if len(uniq) != len(values):
        sums = np.zeros(len(uniq))
        counts = np.zeros(len(uniq))
        np.add.at(sums, inverse, ranks)
        np.add.at(counts, inverse, 1.0)
        ranks = (sums / counts)[inverse]
    return ranks - ranks.mean()


def fit_residuals(scores: np.ndarray, feats: dict[str, np.ndarray]) -> tuple[np.ndarray, float]:
    """Return (residuals, adjusted R^2) for score ~ ranked features + squares."""
    cols = []
    for name in FEATURES:
        r = _ranks(feats[name])
        if np.allclose(r, 0.0):  # constant within this protein (e.g. length)
            continue
        cols.append(r)
        cols.append(r**2 - (r**2).mean())
    if not cols:
        y = scores - scores.mean()
        return y, 0.0
    x = np.column_stack([np.ones(len(scores)), *cols])
    y = _ranks(scores)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    ss_res = float(resid @ resid)
    ss_tot = float(y @ y)
    n, p = len(y), x.shape[1] - 1
    if ss_tot <= 0 or n - p - 1 <= 0:
        return resid, 0.0
    r2 = 1.0 - ss_res / ss_tot
    r2_adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)
    return resid, float(r2_adj)


def gc_stratified_retention(scores: np.ndarray, gc: np.ndarray, bins: int = 5) -> float:
    """Pooled within-GC-bin SD / overall SD.

    ~0 means a pure GC detector (scores collapse once GC is held fixed); ~1 means the
    response is GC-independent. This is a **magnitude**, not a gate: §3.3 records that
    an earlier draft had its direction backwards.
    """
    overall = float(np.std(scores))
    if overall <= 0:
        return 0.0
    edges = np.quantile(gc, np.linspace(0, 1, bins + 1))
    idx = np.clip(np.searchsorted(edges, gc, side="right") - 1, 0, bins - 1)
    within = []
    for b in range(bins):
        sel = scores[idx == b]
        if len(sel) >= 2:
            within.append(float(np.var(sel)) * len(sel))
    if not within:
        return 0.0
    pooled = float(np.sqrt(sum(within) / sum(1 for _ in scores)))
    return pooled / overall


def cluster_bootstrap_ci(
    per_protein: dict[str, float], draws: int = BOOTSTRAP_DRAWS
) -> tuple[float, float, float]:
    """Median across proteins, with a percentile CI from resampling whole proteins."""
    keys = sorted(per_protein)
    vals = np.array([per_protein[k] for k in keys], dtype=float)
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    stats = np.empty(draws)
    for i in range(draws):
        stats[i] = np.median(rng.choice(vals, size=len(vals), replace=True))
    return float(np.median(vals)), float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))


def _med(d: dict[str, list[float]]) -> dict[str, float]:
    return {k: statistics.median(v) for k, v in d.items() if v}


def _cross_utr_rho(
    proteins: list[str],
    contexts: list[str],
    residuals: dict[tuple[str, str], np.ndarray],
) -> dict[str, float]:
    """Median pairwise Spearman of one protein's residuals across UTR contexts."""
    out: dict[str, float] = {}
    for protein in proteins:
        pairs = []
        for a, b in combinations(contexts, 2):
            ra, rb = residuals.get((protein, a)), residuals.get((protein, b))
            if ra is None or rb is None or len(ra) != len(rb):
                continue
            pairs.append(spearman(list(ra), list(rb)))
        if pairs:
            out[protein] = statistics.median(pairs)
    return out


def _stable_non_free(
    resid_fraction: dict[str, float], cross_utr: dict[str, float]
) -> dict[str, float]:
    """residual_fraction x max(0, cross-UTR rho) -- size *and* stability, together."""
    return {
        p: resid_fraction[p] * max(0.0, cross_utr.get(p, 0.0))
        for p in resid_fraction
        if p in cross_utr
    }


def analyze(rows: list[dict[str, Any]]) -> dict[str, Any]:
    proteins = sorted({r["protein"] for r in rows})
    contexts = sorted({r["utr_context"] for r in rows})

    # index[(protein, context)] -> rows in stable variant order
    index: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for r in rows:
        index.setdefault((r["protein"], r["utr_context"]), []).append(r)
    for key in index:
        index[key].sort(key=lambda r: r["variant_index"])

    residuals: dict[tuple[str, str], np.ndarray] = {}
    resid_fraction: dict[str, list[float]] = {}
    resid_gc: dict[str, list[float]] = {}
    retention: dict[str, list[float]] = {}
    within_sd: dict[str, list[float]] = {}
    group_means: list[float] = []

    primary = contexts[0]
    for protein in proteins:
        for ctx in contexts:
            recs = index.get((protein, ctx))
            if not recs or len(recs) < 12:
                continue
            scores = np.array([r["score"] for r in recs], dtype=float)
            feats = {f: np.array([r[f] for r in recs], dtype=float) for f in FEATURES}
            resid, r2_adj = fit_residuals(scores, feats)
            residuals[(protein, ctx)] = resid
            resid_fraction.setdefault(protein, []).append(1.0 - r2_adj)
            resid_gc.setdefault(protein, []).append(abs(spearman(list(resid), list(feats["gc"]))))
            retention.setdefault(protein, []).append(gc_stratified_retention(scores, feats["gc"]))
            if ctx == primary:
                within_sd.setdefault(protein, []).append(float(np.std(scores)))
                group_means.append(float(np.mean(scores)))

    # --- the gate: STABLE NON-FREE FRACTION, vs free-feature baselines ---------
    #
    # Stability alone is not enough, and a synthetic ground-truth check proved it: a
    # pure GC detector leaves only numerical dust after the fit, but that dust is
    # deterministic and identical in every context, so its cross-context correlation
    # is exactly 1.0 and a stability-only gate passes the very model we are trying to
    # rule out. Magnitude alone is not enough either -- per-context noise produces a
    # large residual that means nothing. The gate is therefore their product,
    #
    #     stable_non_free = residual_fraction x max(0, cross-UTR residual rho)
    #
    # compared against the same statistic computed for pure-GC and pure-CAI scorers
    # pushed through this identical pipeline. Those baselines are the honest null: if
    # RiboNN is a GC detector it scores like one, and the paired difference is ~0.
    cross_utr = _cross_utr_rho(proteins, contexts, residuals)
    stable_ribonn = _stable_non_free(_med(resid_fraction), cross_utr)

    # The null is "RiboNN is some blend of the free features", not "RiboNN is GC".
    # Single-feature baselines are too weak: the fit is in rank space, so a monotone
    # transform of a *blend* is not linear in the ranked features and leaves a small
    # residual which -- being deterministic -- is perfectly stable across contexts and
    # therefore reads as signal. A synthetic control proved this passes a pure
    # free-feature model. So the floor is the **worst case over many blends**, which is
    # the misspecification this pipeline cannot distinguish from information.
    rng_b = np.random.default_rng(BOOTSTRAP_SEED)
    blend_weights = [rng_b.normal(size=len(FEATURES)) for _ in range(N_BLEND_BASELINES)]
    baseline_specs: list[tuple[str, Any]] = [("gc", "gc"), ("cai", "cai")]
    baseline_specs += [(f"blend{i}", w) for i, w in enumerate(blend_weights)]

    baselines: dict[str, dict[str, float]] = {}
    for bname, spec in baseline_specs:
        b_resid: dict[tuple[str, str], np.ndarray] = {}
        b_frac: dict[str, list[float]] = {}
        for protein in proteins:
            for ctx in contexts:
                recs = index.get((protein, ctx))
                if not recs or len(recs) < 12:
                    continue
                feats = {f: np.array([r[f] for r in recs], dtype=float) for f in FEATURES}
                if isinstance(spec, str):
                    bscores = feats[spec]
                else:
                    cols = []
                    for f in FEATURES:
                        v = feats[f]
                        sd = float(np.std(v))
                        cols.append((v - v.mean()) / sd if sd > 0 else np.zeros_like(v))
                    bscores = np.column_stack(cols) @ np.asarray(spec, dtype=float)
                r, r2a = fit_residuals(bscores, feats)
                b_resid[(protein, ctx)] = r
                b_frac.setdefault(protein, []).append(1.0 - r2a)
        baselines[bname] = _stable_non_free(
            _med(b_frac), _cross_utr_rho(proteins, contexts, b_resid)
        )

    delta = {
        p: stable_ribonn[p] - max(b.get(p, 0.0) for b in baselines.values())
        for p in stable_ribonn
    }
    gate_med, gate_lo, gate_hi = cluster_bootstrap_ci(delta)
    gate_pass = bool(gate_lo > 0.0)
    stab_med, stab_lo, stab_hi = cluster_bootstrap_ci(cross_utr)

    rf_med, rf_lo, rf_hi = cluster_bootstrap_ci(_med(resid_fraction))
    rg_med, rg_lo, rg_hi = cluster_bootstrap_ci(_med(resid_gc))
    rt_med, rt_lo, rt_hi = cluster_bootstrap_ci(_med(retention))

    void_reasons = []
    if rg_med > VOID_RESIDUAL_GC_MAX:
        void_reasons.append(
            f"median |rho(residual, GC)| = {rg_med:.3f} > {VOID_RESIDUAL_GC_MAX} "
            "-- the free-feature fit did not remove GC"
        )
    wob = None
    if within_sd and len(group_means) > 1:
        med_within = statistics.median(statistics.median(v) for v in within_sd.values())
        between = statistics.pstdev(group_means)
        wob = med_within / between if between > 0 else None
        if wob is not None and wob < SANITY_WITHIN_OVER_BETWEEN_MIN:
            void_reasons.append(
                f"within_over_between = {wob:.3f} < {SANITY_WITHIN_OVER_BETWEEN_MIN} "
                "-- harness sanity floor failed"
            )

    return {
        "pre_registration": "docs/PREREG_ribonn_part3_round2.md",
        "n_proteins": len(proteins),
        "n_contexts": len(contexts),
        "contexts": contexts,
        "void": bool(void_reasons),
        "void_reasons": void_reasons,
        "gate": {
            "name": (
                "stable non-free fraction (residual_fraction x cross-UTR rho), "
                "minus the best of the pure-GC / pure-CAI baselines"
            ),
            "median": gate_med,
            "ci95": [gate_lo, gate_hi],
            "passes": gate_pass and not void_reasons,
            "rule": "CI excludes 0 from above",
            "per_protein": delta,
        },
        "magnitudes": {
            "residual_fraction": {"median": rf_med, "ci95": [rf_lo, rf_hi]},
            "gc_stratified_retention": {"median": rt_med, "ci95": [rt_lo, rt_hi]},
            "residual_gc_correlation": {"median": rg_med, "ci95": [rg_lo, rg_hi]},
            "cross_utr_residual_rho": {"median": stab_med, "ci95": [stab_lo, stab_hi]},
            "stable_non_free_ribonn": {"median": statistics.median(stable_ribonn.values())
                                       if stable_ribonn else float("nan"), "ci95": [float("nan")] * 2},
            "stable_non_free_baseline_floor": {
                "median": max(
                    (statistics.median(b.values()) for b in baselines.values() if b),
                    default=float("nan"),
                ),
                "ci95": [float("nan")] * 2,
                "n_baselines": len(baselines),
            },
            "within_over_between": wob,
        },
        "honesty": (
            "Diagnostics only. Every input score is an uncalibrated RiboNN output in "
            "arbitrary CLR-residual units; no result here can promote a backend, and a "
            "score difference is a statement about the model's output, not about "
            "translation efficiency (CLAUDE.md §6/§10.6)."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    rows = [json.loads(line) for line in args.scores.read_text().splitlines() if line.strip()]
    report = analyze(rows)

    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    g = report["gate"]
    print(f"proteins={report['n_proteins']}  contexts={report['n_contexts']}")
    print()
    if report["void"]:
        print("RUN IS VOID -- neither pass nor fail:")
        for reason in report["void_reasons"]:
            print(f"  - {reason}")
        print()
    print("GATE  cross-UTR stability of the free-feature residual")
    print(f"  median {g['median']:.4f}   95% CI [{g['ci95'][0]:.4f}, {g['ci95'][1]:.4f}]")
    print(f"  rule: {g['rule']}  ->  {'PASS' if g['passes'] else 'FAIL'}")
    print()
    print("MAGNITUDES (for judgement, not gates)")
    for key, val in report["magnitudes"].items():
        if isinstance(val, dict):
            print(f"  {key:28s} {val['median']:.4f}  CI [{val['ci95'][0]:.4f}, {val['ci95'][1]:.4f}]")
        elif val is not None:
            print(f"  {key:28s} {val:.4f}")
    print()
    print(report["honesty"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
