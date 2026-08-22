#!/usr/bin/env python
"""Score the round-3 panel: temperature-ladder variants, in every UTR context.

Implements the data-collection half of `docs/PREREG_ribonn_part3_round3.md`. It makes no
decisions and computes no endpoints -- `scripts/prereg_round2_analyze.py` is reused
**verbatim** for that, so the deciding code is byte-identical to what was frozen before
round 2's data existed.

Two things differ from round 2's scorer, and both come from round 2's void:

1. **A temperature ladder, not a single temperature.** Round 2 drew every variant at
   `temperature=1.0`, which samples the natural codon distribution and clusters tightly:
   within-protein GC spans of 0.020-0.055, against 0.062-0.145 for the optimizer outputs
   that produced round 1's usable signal. RiboNN's response was small because the design
   space was small. Measured, temperature does not *widen* a sample (each rung spans
   ~0.02-0.05) -- it *moves* it, so the **union** across rungs is what reaches
   deployment-relevant amplitude (0.093-0.113 GC, 0.247-0.288 CAI).

2. **Adequacy is checked on the inputs, before scoring.** Round 2's floor
   (`within_over_between`) was computed from the model's own output, so it could not
   distinguish "the panel is too narrow" from "the model is insensitive" -- opposite
   conclusions, same number. The GC/CAI spans below are properties of the sequences, not
   of RiboNN, and they are checked **first**: an inadequate panel costs seconds instead of
   three hours.

The positive control is the one legitimate model-output sanity test here, because it has
a known expected direction independent of the hypothesis: the same sequence scored under
two different UTR contexts must move.

Honesty: every score is an **uncalibrated** model output in arbitrary CLR-residual units.
Nothing here can promote a backend, and a score difference is a statement about RiboNN's
output, not about translation efficiency.

Usage::

    python scripts/prereg_round3_score.py \\
        --panel scripts/data/prereg_round3_panel.json \\
        --out   round3_scores.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

from bt4 import api
from bt4._accel import gc_count
from bt4.biomodels.codon.tables import load_table
from bt4.biomodels.codon.tai import load_tai_table
from bt4.biomodels.expression import RiboNNExpressionModel

# Bound by the pre-registration; not tunables.
TOP_K = 5
BATCH_SIZE = 64
NUM_WORKERS = 0
SPECIES = "human"
CELL_TYPES = ("HEK293T",)
ORGANISM = "homo_sapiens"

# Section 2: six rungs, seven draws each = 42 variants per protein.
LADDER = (0.4, 0.7, 1.0, 1.6, 2.6, 4.2)
PER_RUNG = 7

# Section 3: anchored to the *minimum round 1 actually achieved* (PDE3A: GC 0.0622,
# CAI 0.0905), not chosen for convenience, and verified achievable before freezing.
MIN_GC_SPAN = 0.060
MIN_CAI_SPAN = 0.090
MAX_INADEQUATE = 2
MIN_CONTROL_DELTA = 0.01


def gc3_fraction(dna: str) -> float:
    """GC fraction at codon third positions (identical to rounds 1 and 2)."""
    thirds = dna[2::3]
    return sum(1 for b in thirds if b in "GC") / len(thirds) if thirds else 0.0


def cpg_density(dna: str) -> float:
    """CpG dinucleotides per position."""
    if len(dna) < 2:
        return 0.0
    return sum(1 for i in range(len(dna) - 1) if dna[i : i + 2] == "CG") / (len(dna) - 1)


def build_ladder_variants(panel: dict[str, Any]) -> dict[str, list[str]]:
    """Draw the temperature ladder for every protein, at the manifest's seeds."""
    out: dict[str, list[str]] = {}
    for entry in panel["proteins"]:
        symbol = entry["symbol"]
        base_seed = panel["variant_seeds"][symbol]
        dnas: list[str] = []
        for rung, temp in enumerate(LADDER):
            lib = api.library(
                entry["protein"], n=PER_RUNG, seed=base_seed + rung, temperature=temp
            )
            dnas.extend(res.dna for res in lib.results)
        distinct = len(set(dnas))
        if distinct != len(dnas):
            print(
                f"  WARNING {symbol}: only {distinct}/{len(dnas)} variants are distinct",
                file=sys.stderr,
            )
        out[symbol] = dnas
    return out


def check_adequacy(
    variants: dict[str, list[str]], table: Any
) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Return (inadequate protein symbols, per-protein spans). Model-independent."""
    spans: dict[str, dict[str, float]] = {}
    inadequate: list[str] = []
    for symbol, dnas in variants.items():
        gc = [gc_count(d) / len(d) for d in dnas]
        cai = [table.cai(d) for d in dnas]
        gc_span, cai_span = max(gc) - min(gc), max(cai) - min(cai)
        spans[symbol] = {"gc_span": gc_span, "cai_span": cai_span}
        if gc_span < MIN_GC_SPAN or cai_span < MIN_CAI_SPAN:
            inadequate.append(symbol)
    return inadequate, spans


def positive_control(panel: dict[str, Any], dna: str) -> float:
    """Score one sequence under two different UTR contexts; return |difference|."""
    contexts = panel["utr_contexts"][:2]
    scores = []
    for ctx in contexts:
        model = RiboNNExpressionModel(
            species=SPECIES,
            utr5=ctx["utr5"],
            utr3=ctx["utr3"],
            cell_types=CELL_TYPES,
            top_k=TOP_K,
            batch_size=BATCH_SIZE,
            num_workers=NUM_WORKERS,
        )
        scores.append(model.score_sequence(dna).score)
    return abs(scores[0] - scores[1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--check-only",
        action="store_true",
        help="run adequacy + positive control and stop, without scoring",
    )
    args = ap.parse_args()

    panel = json.loads(args.panel.read_text())
    print(f"panel content_hash={panel.get('content_hash')}", file=sys.stderr)
    table = load_table(ORGANISM)
    tai_table = load_tai_table(ORGANISM)

    print(f"drawing ladder {LADDER} x {PER_RUNG} ...", file=sys.stderr)
    variants = build_ladder_variants(panel)

    inadequate, spans = check_adequacy(variants, table)
    print("\n--- adequacy (inputs only; model not consulted) ---", file=sys.stderr)
    print(f"{'protein':12s} {'GC span':>9s} {'CAI span':>9s}  verdict", file=sys.stderr)
    for symbol in sorted(spans):
        s = spans[symbol]
        ok = symbol not in inadequate
        print(
            f"{symbol:12s} {s['gc_span']:9.4f} {s['cai_span']:9.4f}  "
            f"{'ok' if ok else 'INADEQUATE'}",
            file=sys.stderr,
        )
    print(
        f"thresholds: GC >= {MIN_GC_SPAN}, CAI >= {MIN_CAI_SPAN}; "
        f"{len(inadequate)} inadequate (void if > {MAX_INADEQUATE})",
        file=sys.stderr,
    )
    if len(inadequate) > MAX_INADEQUATE:
        print(
            f"\nVOID: {len(inadequate)} proteins fail the adequacy floor -- "
            "the panel does not exercise the axis under test. Nothing scored.",
            file=sys.stderr,
        )
        return 2
    for symbol in inadequate:
        print(f"dropping inadequate protein: {symbol}", file=sys.stderr)
        variants.pop(symbol)

    probe = next(iter(variants.values()))[0]
    delta = positive_control(panel, probe)
    print(f"\npositive control |delta| across two UTR contexts = {delta:.6f}", file=sys.stderr)
    if delta < MIN_CONTROL_DELTA:
        print(
            f"VOID: positive control below {MIN_CONTROL_DELTA} -- the UTR context is not "
            "reaching the model, so nothing downstream is interpretable.",
            file=sys.stderr,
        )
        return 2
    if args.check_only:
        print("--check-only: adequacy and control pass; stopping before scoring.", file=sys.stderr)
        return 0

    features: dict[str, dict[str, float]] = {}
    for dnas in variants.values():
        for dna in dnas:
            if dna not in features:
                features[dna] = {
                    "gc": gc_count(dna) / len(dna),
                    "gc3": gc3_fraction(dna),
                    "cai": table.cai(dna),
                    "tai": tai_table.tai(dna),
                    "cpg": cpg_density(dna),
                    "length_nt": float(len(dna)),
                }

    done: set[str] = set()
    if args.out.exists():
        for line in args.out.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["utr_context"])
        if done:
            print(f"resuming; already scored: {sorted(done)}", file=sys.stderr)

    flat_symbols = [s for s in variants for _ in variants[s]]
    flat_dnas = [d for s in variants for d in variants[s]]

    with args.out.open("a", encoding="utf-8", newline="\n") as fh:
        for ctx in panel["utr_contexts"]:
            name = ctx["name"]
            if name in done:
                continue
            model = RiboNNExpressionModel(
                species=SPECIES,
                utr5=ctx["utr5"],
                utr3=ctx["utr3"],
                cell_types=CELL_TYPES,
                top_k=TOP_K,
                batch_size=BATCH_SIZE,
                num_workers=NUM_WORKERS,
            )
            print(
                f"scoring {len(flat_dnas)} sequences in context {name} "
                f"(5'={ctx['utr5_len']}, 3'={ctx['utr3_len']}) ...",
                file=sys.stderr,
            )
            t0 = time.time()
            results = model.score_many(flat_dnas)
            elapsed = time.time() - t0
            per = elapsed / max(len(flat_dnas), 1)
            print(f"  {elapsed:.0f}s ({per:.2f}s/seq)", file=sys.stderr)
            for i, (symbol, dna, res) in enumerate(
                zip(flat_symbols, flat_dnas, results, strict=True)
            ):
                fh.write(
                    json.dumps(
                        {
                            "protein": symbol,
                            "variant_index": i,
                            "utr_context": name,
                            "score": res.score,
                            "units": res.units,
                            "calibrated": res.calibrated,
                            **features[dna],
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )
            fh.flush()

    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
