#!/usr/bin/env python
"""Score the frozen round-2 panel: every variant, in every UTR context.

Implements the *data collection* half of `docs/PREREG_ribonn_part3_round2.md`. It makes
no decisions and computes no endpoints -- it emits a long-format table that
``prereg_round2_analyze.py`` consumes. Keeping the two apart is deliberate: the
analysis can then be re-read without re-spending four hours of CPU, and the scoring
code cannot be quietly tuned in response to an endpoint it never sees.

What it does, per the pre-registration:

* draws ``--n-variants`` synonymous variants per protein with ``api.library`` at the
  seed recorded in the panel manifest, so the variant set is reproducible;
* scores **all** proteins' variants for a given UTR context in **one**
  ``score_many`` invocation, because RiboNN's per-invocation overhead is large and
  fixed -- 4 invocations total rather than 64;
* recomputes every free feature (GC, GC3, CAI, tAI, CpG density, length) from the
  variant's own DNA, never from an accumulator (CLAUDE.md invariant #2);
* writes each context's rows as soon as they exist, so a 3.5 h run survives an
  interruption and resumes at the next context.

Honesty: every score is an **uncalibrated** model output in arbitrary CLR-residual
units. Nothing here can promote a backend, and a score difference is a statement about
RiboNN's output, not about translation efficiency.

Usage::

    python scripts/prereg_round2_score.py \\
        --panel scripts/data/prereg_round2_panel.json \\
        --out   round2_scores.jsonl
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
TEMPERATURE = 1.0
ORGANISM = "homo_sapiens"


def gc3_fraction(dna: str) -> float:
    """GC fraction at codon third positions.

    Identical in definition to ``scripts/ribonn_sensitivity.py``'s function of the same
    name, so round 1 and round 2 measure the same quantity.
    """
    thirds = dna[2::3]
    if not thirds:
        return 0.0
    return sum(1 for base in thirds if base in "GC") / len(thirds)


def cpg_density(dna: str) -> float:
    """CpG dinucleotides per position (0.0 for a sequence shorter than 2)."""
    if len(dna) < 2:
        return 0.0
    return sum(1 for i in range(len(dna) - 1) if dna[i : i + 2] == "CG") / (len(dna) - 1)


def build_variants(panel: dict[str, Any], n: int) -> dict[str, list[str]]:
    """Draw ``n`` synonymous variants per protein at the manifest's recorded seeds."""
    out: dict[str, list[str]] = {}
    for entry in panel["proteins"]:
        symbol = entry["symbol"]
        seed = panel["variant_seeds"][symbol]
        lib = api.library(entry["protein"], n=n, seed=seed, temperature=TEMPERATURE)
        dnas = [res.dna for res in lib.results]
        if lib.distinct != len(dnas):
            # Recorded rather than silently tolerated: duplicates would shrink the
            # effective within-protein n without changing the row count. `distinct` is
            # the library's own honest diversity statistic, not a recount here.
            print(
                f"  WARNING {symbol}: only {lib.distinct}/{len(dnas)} variants are distinct",
                file=sys.stderr,
            )
        out[symbol] = dnas
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-variants", type=int, default=40)
    args = ap.parse_args()

    panel = json.loads(args.panel.read_text())
    print(f"panel content_hash={panel.get('content_hash')}", file=sys.stderr)
    print(
        f"proteins={len(panel['proteins'])} contexts={len(panel['utr_contexts'])} "
        f"variants={args.n_variants}",
        file=sys.stderr,
    )

    table = load_table(ORGANISM)
    tai_table = load_tai_table(ORGANISM)

    print("drawing variants ...", file=sys.stderr)
    variants = build_variants(panel, args.n_variants)

    # Free features depend only on the CDS, so compute once and reuse per context.
    features: dict[str, dict[str, float]] = {}
    for dnas in variants.values():
        for dna in dnas:
            if dna in features:
                continue
            features[dna] = {
                "gc": gc_count(dna) / len(dna),
                "gc3": gc3_fraction(dna),
                "cai": table.cai(dna),
                "tai": tai_table.tai(dna),
                "cpg": cpg_density(dna),
                "length_nt": float(len(dna)),
            }

    done_contexts: set[str] = set()
    if args.out.exists():
        for line in args.out.read_text().splitlines():
            if line.strip():
                done_contexts.add(json.loads(line)["utr_context"])
        if done_contexts:
            print(f"resuming; already scored: {sorted(done_contexts)}", file=sys.stderr)

    flat_symbols = [s for s in variants for _ in variants[s]]
    flat_dnas = [d for s in variants for d in variants[s]]

    with args.out.open("a", encoding="utf-8", newline="\n") as fh:
        for ctx in panel["utr_contexts"]:
            name = ctx["name"]
            if name in done_contexts:
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
            print(f"  {elapsed:.0f}s ({elapsed / max(len(flat_dnas), 1):.2f}s/seq)", file=sys.stderr)

            for i, (symbol, dna, res) in enumerate(
                zip(flat_symbols, flat_dnas, results, strict=True)
            ):
                row = {
                    "protein": symbol,
                    "variant_index": i,
                    "utr_context": name,
                    "score": res.score,
                    "units": res.units,
                    "calibrated": res.calibrated,
                    **features[dna],
                }
                fh.write(json.dumps(row, sort_keys=True) + "\n")
            fh.flush()

    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
