"""Run-to-run reproducibility of anonymized codon-optimization algorithms.

A standalone report over Ranaghan et al. (2021), *BMC Biology* 19:36 (DOI
10.1186/s12915-021-00968-8, CC BY 4.0), **Table 4**: three algorithms each
optimized three human genes (KRas4B, Beclin 1, PDE3A) ten independent times. See
the attribution + honesty caveats in ``scripts/data/ranaghan2021_tab4.LICENSE.md``.

What this is (and is not):

* This is a **reproducibility / variability** view, not a named-tool scoreboard.
  The paper anonymizes the algorithms as *Algorithm 1/2/3*; do NOT map them to
  vendors. The named-tool head-to-head lives separately in
  ``scripts/compare_tools.py`` over ``kras_ranaghan2021.fasta``.
* The **ten runs of one algorithm are repeat runs of a single stochastic tool**,
  so their spread measures that tool's *determinism* (does the same input give
  the same output?), not a difference between tools. Ten runs are not ten tools.
* Every metric is **recomputed by BT4's own functions** from each delivered
  nucleotide sequence (CAI/tAI from the bundled ``homo_sapiens`` tables, GC% and
  the longest homopolymer run from ``bt4._accel``, CpG by direct counting).
  Nothing is copied from the paper; nothing is fabricated. BT4 is added as a
  reference row and is never claimed "better" - being an exact deterministic
  solve, its run-to-run spread is zero by construction, which is the honest
  contrast this panel draws.

Run it directly::

    python scripts/compare_reproducibility.py
    python scripts/compare_reproducibility.py --json
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics
from collections.abc import Mapping, Sequence

from bt4 import api
from bt4._accel import gc_count, max_homopolymer_run
from bt4.biomodels.codon.tables import CodonUsageTable, load_table
from bt4.biomodels.codon.tai import TaiTable, load_tai_table
from bt4.domain.genetic_code import STOP, translate

__all__ = [
    "DEFAULT_PANEL_PATH",
    "load_panel",
    "main",
    "parse_header",
    "reproducibility",
]

DEFAULT_PANEL_PATH = pathlib.Path(__file__).parent / "data" / "ranaghan2021_tab4.fasta"

# The compositional/adaptation metrics recomputed per sequence, in report order.
_METRICS: tuple[str, ...] = ("cai", "tai", "gc_pct", "cpg", "max_homo")

_BT4 = "BT4"


def load_panel(path: str | pathlib.Path = DEFAULT_PANEL_PATH) -> list[tuple[str, str]]:
    """Return the panel's ``(header, sequence)`` records via :func:`bt4.api.read_fasta`."""
    return api.read_fasta(path)


def parse_header(header: str) -> tuple[str, str, int | None]:
    """Split a panel header into ``(protein, source, run_or_None)``.

    Headers look like ``KRas4B|Native|acc=...`` or
    ``KRas4B|Algorithm1|run3|acc=...``. ``run`` is ``None`` for the Native record.
    """
    parts = header.split("|")
    protein = parts[0]
    source = parts[1] if len(parts) > 1 else "?"
    run: int | None = None
    if len(parts) > 2 and parts[2].startswith("run"):
        try:
            run = int(parts[2][3:])
        except ValueError:
            run = None
    return protein, source, run


def _protein(dna: str) -> str:
    """Translate ``dna`` and strip a single trailing stop, for grouping/BT4 input."""
    aa = translate(dna)
    return aa[:-1] if aa.endswith(STOP) else aa


def _metrics(dna: str, table: CodonUsageTable, tai_table: TaiTable) -> dict[str, float]:
    """Recompute every reported metric from ``dna`` using BT4's own functions."""
    seq = dna.upper()
    return {
        "cai": table.cai(seq),
        "tai": tai_table.tai(seq),
        "gc_pct": 100.0 * gc_count(seq) / len(seq),
        "cpg": float(seq.count("CG")),
        "max_homo": float(max_homopolymer_run(seq)),
    }


def _aggregate(values: Sequence[float]) -> dict[str, float]:
    """Return mean, span (max - min), and population std of ``values``."""
    if not values:
        return {"mean": 0.0, "span": 0.0, "std": 0.0}
    return {
        "mean": statistics.fmean(values),
        "span": max(values) - min(values),
        "std": statistics.pstdev(values) if len(values) > 1 else 0.0,
    }


def reproducibility(
    records: Sequence[tuple[str, str]], config: api.OptimizeConfig | None = None
) -> list[dict[str, object]]:
    """Summarize per-``(protein, source)`` run-to-run spread of every metric.

    For each protein: one ``Native`` reference row, one row per algorithm
    aggregating its ten runs (mean/span/std of each metric), and one deterministic
    ``BT4`` reference row (BT4 optimizes the Native protein through
    :func:`bt4.api.optimize`; its spread is zero). All metrics are recomputed from
    the sequences by :func:`_metrics`.

    Args:
        records: ``(header, dna)`` panel records (see :func:`load_panel`).
        config: Optimization config for the BT4 rows; defaults to the engine
            default. The panel proteins are human, so the human tables are used.

    Returns:
        One row dict per ``(protein, source)`` with ``protein``, ``source``,
        ``n_runs``, and a ``{mean, span, std}`` block per metric.
    """
    cfg = config if config is not None else api.OptimizeConfig()
    table = load_table(cfg.organism, reference_set=cfg.reference_set)
    tai_table = load_tai_table(cfg.organism)

    # protein -> source -> list of per-run metric dicts; protein -> native protein aa.
    grouped: dict[str, dict[str, list[dict[str, float]]]] = {}
    native_protein: dict[str, str] = {}
    order: list[str] = []
    for header, dna in records:
        protein, source, _run = parse_header(header)
        if protein not in grouped:
            grouped[protein] = {}
            order.append(protein)
        grouped[protein].setdefault(source, []).append(_metrics(dna, table, tai_table))
        if source == "Native":
            native_protein[protein] = _protein(dna)

    rows: list[dict[str, object]] = []
    for protein in order:
        sources = grouped[protein]
        for source in sorted(sources):
            per_run = sources[source]
            row: dict[str, object] = {
                "protein": protein,
                "source": source,
                "n_runs": len(per_run),
            }
            for metric in _METRICS:
                row[metric] = _aggregate([m[metric] for m in per_run])
            rows.append(row)
        # BT4 reference: a single deterministic solve of the native protein.
        if protein in native_protein:
            result = api.optimize(native_protein[protein], cfg)
            metrics = _metrics(result.dna, table, tai_table)
            bt4_row: dict[str, object] = {"protein": protein, "source": _BT4, "n_runs": 1}
            for metric in _METRICS:
                bt4_row[metric] = _aggregate([metrics[metric]])
            rows.append(bt4_row)
    return rows


def board(
    records: Sequence[tuple[str, str]],
    config: api.OptimizeConfig | None = None,
) -> dict[str, object]:
    """The reproducibility rows plus the tables they were scored with.

    Same reason as ``compare_tools.board``: a ``cai`` column without its reference
    set does not say what it measured.
    """
    cfg = config if config is not None else api.OptimizeConfig()
    reference_set = cfg.reference_set or api.default_reference_set(cfg.organism)
    return {
        "organism": cfg.organism,
        "codon_reference_set": reference_set,
        "rows": reproducibility(records, cfg),
    }


def _agg(row: Mapping[str, object], metric: str, key: str) -> float:
    """Read one aggregate value (``mean``/``span``/``std``) of a metric from a row."""
    block = row[metric]
    assert isinstance(block, dict)
    return float(block[key])


def _format_table(rows: Sequence[Mapping[str, object]]) -> str:
    """Render rows as a fixed-width table of ``mean`` and ``span`` (run-to-run)."""
    headers = ["protein", "source", "n", "cai", "cai_span", "tai", "gc%", "gc_span", "cpg_span"]
    lines: list[list[str]] = [
        [
            str(row["protein"]),
            str(row["source"]),
            str(row["n_runs"]),
            f"{_agg(row, 'cai', 'mean'):.3f}",
            f"{_agg(row, 'cai', 'span'):.3f}",
            f"{_agg(row, 'tai', 'mean'):.3f}",
            f"{_agg(row, 'gc_pct', 'mean'):.1f}",
            f"{_agg(row, 'gc_pct', 'span'):.1f}",
            f"{_agg(row, 'cpg', 'span'):.0f}",
        ]
        for row in rows
    ]
    widths = [
        max(len(headers[i]), *(len(line[i]) for line in lines)) if lines else len(headers[i])
        for i in range(len(headers))
    ]
    out = ["  ".join(h.ljust(widths[i]) for i, h in enumerate(headers))]
    out.append("  ".join("-" * widths[i] for i in range(len(headers))))
    out.extend("  ".join(line[i].ljust(widths[i]) for i in range(len(headers))) for line in lines)
    out.append("")
    out.append("span = run-to-run (max - min) across an algorithm's repeat runs;")
    out.append("Native and BT4 are single sequences (n=1, span=0). Tools anonymized; not a")
    out.append("named-tool scoreboard. Metrics recomputed by BT4; no 'better' claim.")
    return "\n".join(out)


def main(argv: Sequence[str] | None = None) -> int:
    """Print the reproducibility table (or ``--json``) over the Tab 4 panel."""
    parser = argparse.ArgumentParser(
        description="Run-to-run reproducibility of anonymized optimizers (Ranaghan 2021 Tab 4)."
    )
    parser.add_argument("--panel", default=str(DEFAULT_PANEL_PATH), help="panel FASTA path")
    parser.add_argument("--organism", default="homo_sapiens", help="codon/tRNA table key")
    parser.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
        help="which of the organism's CAI reference sets to solve against AND "
             "recompute every row with (default: the organism's own default)",
    )
    parser.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    args = parser.parse_args(argv)

    records = load_panel(args.panel)
    config = api.OptimizeConfig(
        organism=args.organism, reference_set=args.reference_set
    )
    payload = board(records, config)
    rows = payload["rows"]
    assert isinstance(rows, list)

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print("BT4 vs anonymized optimizers - run-to-run reproducibility (Ranaghan 2021, Tab 4)")
        print("Recomputed metrics; tools anonymized (Algorithm 1/2/3); not a named-tool board.")
        # Which CAI these rows report is part of the claim, not a detail.
        print(f"Tables: {payload['organism']} / {payload['codon_reference_set']} reference set")
        print()
        print(_format_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
