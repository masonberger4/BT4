"""Uncertainty / sensitivity analysis over the BT4 public API.

CLAUDE.md section 8 asks BT4 to "report sensitivity of CAI/GC to table choice and
solver budget; no point estimate presented as ground truth". This standalone
script does exactly that for one protein, entirely through :mod:`bt4.api`:

* Axis (a) -- codon-table / organism choice. It sweeps every organism returned
  by :func:`bt4.api.available_organisms`, optimizes the protein under each, and
  tabulates the delivered result's CAI, GC%, and (when the organism ships a tRNA
  table) tAI. Entries that are not usable codon-usage tables -- e.g. the bundled
  tRNA-only ``*.trna`` tables surfaced by the organism discovery -- are reported
  as skipped with their reason, never silently dropped.
* Axis (b) -- solver budget. For a single organism it compares the exact DP
  (``beam=None``) against a few beam widths, tabulating the same metrics plus the
  optimality certificate, and flags whether the certificate degrades from
  ``proven_optimal`` to ``beam_truncated``.

The whole point is honesty about non-uniqueness: the numbers reported are a
*spread*, driven by which table and how much solver budget you pick, not a single
ground-truth answer. Every metric is recomputed from the delivered sequence by
the API (CAI and GC% from the result audit; tAI as the additive ``tai_logw``
objective recomputed on the delivered DNA via :func:`bt4.api.validate`). The tAI
figure is the log-sum-of-relative-adaptiveness objective (larger, i.e. closer to
zero, is better), recomputed on the delivered sequence -- it is not the tAI of a
separately re-optimized sequence, and it is deliberately not dressed up as a
calibrated expression prediction.

Determinism: the analysis performs no sampling of its own; it threads the config
seed into every solve and never touches the global RNG or the wall clock, so an
identical invocation yields byte-identical output.

Run it directly::

    python scripts/sensitivity.py --protein MKTAYIAKQR
    python scripts/sensitivity.py --protein MKTAYIAKQR --beams 1,2,5 --json

This is a report script, not part of the ``bt4`` package import graph; it reaches
only into the stable :mod:`bt4.api` surface.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from dataclasses import replace
from typing import cast

from bt4 import api

__all__ = [
    "DEFAULT_BEAMS",
    "DEFAULT_PROTEIN",
    "HONESTY_NOTE",
    "analyze",
    "budget_sensitivity",
    "main",
    "organism_sensitivity",
    "tai_logw",
]

# A small demonstration protein with several genuinely degenerate residues, so
# the codon trellis has real choices to trade off across tables and beam widths.
DEFAULT_PROTEIN = "MKTAYIAKQR"

# Beam widths compared against the exact DP on the solver-budget axis. Small
# widths reliably prune (and so honestly report a truncated certificate) for a
# protein with synonymous codons; the exact solve remains proven-optimal.
DEFAULT_BEAMS: tuple[int, ...] = (1, 2, 5)

# The number of decimal places metrics are rounded to for both table and JSON.
_ROUND = 6

# Certificate status string an exact, fully-explored solve carries.
_PROVEN = "proven_optimal"

HONESTY_NOTE = (
    "These are a SPREAD, not a ground-truth answer: CAI/GC/tAI move with the "
    "codon table you pick and with how much solver budget you spend. tAI is the "
    "additive tai_logw objective (larger = better) recomputed on the delivered "
    "sequence via the API, reported only for organisms shipping a tRNA table; it "
    "is not a calibrated expression prediction. Beam-capped solves are labeled "
    "beam_truncated exactly when the beam actually pruned states."
)


def tai_logw(dna: str, organism: str) -> float | None:
    """Recompute the additive tAI (``tai_logw``) objective for ``dna``, if available.

    tAI is defined only for organisms that ship a tRNA gene-copy-number table.
    Availability is detected purely through the public surface: organism ``X``
    has tAI data iff it appears in :func:`bt4.api.available_tai_organisms`. When
    available, the value is recomputed
    on the exact ``dna`` by :func:`bt4.api.validate` with the tAI term switched on
    (the weight value does not affect the recomputed score, only whether the term
    is active).

    Args:
        dna: A delivered coding sequence (length a multiple of three).
        organism: The organism whose tRNA table defines tAI.

    Returns:
        The ``tai_logw`` objective value (larger, i.e. closer to zero, is
        better), or ``None`` when ``organism`` ships no tRNA table.
    """
    if organism not in api.available_tai_organisms():
        return None
    report = api.validate(dna, api.OptimizeConfig(organism=organism, tai_weight=1.0))
    if "tai_logw" not in report.metrics.objective.terms():
        return None
    return report.metrics.objective.get("tai_logw")


def _delivered_metrics(
    protein: str, config: api.OptimizeConfig
) -> tuple[float, float, float | None, str]:
    """Optimize once and pull the delivered CAI, GC%, tAI, and certificate status.

    tAI is neutralized in the objective (``tai_weight=0``) so the delivered CAI
    and GC% stay comparable across organisms; tAI is then recomputed on that same
    delivered sequence via :func:`tai_logw`.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: The base configuration (its ``organism`` selects the table).

    Returns:
        A tuple ``(cai, gc_percent, tai_logw_or_None, certificate_status)``.
    """
    result = api.optimize(protein, replace(config, tai_weight=0.0))
    cai = float(cast("float", result.audit["cai"]))
    gc_percent = float(cast("float", result.audit["gc_percent"]))
    tlw = tai_logw(result.dna, config.organism)
    return cai, gc_percent, tlw, result.certificate.status.value


def _spread(rows: Sequence[dict[str, object]], key: str) -> dict[str, object]:
    """Summarize the min/max/range of a numeric metric across ``rows``.

    Args:
        rows: Result rows, each possibly carrying ``key`` (``None`` values and
            missing keys are ignored).
        key: The metric column to summarize.

    Returns:
        A dict with ``n`` (count of present values) and, when non-empty,
        ``min``/``max``/``range`` (rounded); ``min``/``max``/``range`` are
        ``None`` when no value is present.
    """
    nums = [float(cast("float", r[key])) for r in rows if r.get(key) is not None]
    if not nums:
        return {"n": 0, "min": None, "max": None, "range": None}
    lo, hi = min(nums), max(nums)
    return {
        "n": len(nums),
        "min": round(lo, _ROUND),
        "max": round(hi, _ROUND),
        "range": round(hi - lo, _ROUND),
    }


def organism_sensitivity(
    protein: str,
    config: api.OptimizeConfig | None = None,
    organisms: Sequence[str] | None = None,
) -> dict[str, object]:
    """Sweep codon-table / organism choice and report the metric spread.

    Each usable organism yields one exact solve; entries that are not usable
    codon-usage tables (for example the tRNA-only ``*.trna`` tables surfaced by
    organism discovery) are recorded under ``skipped`` with the failure reason
    rather than silently dropped.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Base configuration; defaults to :class:`bt4.api.OptimizeConfig`.
            Its ``organism`` is overridden per sweep entry.
        organisms: Explicit organism list to sweep; defaults to the sorted output
            of :func:`bt4.api.available_organisms`.

    Returns:
        A JSON-serializable dict with ``axis``, ``protein``, ``rows`` (one per
        usable organism), ``skipped`` (organism + reason), and ``spread``
        (min/max/range of ``cai``, ``gc_percent``, and ``tai_logw``).
    """
    cfg = config if config is not None else api.OptimizeConfig()
    names = sorted(organisms) if organisms is not None else sorted(api.available_organisms())
    rows: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for name in names:
        try:
            cai, gc_percent, tlw, status = _delivered_metrics(protein, replace(cfg, organism=name))
        except ValueError as exc:
            skipped.append({"organism": name, "reason": str(exc)})
            continue
        rows.append(
            {
                "organism": name,
                "cai": round(cai, _ROUND),
                "gc_percent": round(gc_percent, _ROUND),
                "tai_available": tlw is not None,
                "tai_logw": None if tlw is None else round(tlw, _ROUND),
                "certificate": status,
            }
        )
    spread = {metric: _spread(rows, metric) for metric in ("cai", "gc_percent", "tai_logw")}
    return {
        "axis": "organism",
        "protein": protein,
        "rows": rows,
        "skipped": skipped,
        "spread": spread,
    }


def budget_sensitivity(
    protein: str,
    config: api.OptimizeConfig | None = None,
    beams: Sequence[int] = DEFAULT_BEAMS,
) -> dict[str, object]:
    """Sweep solver budget (exact DP vs beam widths) for one organism.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Base configuration; defaults to :class:`bt4.api.OptimizeConfig`.
            Its ``organism`` fixes the table for the whole sweep.
        beams: Beam widths to compare against the exact DP; each must be ``>= 1``.

    Returns:
        A JSON-serializable dict with ``axis``, ``protein``, ``organism``,
        ``rows`` (the exact solve first, then one per beam width), ``spread``
        (min/max/range of ``cai``, ``gc_percent``, ``tai_logw``), and
        ``certificate_degrades`` (True iff the exact solve is proven-optimal yet
        some beam-capped solve is not).
    """
    cfg = config if config is not None else api.OptimizeConfig()
    budgets: list[int | None] = [None, *beams]
    rows: list[dict[str, object]] = []
    for beam in budgets:
        cai, gc_percent, tlw, status = _delivered_metrics(protein, replace(cfg, beam=beam))
        rows.append(
            {
                "budget": "exact" if beam is None else f"beam={beam}",
                "beam": beam,
                "cai": round(cai, _ROUND),
                "gc_percent": round(gc_percent, _ROUND),
                "tai_logw": None if tlw is None else round(tlw, _ROUND),
                "certificate": status,
            }
        )
    exact_status = rows[0]["certificate"]
    degrades = exact_status == _PROVEN and any(
        r["certificate"] != _PROVEN for r in rows[1:]
    )
    spread = {metric: _spread(rows, metric) for metric in ("cai", "gc_percent", "tai_logw")}
    return {
        "axis": "budget",
        "protein": protein,
        "organism": cfg.organism,
        "rows": rows,
        "spread": spread,
        "certificate_degrades": degrades,
    }


def analyze(
    protein: str,
    config: api.OptimizeConfig | None = None,
    organisms: Sequence[str] | None = None,
    beams: Sequence[int] = DEFAULT_BEAMS,
) -> dict[str, object]:
    """Run both sensitivity axes and bundle them with an honesty note.

    Args:
        protein: A stop-free single-letter amino-acid string.
        config: Base configuration; defaults to :class:`bt4.api.OptimizeConfig`.
        organisms: Organisms to sweep on the table axis (see
            :func:`organism_sensitivity`).
        beams: Beam widths to compare on the budget axis (see
            :func:`budget_sensitivity`).

    Returns:
        A JSON-serializable dict with ``protein``, ``seed``, ``base_organism``,
        ``honesty_note``, ``organism_sensitivity``, and ``budget_sensitivity``.
    """
    cfg = config if config is not None else api.OptimizeConfig()
    return {
        "protein": protein,
        "seed": cfg.seed,
        "base_organism": cfg.organism,
        "honesty_note": HONESTY_NOTE,
        "organism_sensitivity": organism_sensitivity(protein, cfg, organisms),
        "budget_sensitivity": budget_sensitivity(protein, cfg, beams),
    }


def _fmt(value: object) -> str:
    """Format one table cell: ``n/a`` for ``None``, four decimals for floats."""
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _format_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Render fixed-width columns from ``headers`` and pre-ordered cell ``rows``."""
    cells = [[_fmt(v) for v in row] for row in rows]
    widths = [
        max(len(headers[i]), *(len(r[i]) for r in cells)) if cells else len(headers[i])
        for i in range(len(headers))
    ]
    lines = ["  ".join(h.ljust(widths[i]) for i, h in enumerate(headers))]
    lines.append("  ".join("-" * widths[i] for i in range(len(headers))))
    lines.extend(
        "  ".join(r[i].ljust(widths[i]) for i in range(len(headers))) for r in cells
    )
    return "\n".join(lines)


def _render(report: dict[str, object]) -> str:
    """Render the full analysis as a human-readable multi-section report."""
    org = cast("dict[str, object]", report["organism_sensitivity"])
    bud = cast("dict[str, object]", report["budget_sensitivity"])
    org_rows = cast("list[dict[str, object]]", org["rows"])
    bud_rows = cast("list[dict[str, object]]", bud["rows"])
    skipped = cast("list[dict[str, object]]", org["skipped"])

    out: list[str] = []
    out.append(f"BT4 sensitivity analysis  protein={report['protein']}  seed={report['seed']}")
    out.append(f"NOTE: {HONESTY_NOTE}")
    out.append("")
    out.append("(a) codon-table / organism sensitivity")
    out.append(
        _format_table(
            ("organism", "cai", "gc_percent", "tai_logw", "certificate"),
            [
                (r["organism"], r["cai"], r["gc_percent"], r["tai_logw"], r["certificate"])
                for r in org_rows
            ],
        )
    )
    out.append(_spread_line(cast("dict[str, object]", org["spread"])))
    for s in skipped:
        out.append(f"  skipped {s['organism']}: {s['reason']}")
    out.append("")
    out.append(f"(b) solver-budget sensitivity  organism={bud['organism']}")
    out.append(
        _format_table(
            ("budget", "cai", "gc_percent", "tai_logw", "certificate"),
            [
                (r["budget"], r["cai"], r["gc_percent"], r["tai_logw"], r["certificate"])
                for r in bud_rows
            ],
        )
    )
    out.append(_spread_line(cast("dict[str, object]", bud["spread"])))
    out.append(f"  certificate_degrades: {bud['certificate_degrades']}")
    return "\n".join(out)


def _spread_line(spread: dict[str, object]) -> str:
    """Render the per-metric min/max/range summary as one line."""
    parts: list[str] = []
    for metric, summary in spread.items():
        s = cast("dict[str, object]", summary)
        parts.append(
            f"{metric} range={_fmt(s['range'])} "
            f"[min={_fmt(s['min'])}, max={_fmt(s['max'])}]"
        )
    return "  spread: " + "; ".join(parts)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the sensitivity analysis for one protein and print a table (or JSON).

    Args:
        argv: Optional argument vector (defaults to ``sys.argv[1:]``).

    Returns:
        Process exit code (``0`` on success).
    """
    parser = argparse.ArgumentParser(
        description="Report how BT4's delivered CAI/GC/tAI depend on table choice and beam budget."
    )
    parser.add_argument(
        "--protein",
        default=DEFAULT_PROTEIN,
        help=f"Stop-free protein to analyze (default: {DEFAULT_PROTEIN}).",
    )
    parser.add_argument(
        "--organism",
        default="homo_sapiens",
        help="Base organism for the solver-budget axis (default: homo_sapiens).",
    )
    parser.add_argument(
        "--beams",
        default=",".join(str(b) for b in DEFAULT_BEAMS),
        help="Comma-separated beam widths for the budget axis (default: 1,2,5).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Master seed threaded into every solve (default: 0).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the analysis as JSON instead of a formatted report.",
    )
    args = parser.parse_args(argv)

    beams = tuple(int(b) for b in args.beams.split(",") if b.strip())
    if not beams or any(b < 1 for b in beams):
        parser.error("--beams must be a comma-separated list of integers >= 1")

    config = api.OptimizeConfig(organism=args.organism, seed=args.seed)
    report = analyze(args.protein, config, beams=beams)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(_render(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
