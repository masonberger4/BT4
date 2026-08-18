#!/usr/bin/env python3
"""Run the expression acceptance gate on a measured CDS-variant panel, against baselines.

This is the script that decides whether a learned expression head has earned
``calibrated=True`` for BT4's regime. It loads a panel
(:mod:`bt4.biomodels.expression.panel`), scores it with the chosen backend in as few
invocations as the panel's UTR contexts allow, runs
:func:`bt4.biomodels.expression.verify_expression_gate`, and runs **the same gate on a
fixed set of dumb baselines** so the head's number is never read in isolation.

**Why the baselines are not optional.** A within-protein Spearman of 0.3 is worthless if
plain CAI scores 0.35 -- BT4 already computes CAI for free, *inside* the optimizer loop,
so a head that cannot beat it has earned nothing. And because split conformal is valid
for any score function, a **constant predictor** achieves exactly valid coverage; it is
included precisely so that its "pass" on the coverage axis is visible in the same table.
The permutation baseline is the null: the same predictions against shuffled labels.

**The decision rule this script reports** is the one that survives scrutiny: the head's
cluster-bootstrap CI lower bound on the primary metric must exceed **every** baseline's
point estimate. An absolute Spearman threshold is a pre-commitment, not a standard --
there is no community-standard cutoff -- so ``--min-spearman`` is recorded in the output
rather than presented as authoritative.

**Nothing here flips a flag.** The gate returns a report; promotion is a separate,
deliberate step against a recorded attestation. Run the gate **once**, against
thresholds written down beforehand: a few-hundred-row panel does not support a search
over thresholds, and re-running until something passes converts a validation into one.

Run it directly::

    python scripts/run_expression_gate.py --panel panel.tsv \
        --backend ribonn --species human --cell-type HEK293T \
        --within-group --recalibrate --json > gate_result.json

On Windows use ``^`` for line continuations; ``--num-workers 0`` is already the default
and is required there.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import asdict

from bt4 import api
from bt4.pipeline.expression_gate import BASELINES

__all__ = ["BASELINES", "build_report", "main"]


def build_report(comparison: api.GateComparison) -> dict[str, object]:
    """Turn a :class:`~bt4.api.GateComparison` into a JSON-ready report.

    The judgement lives in :mod:`bt4.pipeline.expression_gate`; this only shapes it for
    printing, so the CLI (``bt4 expression-gate``) and this script can never disagree
    about what a result means.
    """
    return {
        "panel_hash": comparison.panel_hash,
        "backend": {
            "name": comparison.backend,
            "calibrated": comparison.backend_calibrated,
        },
        "settings": asdict(comparison.settings),
        "notes": list(comparison.notes),
        "head": asdict(comparison.head),
        "baselines": {name: asdict(report) for name, report in comparison.baselines},
        "verdict": {
            "gate_passed": comparison.head.passed,
            "best_baseline": comparison.best_baseline,
            "best_baseline_spearman": comparison.best_baseline_spearman,
            "head_spearman_ci_low": comparison.head.spearman_ci_low,
            "beats_every_baseline": comparison.beats_every_baseline,
            "interval_is_informative": comparison.interval_is_informative,
            "promotable": comparison.promotable,
        },
        "honesty": (
            "This script flips nothing. 'promotable' means the pre-registered "
            "conditions held on THIS panel; promotion is a separate, recorded step. "
            "min_spearman is a pre-commitment, not a community standard (none exists). "
            "The constant baseline is present because split conformal is valid for any "
            "score function, so its coverage 'pass' must be visible next to the head's."
        ),
    }


def _render(report: Mapping[str, object]) -> str:
    """Render the report as a table (the ``--json`` alternative)."""
    head = report["head"]
    assert isinstance(head, dict)
    baselines = report["baselines"]
    assert isinstance(baselines, dict)
    verdict = report["verdict"]
    assert isinstance(verdict, dict)
    backend = report["backend"]
    assert isinstance(backend, dict)
    settings = report["settings"]
    assert isinstance(settings, dict)

    mode = "within-protein" if settings["within_group"] else "POOLED (not BT4's regime)"
    flag = "calibrated" if backend["calibrated"] else "UNCALIBRATED"
    panel_hash = report["panel_hash"]
    assert isinstance(panel_hash, str)
    lines = [
        f"panel:    sha256 {panel_hash[:16]}...",
        f"backend:  {backend['name']}  [{flag}]",
        f"mode:     {mode}"
        + (", link fitted on calibration fold" if settings["recalibrate"] else ""),
        "",
        f"{'':<14}{'spearman':>10}{'CI low':>9}{'CI high':>9}{'coverage':>10}"
        f"{'width/IQR':>11}",
    ]

    def _row(label: str, data: Mapping[str, object]) -> str:
        return (
            f"{label:<14}{float(data['spearman']):>10.3f}"
            f"{float(data['spearman_ci_low']):>9.3f}"
            f"{float(data['spearman_ci_high']):>9.3f}"
            f"{float(data['empirical_coverage']):>10.3f}"
            f"{float(data['width_over_iqr']):>11.3f}"
        )

    lines.append(_row("HEAD", head))
    for name, data in baselines.items():
        lines.append(_row(f"  {name}", data))

    lines += [
        "",
        f"gate passed (thresholds)   : {verdict['gate_passed']}",
        f"beats every baseline       : {verdict['beats_every_baseline']} "
        f"(best: {verdict['best_baseline']} at "
        f"{float(verdict['best_baseline_spearman']):.3f})",
        f"interval is informative    : {verdict['interval_is_informative']}",
        f"PROMOTABLE on this panel   : {verdict['promotable']}",
        "",
        f"honesty: {report['honesty']}",
    ]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else "",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--panel", required=True, help="panel TSV (see panel.py)")
    parser.add_argument("--backend", default="ribonn")
    parser.add_argument("--species", default="human", choices=("human", "mouse"))
    parser.add_argument(
        "--cell-type", action="append", dest="cell_types", default=None,
        help="restrict RiboNN to this cell type (repeatable); match it to the panel",
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument(
        "--within-group", action="store_true",
        help="score inside each protein (the strict bar; BT4's actual regime)",
    )
    parser.add_argument(
        "--recalibrate", action="store_true",
        help="fit measured ~ a*pred + b on the calibration fold before residuals",
    )
    parser.add_argument("--min-spearman", type=float, default=0.30)
    parser.add_argument("--target-coverage", type=float, default=0.90)
    parser.add_argument("--coverage-tolerance", type=float, default=0.05)
    parser.add_argument("--calibration-fraction", type=float, default=0.50)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--baselines", default=",".join(BASELINES),
        help="comma-separated; the head must beat all of them",
    )
    parser.add_argument("--organism", default="homo_sapiens")
    parser.add_argument("--reference-set", default=None)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    baselines = [name.strip() for name in args.baselines.split(",") if name.strip()]

    try:
        panel = api.read_panel(args.panel)
        comparison = api.expression_gate(
            panel,
            args.backend,
            settings=api.GateSettings(
                within_group=args.within_group,
                recalibrate=args.recalibrate,
                target_coverage=args.target_coverage,
                coverage_tolerance=args.coverage_tolerance,
                min_spearman=args.min_spearman,
                calibration_fraction=args.calibration_fraction,
                bootstrap_resamples=args.bootstrap_resamples,
                seed=args.seed,
            ),
            baselines=baselines,
            species=args.species,
            cell_types=tuple(args.cell_types or ()),
            top_k=args.top_k,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            organism=args.organism,
            reference_set=args.reference_set,
        )
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    report = build_report(comparison)
    report["panel_summary"] = panel.describe()

    if not args.within_group:
        print(
            "warning: pooled mode credits between-protein skill, which is NOT the "
            "regime BT4 deploys in. Pass --within-group for the strict bar.",
            file=sys.stderr,
        )

    print(json.dumps(report, indent=2, sort_keys=True) if args.json else _render(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
