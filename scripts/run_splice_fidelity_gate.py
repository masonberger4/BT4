"""Run the splice integration-fidelity gate and, on a pass, emit an attestation.

Steps **A5** and **A6** of
[`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md).
Takes the panel built by ``scripts/make_splice_panel.py`` and the expectations
captured by ``scripts/capture_pangolin_panel.py`` (which never imports ``bt4``),
runs :func:`~bt4.biomodels.splice.verify_pangolin_fidelity`, and reports.

**What a pass does and does not mean.** Passing proves BT4's adapter reproduces
the published model's own numbers -- *integration fidelity*. It does **not** show
those numbers are calibrated probabilities in BT4's regime; that is a separate,
still-unmet gate (Part B of the runbook), and the honest scope limit is that these
models score median prAUC 0.419 on **exonic** variants versus 0.773 intronic
(Smith & Kitzman 2023), while BT4 designs coding sequence.

**Panel-strength reporting, and why it is here.** A fidelity gate on a panel where
the model scores ~0 everywhere is nearly vacuous: a wrong output channel, a wrong
fold set, or a transposed one-hot would still "match" within tolerance because
both sides are near zero. So this script reports the **spread** of the captured
expectations and warns when that spread is too small to discriminate. A gate that
passes on a flat panel is reported as passing *and* as weak -- never as a clean
bill of health (CLAUDE.md sections 5, 10.6).

Run it::

    python scripts/run_splice_fidelity_gate.py \\
        --panel panel_sequences.json --captured expected_pangolin.json

    python scripts/run_splice_fidelity_gate.py ... \\
        --attest-out src/bt4/biomodels/splice/data/pangolin.attestation.json

Exit status is ``0`` when the gate passes, ``1`` when it fails.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

__all__ = ["load_cases", "main", "spread_report"]

# Below this range of peak scores the panel cannot distinguish a correct wiring
# from several wrong ones, so a pass on it is reported as weak.
MIN_USEFUL_PEAK_SPREAD = 0.05
"""Minimum (max peak - min peak) across the panel for the gate to be discriminating."""

MIN_USEFUL_MAX_PEAK = 0.10
"""At least one case should score meaningfully above zero, or nothing is exercised."""


def spread_report(cases: Sequence[dict]) -> dict[str, float]:
    """Return peak-score statistics across the captured panel.

    Args:
        cases: Captured cases, each with ``expected_site_scores``.

    Returns:
        ``min_peak`` / ``max_peak`` / ``mean_peak`` / ``spread`` over per-case peaks.
    """
    peaks = [max(c["expected_site_scores"]) for c in cases if c["expected_site_scores"]]
    if not peaks:
        return {"min_peak": 0.0, "max_peak": 0.0, "mean_peak": 0.0, "spread": 0.0}
    return {
        "min_peak": min(peaks),
        "max_peak": max(peaks),
        "mean_peak": sum(peaks) / len(peaks),
        "spread": max(peaks) - min(peaks),
    }


def load_cases(captured: dict) -> list:
    """Turn a capture payload into ``FidelityCase`` objects.

    Split out from :func:`main` so the parsing and pairing logic is testable
    without torch or the licensed weights installed (the gate itself needs the
    real adapter, but everything around it should still be covered by CI).
    """
    from bt4.biomodels.splice import FidelityCase

    return [
        FidelityCase(
            sequence=c["sequence"],
            expected_site_scores=tuple(c["expected_site_scores"]),
        )
        for c in captured["cases"]
    ]


def panels_match(captured: dict, panel: dict) -> bool:
    """Return whether ``captured`` was produced from ``panel``.

    Binding the two by content hash is what stops a stale capture from being
    silently compared against a regenerated panel -- the scores would still be
    numbers, they would just be the wrong ones.
    """
    return captured.get("panel_content_hash") == panel.get("content_hash")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the fidelity gate over a captured panel and report honestly."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--captured", required=True, help="Captured scores JSON.")
    parser.add_argument("--panel", default=None, help="Panel JSON, to check the pairing.")
    parser.add_argument(
        "--tolerance", type=float, default=1e-3, help="Max abs deviation (default 1e-3)."
    )
    parser.add_argument("--attest-out", default=None, help="On a pass, write the attestation here.")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a table.")
    args = parser.parse_args(argv)

    import bt4
    from bt4.biomodels.splice import (
        PangolinSplicePredictor,
        attest_backend,
        verify_pangolin_fidelity,
    )
    from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256

    captured = json.loads(Path(args.captured).read_text(encoding="utf-8"))
    cases_raw = captured["cases"]

    # Bind the capture to the panel it came from, so a mismatched pairing is caught
    # rather than silently compared against the wrong sequences.
    if args.panel:
        panel = json.loads(Path(args.panel).read_text(encoding="utf-8"))
        if not panels_match(captured, panel):
            print("ERROR: captured scores were not produced from this panel")
            print(f"  panel    content_hash: {panel.get('content_hash')}")
            print(f"  captured panel_hash  : {captured.get('panel_content_hash')}")
            return 1

    spread = spread_report(cases_raw)
    panel_weak = (
        spread["spread"] < MIN_USEFUL_PEAK_SPREAD or spread["max_peak"] < MIN_USEFUL_MAX_PEAK
    )

    cases = load_cases(captured)
    report = verify_pangolin_fidelity(PangolinSplicePredictor(), cases, tolerance=args.tolerance)

    result = {
        "passed": report.passed,
        "max_abs_deviation": report.max_abs_deviation,
        "n_cases": report.n_cases,
        "tolerance": report.tolerance,
        "panel_weak": panel_weak,
        **spread,
    }

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("\n=== Pangolin integration-fidelity gate ===")
        print(f"  cases              {report.n_cases}")
        print(f"  tolerance          {report.tolerance}")
        print(f"  max abs deviation  {report.max_abs_deviation:.3e}")
        print(f"  PASSED             {report.passed}")
        print("\n  panel peak P(splice):")
        print(f"    min {spread['min_peak']:.4f}   max {spread['max_peak']:.4f}   "
              f"mean {spread['mean_peak']:.4f}   spread {spread['spread']:.4f}")
        if panel_weak:
            print("\n  WARNING: this panel barely discriminates. Every case scores in a")
            print("  narrow band, so a wrong channel/fold/one-hot could still match within")
            print("  tolerance. Treat a pass here as weak evidence and widen the panel.")
        else:
            print("\n  Panel spread is wide enough to exercise the wiring.")

    if not report.passed:
        print("\nGate FAILED. This is a defect in BT4's adapter, not a reason to loosen")
        print("the tolerance -- see docs/DESIGN_splice_cnn_calibration.md step A5 for the")
        print("deviation-size diagnosis table.")
        return 1

    if args.attest_out:
        attestation = attest_backend(
            "pangolin", report, dict(PINNED_WEIGHT_SHA256), bt4_version=bt4.__version__
        )
        out = Path(args.attest_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(attestation.to_dict(), indent=2, sort_keys=True) + "\n"
        out.write_text(text, encoding="utf-8")
        print(f"\nwrote attestation to {out}")
        print(f"content_hash: {attestation.content_hash()}")
        print("\nThis file carries only license-clean scalars plus the public weight")
        print("SHA-256s -- never a raw per-position score. It is safe to commit.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
