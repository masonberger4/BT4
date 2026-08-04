"""Two-backend splice comparison harness -- agreement as an uncertainty signal.

CLAUDE.md section 6 makes running more than one splice backend and **reporting
their agreement** a first-class uncertainty signal, not redundancy. This script
scores a panel of candidate coding sequences against a reference with every
*available* :class:`~bt4.biomodels.splice.base.SplicePredictor` backend -- always
the honestly-labeled PWM baseline, plus the wrapped **Pangolin** and **SpliceAI**
CNNs when they are installed and their weights resolve -- and prints each
backend's Delta-splicing ranking together with the cross-backend agreement report
(:func:`bt4.biomodels.splice.backend_agreement`). With both CNNs installed, this
is agreement between two real, independently-trained splice models -- the §6
"run both, report agreement" uncertainty signal.

Honest framing (read before reading the numbers):

* **The baseline is uncalibrated.** The PWM baseline reports
  ``calibrated is False`` -- its per-position scores are consensus
  pseudo-probabilities, not calibrated splice probabilities. Pangolin, until it
  passes its integration-fidelity gate, also reports ``calibrated is False`` (it
  is a faithfully-wrapped published model, but wrapping is not yet verified). So
  every number here is a *predicted* cryptic-splice-risk prior, never a validated
  expression claim (the CAI-as-weak-proxy caution, CLAUDE.md sections 1 and 6).
* **Delta-splicing is larger-is-better.** A positive value means the candidate
  carries *less* predicted splice risk than the reference; negative means more.
* **This reports, it does not judge.** Where the backends agree on the ranking
  (high Spearman) that is corroboration; where they disagree, that disagreement
  is surfaced as uncertainty, not resolved.

Neither CNN is bundled with BT4 (Pangolin is GPL-3.0; SpliceAI code is PolyForm
Strict 1.0.0 and its weights are CC BY-NC 4.0, noncommercial). Install them
separately (https://github.com/tkzeng/Pangolin, https://github.com/Illumina/SpliceAI,
each of which provides its weights) and BT4 picks them up. Without them the
harness runs with the baseline alone and says so.

Being a standalone script (not part of the ``bt4`` package import graph), it runs
through the stable :mod:`bt4.api` / :mod:`bt4.io` surfaces and reaches into
:mod:`bt4.biomodels.splice` only for the backends the public API does not
surface. Run it directly::

    # demo: an optimized reference vs a sampled library, ranked by splice risk
    python scripts/compare_splice_backends.py

    # your own sequences (first FASTA record = reference, rest = candidates)
    python scripts/compare_splice_backends.py --fasta my_variants.fasta --json
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence

from bt4 import api
from bt4.biomodels.splice import (
    ConsensusPwmSplicePredictor,
    PangolinSplicePredictor,
    SpliceAiSplicePredictor,
    SplicePredictor,
    backend_agreement,
)
from bt4.io import read_fasta

__all__ = ["available_backends", "build_report", "main"]

# A short demo protein used when no --fasta panel is supplied. The reference is
# BT4's deterministic optimum; the candidates are a sampled library (Phase 5),
# so the harness answers "do the backends agree on which library members lower
# splice risk relative to the optimum?".
_DEMO_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"


def available_backends() -> tuple[list[SplicePredictor], list[str]]:
    """Return the runnable splice backends and notes on any that were skipped.

    The PWM baseline is always available. The Pangolin and SpliceAI CNN backends
    are each included only when their ``available()`` reports true (the heavy dep
    installed and their weights resolvable). When BOTH CNNs are present, agreement
    between two real, independently-trained CNNs becomes reachable -- the §6
    "run both, report agreement" uncertainty signal.

    Returns:
        A ``(backends, notes)`` pair: the usable backends and human-readable
        notes about any unavailable ones.
    """
    backends: list[SplicePredictor] = [ConsensusPwmSplicePredictor()]
    notes: list[str] = []

    pangolin = PangolinSplicePredictor()
    if pangolin.available():
        backends.append(pangolin)
        if not pangolin.calibrated:
            notes.append(
                "Pangolin is available but reports calibrated=False (its "
                "integration-fidelity gate has not been run); its scores are a "
                "prior, not a validated result."
            )
    else:
        notes.append(
            "Pangolin backend unavailable (install torch + the GPL 'pangolin' "
            "package from https://github.com/tkzeng/Pangolin)."
        )

    spliceai = SpliceAiSplicePredictor()
    if spliceai.available():
        backends.append(spliceai)
        if not spliceai.calibrated:
            notes.append(
                "SpliceAI is available but reports calibrated=False (its "
                "integration-fidelity gate has not been run); its scores are a "
                "prior, not a validated result. SpliceAI weights are CC BY-NC 4.0 "
                "(noncommercial)."
            )
    else:
        notes.append(
            "SpliceAI backend unavailable (install tensorflow + the CC BY-NC "
            "'spliceai' package from https://github.com/Illumina/SpliceAI)."
        )

    if len(backends) == 1:
        notes.append("Only the PWM baseline is available; showing it alone.")
    return backends, notes


def build_report(
    reference: str,
    candidates: Sequence[str],
) -> dict[str, object]:
    """Build the JSON-ready comparison report for ``candidates`` vs ``reference``.

    Args:
        reference: The reference coding sequence.
        candidates: Candidate coding sequences to rank by splice risk.

    Returns:
        A dictionary with the backends used, per-backend Delta-splicing per
        candidate, the pairwise rank correlations, the sign-agreement fraction,
        and any availability notes.
    """
    backends, notes = available_backends()
    report = backend_agreement(backends, list(candidates), reference)
    return {
        "reference_len_nt": len(reference),
        "n_candidates": report.n_candidates,
        "backends": [
            {"name": b.name, "calibrated": b.calibrated} for b in backends
        ],
        "delta_splicing_by_backend": {
            name: list(vals) for name, vals in report.delta_by_backend.items()
        },
        "rank_correlations": [
            {"backends": list(pair), "spearman": corr}
            for pair, corr in report.rank_correlations.items()
        ],
        "sign_agreement": report.sign_agreement,
        "notes": notes,
    }


def _print_report(report: dict[str, object]) -> None:
    """Print a human-readable rendering of :func:`build_report`'s output."""
    backends = report["backends"]
    assert isinstance(backends, list)
    deltas = report["delta_splicing_by_backend"]
    assert isinstance(deltas, dict)

    print("=== Splice-backend comparison (Delta-splicing; larger = less risk) ===")
    print(f"reference: {report['reference_len_nt']} nt, {report['n_candidates']} candidates")
    print()
    names = [b["name"] for b in backends]
    header = "cand  " + "  ".join(f"{name[:22]:>22}" for name in names)
    print(header)
    print("-" * len(header))
    n = int(report["n_candidates"])  # type: ignore[arg-type]
    for idx in range(n):
        row = f"{idx:>4}  " + "  ".join(f"{deltas[name][idx]:>22.4f}" for name in names)
        print(row)
    print()
    for b in backends:
        flag = "calibrated" if b["calibrated"] else "UNCALIBRATED"
        print(f"  backend {b['name']}: {flag}")
    corrs = report["rank_correlations"]
    assert isinstance(corrs, list)
    if corrs:
        print()
        print("  pairwise Spearman rank agreement of Delta-splicing:")
        for entry in corrs:
            a, b = entry["backends"]
            print(f"    {a} vs {b}: rho = {entry['spearman']:.4f}")
    print()
    print(f"  sign agreement across backends: {report['sign_agreement']:.2%}")
    notes = report["notes"]
    assert isinstance(notes, list)
    for note in notes:
        print(f"  note: {note}")


def _demo_panel() -> tuple[str, list[str]]:
    """Return ``(reference, candidates)`` for the no-argument demo.

    Reference is BT4's deterministic optimum for :data:`_DEMO_PROTEIN`; candidates
    are a small sampled library (deterministic from the seed).
    """
    reference = api.optimize(_DEMO_PROTEIN).dna
    lib = api.library(_DEMO_PROTEIN, n=6, seed=0, temperature=1.5)
    candidates = [member.dna for member in lib.results]
    return reference, candidates


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point.

    Args:
        argv: Optional argument vector (defaults to ``sys.argv``).

    Returns:
        Process exit code (``0`` on success).
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        help="FASTA file; the first record is the reference, the rest are candidates.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a table.")
    args = parser.parse_args(argv)

    if args.fasta:
        records = read_fasta(args.fasta)
        if len(records) < 2:
            parser.error("--fasta must contain at least two records (reference + >=1 candidate)")
        reference = records[0][1]
        candidates = [seq for _, seq in records[1:]]
    else:
        reference, candidates = _demo_panel()

    report = build_report(reference, candidates)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_report(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
