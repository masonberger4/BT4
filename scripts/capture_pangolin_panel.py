"""Capture upstream Pangolin's own per-position scores for a fidelity panel.

Step **A4b** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md).

**This script must never import ``bt4``, and enforces that on itself.** The whole
point of the integration-fidelity gate is to check BT4's adapter against an
*independent* source of truth. Capturing the expected scores with the very adapter
under test would make the gate pass no matter what the adapter did -- a test that
cannot fail. So this file drives the user's installed ``pangolin`` package
directly, and :func:`_assert_bt4_not_imported` fails loudly if ``bt4`` ever
appears in ``sys.modules``.

**It mirrors the Pangolin CLI, not ``custom_usage.py``.** Upstream ships two
disjoint model sets and its own two entry points disagree about which to use:

* ``pangolin/pangolin.py`` (the **CLI**) loads ``final.{1,2,3}.{0,2,4,6}.3.v2`` --
  12 fine-tuned models -- and reads output channels ``[1, 4, 7, 10]``;
* ``scripts/custom_usage.py`` (the **example**) loads ``final.{1..5}.{i}.3`` -- the
  older base models -- via a different channel map.

BT4's adapter tracks the CLI, because that is the model users actually run. A
capture built by copying ``custom_usage.py`` would produce different numbers and
fail the gate for a reason that has nothing to do with BT4. This script therefore
reimplements the **CLI's** path: v2 weights, folds 1-3, channels [1,4,7,10],
averaged across folds and then across the four tissue heads.

Weights are located from ``--model-dir``, then ``$BT4_PANGOLIN_MODEL_DIR``, then
the ``models`` directory inside the installed ``pangolin`` package -- the same
order BT4's adapter uses, so a maintainer who never needed to set the variable
does not have to discover it here. Resolution deliberately imports ``pangolin``
rather than ``pkg_resources`` (removed in setuptools >= 81, which
``pip install torch`` installs). Importing ``pangolin`` is fine; importing ``bt4``
is what would compromise the capture.

The captured scores are **GPL-derived model outputs**. Write them outside the BT4
repository and never commit them; only the license-clean scalars of a passing gate
(a :class:`~bt4.biomodels.splice.FidelityAttestation`) may be committed.

Runtime: each sequence is padded to length + 10,000 and pushed through 12 models
on CPU, so expect on the order of a minute per sequence. Progress is printed.

Run it::

    python scripts/capture_pangolin_panel.py \\
        --panel panel_sequences.json --out expected_pangolin.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

__all__ = ["capture_panel", "main", "resolve_model_dir", "site_scores"]

PANGOLIN_FLANK = 5000
"""``N`` bases padded each side; Pangolin returns scores for the middle N-10000."""

TISSUE_WEIGHT_INDEX = (0, 2, 4, 6)
"""The CLI's P(splice) weight indices: Heart, Liver, Brain, Testis."""

TISSUE_CHANNEL = {0: 1, 2: 4, 4: 7, 6: 10}
"""The CLI's output channel per tissue head (its ``[1, 4, 7, 10]``)."""

CV_FOLDS = (1, 2, 3)
"""The CLI's three production folds (``range(1, 4)``)."""

WEIGHT_TEMPLATE = "final.{fold}.{index}.3.v2"


def _assert_bt4_not_imported() -> None:
    """Fail loudly if ``bt4`` is importable into this process.

    An independence guard, not a style rule: if the expectations were produced
    with BT4 in the loop, the gate they feed proves nothing.
    """
    leaked = sorted(m for m in sys.modules if m == "bt4" or m.startswith("bt4."))
    if leaked:
        raise RuntimeError(
            "bt4 is imported in this process -- the capture must be independent of "
            f"the adapter under test (found: {leaked[:5]})"
        )


def resolve_model_dir(explicit: str | None = None) -> Path | None:
    """Resolve the Pangolin weights directory, or ``None``.

    Order: ``explicit`` -> ``$BT4_PANGOLIN_MODEL_DIR`` -> the ``models`` directory
    of the installed ``pangolin`` package. This mirrors
    :meth:`PangolinSplicePredictor.weights_dir` **without** importing ``bt4``:
    the independence rule forbids the adapter under test, not the upstream package
    whose scores are being captured.

    Never raises -- an unresolvable location is reported by the caller.
    """
    import os

    for candidate in (explicit, os.environ.get("BT4_PANGOLIN_MODEL_DIR")):
        if candidate:
            path = Path(candidate)
            return path if path.is_dir() else None
    try:
        import pangolin
    except ImportError:
        return None
    pkg_file = getattr(pangolin, "__file__", None)
    if not pkg_file:
        return None
    models = Path(pkg_file).resolve().parent / "models"
    return models if models.is_dir() else None


def _one_hot(seq: str):
    """One-hot encode to Pangolin's channel-major layout, using its own IN_MAP."""
    import numpy as np

    in_map = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    s = (
        seq.upper()
        .replace("A", "1")
        .replace("C", "2")
        .replace("G", "3")
        .replace("T", "4")
        .replace("N", "0")
    )
    return in_map[np.asarray(list(map(int, s))).astype("int8")]


def _load_models(model_dir: Path):
    """Load the CLI's 12 fine-tuned P(splice) models, keyed by tissue weight index."""
    import torch
    from pangolin.model import AR, L, Pangolin, W

    models: dict[int, list] = {}
    for index in TISSUE_WEIGHT_INDEX:
        fold_models = []
        for fold in CV_FOLDS:
            path = model_dir / WEIGHT_TEMPLATE.format(fold=fold, index=index)
            if not path.is_file():
                raise FileNotFoundError(f"missing Pangolin weight file: {path}")
            model = Pangolin(L, W, AR)
            state = torch.load(path, map_location=torch.device("cpu"))
            model.load_state_dict(state)
            model.eval()
            fold_models.append(model)
        models[index] = fold_models
    return models


def site_scores(seq: str, models: dict[int, list]) -> list[float]:
    """Return upstream Pangolin's combined per-position ``P(splice)`` for ``seq``.

    Averages the three production folds within each tissue head, then averages the
    four tissue heads -- the tissue-agnostic readout appropriate for a coding
    sequence redesign, and the same reduction BT4's adapter performs.
    """
    import numpy as np
    import torch

    padded = "N" * PANGOLIN_FLANK + seq.upper() + "N" * PANGOLIN_FLANK
    x = torch.from_numpy(np.expand_dims(_one_hot(padded).T, axis=0)).float()

    per_tissue = []
    with torch.no_grad():
        for index in TISSUE_WEIGHT_INDEX:
            channel = TISSUE_CHANNEL[index]
            folds = [m(x)[0, channel, :].numpy() for m in models[index]]
            per_tissue.append(np.mean(folds, axis=0))
    combined = np.mean(per_tissue, axis=0)

    if len(combined) != len(seq):
        raise RuntimeError(
            f"Pangolin returned {len(combined)} scores for a {len(seq)} nt sequence; "
            "the context crop differs from the expected 10,000"
        )
    return [float(v) for v in combined]


def capture_panel(panel_path: Path, model_dir: Path, *, limit: int | None = None) -> dict:
    """Score every sequence in the panel and return the capture payload."""
    _assert_bt4_not_imported()

    payload = json.loads(panel_path.read_text(encoding="utf-8"))
    sequences = payload["sequences"]
    if limit is not None:
        sequences = sequences[:limit]

    print(f"loading 12 Pangolin models from {model_dir} ...", flush=True)
    models = _load_models(model_dir)
    print(f"scoring {len(sequences)} sequences (this is slow on CPU)\n", flush=True)

    captured = []
    for i, entry in enumerate(sequences, start=1):
        seq = entry["sequence"]
        print(f"  [{i}/{len(sequences)}] {entry['id']:22} {len(seq):>5} nt ...", end="", flush=True)
        scores = site_scores(seq, models)
        peak = max(scores) if scores else 0.0
        # Report WHERE the peak is, not just how high. A peak pinned to position 0
        # or the last base across every sequence is the documented N-padding
        # boundary artifact (OpenSpliceAI, eLife 2025) rather than sequence-driven
        # signal -- the adapter fills a ~10 kb window, so a short CDS is mostly
        # padding. That does not invalidate a fidelity gate (both sides reproduce
        # the artifact identically), but it does mean the panel is exercising the
        # boundary rather than the biology, which is worth seeing.
        at = scores.index(peak) if scores else 0
        edge = "" if 3 < at < len(scores) - 4 else "  <-- at sequence edge"
        print(f" peak P(splice)={peak:.4f} at {at}/{len(scores)}{edge}", flush=True)
        captured.append({"id": entry["id"], "sequence": seq, "expected_site_scores": scores})

    return {
        "schema_version": 1,
        "backend": "pangolin",
        "panel_content_hash": payload.get("content_hash"),
        "n_cases": len(captured),
        "cases": captured,
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Capture upstream Pangolin scores for a panel and write them as JSON."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--panel", required=True, help="Panel JSON from make_splice_panel.py.")
    parser.add_argument("--out", required=True, help="Where to write the captured scores.")
    parser.add_argument(
        "--model-dir",
        default=None,
        help="Pangolin weights directory (default: $BT4_PANGOLIN_MODEL_DIR, "
        "then the installed pangolin package).",
    )
    parser.add_argument("--limit", type=int, default=None, help="Only score the first N sequences.")
    args = parser.parse_args(argv)

    model_dir = resolve_model_dir(args.model_dir)
    if model_dir is None:
        parser.error(
            "no Pangolin weights directory resolved: pass --model-dir, set "
            "$BT4_PANGOLIN_MODEL_DIR, or install the 'pangolin' package"
        )

    payload = capture_panel(Path(args.panel), model_dir, limit=args.limit)
    Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    peaks = [max(c["expected_site_scores"]) for c in payload["cases"] if c["expected_site_scores"]]
    print(f"\nwrote {args.out} ({payload['n_cases']} cases)")
    if peaks:
        print(f"peak P(splice) across panel: min={min(peaks):.4f} max={max(peaks):.4f}")
    print("\nThese scores are GPL-derived model outputs: keep them OUT of the BT4 repo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
