"""Capture upstream SpliceAI's own per-position scores for a fidelity panel.

Step **A4b** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md),
the SpliceAI sibling of ``scripts/capture_pangolin_panel.py``.

**This script must never import ``bt4``, and enforces that on itself.** The whole
point of the integration-fidelity gate is to check BT4's adapter against an
*independent* source of truth. Capturing the expected scores with the very adapter
under test would make the gate pass no matter what the adapter did -- a test that
cannot fail. So this file drives the user's installed ``spliceai`` package and
Keras directly, and :func:`_assert_bt4_not_imported` fails loudly if ``bt4`` ever
appears in ``sys.modules``.

**It uses upstream's own one-hot encoder, deliberately.** Unlike Pangolin -- whose
CLI encodes inline, forcing the Pangolin capture to reimplement it -- SpliceAI
ships ``spliceai.utils.one_hot_encode`` as a reusable function, so this script
imports it. That is strictly stronger evidence: BT4's adapter re-derives the
encoding in pure Python (``_one_hot_rows``, position-major ``[L][4]``), and a
transposed layout, a wrong base order, or a mishandled ``N`` would show up as a
gate failure rather than being reproduced identically on both sides.

There is **no fallback encoder**. If ``one_hot_encode`` cannot run, the script
refuses instead of substituting its own: a capture that quietly stops being
independent is worse than one that does not run. The usual cause is NumPy 2 --
``spliceai/utils.py`` still calls the long-removed ``np.fromstring`` -- so pin
``numpy<2`` (see the runbook's A1 install table).

Weights are located from ``--model-dir``, then ``$BT4_SPLICEAI_MODEL_DIR``, then
the ``models`` directory inside the installed ``spliceai`` package -- the same
order BT4's adapter uses, so a maintainer who never needed to set the variable
does not have to discover it here.

The captured scores are **CC BY-NC 4.0 model outputs**. Write them outside the BT4
repository and never commit them; only the license-clean scalars of a passing gate
(a :class:`~bt4.biomodels.splice.FidelityAttestation`) may be committed. That
noncommercial term is also why SpliceAI, once attested, stays usable only within
BT4's noncommercial scope.

Runtime: each sequence is padded to length + 10,000 and pushed through 5 Keras
models on CPU, so expect roughly a minute per sequence. Progress is printed.

Run it (in the SpliceAI environment -- TensorFlow 2.15 and PyTorch do not
coexist, so this is a different virtualenv from the Pangolin capture)::

    python scripts/capture_spliceai_panel.py \\
        --panel panel_sequences.json --out expected_spliceai.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

__all__ = ["acceptor_donor_scores", "capture_panel", "main", "resolve_model_dir"]

SPLICEAI_FLANK = 5000
"""``N`` bases padded each side; SpliceAI returns scores for the middle N-10000."""

ENSEMBLE_SIZE = 5
"""SpliceAI's published ensemble: ``spliceai1.h5`` .. ``spliceai5.h5``, averaged."""

ACCEPTOR_CHANNEL = 1
"""Output channel holding the per-position acceptor probability (``y[:, 1]``)."""

DONOR_CHANNEL = 2
"""Output channel holding the per-position donor probability (``y[:, 2]``).

Channel 0 is the *null* class of SpliceAI's 3-way softmax and is not used: unlike
Pangolin's single combined ``P(splice)``, SpliceAI genuinely separates acceptor
from donor, which is why this capture records two tracks and Pangolin's records one.
"""

WEIGHT_TEMPLATE = "spliceai{index}.h5"


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
    """Resolve the SpliceAI weights directory, or ``None``.

    Order: ``explicit`` -> ``$BT4_SPLICEAI_MODEL_DIR`` -> the ``models`` directory
    of the installed ``spliceai`` package. Mirrors
    :meth:`SpliceAiSplicePredictor.weights_dir` **without** importing ``bt4``:
    the independence rule forbids the adapter under test, not the upstream package
    whose scores are being captured.

    Uses ``find_spec`` rather than importing ``spliceai`` here, so merely locating
    the weights does not run the package's ``__init__`` (which installs a SIGINT
    handler). The import happens later, when scoring actually needs it.

    Never raises -- an unresolvable location is reported by the caller.
    """
    import importlib.util
    import os

    for candidate in (explicit, os.environ.get("BT4_SPLICEAI_MODEL_DIR")):
        if candidate:
            path = Path(candidate)
            return path if path.is_dir() else None
    try:
        spec = importlib.util.find_spec("spliceai")
    except (ImportError, ValueError):
        return None
    if spec is None or spec.origin is None:
        return None
    models = Path(spec.origin).resolve().parent / "models"
    return models if models.is_dir() else None


def _upstream_one_hot() -> Any:
    """Return upstream SpliceAI's own ``one_hot_encode``, or refuse.

    Deliberately has no fallback. Re-deriving the encoding here would make the
    capture agree with BT4's adapter by construction on exactly the axis the gate
    is meant to test.

    Raises:
        RuntimeError: If ``spliceai.utils`` cannot be imported, naming the usual
            cause.
    """
    try:
        from spliceai.utils import one_hot_encode
    except ImportError as exc:  # pragma: no cover - needs the licensed install
        raise RuntimeError(
            "cannot import spliceai.utils.one_hot_encode. Install the 'spliceai' "
            "package (pip install spliceai==1.3.1) and pin 'numpy<2' -- "
            "spliceai/utils.py still calls np.fromstring, removed in NumPy 2. This "
            "script will NOT substitute its own encoder: capturing with a "
            "re-derived encoding would make the gate agree with BT4 by "
            f"construction. Underlying error: {exc}"
        ) from exc
    return one_hot_encode


def _load_ensemble(model_dir: Path) -> list[Any]:
    """Load SpliceAI's five published Keras models from ``model_dir``."""
    for module_name in ("tensorflow.keras.models", "keras.models", "tf_keras.models"):
        try:
            module = __import__(module_name, fromlist=["load_model"])
        except ImportError:
            continue
        load_model = module.load_model
        break
    else:  # pragma: no cover - needs the licensed install
        raise RuntimeError(
            "no Keras entry point importable. Install tensorflow==2.15.* (TF >= 2.16 "
            "defaults to Keras 3, which cannot load these 2019 Keras-2 .h5 graphs), "
            "or add tf_keras and set TF_USE_LEGACY_KERAS=1"
        )

    models = []
    for index in range(1, ENSEMBLE_SIZE + 1):
        path = model_dir / WEIGHT_TEMPLATE.format(index=index)
        if not path.is_file():
            raise FileNotFoundError(f"missing SpliceAI weight file: {path}")
        models.append(load_model(str(path)))
    return models


def acceptor_donor_scores(
    seq: str, models: list[Any], one_hot_encode: Any
) -> tuple[list[float], list[float]]:
    """Return upstream SpliceAI's per-position ``(acceptor, donor)`` for ``seq``.

    Follows the official README's custom-sequence recipe exactly: pad 5,000 ``N``
    each side, encode with upstream's own encoder, average the five models'
    predictions, and read channels 1 and 2.

    ``verbose=0`` only suppresses Keras' progress bar and cannot change a score.
    """
    import numpy as np

    padded = "N" * SPLICEAI_FLANK + seq.upper() + "N" * SPLICEAI_FLANK
    x = one_hot_encode(padded)[None, :]
    y = np.mean([model.predict(x, verbose=0) for model in models], axis=0)

    acceptor = [float(v) for v in y[0, :, ACCEPTOR_CHANNEL]]
    donor = [float(v) for v in y[0, :, DONOR_CHANNEL]]
    if len(acceptor) != len(seq) or len(donor) != len(seq):
        raise RuntimeError(
            f"SpliceAI returned {len(acceptor)} scores for a {len(seq)} nt sequence; "
            "the context crop differs from the expected 10,000"
        )
    return acceptor, donor


def capture_panel(panel_path: Path, model_dir: Path, *, limit: int | None = None) -> dict:
    """Score every sequence in the panel and return the capture payload."""
    _assert_bt4_not_imported()

    payload = json.loads(panel_path.read_text(encoding="utf-8"))
    sequences = payload["sequences"]
    if limit is not None:
        sequences = sequences[:limit]

    one_hot_encode = _upstream_one_hot()
    print(f"loading 5 SpliceAI models from {model_dir} ...", flush=True)
    models = _load_ensemble(model_dir)
    print(f"scoring {len(sequences)} sequences (this is slow on CPU)\n", flush=True)

    captured = []
    for i, entry in enumerate(sequences, start=1):
        seq = entry["sequence"]
        print(f"  [{i}/{len(sequences)}] {entry['id']:22} {len(seq):>5} nt ...", end="", flush=True)
        acceptor, donor = acceptor_donor_scores(seq, models, one_hot_encode)
        # Report WHERE each peak is, not just how high. A peak pinned to the first
        # or last base across every sequence is the documented N-padding boundary
        # artifact (OpenSpliceAI, eLife 2025) rather than sequence-driven signal --
        # the adapter fills a ~10 kb window, so a short CDS is mostly padding. That
        # does not invalidate a fidelity gate (both sides reproduce the artifact
        # identically), but it does mean the panel is exercising the boundary
        # rather than the biology, which is worth seeing.
        parts = []
        for label, track in (("A", acceptor), ("D", donor)):
            peak = max(track) if track else 0.0
            at = track.index(peak) if track else 0
            edge = "*" if track and not 3 < at < len(track) - 4 else ""
            parts.append(f"{label}={peak:.4f}@{at}{edge}")
        print(f" {' '.join(parts)}  (*=at edge)", flush=True)
        captured.append(
            {
                "id": entry["id"],
                "sequence": seq,
                "expected_acceptor": acceptor,
                "expected_donor": donor,
            }
        )

    return {
        "schema_version": 1,
        "backend": "spliceai",
        "panel_content_hash": payload.get("content_hash"),
        "n_cases": len(captured),
        "cases": captured,
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Capture upstream SpliceAI scores for a panel and write them as JSON."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--panel", required=True, help="Panel JSON from make_splice_panel.py.")
    parser.add_argument("--out", required=True, help="Where to write the captured scores.")
    parser.add_argument(
        "--model-dir",
        default=None,
        help="SpliceAI weights directory (default: $BT4_SPLICEAI_MODEL_DIR, "
        "then the installed spliceai package).",
    )
    parser.add_argument("--limit", type=int, default=None, help="Only score the first N sequences.")
    args = parser.parse_args(argv)

    model_dir = resolve_model_dir(args.model_dir)
    if model_dir is None:
        parser.error(
            "no SpliceAI weights directory resolved: pass --model-dir, set "
            "$BT4_SPLICEAI_MODEL_DIR, or install the 'spliceai' package"
        )

    payload = capture_panel(Path(args.panel), model_dir, limit=args.limit)
    Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    peaks = [
        max((*case["expected_acceptor"], *case["expected_donor"]))
        for case in payload["cases"]
        if case["expected_acceptor"] or case["expected_donor"]
    ]
    print(f"\nwrote {args.out} ({payload['n_cases']} cases)")
    if peaks:
        print(f"peak site probability across panel: min={min(peaks):.4f} max={max(peaks):.4f}")
    print("\nThese scores are CC BY-NC 4.0 model outputs: keep them OUT of the BT4 repo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
