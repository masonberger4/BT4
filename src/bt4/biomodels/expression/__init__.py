"""Learned-expression models behind the ``ExpressionPredictor`` contract.

This package is the **scaffold** for BT4's Phase 4 learned-expression head
(CLAUDE.md §6, §9): a swappable :class:`ExpressionPredictor` contract and an
honestly-labeled neutral placeholder (:class:`NullExpressionModel`). No
calibrated model ships yet, so :func:`default` returns the placeholder -- which
scores every sequence ``0.0`` and reports ``calibrated is False``.

A future validated head (trained on real MPRA / ribosome-load data, hash-pinned,
passing a held-out gate reporting calibration and uncertainty per CLAUDE.md §8)
would live beside this module behind a lazy import and be selected by
:func:`default` ahead of the placeholder -- exactly as
:class:`~bt4.biomodels.folding.ViennaFoldingModel` is preferred over the folding
baseline. Until then no such stub exists: a fake-weight expression model would be
the dishonest placeholder CLAUDE.md §10.6 forbids.

This package depends only on the standard library.
"""

from __future__ import annotations

from bt4.biomodels.expression.base import (
    BatchExpressionPredictor,
    ExpressionPredictor,
    ExpressionResult,
)
from bt4.biomodels.expression.baseline import NullExpressionModel
from bt4.biomodels.expression.gate import (
    ExpressionEvalCase,
    ExpressionGateReport,
    run_expression_gate,
    verify_expression_gate,
)
from bt4.biomodels.expression.ribonn import (
    PINNED_WEIGHT_SHA256,
    RiboNNExpressionModel,
    RiboNNFoldPrediction,
    load_pinned_sha256,
)

__all__ = [
    "PINNED_WEIGHT_SHA256",
    "BatchExpressionPredictor",
    "ExpressionEvalCase",
    "ExpressionGateReport",
    "ExpressionPredictor",
    "ExpressionResult",
    "NullExpressionModel",
    "RiboNNExpressionModel",
    "RiboNNFoldPrediction",
    "available_backends",
    "default",
    "load_pinned_sha256",
    "resolve_backend",
    "run_expression_gate",
    "verify_expression_gate",
]

# Public backend registry (CLAUDE.md §10.9: registries are public, never private
# symbols crossing a layer). Aliases map to a canonical key.
_BACKENDS = {
    "null": "null",
    "placeholder": "null",
    "none": "null",
    "ribonn": "ribonn",
}


def default() -> ExpressionPredictor:
    """Return the best available expression predictor, never crashing.

    Returns:
        The neutral :class:`NullExpressionModel` placeholder (``calibrated is
        False``) until a validated, hash-pinned expression head ships and is
        selected here ahead of it. The placeholder claims nothing and must not be
        read as an expression prediction (CLAUDE.md §6, §10.6).
    """
    return NullExpressionModel()


def available_backends() -> tuple[str, ...]:
    """Return the expression backends that can actually run on this machine.

    ``"null"`` (the neutral placeholder) is always available. ``"ribonn"`` is
    listed only when :meth:`RiboNNExpressionModel.available` reports that the
    user's own RiboNN checkout, its ``<species>`` weight directory, and the heavy
    deps all resolve -- so a frontend can offer the wrapped head *only* when
    selecting it would work, and explain its absence otherwise.

    Availability is emphatically **not** calibration: a listed RiboNN is still
    ``calibrated is False`` until it passes the CDS-variant acceptance gate
    (CLAUDE.md §6/§10.6), so an uncalibrated head must never be shown as a
    ranking. Never raises.

    Returns:
        Backend names, always beginning with ``"null"``.
    """
    names = ["null"]
    try:
        if RiboNNExpressionModel().available():
            names.append("ribonn")
    except (OSError, ValueError, ImportError):
        # Probing a user-supplied path / optional dependency must never break a
        # caller that is only listing its options: an unreadable $BT4_RIBONN_DIR
        # or a broken torch install means "not available here", not an error.
        pass
    return tuple(names)


def resolve_backend(
    name: str,
    *,
    species: str = "human",
    utr5: str = "",
    utr3: str = "",
    top_k: int = 5,
    batch_size: int = 64,
    num_workers: int = 0,
    cell_types: tuple[str, ...] = (),
) -> ExpressionPredictor:
    """Construct an expression backend by name (the mirror of the splice resolver).

    Args:
        name: ``"null"`` (aliases ``"placeholder"`` / ``"none"``) for the neutral
            placeholder, or ``"ribonn"`` for the wrapped RiboNN head.
        species: RiboNN weight set -- ``"human"`` or ``"mouse"``. Ignored by
            ``"null"``.
        utr5: Fixed 5' UTR context held constant while the CDS varies. RiboNN
            **requires** it non-empty to score (its loader reads an all-empty UTR
            column as NaN, and the UTRs carry most of its signal).
        utr3: Fixed 3' UTR context, as ``utr5``.
        top_k: Number of RiboNN cross-validation runs to ensemble.
        batch_size: RiboNN inference batch size (memory/speed only -- it cannot
            change a score, since RiboNN pads to a fixed width and does not shuffle
            when predicting). Defaults below RiboNN's own 1024, which OOMs a CPU box.
        cell_types: Which RiboNN per-cell-type outputs to average. Empty (default)
            averages all of them; naming one (e.g. ``("HEK293T",)``) is required to
            compare honestly against a single-cell-line measurement. An unmatched name
            raises when scoring.
        num_workers: RiboNN dataloader worker count. Defaults to ``0``, which is
            required wherever the multiprocessing start method is *spawn* (Windows,
            macOS) because the adapter scores from a mutated ``sys.path`` and a
            temporary working directory that a spawned worker does not inherit.

    Returns:
        An :class:`ExpressionPredictor`. Constructing a RiboNN backend does not
        load weights or import torch -- those happen lazily on the first score --
        so this stays cheap and cannot fail on a missing checkout.

    Raises:
        ValueError: If ``name`` is not a known backend, or RiboNN rejects
            ``species`` / ``top_k``.
    """
    key = _BACKENDS.get(name.strip().lower())
    if key is None:
        raise ValueError(
            f"unknown expression backend {name!r}; choose from {sorted(set(_BACKENDS))}"
        )
    if key == "ribonn":
        return RiboNNExpressionModel(
            species=species,
            top_k=top_k,
            utr5=utr5,
            utr3=utr3,
            batch_size=batch_size,
            num_workers=num_workers,
            cell_types=cell_types,
        )
    return NullExpressionModel()
