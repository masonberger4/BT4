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

from bt4.biomodels.expression.base import ExpressionPredictor, ExpressionResult
from bt4.biomodels.expression.baseline import NullExpressionModel
from bt4.biomodels.expression.gate import (
    ExpressionEvalCase,
    ExpressionGateReport,
    run_expression_gate,
    verify_expression_gate,
)

__all__ = [
    "ExpressionEvalCase",
    "ExpressionGateReport",
    "ExpressionPredictor",
    "ExpressionResult",
    "NullExpressionModel",
    "default",
    "run_expression_gate",
    "verify_expression_gate",
]


def default() -> ExpressionPredictor:
    """Return the best available expression predictor, never crashing.

    Returns:
        The neutral :class:`NullExpressionModel` placeholder (``calibrated is
        False``) until a validated, hash-pinned expression head ships and is
        selected here ahead of it. The placeholder claims nothing and must not be
        read as an expression prediction (CLAUDE.md §6, §10.6).
    """
    return NullExpressionModel()
