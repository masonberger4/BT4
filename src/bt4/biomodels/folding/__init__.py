"""5' mRNA folding free-energy models behind the ``FoldingModel`` contract.

This package provides the folding half of BT4's non-local biology (CLAUDE.md
sections 6 and 4.3): a swappable :class:`FoldingModel` contract, a calibrated
ViennaRNA backend (:class:`ViennaFoldingModel`), and an honestly-labeled,
dependency-free fallback (:class:`BaselinePairingProxyModel`).

Use :func:`default` to obtain a working model: it returns the calibrated
ViennaRNA backend when its bindings import, otherwise the labeled baseline. It
never crashes and never imports ViennaRNA at package load.

This package depends only on :mod:`bt4.domain`, the standard library, and --
lazily, inside methods -- the optional ViennaRNA binding.
"""

from __future__ import annotations

from bt4.biomodels.folding.base import (
    DEFAULT_FIVE_PRIME_WINDOW,
    DEFAULT_LEADER_WINDOW,
    FoldingModel,
    FoldingResult,
    five_prime_window,
    junction_window,
)
from bt4.biomodels.folding.baseline import BaselinePairingProxyModel
from bt4.biomodels.folding.vienna import ViennaFoldingModel

__all__ = [
    "DEFAULT_FIVE_PRIME_WINDOW",
    "DEFAULT_LEADER_WINDOW",
    "BaselinePairingProxyModel",
    "FoldingModel",
    "FoldingResult",
    "ViennaFoldingModel",
    "default",
    "five_prime_window",
    "junction_window",
]


def default() -> FoldingModel:
    """Return the best available folding model, never crashing.

    Returns:
        A calibrated :class:`ViennaFoldingModel` when the ViennaRNA bindings can
        be imported, otherwise an uncalibrated
        :class:`BaselinePairingProxyModel`. The baseline is honestly labeled
        (``calibrated is False``) so its numbers are never mistaken for real
        deltaG (CLAUDE.md sections 4.3 and 10.6).
    """
    try:
        if ViennaFoldingModel.available():
            return ViennaFoldingModel()
    except Exception:
        pass
    return BaselinePairingProxyModel()
