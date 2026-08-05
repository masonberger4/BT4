"""The ``ExpressionPredictor`` contract for a learned expression head.

This is the *scaffold* for BT4's Phase 4 learned-expression model (CLAUDE.md §6,
§9). It defines the contract and the honesty flags; it deliberately ships **no**
calibrated model, because none can be honestly built yet.

Why a scaffold and not a shipped objective. Expression is the ultimate target a
back-translation tool cares about, but it is **non-local and learned**: it does
not decompose over a bounded codon-trellis context, so it is not an
:class:`~bt4.domain.contracts.ObjectiveTerm` (which must satisfy
``delta == score``). A validated expression head instead **reranks the Pareto
frontier** as a post-solve pass (see :mod:`bt4.pipeline.rerank`), never inside the
per-move optimizer loop.

The honesty rules (identical in spirit to :mod:`bt4.biomodels.splice` and
:mod:`bt4.biomodels.folding`):

* **``calibrated`` is a first-class flag, and it is earned, not assigned.** A
  backend may return ``calibrated is True`` **only** when its score is a
  validated prediction from a hash-pinned model that passed a held-out,
  homology/chromosome-grouped acceptance gate reporting PR-AUC/MCC/ECE/Brier and
  conformal coverage, *measured on data from the regime it serves* (CLAUDE.md
  §8). "Calibrated" is a property proven on a distribution, not a name.
* **No placeholder is presented as a feature (CLAUDE.md §10.6), and no
  hand-weighted composite is relabeled "expression" (§10.5).** Expression has no
  structural ground-truth anchor the way folding (ΔG) and splice (site motifs)
  do, so a hand-weighted CAI+GC+ΔG scalar dressed up as an expression prediction
  is the exact §10.5 magic-scalar trap. Until a validated head exists,
  :func:`bt4.biomodels.expression.default` returns a **neutral placeholder**
  (:class:`~bt4.biomodels.expression.baseline.NullExpressionModel`) that scores
  every sequence ``0.0`` and claims nothing.
* **An uncalibrated score never steers delivery.** The rerank hook may *annotate*
  each frontier point with an uncalibrated score for reporting, but it only
  re-picks the delivered point when the predictor is calibrated (see
  :func:`bt4.pipeline.rerank.rerank_by_expression`).

**Orientation is fixed: larger is better** -- a higher score means higher
predicted expression, matching every other BT4 objective.

This module depends only on the standard library.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

__all__ = ["BatchExpressionPredictor", "ExpressionPredictor", "ExpressionResult"]


@dataclass(frozen=True, slots=True)
class ExpressionResult:
    """A single expression score for one coding sequence (immutable).

    Attributes:
        score: The predicted-expression score, oriented **larger is better**.
            A calibrated prediction when ``calibrated`` is ``True``; an
            uncalibrated placeholder value (do not interpret) otherwise.
        model_name: The :attr:`ExpressionPredictor.name` of the producing backend.
        calibrated: Mirror of :attr:`ExpressionPredictor.calibrated` -- ``False``
            marks the score as an uncalibrated placeholder that must not be read
            as an expression prediction or used to steer delivery.
        units: A short label for the score's units (e.g. ``"log-TE"`` for a
            calibrated head, ``"none (placeholder)"`` for the null model).
    """

    score: float
    model_name: str
    calibrated: bool
    units: str


@runtime_checkable
class ExpressionPredictor(Protocol):
    """A backend that scores a coding sequence for predicted expression.

    Backends are swappable behind this contract (CLAUDE.md §4.3): consumers depend
    only on the protocol and never on a concrete backend. Every implementation
    must obey the orientation and honesty rules documented at the module level.
    """

    @property
    def name(self) -> str:
        """Stable identifier for the backend (read-only).

        Declared as a read-only property so concrete backends may be frozen
        dataclasses exposing ``name`` as a property.
        """
        ...

    @property
    def calibrated(self) -> bool:
        """Honesty flag: ``True`` only for a validated, hash-pinned, gated model.

        Returns ``True`` **only** when the score is a validated expression
        prediction whose calibration was measured on held-out data from the
        deployment regime (CLAUDE.md §6, §8, §10.6). The placeholder returns
        ``False`` so its value is never mistaken for a prediction.
        """
        ...

    def score_sequence(self, dna: str) -> ExpressionResult:
        """Return the predicted-expression score for ``dna``.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            An :class:`ExpressionResult` whose ``score`` is a calibrated
            prediction when :attr:`calibrated` is ``True`` and an uncalibrated
            placeholder otherwise.
        """
        ...


@runtime_checkable
class BatchExpressionPredictor(ExpressionPredictor, Protocol):
    """An :class:`ExpressionPredictor` that can score a whole set in one call.

    Some backends (notably :class:`~bt4.biomodels.expression.RiboNNExpressionModel`)
    have a large fixed *per-invocation* cost, so scoring a candidate set one
    sequence at a time pays it N times. Such backends additionally implement
    :meth:`score_many`, which scores the whole set in a single invocation.

    Consumers detect the capability with ``isinstance(backend,
    BatchExpressionPredictor)`` (this Protocol is ``runtime_checkable``) and fall
    back to :meth:`~ExpressionPredictor.score_sequence` otherwise. The two paths
    **must** agree: ``score_many`` returns one result per input, in input order,
    equal to calling ``score_sequence`` on each (only cheaper).
    """

    def score_many(self, dnas: list[str]) -> list[ExpressionResult]:
        """Score every DNA in ``dnas`` in one invocation, in input order."""
        ...
