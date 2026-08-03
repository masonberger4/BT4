"""The neutral placeholder expression model (honestly claims nothing).

:class:`NullExpressionModel` is what :func:`bt4.biomodels.expression.default`
returns until a *validated* expression head exists. It scores every sequence
``0.0`` and reports ``calibrated is False``.

Why a null model rather than a computed baseline. Folding and splice ship
uncalibrated *baselines* that still compute something (a Nussinov base-pairing
proxy; a consensus/PWM site score) because each has a defensible structural
anchor. Expression has **no** such anchor: the honest options are a validated
learned head (Phase 4) or nothing. A hand-weighted CAI+GC+ΔG composite dressed up
as "expression" would be the §10.5 magic-scalar / §10.6 placeholder-as-feature
trap. So the placeholder returns a neutral, information-free score: it exercises
the contract and lets the rerank plumbing be written and tested, without ever
implying an expression prediction or steering which sequence ships.

This module depends only on the standard library.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from bt4.biomodels.expression.base import ExpressionResult

__all__ = ["NullExpressionModel"]


@dataclass(frozen=True, slots=True)
class NullExpressionModel:
    """A neutral, uncalibrated placeholder: every sequence scores ``0.0``.

    Reports ``calibrated is False`` and a score of ``0.0`` for every input, so it
    is information-free by construction. It exists only to satisfy the
    :class:`~bt4.biomodels.expression.base.ExpressionPredictor` contract and give
    :func:`bt4.biomodels.expression.default` a non-crashing return until a
    validated head is built; it must never be read as an expression prediction.
    """

    name: str = field(default="null_expression", init=False)

    @property
    def calibrated(self) -> bool:
        """Return ``False`` -- this placeholder is never calibrated."""
        return False

    def score_sequence(self, dna: str) -> ExpressionResult:
        """Return a neutral ``0.0`` score for any ``dna``.

        Args:
            dna: A coding sequence (ignored -- the placeholder is information-free).

        Returns:
            An :class:`~bt4.biomodels.expression.base.ExpressionResult` with score
            ``0.0``, ``calibrated=False``, and placeholder units.
        """
        return ExpressionResult(
            score=0.0,
            model_name=self.name,
            calibrated=False,
            units="none (placeholder)",
        )
