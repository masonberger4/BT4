"""An honestly-labeled, uncalibrated folding fallback -- NOT thermodynamics.

:class:`BaselinePairingProxyModel` exists so that
:func:`bt4.biomodels.folding.default` can always return a working
:class:`~bt4.biomodels.folding.base.FoldingModel` even when the ViennaRNA
bindings are absent (CLAUDE.md section 4.3: "else a safe baseline; NEVER
crashes"). It is a crude self-complementarity proxy, **not** a free-energy
model.

Honesty (CLAUDE.md sections 6 and 10.6, "no placeholder model presented as a
feature"):

* :attr:`BaselinePairingProxyModel.calibrated` is ``False`` and its
  :attr:`~BaselinePairingProxyModel.name` is ``"baseline-pairing-proxy"`` -- the
  name screams baseline so no consumer can mistake it for real thermodynamics.
* The returned ``dg`` is in **arbitrary units, not kcal/mol**. It must never be
  presented as a calibrated deltaG. It exists only to give the refinement layer
  a smooth, dependency-free structural signal (more self-complementary 5' end
  => more negative proxy) and to keep :func:`default` non-crashing.

The proxy maximizes weighted Watson-Crick/wobble base pairs via a small
Nussinov dynamic program (GC pairs weighted above AT above GT wobble, with a
minimum hairpin loop). It captures self-complementarity honestly but is
thermodynamically meaningless: no stacking, no loop penalties, no temperature,
no calibrated parameters. It is ``O(n^3)`` in the window length, so it is only
ever run over a bounded 5' window for :meth:`score_sequence`; folding a whole
long sequence with it is discouraged.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt4.biomodels.folding.base import (
    DEFAULT_FIVE_PRIME_WINDOW,
    FoldingResult,
    five_prime_window,
)

__all__ = ["BaselinePairingProxyModel"]

_MIN_LOOP: int = 3
"""Minimum number of unpaired nucleotides enclosed by a base pair (hairpin loop)."""

_PAIR_WEIGHT: dict[tuple[str, str], float] = {
    ("G", "C"): 3.0,
    ("C", "G"): 3.0,
    ("A", "T"): 2.0,
    ("T", "A"): 2.0,
    ("G", "T"): 1.0,
    ("T", "G"): 1.0,
}
"""Crude, non-thermodynamic pair weights (GC > AT > GT wobble)."""


def _max_pair_weight(seq: str) -> float:
    """Return the maximum weighted base-pair sum over ``seq`` (Nussinov DP).

    Args:
        seq: An upper-cased ``{A,C,G,T}`` sequence.

    Returns:
        The best achievable sum of :data:`_PAIR_WEIGHT` over a set of
        non-crossing base pairs respecting the :data:`_MIN_LOOP` hairpin loop.
        ``0.0`` when the sequence is too short to form any pair.
    """
    n = len(seq)
    if n < _MIN_LOOP + 2:
        return 0.0
    dp: list[list[float]] = [[0.0] * n for _ in range(n)]
    for span in range(_MIN_LOOP + 1, n):
        for i in range(n - span):
            j = i + span
            best = dp[i + 1][j]  # leave i unpaired
            for k in range(i + _MIN_LOOP + 1, j + 1):
                weight = _PAIR_WEIGHT.get((seq[i], seq[k]))
                if weight is None:
                    continue
                left = dp[i + 1][k - 1] if k - 1 >= i + 1 else 0.0
                right = dp[k + 1][j] if k + 1 <= j else 0.0
                cand = weight + left + right
                if cand > best:
                    best = cand
            dp[i][j] = best
    return dp[0][n - 1]


@dataclass(frozen=True, slots=True)
class BaselinePairingProxyModel:
    """Uncalibrated self-complementarity proxy for 5' folding (NOT deltaG).

    A safe, dependency-free fallback :class:`~bt4.biomodels.folding.base.FoldingModel`.
    Its ``dg`` is ``-1 x`` the maximum weighted base-pair sum of the window, so a
    more self-complementary 5' end yields a more negative (more "stable") proxy.
    The value is in arbitrary units and is **never** a calibrated deltaG; see the
    module docstring.

    Attributes:
        five_prime_window: Number of 5' nucleotides scored by
            :meth:`score_sequence`. Defaults to
            :data:`~bt4.biomodels.folding.base.DEFAULT_FIVE_PRIME_WINDOW`.
    """

    five_prime_window: int = DEFAULT_FIVE_PRIME_WINDOW

    def __post_init__(self) -> None:
        """Validate the configured 5' window.

        Raises:
            ValueError: If ``five_prime_window`` is not a positive integer.
        """
        if self.five_prime_window <= 0:
            raise ValueError(
                f"five_prime_window must be positive, got {self.five_prime_window}"
            )

    @property
    def name(self) -> str:
        """Backend identifier -- deliberately screams that this is a baseline."""
        return "baseline-pairing-proxy"

    @property
    def calibrated(self) -> bool:
        """Always ``False``: this proxy is not thermodynamics (honesty flag)."""
        return False

    def fold(self, dna: str, window: int | None = None) -> FoldingResult:
        """Fold the 5' window of ``dna`` with the pairing proxy.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
            window: Number of 5' nucleotides to score, or ``None`` for the whole
                sequence (discouraged for long input; the proxy is ``O(n^3)``).

        Returns:
            A :class:`~bt4.biomodels.folding.base.FoldingResult` whose ``dg`` is
            an uncalibrated proxy in arbitrary units and whose ``structure`` is
            ``None``.
        """
        seq = five_prime_window(dna, window)
        weight = _max_pair_weight(seq)
        dg = -weight if weight else 0.0
        return FoldingResult(
            dg=dg,
            structure=None,
            window=window,
            model_name=self.name,
            calibrated=False,
        )

    def five_prime_dg(self, dna: str, window: int | None = None) -> float:
        """Return the uncalibrated proxy free energy of the 5' window.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
            window: Number of 5' nucleotides to score, or ``None`` for the whole
                sequence.

        Returns:
            The proxy ``dg`` in arbitrary units (NOT kcal/mol). More negative
            means a more self-complementary window.
        """
        return self.fold(dna, window).dg

    def score_sequence(self, dna: str) -> float:
        """Return the 5' window proxy deltaG directly (larger is better).

        See :meth:`bt4.biomodels.folding.base.FoldingModel.score_sequence` for
        the fixed orientation: a weakly self-complementary 5' end (proxy near
        zero) scores higher than a strongly self-complementary one (very
        negative proxy).
        """
        return self.five_prime_dg(dna, self.five_prime_window)
