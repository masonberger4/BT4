"""ViennaRNA-backed folding model -- real thermodynamics, lazily imported.

:class:`ViennaFoldingModel` computes minimum-free-energy deltaG in kcal/mol via
the ViennaRNA Python bindings (the ``RNA`` module). It is the only
:class:`~bt4.biomodels.folding.base.FoldingModel` allowed to report
``calibrated is True``, because its numbers are real thermodynamics from
ViennaRNA's validated, versioned energy parameters (CLAUDE.md sections 6 and
4.3).

The ViennaRNA binding is a heavy, optional dependency, so this module keeps
``import bt4`` lightweight (CLAUDE.md section 3): the ``RNA`` module is imported
**only inside methods**, never at module load. Importing this file therefore
never requires ViennaRNA. :meth:`ViennaFoldingModel.available` reports whether
the binding can be imported without raising, and :func:`bt4.biomodels.folding.default`
uses it to fall back to the baseline when ViennaRNA is absent.

This module depends only on :mod:`bt4.domain`, the standard library, and -- lazily,
inside methods -- the optional ViennaRNA binding.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from bt4.biomodels.folding.base import (
    DEFAULT_FIVE_PRIME_WINDOW,
    FoldingResult,
    five_prime_window,
)

__all__ = ["ViennaFoldingModel"]


def _import_rna() -> Any:
    """Import and return the ViennaRNA ``RNA`` module (lazy, guarded).

    Returns:
        The imported ViennaRNA bindings module.

    Raises:
        ModuleNotFoundError: If neither ``RNA`` nor ``ViennaRNA`` can be
            imported. Install the ``bt4[fold]`` extra (ViennaRNA) or use the
            baseline model.
    """
    try:
        import RNA  # type: ignore[import-not-found]

        return RNA
    except ImportError:
        pass
    try:
        import ViennaRNA  # type: ignore[import-not-found]

        return ViennaRNA
    except ImportError as exc:
        raise ModuleNotFoundError(
            "ViennaRNA Python bindings ('RNA') are not installed; "
            "install the 'bt4[fold]' extra or use the baseline folding model"
        ) from exc


@dataclass(frozen=True, slots=True)
class ViennaFoldingModel:
    """Calibrated 5' folding deltaG via ViennaRNA (real thermodynamics).

    Computes minimum-free-energy deltaG in kcal/mol for the 5' window using
    ViennaRNA's ``fold_compound`` MFE routine at the configured temperature.
    ``T`` is folded as ``U`` so DNA coding sequences map onto RNA parameters.

    Attributes:
        five_prime_window: Number of 5' nucleotides scored by
            :meth:`score_sequence`. Defaults to
            :data:`~bt4.biomodels.folding.base.DEFAULT_FIVE_PRIME_WINDOW`.
        temperature_c: Folding temperature in degrees Celsius (default 37.0).
    """

    five_prime_window: int = DEFAULT_FIVE_PRIME_WINDOW
    temperature_c: float = 37.0

    def __post_init__(self) -> None:
        """Validate the configured 5' window.

        Raises:
            ValueError: If ``five_prime_window`` is not a positive integer.
        """
        if self.five_prime_window <= 0:
            raise ValueError(
                f"five_prime_window must be positive, got {self.five_prime_window}"
            )

    @staticmethod
    def available() -> bool:
        """Return whether the ViennaRNA binding can be imported (never raises)."""
        try:
            _import_rna()
        except ImportError:
            return False
        return True

    @property
    def name(self) -> str:
        """Backend identifier."""
        return "viennarna-mfe"

    @property
    def calibrated(self) -> bool:
        """Always ``True``: ViennaRNA is real, versioned thermodynamics."""
        return True

    def fold(self, dna: str, window: int | None = None) -> FoldingResult:
        """Fold the 5' window of ``dna`` with ViennaRNA (MFE).

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
            window: Number of 5' nucleotides to fold, or ``None`` for the whole
                sequence.

        Returns:
            A :class:`~bt4.biomodels.folding.base.FoldingResult` with the MFE
            deltaG in kcal/mol and the dot-bracket ``structure``.

        Raises:
            ModuleNotFoundError: If the ViennaRNA binding is not importable.
        """
        seq = five_prime_window(dna, window)
        rna = _import_rna()
        rna_seq = seq.replace("T", "U")
        model_details = rna.md()
        model_details.temperature = self.temperature_c
        fold_compound = rna.fold_compound(rna_seq, model_details)
        structure, mfe = fold_compound.mfe()
        return FoldingResult(
            dg=float(mfe),
            structure=str(structure),
            window=window,
            model_name=self.name,
            calibrated=True,
        )

    def five_prime_dg(self, dna: str, window: int | None = None) -> float:
        """Return the ViennaRNA MFE deltaG (kcal/mol) of the 5' window.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
            window: Number of 5' nucleotides to fold, or ``None`` for the whole
                sequence.

        Returns:
            The minimum-free-energy deltaG in kcal/mol. More negative means a
            more stable base-paired structure.

        Raises:
            ModuleNotFoundError: If the ViennaRNA binding is not importable.
        """
        return self.fold(dna, window).dg

    def score_sequence(self, dna: str) -> float:
        """Return the 5' window MFE deltaG directly (larger is better).

        See :meth:`bt4.biomodels.folding.base.FoldingModel.score_sequence` for
        the fixed orientation: a weakly structured 5' end (deltaG near zero)
        scores higher than a strongly structured one (deltaG very negative).

        Raises:
            ModuleNotFoundError: If the ViennaRNA binding is not importable.
        """
        return self.five_prime_dg(dna, self.five_prime_window)
