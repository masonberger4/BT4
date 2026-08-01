"""The ``FoldingModel`` contract for 5' mRNA secondary-structure free energy.

BT4 treats 5' folding free energy (deltaG) as a genuinely non-local objective
that lives in the refinement layer (CLAUDE.md sections 6 and 4.3). A folding
backend reduces a coding sequence to a scalar deltaG for its 5' window: the
more negative the deltaG, the more stable the base-paired structure. Strong
structure near the start codon occludes ribosome loading, so BT4's design goal
is to *avoid* RBS/Kozak-occluding hairpins (CLAUDE.md section 6).

Two honesty rules shape this module:

* **``calibrated`` is a first-class flag.** A backend may only claim
  ``calibrated is True`` when it computes real thermodynamics from a validated,
  hash-pinned parameter set (the ViennaRNA backend). Any proxy or fallback must
  report ``calibrated is False`` and never be presented as a real deltaG
  (CLAUDE.md sections 4.3 and 10.6 -- "no placeholder model presented as a
  feature").
* **Orientation is fixed and documented.** Every :meth:`FoldingModel.score_sequence`
  follows BT4's convention that **larger is better**: it returns the 5' window
  deltaG *directly*, so a weakly structured 5' end (deltaG near zero) scores
  higher than a strongly structured one (deltaG very negative). Maximizing the
  score therefore opens up the 5' end, matching the biology above.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from bt4.domain.sequence import validate_dna

__all__ = [
    "DEFAULT_FIVE_PRIME_WINDOW",
    "FoldingModel",
    "FoldingResult",
    "five_prime_window",
]

DEFAULT_FIVE_PRIME_WINDOW: int = 48
"""Default 5' window (nucleotides) scored by :meth:`FoldingModel.score_sequence`.

48 nt is 16 codons -- a standard span for the translation-initiation region
where 5' secondary structure most affects ribosome loading.
"""


def five_prime_window(dna: str, window: int | None) -> str:
    """Return the validated, upper-cased 5' window of ``dna``.

    Args:
        dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
        window: Number of 5' nucleotides to keep, or ``None`` for the whole
            sequence. A window longer than ``dna`` yields the whole sequence.

    Returns:
        The upper-cased 5' slice (``dna[:window]``), or the whole upper-cased
        sequence when ``window is None``.

    Raises:
        ValueError: If ``dna`` is empty/non-ACGT, or ``window`` is not positive.
    """
    seq = validate_dna(dna)
    if window is None:
        return seq
    if window <= 0:
        raise ValueError(f"window must be a positive number of nucleotides, got {window}")
    return seq[:window]


@dataclass(frozen=True, slots=True)
class FoldingResult:
    """The outcome of folding one 5' window (immutable).

    Attributes:
        dg: The window's free energy. Real minimum-free-energy deltaG in
            kcal/mol when ``calibrated`` is ``True``; an uncalibrated proxy in
            arbitrary units otherwise. More negative means more stable
            structure.
        structure: A dot-bracket secondary structure if the backend produced
            one, else ``None`` (proxies do not).
        window: The *requested* 5' window in nucleotides (``None`` means the
            whole sequence). A request longer than ``dna`` still folds only the
            whole sequence, so this records the request, not the folded length.
        model_name: The :attr:`FoldingModel.name` of the producing backend.
        calibrated: Mirror of :attr:`FoldingModel.calibrated` for the producing
            backend -- ``False`` marks ``dg`` as an uncalibrated proxy.
    """

    dg: float
    structure: str | None
    window: int | None
    model_name: str
    calibrated: bool


@runtime_checkable
class FoldingModel(Protocol):
    """A backend that scores 5' mRNA secondary-structure free energy.

    Implementations are swappable behind this contract (CLAUDE.md section 4.3):
    consumers depend only on the protocol and never on a concrete backend. Every
    implementation must satisfy the orientation and honesty rules documented at
    the module level.
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
        """Honesty flag: ``True`` only for validated, hash-pinned thermodynamics.

        A backend returns ``True`` here **only** when its deltaG is a real,
        hash-pinned thermodynamic computation. Proxies and safe fallbacks return
        ``False`` so their numbers are never mistaken for calibrated deltaG.
        """
        ...

    def five_prime_dg(self, dna: str, window: int | None = None) -> float:
        """Return the minimum-free-energy deltaG of the 5' window of ``dna``.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).
            window: Number of 5' nucleotides to fold, or ``None`` for the whole
                sequence.

        Returns:
            The deltaG in kcal/mol (calibrated backends) or an uncalibrated
            proxy in arbitrary units (fallback backends). More negative means a
            more stable base-paired structure.
        """
        ...

    def score_sequence(self, dna: str) -> float:
        """Return the objective score of ``dna`` (larger is better).

        Orientation (fixed for every backend): the 5' window deltaG is returned
        *directly*, so a weakly structured 5' end (deltaG near zero) scores
        higher than a strongly structured one (deltaG very negative). Maximizing
        the score opens up the 5' end near the start codon, which aids ribosome
        loading (CLAUDE.md section 6, "avoid RBS/Kozak-occluding hairpins").
        """
        ...
