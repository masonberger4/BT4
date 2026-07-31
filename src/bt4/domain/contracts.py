"""The structural contracts that hold BT4 together (CLAUDE.md §4).

The ``ObjectiveTerm`` and ``Constraint`` protocols are pure vocabulary: they are
written entirely in terms of :class:`~bt4.domain.scope.Scope` and
:class:`~bt4.domain.result.Violation`, so they live in ``domain`` where *every*
layer can speak them. That placement is what lets the optimizer consume the
``Constraint`` contract and the concrete objective/constraint packages implement
it while each still imports **only** ``domain`` (the strict-layering rule of §3).

For ergonomics the protocols are re-exported from :mod:`bt4.objectives` and
:mod:`bt4.constraints`, so ``from bt4.objectives import ObjectiveTerm`` keeps
working -- but the single source of truth is here.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Protocol, runtime_checkable

from bt4.domain.result import Violation
from bt4.domain.scope import Scope

__all__ = ["Constraint", "ObjectiveTerm"]


@runtime_checkable
class ObjectiveTerm(Protocol):
    """One additive-in-its-own-deltas component of the objective vector.

    Implementations are typically frozen dataclasses. The load-bearing contract:

    * ``delta(prefix, codon, pos)`` is the term's contribution from placing
      ``codon`` at codon index ``pos`` after DNA ``prefix``. It must depend only
      on the last ``context_len`` characters of ``prefix`` (plus ``pos`` for
      ``POSITIONAL`` terms).
    * ``score(dna)`` recomputes the whole-sequence value directly from ``dna``
      and MUST equal the running sum of ``delta`` over the codons of ``dna``
      (invariant #4, "delta == score").

    Convention: every term is oriented so that **larger is better**.
    """

    @property
    def name(self) -> str:
        """Stable identifier used as the key in an ``ObjectiveVector``.

        Declared as a read-only property so concrete terms may be frozen
        dataclasses (a read-only ``name`` field satisfies this contract).
        """
        ...

    def scope(self) -> Scope:
        """The locality class of this term."""
        ...

    def context_len(self) -> int:
        """Number of trailing DNA characters ``delta`` may inspect in ``prefix``."""
        ...

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Incremental contribution of ``codon`` at codon index ``pos``."""
        ...

    def score(self, dna: str) -> float:
        """Whole-sequence value of ``dna``; equals the sum of this term's deltas."""
        ...


@runtime_checkable
class Constraint(Protocol):
    """A hard/soft feasibility rule over a coding sequence.

    Implementations are typically frozen dataclasses. Contract:

    * ``ok_suffix(prefix, next_codon)`` is a hard veto: it returns ``False`` iff
      appending ``next_codon`` to ``prefix`` introduces a new hard violation. It
      must depend only on the last ``context_len`` characters of ``prefix``.
    * ``penalty(prefix, next_codon)`` is a soft, non-negative cost for the same
      extension (``0.0`` for a purely hard constraint).
    * ``validate(dna)`` is the whole-sequence audit and the source of truth for
      reported violations. Invariant #3 binds ``ok_suffix`` and ``validate``
      together: a sequence built respecting ``ok_suffix`` has zero hard
      violations, and ``context_len`` must actually suffice for ``ok_suffix``.
    """

    @property
    def name(self) -> str:
        """Stable identifier used in ``Violation.constraint``.

        Declared as a read-only property so concrete constraints may be frozen
        dataclasses (a read-only ``name`` field satisfies this contract).
        """
        ...

    def scope(self) -> Scope:
        """The locality class of this constraint."""
        ...

    def context_len(self) -> int:
        """Number of trailing DNA characters ``ok_suffix`` may inspect in ``prefix``."""
        ...

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` adds a hard violation."""
        ...

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Non-negative soft cost of the extension (``0.0`` when purely hard)."""
        ...

    def validate(self, dna: str) -> Iterable[Violation]:
        """Audit a whole sequence, yielding every violation found."""
        ...
