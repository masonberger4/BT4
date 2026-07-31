"""Concrete local constraints (CLAUDE.md §6, invariant #3).

Two ``LOCAL`` hard constraints whose ``ok_suffix`` veto and ``validate`` audit
are proven to agree (invariant #3, "ok_suffix <=> validate"), with a
``context_len`` that provably suffices for the veto:

* :class:`HomopolymerConstraint` bans any run of one base longer than
  ``max_run``. Its context is ``max_run`` trailing bases - enough to see a run
  the incoming codon would extend past the limit.
* :class:`ForbiddenMotifConstraint` bans a set of literal DNA motifs (optionally
  together with their reverse complements). Its context is ``maxlen - 1`` bases,
  enough for any motif to straddle the boundary into the incoming codon.

Both are purely hard (``penalty`` is ``0.0``); their soft cost is zero.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4._accel import max_homopolymer_run, reverse_complement
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope
from bt4.domain.sequence import validate_dna

__all__ = ["ForbiddenMotifConstraint", "HomopolymerConstraint"]


@dataclass(frozen=True, slots=True)
class HomopolymerConstraint:
    """Ban any homopolymer run longer than ``max_run`` bases.

    Attributes:
        max_run: Maximum allowed run length of a single repeated base (``>= 1``).
    """

    max_run: int
    name: str = field(default="homopolymer", init=False)

    def __post_init__(self) -> None:
        """Validate the run bound.

        Raises:
            ValueError: If ``max_run`` is less than one.
        """
        if self.max_run < 1:
            raise ValueError(f"max_run must be >= 1, got {self.max_run!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``max_run`` - the trailing context ``ok_suffix`` needs."""
        return self.max_run

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` creates an over-long run.

        Only the last ``max_run`` bases of ``prefix`` are inspected. ``prefix``
        is assumed already feasible, so any run exceeding ``max_run`` in the
        window must involve ``next_codon``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        window = prefix[-self.max_run :] + next_codon
        return max_homopolymer_run(window) <= self.max_run

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per maximal run longer than ``max_run``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for each maximal run whose
            length exceeds ``max_run``.
        """
        seq = dna.upper()
        n = len(seq)
        i = 0
        while i < n:
            j = i + 1
            while j < n and seq[j] == seq[i]:
                j += 1
            run_len = j - i
            if run_len > self.max_run:
                yield Violation(
                    constraint="homopolymer",
                    severity=Severity.HARD,
                    start=i,
                    end=j,
                    detail=f"homopolymer run {seq[i]}x{run_len} exceeds max_run {self.max_run}",
                )
            i = j


@dataclass(frozen=True, slots=True)
class ForbiddenMotifConstraint:
    """Ban a set of literal DNA motifs (optionally their reverse complements).

    Attributes:
        motifs: The forbidden motifs (case-insensitive, ACGT only). An empty
            tuple makes the constraint an inert no-op.
        reverse_complement: When ``True``, each motif's reverse complement is
            also forbidden (strand-agnostic banning).
    """

    motifs: tuple[str, ...]
    reverse_complement: bool = False
    name: str = field(default="forbidden_motif", init=False)
    _motifs: frozenset[str] = field(init=False, repr=False, compare=False)
    _maxlen: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Normalize, validate, and (optionally) reverse-complement the motifs.

        Raises:
            ValueError: If any motif is empty or contains a non-ACGT character.
        """
        collected: set[str] = set()
        for motif in self.motifs:
            normalized = validate_dna(motif)
            collected.add(normalized)
            if self.reverse_complement:
                collected.add(reverse_complement(normalized))
        stored = frozenset(collected)
        maxlen = max((len(m) for m in stored), default=0)
        object.__setattr__(self, "_motifs", stored)
        object.__setattr__(self, "_maxlen", maxlen)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``max(0, maxlen - 1)`` - enough for a motif to cross the seam."""
        return max(0, self._maxlen - 1)

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff a motif occurrence overlaps ``next_codon``.

        The window is ``prefix[-(maxlen-1):] + next_codon``; an occurrence counts
        only when its exclusive end index lies beyond the start of ``next_codon``
        (so it actually involves the codon being appended). With no motifs the
        constraint is inert and this always returns ``True``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        if not self._motifs:
            return True
        tail = self._maxlen - 1
        window = (prefix[-tail:] if tail > 0 else "") + next_codon
        seam = len(window) - len(next_codon)
        for motif in self._motifs:
            idx = window.find(motif)
            while idx != -1:
                if idx + len(motif) > seam:
                    return False
                idx = window.find(motif, idx + 1)
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per motif occurrence in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for every occurrence of every
            forbidden motif, in a deterministic (motif-sorted) order.
        """
        seq = dna.upper()
        for motif in sorted(self._motifs):
            idx = seq.find(motif)
            while idx != -1:
                yield Violation(
                    constraint="forbidden_motif",
                    severity=Severity.HARD,
                    start=idx,
                    end=idx + len(motif),
                    detail=f"forbidden motif {motif} at {idx}",
                )
                idx = seq.find(motif, idx + 1)
