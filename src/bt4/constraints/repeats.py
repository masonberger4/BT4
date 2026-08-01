"""Tandem- and inverted-repeat constraints (CLAUDE.md §6, invariant #3).

Two ``LOCAL`` hard constraints that ban repetitive DNA structures whose
``ok_suffix`` veto and ``validate`` audit provably agree (invariant #3,
"ok_suffix <=> validate"), each with a ``context_len`` that provably suffices
for the veto:

* :class:`TandemRepeatConstraint` bans a *tandem repeat*: a ``unit_len``-mer
  repeated ``copies`` times back-to-back -- i.e. a substring of length
  ``span = unit_len * copies`` whose period is ``unit_len`` (``seq[k] ==
  seq[k + unit_len]`` across the whole span). Such runs slip during replication
  and are hard to synthesize. Its context is ``span - 1`` trailing bases: one
  base short of a full repeat, so the incoming codon that completes the final
  copy is always visible and vetoed. A run of *more* than ``copies`` copies is
  banned transitively, since it contains a length-``span`` occurrence.
* :class:`InvertedRepeatConstraint` bans a *hairpin* (stem-loop): a
  ``stem``-length arm ``X`` followed, after a loop of ``0..loop_max`` bases, by
  its own reverse complement -- the pattern ``X . gap . revcomp(X)``. Such stems
  fold into secondary structure that occludes ribosome loading. Its context is
  ``2 * stem + loop_max - 1`` trailing bases: one base short of the longest
  hairpin, so the codon completing the 3' arm is always visible and vetoed.

Both are purely hard (``penalty`` is ``0.0``). Each ``validate`` scans the whole
sequence and reports *every* occurrence, while ``ok_suffix`` -- mirroring
:class:`~bt4.constraints.rules.ForbiddenMotifConstraint` -- only vetoes a repeat
whose 3' end lies past the seam into the incoming codon. A repeat lying wholly
inside the already-feasible prefix is not the new codon's fault (and, by
induction over a feasible build, no such repeat exists), so the two agree.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4._accel import reverse_complement
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["InvertedRepeatConstraint", "TandemRepeatConstraint"]


@dataclass(frozen=True, slots=True)
class TandemRepeatConstraint:
    """Ban a tandem repeat: a ``unit_len``-mer repeated ``copies`` times.

    A violation is any substring of length ``span = unit_len * copies`` whose
    period is ``unit_len`` -- that is, the first ``unit_len`` bases repeated
    ``copies`` times back-to-back. Because a longer run contains a length-``span``
    window, banning exactly ``copies`` copies also bans any greater number.

    Attributes:
        unit_len: Length of the repeated unit (``>= 1``).
        copies: Number of consecutive copies that constitutes a repeat
            (``>= 2``; a single copy is not a repeat).
    """

    unit_len: int
    copies: int = 3
    name: str = field(default="tandem_repeat", init=False)

    def __post_init__(self) -> None:
        """Validate the unit length and copy count.

        Raises:
            ValueError: If ``unit_len`` is less than one or ``copies`` is less
                than two.
        """
        if self.unit_len < 1:
            raise ValueError(f"unit_len must be >= 1, got {self.unit_len!r}")
        if self.copies < 2:
            raise ValueError(f"copies must be >= 2, got {self.copies!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``unit_len * copies - 1`` -- the trailing context needed.

        One base short of a full repeat: enough trailing context that the
        incoming codon completing the final copy is seen across the seam.
        """
        return self.unit_len * self.copies - 1

    def _is_tandem(self, s: str, i: int) -> bool:
        """Return ``True`` iff ``s[i : i + span]`` has period ``unit_len``."""
        u = self.unit_len
        span = u * self.copies
        return all(s[i + k] == s[i + k + u] for k in range(span - u))

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` completes a tandem repeat.

        The window is ``prefix[-(span - 1):] + next_codon``; an occurrence counts
        only when its exclusive end index lies beyond the seam (so it actually
        involves the incoming codon). ``prefix`` is assumed already feasible.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        span = self.unit_len * self.copies
        window = prefix[-(span - 1) :] + next_codon
        seam = len(window) - len(next_codon)
        for i in range(len(window) - span + 1):
            if i + span > seam and self._is_tandem(window, i):
                return False
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per length-``span`` tandem repeat in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for every substring of length
            ``unit_len * copies`` whose period is ``unit_len``, in left-to-right
            order.
        """
        seq = dna.upper()
        span = self.unit_len * self.copies
        for i in range(len(seq) - span + 1):
            if self._is_tandem(seq, i):
                yield Violation(
                    constraint="tandem_repeat",
                    severity=Severity.HARD,
                    start=i,
                    end=i + span,
                    detail=f"tandem repeat {seq[i : i + span]} at {i}",
                )


@dataclass(frozen=True, slots=True)
class InvertedRepeatConstraint:
    """Ban a hairpin: an arm ``X`` followed by ``revcomp(X)`` within a loop.

    A violation is any substring of the form ``X . gap . revcomp(X)`` where
    ``len(X) == stem`` and the loop ``gap`` has length ``0 <= len(gap) <=
    loop_max`` -- i.e. the two ``stem``-length arms are reverse complements of
    each other and so can base-pair into a stem-loop.

    Attributes:
        stem: Length of each arm of the stem (``>= 1``).
        loop_max: Maximum loop length between the arms (``>= 0``; ``0`` is a
            perfect adjacent inverted repeat / palindrome).
    """

    stem: int
    loop_max: int = 0
    name: str = field(default="inverted_repeat", init=False)

    def __post_init__(self) -> None:
        """Validate the stem length and loop bound.

        Raises:
            ValueError: If ``stem`` is less than one or ``loop_max`` is negative.
        """
        if self.stem < 1:
            raise ValueError(f"stem must be >= 1, got {self.stem!r}")
        if self.loop_max < 0:
            raise ValueError(f"loop_max must be >= 0, got {self.loop_max!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``2 * stem + loop_max - 1`` -- the trailing context needed.

        One base short of the longest hairpin (``2 * stem + loop_max``): enough
        trailing context that the incoming codon completing the 3' arm is seen.
        """
        return 2 * self.stem + self.loop_max - 1

    def _is_hairpin(self, s: str, i: int, gap: int) -> bool:
        """Return ``True`` iff ``s[i:]`` opens a hairpin with loop length ``gap``.

        Assumes ``i + 2 * stem + gap <= len(s)`` so both arms are full length.
        """
        stem = self.stem
        left = s[i : i + stem]
        right = s[i + stem + gap : i + 2 * stem + gap]
        return reverse_complement(left) == right

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` completes a hairpin.

        The window is ``prefix[-context_len:] + next_codon``; a hairpin counts
        only when its exclusive end index lies beyond the seam (so it actually
        involves the incoming codon). ``prefix`` is assumed already feasible.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        window = prefix[-self.context_len() :] + next_codon
        seam = len(window) - len(next_codon)
        n = len(window)
        for gap in range(self.loop_max + 1):
            span = 2 * self.stem + gap
            for i in range(n - span + 1):
                if i + span > seam and self._is_hairpin(window, i, gap):
                    return False
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per hairpin in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for every ``X . gap .
            revcomp(X)`` occurrence (``len(X) == stem``, ``0 <= gap <=
            loop_max``), ordered by start then loop length.
        """
        seq = dna.upper()
        n = len(seq)
        for i in range(n):
            for gap in range(self.loop_max + 1):
                span = 2 * self.stem + gap
                if i + span <= n and self._is_hairpin(seq, i, gap):
                    yield Violation(
                        constraint="inverted_repeat",
                        severity=Severity.HARD,
                        start=i,
                        end=i + span,
                        detail=f"inverted repeat {seq[i : i + span]} at {i}",
                    )
