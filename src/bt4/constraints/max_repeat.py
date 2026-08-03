"""Whole-sequence maximum-repeat constraint (CLAUDE.md invariant #3, section 10.1).

:class:`MaxRepeatConstraint` bans *any* repeated substring longer than
``max_length`` anywhere in the sequence -- direct repeats, inverted repeats, and
palindromes -- and is reverse-complement aware. This mirrors BT3's
``almost-there`` repeat ban, but computed honestly and efficiently rather than by
its O(L^3) all-pairs scan.

Why this is a **GLOBAL** constraint. The two copies of a repeat can sit
arbitrarily far apart, so no fixed trailing window sees both. Unlike the LOCAL
:class:`~bt4.constraints.repeats.TandemRepeatConstraint` /
:class:`~bt4.constraints.repeats.InvertedRepeatConstraint` (which bound the two
arms to a fixed span and so declare a finite ``context_len``), this constraint's
``ok_suffix`` must read the *whole* prefix. It is therefore **unsafe to merge
into a bounded-context exact-DP trellis** (doing so is exactly the BT3 mistake of
section 10.1: a global context cap silently over-merges prefixes that differ far
back). The planner keeps it out of the merged DP and enforces it in the
refinement / whole-sequence-validation layer; ``scope()`` returns
:attr:`~bt4.domain.scope.Scope.GLOBAL` to say so, and ``context_len()`` returns
an unbounded sentinel.

Why checking a single k-mer length is sufficient and correct. Let
``k = max_length + 1``. A sequence contains a repeat longer than ``max_length``
**iff** it contains one of length exactly ``k``:

* *Direct.* If a substring of length ``L > max_length`` occurs at two distinct
  starts, its first ``k`` bases form a ``k``-mer that also occurs at both starts.
  Conversely a ``k``-mer occurring twice is itself a direct repeat of length
  ``k > max_length``.
* *Inverted / palindrome.* If a length-``L`` substring ``X`` has
  ``reverse_complement(X)`` occurring elsewhere, then the ``k``-mer
  ``Y = X[L - k:]`` has ``reverse_complement(Y)`` equal to the first ``k`` bases
  of ``reverse_complement(X)``, which also occurs -- a length-``k`` inverted
  repeat. A palindrome is the special case ``reverse_complement(Y) == Y``.

So it suffices to test, for ``k = max_length + 1``, whether any ``k``-mer occurs
at two positions (direct) or has its reverse complement occurring anywhere
(inverted, with the self-match being a palindrome). Any longer repeat necessarily
contains such a ``k``-mer, so nothing is missed and the scan is O(L * k).

The constraint is purely hard (``penalty`` is ``0.0``). ``ok_suffix`` and
``validate`` provably agree (invariant #3): a ``k``-mer is *offending* iff it is
duplicated or its reverse complement occurs anywhere, and ``ok_suffix`` -- reading
the whole prefix -- vetoes any codon that would complete a new offending ``k``-mer
(one ending inside the codon). Since both the duplicate and reverse-complement
relations are symmetric, checking only the newly introduced ``k``-mers is enough
to keep the full sequence free of offending ``k``-mers, so a build respecting
``ok_suffix`` yields zero hard violations under ``validate``.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4._accel import longest_repeat, reverse_complement
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["MaxRepeatConstraint"]


@dataclass(frozen=True, slots=True)
class MaxRepeatConstraint:
    """Ban any repeated substring longer than ``max_length`` (RC-aware, GLOBAL).

    A violation is any ``k``-mer (``k = max_length + 1``) that occurs at two
    positions (a *direct* repeat) or whose reverse complement occurs anywhere in
    the sequence (an *inverted* repeat; a *palindrome* when the ``k``-mer equals
    its own reverse complement). This bans every longer repeat transitively (see
    the module docstring for the sufficiency proof).

    Attributes:
        max_length: The longest repeated substring that is still allowed
            (``>= 1``); anything strictly longer is a violation.
    """

    max_length: int
    name: str = field(default="max_repeat", init=False)

    def __post_init__(self) -> None:
        """Validate the length bound.

        Raises:
            ValueError: If ``max_length`` is less than one.
        """
        if self.max_length < 1:
            raise ValueError(f"max_length must be >= 1, got {self.max_length!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.GLOBAL`.

        The two copies of a repeat may be arbitrarily far apart, so this
        constraint is non-local and must never be merged into a bounded-context
        exact-DP trellis (CLAUDE.md section 10.1); the planner enforces it in the
        whole-sequence validation / refinement layer instead.
        """
        return Scope.GLOBAL

    def context_len(self) -> int:
        """Return ``sys.maxsize`` -- an unbounded sentinel.

        ``ok_suffix`` reads the *entire* prefix (a repeat's earlier copy can lie
        anywhere), so no finite trailing window suffices. Returning ``sys.maxsize``
        states this honestly, and any naive ``prefix[-context_len():]`` slice
        degrades to the whole prefix -- the correct semantics. This is why the
        constraint is ``GLOBAL`` and unsafe to merge into a bounded trellis.
        """
        return sys.maxsize

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` introduces a new repeat.

        Reads the **whole** prefix (necessary and correct for a global
        constraint). Only ``k``-mers (``k = max_length + 1``) whose last base lies
        inside ``next_codon`` are new; for each such ``k``-mer this vetoes the
        codon if the ``k``-mer occurs at another position in ``prefix + next_codon``
        (a direct repeat) or its reverse complement occurs anywhere in
        ``prefix + next_codon`` (an inverted repeat, or a palindrome when the
        ``k``-mer is its own reverse complement).

        Because both relations are symmetric, checking only the newly introduced
        ``k``-mers also catches an *old* ``k``-mer that this codon turns into a
        repeat (its new duplicate or reverse-complement partner is itself a new
        ``k``-mer). ``prefix`` is assumed already feasible.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        k = self.max_length + 1
        combined = prefix + next_codon
        n = len(combined)
        if n < k:
            return True
        first = max(0, len(prefix) - k + 1)
        for s in range(first, n - k + 1):
            kmer = combined[s : s + k]
            # Direct: the same k-mer at any other start position.
            if combined.find(kmer) != s or combined.find(kmer, s + 1) != -1:
                return False
            # Inverted / palindrome: the k-mer's reverse complement occurs anywhere.
            if reverse_complement(kmer) in combined:
                return False
        return True

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per offending ``k``-mer occurrence.

        A ``k``-mer occurrence (``k = max_length + 1``) is offending when it is
        duplicated (*direct*), equals its own reverse complement (*palindrome*),
        or has its reverse complement present elsewhere (*inverted*). One
        violation is emitted per offending occurrence, in ascending start order;
        a position that is both direct and (palindrome/inverted) is reported once,
        preferring ``direct`` then ``palindrome`` then ``inverted``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for every offending ``k``-mer
            occurrence, with a ``detail`` naming the repeat type and positions.
        """
        seq = dna.upper()
        k = self.max_length + 1
        n = len(seq)
        if n < k:
            return
        # Fast path (§7, invariant #3): ``longest_repeat(seq) > max_length`` iff a
        # hard violation exists, so when the longest RC-aware repeat is within the
        # limit there is nothing to report -- skip the O(n*k) k-mer position scan
        # entirely. This is the hot path re-run per SA refinement move; the
        # (optionally Rust-accelerated) ``longest_repeat`` shares this constraint's
        # exact direct/inverted/palindrome notion (cross-checked in the tests).
        if longest_repeat(seq) <= self.max_length:
            return
        positions: dict[str, list[int]] = {}
        for i in range(n - k + 1):
            positions.setdefault(seq[i : i + k], []).append(i)
        for i in range(n - k + 1):
            kmer = seq[i : i + k]
            occurrences = positions[kmer]
            rc = reverse_complement(kmer)
            if len(occurrences) > 1:
                others = [p for p in occurrences if p != i]
                detail = f"direct repeat {kmer} at {i}; also at {others}"
            elif rc == kmer:
                detail = f"palindrome {kmer} at {i}"
            elif rc in positions:
                detail = f"inverted repeat {kmer} at {i}; revcomp at {positions[rc]}"
            else:
                continue
            yield Violation(
                constraint="max_repeat",
                severity=Severity.HARD,
                start=i,
                end=i + k,
                detail=detail,
            )
