"""GC-run constraint (CLAUDE.md §6, invariant #3).

A single ``LOCAL`` hard constraint whose ``ok_suffix`` veto and ``validate``
audit provably agree (invariant #3, "ok_suffix <=> validate"), with a
``context_len`` that provably suffices for the veto:

* :class:`GcRunConstraint` bans any run of *more than* ``max_run`` consecutive
  bases each drawn from the set ``{G, C}``. It is the direct analogue of
  :class:`~bt4.constraints.rules.HomopolymerConstraint`, but over the strong-base
  *set* ``{G, C}`` rather than a single repeated base -- so a mixed run such as
  ``GCGCGC`` counts as one length-6 run. Long GC stretches are GC-rich, form
  stable secondary structure, and are hard to synthesize. Its context is
  ``max_run`` trailing bases -- enough to see a trailing GC run the incoming
  codon would extend past the limit (``prefix`` is already feasible, so its
  trailing GC run is at most ``max_run`` bases long and fits in the window).

This constraint is purely hard (``penalty`` is ``0.0``); its soft cost is zero.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["GcRunConstraint"]


def _max_gc_run(seq: str) -> int:
    """Return the length of the longest run of consecutive ``{G, C}`` bases.

    A "GC run" is a maximal stretch of positions whose base is ``G`` or ``C``;
    the bases may be mixed (``GCGC`` is a run of four). ``seq`` is assumed to be
    upper-cased already, matching the trellis convention.

    Args:
        seq: The (upper-cased) DNA window to scan.
    """
    best = run = 0
    for ch in seq:
        if ch in ("G", "C"):
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


@dataclass(frozen=True, slots=True)
class GcRunConstraint:
    """Ban any run of more than ``max_run`` consecutive ``{G, C}`` bases.

    Consecutive positions whose base is ``G`` or ``C`` -- possibly mixed -- form
    a single GC run (``GCGCGCGC`` is a run of eight). Any such run longer than
    ``max_run`` is a hard violation.

    Attributes:
        max_run: Maximum allowed length of a consecutive ``{G, C}`` run
            (``>= 1``).
    """

    max_run: int
    name: str = field(default="gc_run", init=False)

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
        """Return ``False`` iff appending ``next_codon`` creates an over-long GC run.

        Only the last ``max_run`` bases of ``prefix`` are inspected. ``prefix``
        is assumed already feasible, so its trailing GC run is at most ``max_run``
        bases long and any run exceeding ``max_run`` in the window must involve
        ``next_codon``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        window = prefix[-self.max_run :] + next_codon
        return _max_gc_run(window) <= self.max_run

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per maximal ``{G, C}`` run longer than ``max_run``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for each maximal run of
            consecutive ``{G, C}`` bases whose length exceeds ``max_run``.
        """
        seq = dna.upper()
        n = len(seq)
        i = 0
        while i < n:
            if seq[i] in ("G", "C"):
                j = i + 1
                while j < n and seq[j] in ("G", "C"):
                    j += 1
                run_len = j - i
                if run_len > self.max_run:
                    run = seq[i:j]
                    yield Violation(
                        constraint="gc_run",
                        severity=Severity.HARD,
                        start=i,
                        end=j,
                        detail=f"GC run {run} (length {run_len}) exceeds max_run {self.max_run}",
                    )
                i = j
            else:
                i += 1
