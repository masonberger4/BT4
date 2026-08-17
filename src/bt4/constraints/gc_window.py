"""Windowed GC-content constraint (CLAUDE.md §6, invariant #3).

A whole-sequence GC *count* budget (``gc_min``/``gc_max``) pins the total but says
nothing about how GC is **distributed**: a sequence can hit 50% overall while
carrying a 74%-GC stretch that fails synthesis QC. The commercial synthesis
vendors specify GC the other way round -- as a bound on **every sliding window**
(e.g. 25-65% GC in any 50 bp) -- because a local GC extreme is what actually
breaks oligo synthesis, PCR assembly and sequencing.

:class:`GcWindowConstraint` is that rule, and it is genuinely ``LOCAL``: a bound
on a fixed-length window has *bounded context*, so it belongs in the exact codon
trellis alongside :class:`~bt4.constraints.gc_run.GcRunConstraint` -- not in the
refinement layer. A run under it therefore keeps a ``PROVEN_OPTIMAL`` certificate
rather than degrading to a heuristic (contrast the genuinely non-local
:class:`~bt4.constraints.max_repeat.MaxRepeatConstraint`).

Why ``context_len == window - 1`` suffices (invariant #3). Only *full-length*
windows are constrained. A window whose rightmost base lands inside the incoming
codon is the only kind this extension can newly complete; the furthest such window
back is the one whose last base is the codon's **first** base, and it reaches
``window - 1`` bases into the prefix. Every window that ends inside the prefix was
already checked when it completed (``prefix`` is assumed feasible), so
``prefix[-(window-1):] + next_codon`` shows ``ok_suffix`` every window it must
judge. A sequence built respecting ``ok_suffix`` therefore has zero hard violations
under ``validate``, which checks exactly the same set of full-length windows.

Bounds are given as **fractions** (matching how vendors state them) and converted
once to integer count thresholds, so the comparison is exact integer arithmetic and
cannot drift with floating-point rounding.

This constraint is purely hard (``penalty`` is ``0.0``).
"""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4._accel import gc_count
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["GcWindowConstraint"]

# Tolerance applied when converting a fractional bound to an integer count
# threshold, so a bound that is exactly representable (0.65 * 20 == 13.0) is not
# pushed off by one by binary floating-point representation error.
_EPS = 1e-9


@dataclass(frozen=True, slots=True)
class GcWindowConstraint:
    """Bound the GC fraction of every ``window``-nt sliding window.

    Any full-length window whose GC fraction falls below ``gc_min`` or above
    ``gc_max`` is a hard violation. Windows are evaluated at every start offset
    (step 1); a sequence shorter than ``window`` contains no full window and is
    therefore unconstrained.

    Attributes:
        window: Window length in nucleotides (``>= 1``).
        gc_min: Minimum allowed GC fraction in ``[0, 1]``.
        gc_max: Maximum allowed GC fraction in ``[0, 1]``, ``>= gc_min``.
    """

    window: int
    gc_min: float = 0.0
    gc_max: float = 1.0
    name: str = field(default="gc_window", init=False)

    def __post_init__(self) -> None:
        """Validate the window length and the fractional bounds.

        Raises:
            ValueError: If ``window`` is less than one, a bound lies outside
                ``[0, 1]``, or ``gc_min`` exceeds ``gc_max``.
        """
        if self.window < 1:
            raise ValueError(f"window must be >= 1, got {self.window!r}")
        for label, value in (("gc_min", self.gc_min), ("gc_max", self.gc_max)):
            if not (0.0 <= value <= 1.0):
                raise ValueError(f"{label} must be a fraction in [0, 1], got {value!r}")
        if self.gc_min > self.gc_max:
            raise ValueError(
                f"gc_min ({self.gc_min}) must not exceed gc_max ({self.gc_max})"
            )

    @property
    def min_count(self) -> int:
        """Fewest GC bases a full window may contain (derived from ``gc_min``)."""
        return math.ceil(self.gc_min * self.window - _EPS)

    @property
    def max_count(self) -> int:
        """Most GC bases a full window may contain (derived from ``gc_max``)."""
        return math.floor(self.gc_max * self.window + _EPS)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``window - 1`` -- the trailing context ``ok_suffix`` needs.

        See the module docstring for why this provably suffices: the earliest
        window this extension can complete starts ``window - 1`` bases before the
        incoming codon's first base.
        """
        return self.window - 1

    def _bad_count(self, count: int) -> bool:
        """Return ``True`` iff a full window with ``count`` GC bases is out of bounds."""
        return count < self.min_count or count > self.max_count

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` completes an out-of-bounds window.

        Only windows whose rightmost base lies inside ``next_codon`` are judged --
        a window ending inside the already-feasible ``prefix`` was checked when it
        completed and is not this codon's fault.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        tail = prefix[-(self.window - 1) :] if self.window > 1 else ""
        combined = tail + next_codon
        n = len(combined)
        # Windows ending inside next_codon end at combined indices
        # len(tail) .. n-1; each starts `window` bases earlier. A window that would
        # start before the sequence begins is not full-length and is not judged.
        for end in range(len(tail), n):
            start = end - self.window + 1
            if start < 0:
                continue
            if self._bad_count(gc_count(combined[start : end + 1])):
                return False
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per out-of-bounds window.

        The unit is **one window**, not one merged region, so the hard-violation
        count is literally "how many windows are out of bounds". That granularity
        is load-bearing, not cosmetic: when this rule is refinement-enforced (a
        window too wide for the trellis, see
        :data:`~bt4.pipeline.optimize._GC_WINDOW_TRELLIS_MAX_NT`), the annealer
        drives the *count* down, and a merged count would be a plateau -- easing a
        region from 74% to 66% GC would not change it, so the search would see no
        gradient and stall above the bound even though a compliant sequence exists.
        Overlapping windows do mean one bad stretch yields several violations; they
        share a constraint name and adjacent spans, so a viewer renders them as one
        band.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` per offending window, in
            ascending window-start order.
        """
        seq = dna.upper()
        n = len(seq)
        if n < self.window:
            return
        for start in range(n - self.window + 1):
            count = gc_count(seq[start : start + self.window])
            if count < self.min_count:
                yield self._violation(start, count, "low")
            elif count > self.max_count:
                yield self._violation(start, count, "high")

    def _violation(self, start: int, count: int, direction: str) -> Violation:
        """Build the violation for the offending window starting at ``start``."""
        bound = self.gc_min if direction == "low" else self.gc_max
        relation = "below min" if direction == "low" else "above max"
        return Violation(
            constraint="gc_window",
            severity=Severity.HARD,
            start=start,
            end=start + self.window,
            detail=(
                f"{self.window}-nt window at {start} has GC {count / self.window:.1%}, "
                f"{relation} {bound:.1%}"
            ),
        )
