"""Functional polyadenylation-signal constraint (CLAUDE.md §6).

A premature poly(A) site inside a coding sequence truncates the transcript, so
BT4 already ships a ``poly_a_signal`` forbidden-motif preset that bans the
hexamers ``AATAAA`` and ``ATTAAA``. That preset is *strict but blunt*: a bare
hexamer is *not* a functional poly(A) signal, and ``AATAAA`` occurs by chance
roughly once per 4 kb of random AT-balanced sequence, so banning every occurrence
throws away synonymous choices to avoid sites that would never be used.

The real signal is **bipartite**, and both halves are bound by different proteins:

* the **AAUAAA hexamer**, recognised by CPSF, sitting 10-30 nt upstream of the
  cleavage site, and
* a **downstream sequence element (DSE)** -- a poorly-conserved U-rich or GU-rich
  stretch after the cleavage site -- recognised by CstF.

Cleavage happens between them, and a hexamer with no downstream partner is not
processed. :class:`FunctionalPolyASignalConstraint` encodes that pairing: it
forbids a hexamer **only when a U/GU-rich element follows it** in the window where
a DSE would have to be. That is strictly more permissive than the blunt preset and
strictly closer to the biology -- so the two are offered as a *choice*, not a
replacement (see ``polya_mode`` in :class:`~bt4.pipeline.optimize.OptimizeConfig`).

Honest scope. This is a **structural** rule built from the consensus architecture
above, not a calibrated poly(A)-site predictor: the DSE has no strong consensus,
so the U/GU-rich test here is a documented heuristic with tunable thresholds, and
a run under it makes no claim about measured cleavage efficiency. Real poly(A)-site
prediction is a learned-model problem BT4 does not pretend to solve.

Tractability. The footprint spans the hexamer plus the whole DSE search window
(~45 nt by default), so ``context_len`` is far too wide for the exact trellis --
the same wall the windowed-GC rule hits. The pipeline therefore routes this rule
to the refinement layer and reports its residuals honestly, rather than silently
capping the context (CLAUDE.md §10.1).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["DEFAULT_POLYA_HEXAMERS", "FunctionalPolyASignalConstraint"]

DEFAULT_POLYA_HEXAMERS: tuple[str, ...] = ("AATAAA", "ATTAAA")
"""The canonical hexamer and its commonest single-base variant.

``AAUAAA`` accounts for the majority of mammalian poly(A) sites and ``AUUAAA`` for
most of the remainder; rarer variants exist but are progressively weaker, and
banning them all is what makes the blunt preset over-restrictive.
"""


@dataclass(frozen=True, slots=True)
class FunctionalPolyASignalConstraint:
    """Forbid a poly(A) hexamer only when a downstream U/GU-rich element follows.

    Attributes:
        hexamers: Upstream hexamers to look for (upper-case ACGT).
        dse_start: First offset **after the hexamer's last base** at which the DSE
            search begins. Defaults to 10, since cleavage occurs 10-30 nt after
            the hexamer and the DSE lies beyond the cleavage site.
        dse_end: Offset after the hexamer's last base at which the DSE search
            stops. Defaults to 40.
        dse_window: Length of the sliding window tested for U/GU richness.
        min_u: Minimum ``T`` count in a ``dse_window`` slice for it to count as
            U-rich.
        min_gu: Minimum ``G``+``T`` count in a ``dse_window`` slice for it to count
            as GU-rich (it must also contain at least :attr:`min_gu_u` ``T`` bases,
            so a pure G-run is not mistaken for a GU-rich element).
        min_gu_u: Minimum ``T`` count within a GU-rich slice.
    """

    hexamers: tuple[str, ...] = DEFAULT_POLYA_HEXAMERS
    dse_start: int = 10
    dse_end: int = 40
    dse_window: int = 10
    min_u: int = 6
    min_gu: int = 8
    min_gu_u: int = 3
    name: str = field(default="polya_signal", init=False)

    def __post_init__(self) -> None:
        """Validate the window geometry and thresholds.

        Raises:
            ValueError: If a hexamer is empty, the DSE window is not positive, the
                search range is empty, or a threshold exceeds the window length.
        """
        if not self.hexamers or any(not h for h in self.hexamers):
            raise ValueError("hexamers must be a non-empty tuple of non-empty motifs")
        if self.dse_window < 1:
            raise ValueError(f"dse_window must be >= 1, got {self.dse_window!r}")
        if self.dse_end <= self.dse_start:
            raise ValueError(
                f"dse_end ({self.dse_end}) must exceed dse_start ({self.dse_start})"
            )
        for label, value in (
            ("min_u", self.min_u),
            ("min_gu", self.min_gu),
            ("min_gu_u", self.min_gu_u),
        ):
            if value > self.dse_window:
                raise ValueError(
                    f"{label} ({value}) cannot exceed dse_window ({self.dse_window})"
                )

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`.

        The rule genuinely has bounded context (see :meth:`context_len`); it is the
        *width* of that context, not its boundedness, that keeps it out of the
        trellis. The pipeline routes it accordingly.
        """
        return Scope.LOCAL

    def context_len(self) -> int:
        """Trailing context needed: the hexamer plus the whole DSE search window."""
        return max(len(h) for h in self.hexamers) + self.dse_end - 1

    def _dse_hit(self, window: str) -> int | None:
        """Return the offset of the first U/GU-rich slice in ``window``, else ``None``."""
        for i in range(0, len(window) - self.dse_window + 1):
            slice_ = window[i : i + self.dse_window]
            t = slice_.count("T")
            if t >= self.min_u:
                return i
            if t >= self.min_gu_u and t + slice_.count("G") >= self.min_gu:
                return i
        return None

    def _findings(self, seq: str) -> Iterator[tuple[int, int, str, str]]:
        """Yield ``(hexamer_start, dse_end_index, hexamer, kind)`` for each real signal."""
        n = len(seq)
        for hexamer in self.hexamers:
            start = seq.find(hexamer)
            while start != -1:
                after = start + len(hexamer)
                lo, hi = after + self.dse_start, min(after + self.dse_end, n)
                if hi - lo >= self.dse_window:
                    offset = self._dse_hit(seq[lo:hi])
                    if offset is not None:
                        dse_lo = lo + offset
                        yield (start, dse_lo + self.dse_window, hexamer, "U/GU-rich")
                start = seq.find(hexamer, start + 1)

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` completes a functional signal.

        A signal is charged to this extension only when the DSE that completes it
        ends inside ``next_codon`` -- a signal already finished inside the feasible
        prefix is not this codon's fault.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        cl = self.context_len()
        tail = prefix[-cl:] if cl > 0 else ""
        window = tail + next_codon
        boundary = len(tail)
        return all(dse_end <= boundary for _s, dse_end, _h, _k in self._findings(window))

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per hexamer that has a downstream U/GU element.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` spanning the hexamer through
            the end of the downstream element that completes it.
        """
        seq = dna.upper()
        for start, dse_end, hexamer, kind in self._findings(seq):
            yield Violation(
                constraint="polya_signal",
                severity=Severity.HARD,
                start=start,
                end=dse_end,
                detail=(
                    f"functional poly(A) signal: {hexamer} at {start} with a "
                    f"{kind} downstream element ending at {dse_end}"
                ),
            )
