"""Internal strong-Kozak ATG constraint (CLAUDE.md §6, invariant #3).

An internal ATG in a *strong* translation-initiation (Kozak) context can seed
spurious ribosomal re-initiation. The classic strength signal is LOCAL and
tightly bounded: relative to the A of the ATG (numbered +1), a **purine (A or G)
at position -3** and a **G at position +4** (the base immediately after the ATG's
G) form the strong consensus ``gccRccATGG``. Only those two flanks matter; the
-2 and -1 bases are don't-cares. Because the whole signal fits in the fixed
window ``[-3 .. +4]`` it is a genuine :attr:`~bt4.domain.scope.Scope.LOCAL`
constraint, solvable exactly in the codon trellis.

:class:`InternalStartConstraint` forbids every *internal* ATG (an ``ATG``
substring whose A-index is ``>= min_start``, so the real start codon at index 0
is skipped) that sits in a strong context. "Strong" is configurable; by default
it means *both* a purine at -3 **and** a G at +4. The ATG is the literal
substring ``"ATG"`` anywhere in frame or across codon seams -- the honest
definition treats all ``ATG`` substrings uniformly, whether they are Met codons
or arise from other codons abutting.

Why ``context_len`` suffices (invariant #3). The forbidden footprint spans the 7
bases ``[-3 .. +4]`` (indices ``a-3 .. a+3`` around the A at ``a``), with the
rightmost constrained base at +4 (``a+3``). Exactly like a 7-mer forbidden
motif, the veto must fire when that rightmost base is the *first* base of the
incoming codon; the six preceding bases then live in the prefix, so
``context_len == 6``. Any new forbidden occurrence whose rightmost constrained
base lands inside the incoming codon has its whole footprint visible in
``prefix[-6:] + next_codon``, and ``ok_suffix`` vetoes precisely those -- so a
sequence built respecting ``ok_suffix`` carries zero hard violations under
``validate``.

Scope honesty (CLAUDE.md §4.4, §10.4). The Kozak-strength ATG signal above is
the *only* thing this constraint enforces. The companion idea of **uORF pairing**
-- an out-of-frame internal ATG paired with a downstream in-frame stop -- is a
genuinely NON-LOCAL relationship: the ATG and its stop can be arbitrarily far
apart, so it does not fit a bounded-context ``ok_suffix`` and is deliberately
NOT modelled here. uORF pairing is deferred to the Phase 3 refinement layer; it
is not shipped by this file and is not faked with a padded window.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["InternalStartConstraint"]

_PURINES = frozenset("AG")
_CONTEXT_LEN = 6


@dataclass(frozen=True, slots=True)
class InternalStartConstraint:
    """Forbid internal ATG triplets that sit in a strong Kozak context.

    A strong-Kozak ATG is forbidden when *all enabled* context conditions hold
    for an ATG whose A-index is ``>= min_start``. With both flags ``False`` the
    context is ignored and every internal ATG is forbidden.

    Attributes:
        require_purine_m3: When ``True`` (default), the strong context requires a
            purine (A or G) at position -3. A missing -3 base (the ATG is too
            close to the 5' end for a -3 flank to exist) counts as *not* a
            purine, so such a weak-context ATG is allowed.
        require_g_p4: When ``True`` (default), the strong context requires a G at
            position +4 (the base immediately after the ATG). A missing +4 base
            (the ATG is at the very 3' end) counts as *not* a G, so such an ATG
            is allowed.
        min_start: Only ATGs whose A-index is ``>= min_start`` are "internal".
            The default of ``3`` skips the first codon, which is typically the
            real start.
    """

    require_purine_m3: bool = True
    require_g_p4: bool = True
    min_start: int = 3
    name: str = field(default="internal_start", init=False)

    def __post_init__(self) -> None:
        """Validate the start bound.

        Raises:
            ValueError: If ``min_start`` is negative.
        """
        if self.min_start < 0:
            raise ValueError(f"min_start must be >= 0, got {self.min_start!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``6`` -- the trailing context ``ok_suffix`` needs.

        The forbidden footprint is 7 bases (-3 through +4); catching one whose
        +4 base is the first base of the incoming codon needs the six preceding
        bases from the prefix.
        """
        return _CONTEXT_LEN

    def _is_strong(self, seq: str, w: int, n: int) -> bool:
        """Return ``True`` iff the ATG at index ``w`` in ``seq`` is strong.

        Only the enabled conditions are checked, using the bases of ``seq``. An
        out-of-range flank (``w-3 < 0`` or ``w+3 >= n``) is treated as absent and
        so fails its condition. This reads exactly the flanks the config needs;
        the caller is responsible for the ``min_start`` (internal) test.

        Args:
            seq: The string to inspect (whole DNA, or an ``ok_suffix`` window).
            w: Index within ``seq`` of the ATG's A.
            n: ``len(seq)``.
        """
        if self.require_purine_m3 and (w - 3 < 0 or seq[w - 3] not in _PURINES):
            return False
        if self.require_g_p4 and (w + 3 >= n or seq[w + 3] != "G"):  # noqa: SIM103
            return False
        return True

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` adds a strong internal ATG.

        Only the last ``context_len`` (6) chars of ``prefix`` are inspected. An
        occurrence is charged to this extension only when its rightmost
        constrained base (the +4 base, or the ATG's own G when ``require_g_p4`` is
        off) lands inside ``next_codon`` -- an ATG whose footprint ends inside the
        already-feasible prefix is not this codon's fault. ``prefix`` is assumed
        feasible, so any new strong internal ATG must involve ``next_codon``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        p = len(prefix)
        cl = self.context_len()
        tail = prefix[-cl:]
        base0 = p - len(tail)
        window = tail + next_codon
        n = len(window)
        right = 3 if self.require_g_p4 else 2
        w = window.find("ATG")
        while w != -1:
            abs_a = base0 + w
            if abs_a >= self.min_start and self._is_strong(window, w, n) and abs_a + right >= p:
                return False
            w = window.find("ATG", w + 1)
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per strong internal ATG in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` (anchored on the ATG triplet)
            for every internal ATG in a strong Kozak context.
        """
        seq = dna.upper()
        n = len(seq)
        a = seq.find("ATG")
        while a != -1:
            if a >= self.min_start and self._is_strong(seq, a, n):
                m3 = seq[a - 3] if a - 3 >= 0 else "-"
                p4 = seq[a + 3] if a + 3 < n else "-"
                yield Violation(
                    constraint="internal_start",
                    severity=Severity.HARD,
                    start=a,
                    end=a + 3,
                    detail=f"internal ATG at {a} in strong Kozak context (-3={m3}, +4={p4})",
                )
            a = seq.find("ATG", a + 1)
