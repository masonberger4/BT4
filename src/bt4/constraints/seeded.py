"""Evaluate a LOCAL constraint across the 5' junction (CLAUDE.md §4.2, §5 #3).

The exact-DP trellis starts from an empty prefix, so at codon 0 every LOCAL
constraint is evaluated against *nothing*. A forbidden motif, restriction site,
homopolymer, GC run or GC window formed across the boundary between the user's
5'UTR/backbone and the first codon is therefore structurally unreachable: the
rule is on, the user believes it applies, and the one place a junction defect can
occur is the one place it is not checked.

:class:`SeededConstraint` closes that hole **without changing the ``Constraint``
protocol**. It wraps a constraint and, when the growing prefix is still shorter
than that constraint's own ``context_len``, tops it up from the known upstream
flank. The inner rule is unmodified and unaware; it simply sees the prefix it
would have had if the CDS had not started at index 0.

Two properties make the wrapper safe to drop into the trellis:

* **The trellis accumulator stays CDS-only.** The seed is injected inside
  ``ok_suffix``; nothing is prepended to the sequence being built, so the DP's
  state key, its tie-break, and the delivered DNA are untouched. With no context
  the wrapper is never applied at all, so output stays byte-identical
  (invariant #7).
* **State merging stays exact.** Every prefix in a given trellis layer has the
  same length, so the number of bases borrowed from upstream is the same for all
  of them. Two states sharing a trailing context therefore still share their
  entire future feasibility behaviour, which is what makes merging them sound
  (CLAUDE.md §10.1).

``validate`` is seeded to match, so ``ok_suffix <=> validate`` still holds across
the junction (invariant #3). Violations lying **entirely** inside the upstream
flank are dropped: they are not in the sequence BT4 is designing and no codon
choice can fix them. A violation that straddles the junction is reported with its
span clamped into CDS coordinates, so a reported span always indexes the delivered
sequence.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

from bt4.domain.contracts import Constraint
from bt4.domain.result import Violation
from bt4.domain.scope import Scope

__all__ = ["SeededConstraint", "seed_constraints"]


@dataclass(frozen=True, slots=True)
class SeededConstraint:
    """A LOCAL constraint that also sees the known sequence before the CDS.

    Attributes:
        inner: The constraint being seeded; used unmodified.
        upstream: Known sequence immediately 5' of the CDS. Only the trailing
            ``inner.context_len()`` bases can ever matter, so a long backbone costs
            nothing extra per call.
    """

    inner: Constraint
    upstream: str

    @property
    def name(self) -> str:
        """The wrapped constraint's stable identifier (unchanged)."""
        return self.inner.name

    def scope(self) -> Scope:
        """The wrapped constraint's locality class (unchanged)."""
        return self.inner.scope()

    def context_len(self) -> int:
        """The wrapped constraint's trailing-context need (unchanged).

        The seed tops the prefix up to this length; it never widens the trellis
        state, so seeding costs no state-space growth.
        """
        return self.inner.context_len()

    def _seed(self, prefix: str) -> str:
        """Return ``prefix`` topped up from ``upstream`` to the inner context length."""
        needed = self.inner.context_len() - len(prefix)
        if needed <= 0:
            return prefix  # the prefix already covers everything the rule can read
        return self.upstream[-needed:] + prefix

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Veto exactly as the inner rule would, with the junction visible."""
        return self.inner.ok_suffix(self._seed(prefix), next_codon)

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Soft cost of the inner rule, with the junction visible."""
        return self.inner.penalty(self._seed(prefix), next_codon)

    def validate(self, dna: str) -> Iterator[Violation]:
        """Audit ``dna`` with the junction visible, in CDS coordinates.

        Violations entirely inside the upstream flank are dropped (nothing the
        design can change); a violation straddling the junction is kept with its
        span clamped to the CDS, so every reported span indexes ``dna``.
        """
        seed = self.upstream[-self.inner.context_len() :] if self.inner.context_len() else ""
        if not seed:
            yield from self.inner.validate(dna)
            return
        offset = len(seed)
        for violation in self.inner.validate(seed + dna):
            if violation.end <= offset:
                continue  # wholly upstream: not part of the designed sequence
            yield Violation(
                constraint=violation.constraint,
                severity=violation.severity,
                start=max(0, violation.start - offset),
                end=violation.end - offset,
                detail=(
                    violation.detail
                    if violation.start >= offset
                    else f"{violation.detail} (spans the 5' junction)"
                ),
            )


def seed_constraints(
    constraints: list[Constraint], upstream: str
) -> list[Constraint]:
    """Wrap every constraint so it sees ``upstream``, or return them unchanged.

    Returns the input list untouched when there is no upstream sequence, which is
    what keeps a context-free run byte-identical to before construct context
    existed (invariant #7).
    """
    if not upstream:
        return constraints
    return [SeededConstraint(c, upstream) for c in constraints]
