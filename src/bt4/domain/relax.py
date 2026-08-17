"""Graceful constraint relaxation (CLAUDE.md §4.2, invariant #6).

A hard constraint can make a back-translation instance infeasible. The classic
case is ``avoid_internal_start`` on a protein whose internal Met sits in a
synonymously-forced strong-Kozak context: Met is a single codon, so if the
flanking residues force a purine at -3 and a G at +4, *no* synonymous choice
escapes the strong context. §4.2 requires such a rule to **degrade gracefully**
("a graceful degradation, not a dead-end") rather than abort.

This module is that path. :class:`SoftConstraint` turns a hard veto into a
visible-but-non-binding soft rule: the trellis never rejects a codon on the
relaxed constraint's account (:meth:`SoftConstraint.ok_suffix` is always ``True``),
while :meth:`SoftConstraint.validate` still surfaces the underlying occurrences --
downgraded to :attr:`~bt4.domain.result.Severity.SOFT` -- so nothing is hidden. A
solver that relaxes a constraint reports it in the ``OptimalityCertificate``
(``RELAXED``) and audits the delivered sequence against the *original* hard rule,
so residuals stay honestly reported.

Relaxation is **opt-in per constraint**: a constraint declares itself gracefully
degradable by defining a ``relax()`` method (see
:class:`~bt4.constraints.kozak.InternalStartConstraint`). A constraint that does
not is never silently dropped -- an instance infeasible under it raises, naming
the culprit. :func:`is_relaxable` tests the marker; :func:`relax_constraint`
produces the soft form (using the constraint's own ``relax()`` when present, else
the generic :class:`SoftConstraint` wrapper).

Depends only on :mod:`bt4.domain`.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import cast

from bt4.domain.contracts import Constraint
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["SoftConstraint", "is_relaxable", "relax_constraint"]


@dataclass(frozen=True, slots=True)
class SoftConstraint:
    """A hard :class:`~bt4.domain.contracts.Constraint` relaxed to a soft one.

    Wraps ``inner`` so that its hard veto no longer binds the solver, while its
    occurrences remain visible in an audit. See the module docstring for the full
    contract; the load-bearing points are that ``ok_suffix`` never vetoes and
    ``validate`` downgrades HARD occurrences to SOFT rather than dropping them.
    """

    inner: Constraint

    @property
    def name(self) -> str:
        """The wrapped constraint's stable identifier (unchanged)."""
        return self.inner.name

    def scope(self) -> Scope:
        """The wrapped constraint's locality class (unchanged)."""
        return self.inner.scope()

    def context_len(self) -> int:
        """The wrapped constraint's trailing-context need (unchanged)."""
        return self.inner.context_len()

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Never veto: a relaxed constraint imposes no hard feasibility bar."""
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """A nominal soft cost where the wrapped hard rule *would* have vetoed.

        The exact DP ignores ``penalty`` (it uses only ``ok_suffix``), so this is
        a well-defined soft cost for any penalty-aware consumer, not a solver hook.
        """
        return 0.0 if self.inner.ok_suffix(prefix, next_codon) else 1.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield the wrapped constraint's occurrences, HARD downgraded to SOFT."""
        for v in self.inner.validate(dna):
            if v.severity is Severity.HARD:
                yield Violation(
                    constraint=v.constraint,
                    severity=Severity.SOFT,
                    start=v.start,
                    end=v.end,
                    detail=v.detail,
                )
            else:
                yield v


def is_relaxable(constraint: Constraint) -> bool:
    """Return ``True`` iff ``constraint`` opts in to graceful relaxation.

    A constraint opts in by defining a callable ``relax()`` method. Constraints
    without one are never auto-relaxed (an infeasible instance under them raises,
    naming the culprit) -- so, e.g., a restriction-site ban is never silently
    dropped, while ``avoid_internal_start`` is.
    """
    return callable(getattr(constraint, "relax", None))


def relax_constraint(constraint: Constraint) -> Constraint:
    """Return a soft version of ``constraint`` for graceful degradation (§4.2).

    Uses the constraint's own ``relax()`` when it defines one (a constraint may
    tailor how it degrades), else wraps it in a generic :class:`SoftConstraint`.
    """
    relax = getattr(constraint, "relax", None)
    if callable(relax):
        return cast(Constraint, relax())
    return SoftConstraint(constraint)
