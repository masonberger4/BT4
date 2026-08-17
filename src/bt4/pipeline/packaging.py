"""Viral-vector packaging headroom -- a report, never a lever (CLAUDE.md §6).

An AAV capsid packages a genome of roughly its wild-type size, and a transgene
cassette that overshoots is packaged **truncated**, which is a silent and
expensive failure: the prep looks fine by titre and the construct simply does not
work. Lentiviral vectors are far more forgiving but lose titre as the genome grows.

BT4 cannot fix this. The CDS length is fixed by the protein, and the promoter,
UTRs, poly(A) and ITRs are the user's design -- there is no synonymous choice that
makes a cassette shorter. So this module **reports and never optimizes**: given
what the run actually knows (the designed CDS, plus any construct context the user
supplied), it states the size accounted for, the headroom against a stated limit,
and -- importantly -- **what it could not see**.

That last part is the honest core. A CDS-only run knows nothing about the ITRs,
promoter or poly(A) signal, so its "headroom" would be a flattering fiction if it
did not say so. :attr:`PackagingReport.counted` names exactly which parts were
included, and :attr:`PackagingReport.complete` is ``False`` whenever the construct
is only partly known, so a caller can never present a partial measurement as a
verdict.

Limits are conventions from vector-design practice, not measurements BT4 made:
the AAV figure is the classic ~4.7 kb wild-type genome size (ITR-to-ITR, ITRs
included), and packaging efficiency degrades before the limit rather than at it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from bt4.domain.context import ConstructContext

__all__ = [
    "PACKAGING_LIMITS",
    "PackagingReport",
    "packaging_report",
]

PACKAGING_LIMITS: Mapping[str, int] = MappingProxyType(
    {
        # ITR-to-ITR, ITRs included. Oversized genomes are packaged truncated.
        "aav": 4700,
        # Lentiviral transfer genomes tolerate more but lose titre as they grow;
        # ~8-10 kb is the usual practical ceiling quoted for the transfer plasmid.
        "lvv": 9000,
    }
)
"""Conventional genome-size ceilings by vector system, in nucleotides."""


@dataclass(frozen=True, slots=True)
class PackagingReport:
    """What a run can honestly say about vector packaging size.

    Attributes:
        system: The vector system the limit came from (e.g. ``"aav"``).
        limit_nt: The conventional ceiling used.
        counted_nt: Total length BT4 could actually account for.
        cds_nt: Length of the designed coding sequence.
        context_nt: Length of supplied flanking sequence (0 when none).
        counted: Human-readable names of the parts included in ``counted_nt``.
        complete: ``True`` only when the caller asserted the supplied context spans
            the whole construct. When ``False`` the headroom is an **upper bound**:
            unseen elements (ITRs, promoter, poly(A)) only make it smaller.
    """

    system: str
    limit_nt: int
    counted_nt: int
    cds_nt: int
    context_nt: int
    counted: tuple[str, ...]
    complete: bool

    @property
    def headroom_nt(self) -> int:
        """Nucleotides remaining against the limit (negative when over)."""
        return self.limit_nt - self.counted_nt

    @property
    def over_limit(self) -> bool:
        """True when what BT4 *can already see* exceeds the limit.

        A ``False`` here is not a pass: unseen elements are not counted unless
        :attr:`complete` is ``True``.
        """
        return self.counted_nt > self.limit_nt

    def summary(self) -> str:
        """One honest sentence, including what was not measured."""
        parts = " + ".join(self.counted)
        head = (
            f"{self.system.upper()} packaging: {self.counted_nt} nt counted "
            f"({parts}) against a ~{self.limit_nt} nt limit"
        )
        if self.over_limit:
            head += f" -- OVER by {-self.headroom_nt} nt"
        else:
            head += f" -- {self.headroom_nt} nt headroom"
        if not self.complete:
            head += (
                ". This counts only what BT4 was given; ITRs, promoter and poly(A) "
                "are not included, so the real headroom is smaller"
            )
        return head + "."


def packaging_report(
    cds: str,
    *,
    system: str = "aav",
    context: ConstructContext | None = None,
    complete: bool = False,
    limit_nt: int | None = None,
) -> PackagingReport:
    """Report packaging headroom for ``cds`` (and any supplied context).

    Args:
        cds: The designed coding sequence.
        system: Key into :data:`PACKAGING_LIMITS` (``"aav"``, ``"lvv"``).
        context: Construct context, when the user supplied flanking sequence.
        complete: Set ``True`` only if ``context`` spans the *entire* construct
            (ITR to ITR). Left ``False``, the report says plainly that it is an
            upper bound.
        limit_nt: Override the conventional limit.

    Returns:
        The :class:`PackagingReport`.

    Raises:
        ValueError: On an unknown ``system`` with no explicit ``limit_nt``.
    """
    if limit_nt is None:
        if system not in PACKAGING_LIMITS:
            known = ", ".join(sorted(PACKAGING_LIMITS))
            raise ValueError(f"unknown vector system {system!r}; known: {known}")
        limit_nt = PACKAGING_LIMITS[system]
    cds_nt = len(cds)
    context_nt = 0
    counted = ["CDS"]
    if context is not None and not context.is_empty:
        context_nt = len(context.upstream) + len(context.downstream)
        counted.append("supplied context")
    return PackagingReport(
        system=system,
        limit_nt=limit_nt,
        counted_nt=cds_nt + context_nt,
        cds_nt=cds_nt,
        context_nt=context_nt,
        counted=tuple(counted),
        complete=complete,
    )
