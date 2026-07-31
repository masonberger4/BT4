"""Restriction-site avoidance as a LOCAL hard constraint (CLAUDE.md §6).

A :class:`RestrictionSiteConstraint` bans the recognition sequences of one or
more restriction enzymes (and any explicitly supplied extra sites) from the
coding sequence. Sites are double-stranded, so each site's IUPAC-aware reverse
complement is *always* banned too -- a palindromic site (the common case) simply
folds back onto itself, while an asymmetric site contributes both strands.

Sites may carry IUPAC ambiguity codes (e.g. HinfI ``GANTC``, DraIII
``CACNNNGTG``); matching goes through :mod:`bt4.constraints.iupac`. The
constraint is purely hard (``penalty`` is ``0.0``) and its ``context_len`` is
``maxlen - 1`` -- exactly enough for any banned site to straddle the seam into
the incoming codon, so ``ok_suffix`` and ``validate`` agree (invariant #3).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.constraints.iupac import find_iupac, is_iupac, reverse_complement_iupac
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope
from bt4.domain.sequence import validate_dna

__all__ = ["ENZYMES", "RestrictionSiteConstraint", "available_enzymes"]

ENZYMES: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "NdeI": "CATATG",
    "NcoI": "CCATGG",
    "EcoRV": "GATATC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "HinfI": "GANTC",
    "DraIII": "CACNNNGTG",
}
"""Catalog of common restriction enzymes to their (IUPAC) recognition sites.

Sites are textbook-correct 5'->3' recognition sequences. ``HinfI`` (``GANTC``)
and ``DraIII`` (``CACNNNGTG``) are degenerate, exercising the IUPAC matcher.
"""


def available_enzymes() -> tuple[str, ...]:
    """Return the catalog's enzyme names in sorted order.

    Returns:
        A deterministic (alphabetically sorted) tuple of enzyme names.
    """
    return tuple(sorted(ENZYMES))


@dataclass(frozen=True, slots=True)
class RestrictionSiteConstraint:
    """Ban restriction recognition sites (both strands) from a coding sequence.

    Built from named enzymes, explicit extra sites, or both. Every resolved site
    contributes its IUPAC reverse complement as well, because restriction sites
    are double-stranded.

    Attributes:
        enzymes: Names resolved against :data:`ENZYMES`.
        extra_sites: Additional IUPAC recognition sequences (upper-cased and
            validated at construction).
    """

    enzymes: tuple[str, ...] = ()
    extra_sites: tuple[str, ...] = ()
    name: str = field(default="restriction_site", init=False)
    _sites: tuple[str, ...] = field(init=False, repr=False, compare=False)
    _maxlen: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Resolve, validate, and reverse-complement the banned site set.

        Raises:
            ValueError: If an enzyme name is unknown or an extra site is not a
                valid non-empty IUPAC string.
        """
        collected: set[str] = set()
        for enzyme in self.enzymes:
            if enzyme not in ENZYMES:
                raise ValueError(
                    f"unknown enzyme {enzyme!r}; known enzymes: {available_enzymes()}"
                )
            collected.add(ENZYMES[enzyme])
        for raw in self.extra_sites:
            site = raw.strip().upper()
            if not is_iupac(site):
                raise ValueError(f"invalid IUPAC restriction site: {raw!r}")
            collected.add(site)
        # Restriction sites are double-stranded: always ban each reverse complement.
        with_rc = set(collected)
        for site in collected:
            with_rc.add(reverse_complement_iupac(site))
        stored = tuple(sorted(with_rc))
        maxlen = max((len(site) for site in stored), default=0)
        object.__setattr__(self, "_sites", stored)
        object.__setattr__(self, "_maxlen", maxlen)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``max(0, maxlen - 1)`` - enough for a site to cross the seam."""
        return max(0, self._maxlen - 1)

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff a banned site occurrence overlaps ``next_codon``.

        The window is ``prefix[-(maxlen-1):] + next_codon``; an occurrence counts
        only when its exclusive end index lies beyond the start of ``next_codon``
        (so it actually involves the codon being appended). With no sites the
        constraint is inert and this always returns ``True``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.

        Returns:
            ``True`` when appending ``next_codon`` introduces no banned site.
        """
        if not self._sites:
            return True
        tail = self._maxlen - 1
        window = (prefix[-tail:] if tail > 0 else "") + next_codon
        seam = len(window) - len(next_codon)
        for site in self._sites:
            for idx in find_iupac(window, site):
                if idx + len(site) > seam:
                    return False
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per banned site occurrence in ``dna``.

        Args:
            dna: The whole coding sequence to audit (ACGT, case-insensitive).

        Yields:
            A :class:`~bt4.domain.result.Violation` for every occurrence of every
            banned site, in a deterministic (site-sorted, then position) order.
        """
        if not dna:
            return
        seq = validate_dna(dna)
        for site in self._sites:
            for idx in find_iupac(seq, site):
                yield Violation(
                    constraint="restriction_site",
                    severity=Severity.HARD,
                    start=idx,
                    end=idx + len(site),
                    detail=f"restriction site {site} at {idx}",
                )
