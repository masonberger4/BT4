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

import difflib
import json
from collections.abc import Iterator
from dataclasses import dataclass, field
from importlib.resources import files
from types import MappingProxyType

from bt4.constraints.iupac import find_iupac, is_iupac, reverse_complement_iupac
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope
from bt4.domain.sequence import validate_dna

__all__ = [
    "ENZYMES",
    "RestrictionSiteConstraint",
    "available_enzymes",
    "enzyme_provenance",
    "enzyme_suggestions",
    "resolve_enzyme",
    "unknown_enzyme_message",
]

_DATA_PACKAGE = "bt4.constraints.data"
_CATALOG_FILE = "rebase_enzymes.tsv"
_PROVENANCE_FILE = "rebase_enzymes.provenance.json"

# How many near-miss suggestions an unknown-enzyme error offers. The catalog has
# hundreds of entries, so listing them all would bury the answer in noise.
_SUGGESTIONS = 5


def _load_catalog() -> dict[str, str]:
    """Parse the bundled REBASE-derived ``enzyme<TAB>site`` catalog."""
    text = files(_DATA_PACKAGE).joinpath(_CATALOG_FILE).read_text(encoding="utf-8")
    catalog: dict[str, str] = {}
    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("enzyme\t"):
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError(f"{_CATALOG_FILE}:{line_no}: expected 2 columns")
        name, site = parts[0].strip(), parts[1].strip().upper()
        if not name or not is_iupac(site):
            raise ValueError(f"{_CATALOG_FILE}:{line_no}: bad entry {line!r}")
        catalog[name] = site
    if not catalog:
        raise ValueError(f"{_CATALOG_FILE} is empty")
    return catalog


ENZYMES: MappingProxyType[str, str] = MappingProxyType(_load_catalog())
"""Catalog mapping restriction-enzyme name to its IUPAC recognition site.

Derived from a version-pinned **REBASE** release (Roberts et al.) by
``scripts/build_enzyme_catalog.py`` -- every commercially available restriction
enzyme with a single fully-specified recognition sequence -- rather than
hand-typed, so every site is checkable against a cited source and the shipped bytes are
content-hashed in :func:`enzyme_provenance`. Isoschizomers are included, so a
user can name the enzyme they actually own (``KpnI`` and ``Acc65I`` both resolve
to ``GGTACC``).

Read-only: the catalog is shipped data, not a mutable registry. BT4 models the
*recognition sequence* only -- not cut position, star activity, methylation
sensitivity, or buffer conditions.

**What the catalog covers.** Every commercially available *restriction* enzyme
REBASE lists with a single fully-specified site: Type II, Type II
bifunctional/IIG, Type II modification-dependent/IIM (``DpnI`` among them), and
one Type III. Methyltransferases and homing endonucleases are excluded -- they do
not define a site a coding sequence must avoid. The modification-dependent ones
are kept deliberately: avoiding ``DpnI``'s ``GATC`` is a mainstream design goal
precisely *because* a plasmid grown in a standard dam+ strain is Dam-methylated
and is cut by it.

**Some genuine sites are highly degenerate.** ``MspJI`` is ``CNNR`` and ``BslI``
is ``CCNNNNNNNGG`` -- their true recognition sequences, not placeholders. Banning
a near-wildcard inside a coding sequence can be genuinely unsatisfiable, and BT4
reports that honestly: an infeasible run raises
:class:`~bt4.optimize.InfeasibleError` naming ``restriction_site``, and a
feasible one really does contain zero occurrences. What never happens is the
third option BT3 would have taken -- quietly accepting the request and returning
a sequence that still contains the site.
"""


def available_enzymes() -> tuple[str, ...]:
    """Return the catalog's enzyme names in sorted order.

    Returns:
        A deterministic (alphabetically sorted) tuple of enzyme names.
    """
    return tuple(sorted(ENZYMES))


def enzyme_provenance() -> dict[str, object]:
    """Return the enzyme catalog's provenance sidecar as a plain mapping.

    Records the REBASE version and URL, the source file's SHA-256, the selection
    rule, the shipped catalog's own digest, and REBASE's citation/licensing
    terms -- so a third party can re-derive and re-verify the catalog rather than
    taking its numbers on trust (CLAUDE.md §8).
    """
    text = files(_DATA_PACKAGE).joinpath(_PROVENANCE_FILE).read_text(encoding="utf-8")
    data: dict[str, object] = json.loads(text)
    return data


def resolve_enzyme(name: str) -> str:
    """Resolve an enzyme name to its recognition site, case-insensitively.

    Args:
        name: An enzyme name; matched exactly first, then case-insensitively
            (so ``ecori`` finds ``EcoRI``).

    Returns:
        The IUPAC recognition sequence.

    Raises:
        ValueError: If no catalog entry matches. The message offers the closest
            names rather than dumping the whole catalog -- with hundreds of
            entries, a full listing hides the answer instead of giving it -- and
            it shows each suggestion's **site**, because the suggestions are
            matched on *spelling*, not on recognition sequence. A name-similar
            enzyme usually cuts something completely different (``NotI`` is
            ``GCGGCCGC``, ``NcoI`` is ``CCATGG``), so presenting a bare list
            would invite a user to accept a substitute that does not ban the
            site they care about -- and the run would then report
            proven-optimal with zero violations while their real site remained.
    """
    if name in ENZYMES:
        return ENZYMES[name]
    folded = {key.lower(): key for key in ENZYMES}
    hit = folded.get(name.strip().lower())
    if hit is not None:
        return ENZYMES[hit]
    raise ValueError(unknown_enzyme_message(name))


def enzyme_suggestions(name: str, limit: int = _SUGGESTIONS) -> tuple[str, ...]:
    """Return catalog names spelled similarly to ``name`` (NOT site-similar).

    A pure fuzzy *name* match, exposed so a frontend can show the same list the
    engine would -- and, like the engine, must label it as spelling-matched and
    show each site, since these enzymes generally do not share a recognition
    sequence with each other or with what the user typed.
    """
    return tuple(
        difflib.get_close_matches(name, available_enzymes(), n=limit, cutoff=0.6)
    )


def unknown_enzyme_message(name: str) -> str:
    """Build the shared "unknown enzyme" explanation used by every surface.

    Kept in one place so the CLI, the API and BT4 Studio cannot drift into
    telling the user different things about the same miss.
    """
    close = enzyme_suggestions(name)
    parts = [f"unknown enzyme {name!r}"]
    if close:
        shown = ", ".join(f"{hit} ({ENZYMES[hit]})" for hit in close)
        parts.append(
            f"closest catalog names by SPELLING, not by recognition site: {shown}"
        )
    parts.append(
        f"{len(ENZYMES)} enzymes in the catalog (see `bt4 enzymes`). If BT4 does "
        "not carry yours, ban its recognition sequence directly instead of "
        "substituting a similarly-named enzyme"
    )
    return "; ".join(parts)


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
            collected.add(resolve_enzyme(enzyme))
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
