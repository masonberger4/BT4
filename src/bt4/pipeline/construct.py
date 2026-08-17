"""Audit the whole construct: insert plus backbone (CLAUDE.md §1, §6).

Every other audit in BT4 answers "is the *coding sequence* clean?". This one
answers the question a person actually has before they order DNA: **is the thing I
am about to build clean?** Those differ in ways that only appear once the CDS is
placed in its construct:

* a forbidden motif, restriction site, homopolymer or GC extreme formed **across a
  junction** between the insert and the sequence flanking it;
* a restriction site that is fine in the insert but **no longer unique in the
  plasmid**, so the enzyme a user planned to linearize with now cuts twice;
* a repeat that only exists because the insert happens to echo the backbone.

Two things keep the report honest rather than noisy:

* **Coordinates are construct coordinates**, and every finding says whether it lies
  in the insert, in the flank, or across the junction. A user can act on the first,
  must redesign around the second, and needs to know which they are looking at.
* **Repeats by construction are maskable.** AAV ITRs are a 145-nt palindrome and
  lentiviral LTRs are duplicated by design, so an unmasked repeat audit of exactly
  the two vector systems that most need one would be dominated by findings the user
  cannot change. ``ConstructContext.masked_spans`` excludes them, and the count of
  what was excluded is reported rather than silently dropped.

This is **reporting only**. It never edits a sequence and never changes what the
optimizer delivered; the CDS-level rules remain the things the solver enforces.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt4.constraints.iupac import find_iupac, reverse_complement_iupac
from bt4.constraints.restriction import resolve_enzyme
from bt4.domain import ConstructContext, Severity, Violation, validate_dna
from bt4.pipeline.optimize import (
    OptimizeConfig,
    _build_constraints,
    _build_global_constraints,
)

__all__ = ["ConstructAudit", "EnzymeOccurrence", "audit_construct"]


@dataclass(frozen=True, slots=True)
class EnzymeOccurrence:
    """Where one requested enzyme's recognition site occurs in the construct.

    Attributes:
        enzyme: The enzyme name as the user asked for it.
        site: Its recognition sequence (IUPAC).
        positions: Every occurrence start in construct coordinates, both strands.
        in_cds: How many of those lie inside the designed coding sequence.
        unique: ``True`` when the site occurs exactly once in the whole construct --
            the property that matters when the enzyme is meant to linearize it.
    """

    enzyme: str
    site: str
    positions: tuple[int, ...]
    in_cds: int

    @property
    def count(self) -> int:
        """Total occurrences in the construct."""
        return len(self.positions)

    @property
    def unique(self) -> bool:
        """True when the enzyme cuts the construct exactly once."""
        return len(self.positions) == 1


@dataclass(frozen=True, slots=True)
class ConstructAudit:
    """The result of auditing an assembled construct.

    Attributes:
        construct: The assembled sequence (upstream + CDS + downstream).
        cds_start: Index of the CDS's first base in ``construct``.
        cds_end: Index one past the CDS's last base.
        violations: Findings in **construct** coordinates, masked spans removed.
        junction_violations: How many of those cross a junction -- the findings a
            CDS-only audit structurally could not see.
        masked_violations: How many findings were dropped because they lie in a
            user-declared masked span (ITRs, LTRs).
        enzymes: Per-enzyme occurrence report (uniqueness).
    """

    construct: str
    cds_start: int
    cds_end: int
    violations: tuple[Violation, ...]
    junction_violations: int
    masked_violations: int
    enzymes: tuple[EnzymeOccurrence, ...]

    @property
    def is_clean(self) -> bool:
        """True when no hard violation survived masking."""
        return not any(v.severity is Severity.HARD for v in self.violations)


def _crosses(violation: Violation, cds_start: int, cds_end: int) -> bool:
    """True when ``violation`` spans a boundary between the CDS and a flank."""
    starts_outside = violation.start < cds_start or violation.start >= cds_end
    ends_inside = cds_start < violation.end <= cds_end
    starts_inside = cds_start <= violation.start < cds_end
    ends_outside = violation.end > cds_end
    return (starts_outside and ends_inside) or (starts_inside and ends_outside)


def _region(violation: Violation, cds_start: int, cds_end: int) -> str:
    if _crosses(violation, cds_start, cds_end):
        return "junction"
    return "insert" if violation.start >= cds_start and violation.end <= cds_end else "flank"


def audit_construct(cds: str, config: OptimizeConfig) -> ConstructAudit:
    """Audit ``cds`` inside the construct described by ``config.context``.

    Runs every rule the config sets over the **assembled** construct, so findings
    that only exist once the insert is placed in its flanks are surfaced. With no
    context this degrades to auditing the bare CDS, which is exactly what the
    CDS-level report already says -- it is never wrong, just not informative.

    Args:
        cds: The delivered coding sequence.
        config: The run's configuration; its ``context`` supplies the flanks, its
            constraints supply the rules, and its enzymes drive the uniqueness
            report.

    Returns:
        A :class:`ConstructAudit` in construct coordinates.

    Raises:
        ValueError: On a non-ACGT coding sequence.
    """
    sequence = validate_dna(cds)
    context = config.context if config.context is not None else ConstructContext()
    construct = context.assemble(sequence)
    cds_start = context.cds_offset
    cds_end = cds_start + len(sequence)

    # Rebuild the rules *unseeded*: we are validating the assembled construct
    # directly, so a seeded wrapper would prepend the leader a second time.
    bare = OptimizeConfig(
        **{
            **{f: getattr(config, f) for f in OptimizeConfig.__dataclass_fields__},
            "context": None,
        }
    )
    rules = [*_build_constraints(bare), *_build_global_constraints(bare)]

    kept: list[Violation] = []
    masked = 0
    junction = 0
    for rule in rules:
        for violation in rule.validate(construct):
            if context.is_masked(violation.start, violation.end):
                masked += 1
                continue
            region = _region(violation, cds_start, cds_end)
            if region == "junction":
                junction += 1
            kept.append(
                Violation(
                    constraint=violation.constraint,
                    severity=violation.severity,
                    start=violation.start,
                    end=violation.end,
                    detail=f"[{region}] {violation.detail}",
                )
            )

    enzymes = tuple(
        _enzyme_occurrence(name, construct, cds_start, cds_end)
        for name in config.restriction_enzymes
    )
    return ConstructAudit(
        construct=construct,
        cds_start=cds_start,
        cds_end=cds_end,
        violations=tuple(kept),
        junction_violations=junction,
        masked_violations=masked,
        enzymes=enzymes,
    )


def _enzyme_occurrence(
    name: str, construct: str, cds_start: int, cds_end: int
) -> EnzymeOccurrence:
    """Count one enzyme's sites in ``construct`` on both strands."""
    site = resolve_enzyme(name)
    hits = set(find_iupac(construct, site))
    rc = reverse_complement_iupac(site)
    if rc != site:  # a palindromic site would otherwise be counted twice
        hits.update(find_iupac(construct, rc))
    positions = tuple(sorted(hits))
    in_cds = sum(1 for p in positions if cds_start <= p < cds_end)
    return EnzymeOccurrence(enzyme=name, site=site, positions=positions, in_cds=in_cds)
