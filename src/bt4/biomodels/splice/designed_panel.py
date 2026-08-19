"""A panel of **designed synonymous CDS variants** — BT4's actual regime.

Every splice measurement BT4 has made so far is *recall on natural splice sites in
natural genes*: the GENCODE site panel asks whether a model finds annotated sites, and
the ``splicebench2023`` variant panel asks whether it ranks natural single-nucleotide
changes. BT4 emits neither. It emits a **synonymous re-encoding of one protein's coding
sequence**, and the question there is *specificity*: does a model stay correctly silent
on a clean design, and correctly fire on one that creates a cryptic site?

**This panel cannot answer that question, and it does not pretend to.**

There is no splice ground truth for a designed CDS. Nobody has assayed whether these
sequences splice, they appear in no annotation, and a motif is not a site -- which is the
entire distinction between a CNN and the ``gt_ag`` baseline. A panel claiming labels here
would be manufacturing them. So this format carries **no label column at all**, and the
reader refuses to invent one: :class:`DesignedCdsPanel` has no notion of a positive, and
nothing built on it may report ``passed`` or ``promotable``.

**What it can measure, without any labels:**

1. **Response to synonymous change.** Every member of a group encodes the *same protein*
   and differs only in codon choice. If a model's splice scores barely move across them,
   then routing it into BT4's design loop cannot help, because synonymous positions are
   the only thing BT4 changes. This is the most decision-relevant fact available about
   these models in this regime, and it needs no ground truth.
2. **Cross-backend agreement on designed sequence.** Agreement is label-free, so it is
   measurable exactly where truth is not. Set against the same figure on natural genomic
   sequence, a large drop is direct evidence that the models are out of distribution on
   designed CDS -- again with nothing assumed.
3. **A native reference.** Each group carries the real human CDS, so a design can be
   scored as a Δ against the sequence evolution actually uses. That is BT4's own
   ``delta_splicing`` framing and introduces no external claim.

**The within-group structure is the point**, and it mirrors what the expression gate
learned the hard way (``verify_expression_gate``'s ``within_group``): pooling across
proteins credits a model for telling *different proteins* apart, which is exactly the
skill a natural-gene-trained model has and exactly the one BT4 cannot use. The question
is always asked inside a protein.

The reader **verifies** the defining property rather than trusting the file: every
member of a group must translate to the same protein. A "synonymous" panel whose members
are not synonymous is not a weaker panel, it is a different experiment.
"""

from __future__ import annotations

import csv
import hashlib
import io
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from bt4.biomodels._csv import relaxed_field_size
from bt4.domain.genetic_code import translate
from bt4.domain.sequence import validate_dna

__all__ = [
    "PANEL_COLUMNS",
    "ROLES",
    "DesignedCdsPanel",
    "DesignedMember",
    "panel_from_members",
    "read_designed_cds_panel",
]

ROLES: tuple[str, ...] = ("native", "designed")
"""What a member is.

``native`` is the real coding sequence the organism uses -- exactly one per group, and
the reference a Δ is measured against. ``designed`` is any synonymous re-encoding,
whether by a commercial optimizer, an academic tool, or BT4 itself. The distinction is
deliberately only two-valued: *who* designed a sequence is provenance for the ``note``
column, not a property the analysis may branch on, because a panel that treats BT4's own
output as a special category invites exactly the self-flattering comparison it exists to
avoid."""

PANEL_COLUMNS: tuple[str, ...] = ("group", "member_id", "role", "cds", "note")
"""The format's columns. There is **no label column, by design** -- see the module
docstring. Adding one would require splice ground truth for designed sequence, which does
not exist."""

_REQUIRED: tuple[str, ...] = ("group", "member_id", "role", "cds")


@dataclass(frozen=True, slots=True)
class DesignedMember:
    """One coding sequence in the panel.

    Attributes:
        group: The protein. All members of a group encode it identically.
        member_id: Unique within the panel.
        role: One of :data:`ROLES`.
        cds: The coding sequence, validated ACGT.
        note: Free-text provenance -- the designing tool, the run, the citation.
    """

    group: str
    member_id: str
    role: str
    cds: str
    note: str = ""

    @property
    def protein(self) -> str:
        """The translated protein, including any trailing stop."""
        return translate(self.cds)


@dataclass(frozen=True, slots=True)
class DesignedCdsPanel:
    """Designed synonymous CDS variants, grouped by protein.

    Attributes:
        members: Every member, in file order.
        source: Where it was read from.
        provenance: How the sequences were obtained, verbatim and required -- a
            designed-CDS panel is only interpretable next to who designed it and under
            what objective, and unlike a label this *is* knowable.
    """

    members: tuple[DesignedMember, ...]
    source: str
    provenance: str

    def __len__(self) -> int:
        return len(self.members)

    @property
    def groups(self) -> tuple[str, ...]:
        """The protein groups, sorted."""
        return tuple(sorted({member.group for member in self.members}))

    def group_members(self, group: str) -> tuple[DesignedMember, ...]:
        """Every member of ``group``, in file order."""
        return tuple(member for member in self.members if member.group == group)

    def native(self, group: str) -> DesignedMember:
        """The group's native reference.

        Raises:
            KeyError: If the group has no native member -- which
                :func:`panel_from_members` already refuses, so this only fires on a
                hand-built panel.
        """
        for member in self.group_members(group):
            if member.role == "native":
                return member
        raise KeyError(f"group {group!r} has no native member")

    def designed(self, group: str) -> tuple[DesignedMember, ...]:
        """Every designed (non-native) member of ``group``."""
        return tuple(m for m in self.group_members(group) if m.role != "native")

    def content_hash(self) -> str:
        """A stable, order-independent hash of the panel's content (invariant #7)."""
        digest = hashlib.sha256()
        for member in sorted(self.members, key=lambda m: (m.group, m.member_id)):
            digest.update(
                f"{member.group}\x1f{member.member_id}\x1f{member.role}\x1f"
                f"{member.cds}\x1e".encode()
            )
        return digest.hexdigest()

    def describe(self) -> dict[str, object]:
        """A summary suitable for printing or embedding in a report."""
        return {
            "source": self.source,
            "provenance": self.provenance,
            "content_hash": self.content_hash(),
            "n_members": len(self.members),
            "groups": list(self.groups),
            "group_sizes": {g: len(self.group_members(g)) for g in self.groups},
            "designed_per_group": {g: len(self.designed(g)) for g in self.groups},
            "cds_lengths": {g: len(self.native(g).cds) for g in self.groups},
            "carries_splice_labels": False,
        }


def panel_from_members(
    members: Iterable[DesignedMember], *, source: str = "", provenance: str
) -> DesignedCdsPanel:
    """Validate ``members`` into a panel.

    Enforces the properties that make it *this* experiment rather than another:

    * every member's CDS is ACGT and a whole number of codons;
    * ``member_id`` is unique;
    * ``role`` is one of :data:`ROLES`;
    * exactly **one** native per group, since a Δ needs a single reference;
    * **every member of a group translates to the same protein.** This is the defining
      property and is verified rather than trusted -- a panel whose members are not
      synonymous is not a weaker synonymous panel, it is a different experiment, and
      every conclusion drawn from it would be about the wrong thing.

    Args:
        members: The members.
        source: Where they came from, for the report.
        provenance: How the sequences were obtained. Required -- see
            :attr:`DesignedCdsPanel.provenance`.

    Returns:
        The validated panel.

    Raises:
        ValueError: On any violation above, naming the member at fault.
    """
    collected = tuple(members)
    if not collected:
        raise ValueError("a designed-CDS panel needs at least one member")
    if not provenance.strip():
        raise ValueError(
            "provenance is required: a designed-CDS panel is only interpretable next to "
            "who designed the sequences and under what objective"
        )

    seen: set[str] = set()
    for member in collected:
        if member.member_id in seen:
            raise ValueError(f"duplicate member_id {member.member_id!r}")
        seen.add(member.member_id)
        if member.role not in ROLES:
            raise ValueError(
                f"member {member.member_id!r} has role {member.role!r}; "
                f"expected one of {list(ROLES)}"
            )
        validate_dna(member.cds)
        if len(member.cds) % 3:
            raise ValueError(
                f"member {member.member_id!r} is {len(member.cds)} nt, not a whole "
                "number of codons -- it cannot be a coding sequence"
            )

    by_group: dict[str, list[DesignedMember]] = {}
    for member in collected:
        by_group.setdefault(member.group, []).append(member)

    for group, group_members in by_group.items():
        natives = [m for m in group_members if m.role == "native"]
        if len(natives) != 1:
            raise ValueError(
                f"group {group!r} has {len(natives)} native members; exactly one is "
                "required, because it is the reference every Δ in the group is measured "
                "against"
            )
        proteins = {m.protein for m in group_members}
        if len(proteins) != 1:
            offenders = sorted(
                m.member_id for m in group_members if m.protein != natives[0].protein
            )
            raise ValueError(
                f"group {group!r} is not synonymous: its members encode "
                f"{len(proteins)} distinct proteins. Offending member(s) relative to the "
                f"native: {offenders[:5]}. Synonymy is this panel's defining property, "
                "so a group that lacks it is a different experiment, not a weaker one"
            )

    return DesignedCdsPanel(members=collected, source=source, provenance=provenance)


def _member_from_mapping(row: dict[str, str | None], line_number: int) -> DesignedMember:
    """Build one member from a parsed row, naming the line on failure."""
    values: dict[str, str] = {k: (v or "").strip() for k, v in row.items()}
    missing = [column for column in _REQUIRED if not values.get(column)]
    if missing:
        raise ValueError(f"line {line_number}: missing value(s) for {missing}")
    return DesignedMember(
        group=values["group"],
        member_id=values["member_id"],
        role=values["role"].lower(),
        cds=values["cds"].upper(),
        note=values.get("note", ""),
    )


def _parse(handle: Iterable[str], source: str, provenance: str) -> DesignedCdsPanel:
    """Parse an open tab-separated panel."""
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"{source or 'panel'} is empty (no header row)")
    header = [(name or "").strip() for name in reader.fieldnames]
    missing = [column for column in _REQUIRED if column not in header]
    if missing:
        raise ValueError(
            f"{source or 'panel'} header is missing required column(s) {missing}; "
            f"expected at least {list(_REQUIRED)}"
        )
    unknown = [name for name in header if name and name not in PANEL_COLUMNS]
    if unknown:
        # A `label` column is the one worth refusing loudly: it would mean someone
        # believed this panel carries splice ground truth, and everything downstream
        # would inherit that belief.
        raise ValueError(
            f"{source or 'panel'} has unrecognised column(s) {unknown}; known columns "
            f"are {list(PANEL_COLUMNS)}. Note there is deliberately NO label column: "
            "designed coding sequence has no splice ground truth, and a panel that "
            "claims one is manufacturing it"
        )
    # A CDS can exceed csv's 131,072-char default only for a very long protein, but the
    # cap is process-global state and this reader should behave like the others.
    with relaxed_field_size():
        members = [
            _member_from_mapping({k: v for k, v in row.items() if k}, line_number)
            for line_number, row in enumerate(reader, start=2)
        ]
    return panel_from_members(members, source=source, provenance=provenance)


def read_designed_cds_panel(path: str | Path, *, provenance: str) -> DesignedCdsPanel:
    """Read a tab-separated designed-CDS panel from ``path``.

    Args:
        path: The panel file.
        provenance: How the sequences were obtained, verbatim. Required.

    Returns:
        A validated :class:`DesignedCdsPanel`.

    Raises:
        ValueError: On a malformed header, a bad row, or a failed validation.
        OSError: If ``path`` cannot be read.
    """
    text = Path(path).read_text(encoding="utf-8")
    return _parse(io.StringIO(text), str(path), provenance)


def write_designed_cds_panel(
    members: Sequence[DesignedMember], path: str | Path
) -> None:
    """Write ``members`` as a tab-separated panel, deterministically."""
    lines: Iterator[str] = iter(
        ["\t".join(PANEL_COLUMNS)]
        + [
            "\t".join(
                (m.group, m.member_id, m.role, m.cds, m.note.replace("\t", " "))
            )
            for m in members
        ]
    )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
