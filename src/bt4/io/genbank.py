"""GenBank flat-file I/O: the annotated construct a user actually opens.

BT4 already reports every residual violation in its JSON audit and in BT4 Studio.
But the file a molecular biologist opens is a **map** -- SnapGene, Benchling, ApE
-- and a defect that is only in a JSON blob is a defect nobody sees at the bench.
This module closes that gap: :func:`write_genbank` renders a designed sequence as
a GenBank record whose ``misc_feature`` annotations carry the **residual
violations the optimizer could not remove**, positioned exactly where they occur.
A dispersed repeat the refinement layer could not clear stops being a number in a
table and becomes a labelled span on the map.

It reads, too. :func:`parse_genbank` recovers the sequence, topology and features
from a GenBank file, which is how a user turns the vector they already have into
:class:`~bt4.domain.context.ConstructContext`: open the backbone, take the
sequence flanking the insertion point, design the CDS *in that context*.
:func:`context_from_genbank` does exactly that in one call.

Two properties are load-bearing.

* **Deterministic.** The record carries no timestamp. The LOCUS date field is
  conventional in GenBank but a clock reading would break invariant #7 (identical
  input must give byte-identical output), and BT4 already stamps something
  strictly more useful: the run's config hash and git commit, written into the
  COMMENT block. Readers tolerate the missing date; a reproducible file is worth
  more than a decorative one.
* **Honest about scope.** The COMMENT block states the optimality certificate and
  the *counts* of residual hard/soft violations, so a map that carries annotations
  never implies the design is clean. When a violation could not be located it is
  still listed in the comment rather than silently dropped.

Pure stdlib; depends only on :mod:`bt4.domain`.
"""

from __future__ import annotations

import re
import textwrap
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from bt4.domain.context import CIRCULAR, LINEAR, ConstructContext
from bt4.domain.result import Violation

if TYPE_CHECKING:
    from bt4.domain import Result

__all__ = [
    "GenBankFeature",
    "GenBankRecord",
    "context_from_genbank",
    "parse_genbank",
    "read_genbank",
    "write_genbank",
]

_SEQ_LINE = 60
_SEQ_GROUP = 10
_QUALIFIER_INDENT = " " * 21
_FEATURE_KEY_WIDTH = 16
_LINE_WIDTH = 79

# GenBank's LOCUS name field is fixed-width and space-delimited, so a name with a
# space (or an over-long one) would corrupt the line for every downstream parser.
_LOCUS_SAFE = re.compile(r"[^A-Za-z0-9_.-]")
_MAX_LOCUS = 16


@dataclass(frozen=True, slots=True)
class GenBankFeature:
    """One feature table entry.

    Attributes:
        key: The feature key (``"CDS"``, ``"misc_feature"``, ``"source"``, ...).
        start: Zero-based inclusive start in the record's sequence.
        end: Zero-based exclusive end (so ``[start, end)``, matching
            :class:`~bt4.domain.result.Violation`).
        qualifiers: Ordered ``(name, value)`` pairs. A value of ``None`` renders a
            bare qualifier (``/pseudo``); everything else is rendered quoted unless
            it is an integer-valued string such as ``codon_start``.
    """

    key: str
    start: int
    end: int
    qualifiers: tuple[tuple[str, str | None], ...] = ()

    @property
    def location(self) -> str:
        """The GenBank location string (1-based inclusive)."""
        return f"{self.start + 1}..{self.end}"


@dataclass(frozen=True, slots=True)
class GenBankRecord:
    """A parsed GenBank record -- enough of one for BT4's purposes.

    Attributes:
        locus: The LOCUS name.
        sequence: The upper-cased nucleotide sequence.
        topology: ``"linear"`` or ``"circular"``.
        definition: The DEFINITION text (may be empty).
        features: The feature table, in file order.
    """

    locus: str
    sequence: str
    topology: str = LINEAR
    definition: str = ""
    features: tuple[GenBankFeature, ...] = field(default=())

    def feature(self, key: str) -> GenBankFeature | None:
        """Return the first feature with ``key``, or ``None``."""
        return next((f for f in self.features if f.key == key), None)


def _safe_locus(name: str) -> str:
    """Return ``name`` reduced to a legal, fixed-width LOCUS token."""
    cleaned = _LOCUS_SAFE.sub("_", name.strip()) or "BT4"
    return cleaned[:_MAX_LOCUS]


def _render_qualifier(name: str, value: str | None) -> list[str]:
    """Render one qualifier, wrapped to the GenBank column width."""
    if value is None:
        return [f"{_QUALIFIER_INDENT}/{name}"]
    # codon_start and friends are numeric and conventionally unquoted.
    body = value if value.isdigit() else '"' + value.replace('"', "''") + '"'
    text = f"/{name}={body}"
    return textwrap.wrap(
        text,
        width=_LINE_WIDTH,
        initial_indent=_QUALIFIER_INDENT,
        subsequent_indent=_QUALIFIER_INDENT,
        break_long_words=True,
        break_on_hyphens=False,
    ) or [f"{_QUALIFIER_INDENT}{text}"]


def _render_feature(feature: GenBankFeature) -> list[str]:
    """Render one feature block (key line plus its qualifiers)."""
    key = feature.key.ljust(_FEATURE_KEY_WIDTH)
    lines = [f"     {key}{feature.location}"]
    for name, value in feature.qualifiers:
        lines.extend(_render_qualifier(name, value))
    return lines


def _render_origin(sequence: str) -> list[str]:
    """Render the ORIGIN block: 60 lower-cased bases per line in groups of ten."""
    lines = ["ORIGIN"]
    seq = sequence.lower()
    for offset in range(0, len(seq), _SEQ_LINE):
        chunk = seq[offset : offset + _SEQ_LINE]
        groups = [chunk[i : i + _SEQ_GROUP] for i in range(0, len(chunk), _SEQ_GROUP)]
        lines.append(f"{offset + 1:>9} " + " ".join(groups))
    return lines


def _violation_features(result: Result, cds_offset: int) -> list[GenBankFeature]:
    """``misc_feature`` annotations for the violations on the delivered sequence.

    Positions are shifted by ``cds_offset`` so they land correctly when the record
    is the whole construct rather than the bare CDS. This is the point of the
    module: a residual the optimizer could not remove appears **on the map**, at
    the base where it occurs, rather than only in a JSON audit.

    Overlapping or touching violations of the *same* constraint and severity are
    merged into one annotated span. That is a rendering decision, not an accounting
    one: a sliding rule reports one violation per offending window (19 overlapping
    windows for a single unremovable repeat is a real count the engine needs, and
    the COMMENT block still states it), but 19 stacked features over one region
    make a map unreadable. The merged feature names how many findings it covers, so
    nothing is hidden -- the count is reported, just not drawn nineteen times.
    """
    ordered = sorted(
        result.violations, key=lambda v: (v.constraint, v.severity.value, v.start, v.end)
    )
    features: list[GenBankFeature] = []
    for constraint, severity, start, end, members, detail in _merge_violations(ordered):
        note = f"BT4 {constraint} [{severity}]"
        if members > 1:
            note += f" x{members} overlapping findings"
        if detail:
            note += f": {detail}"
        label = f"{constraint} ({severity})"
        if members > 1:
            label += f" x{members}"
        features.append(
            GenBankFeature(
                key="misc_feature",
                start=start + cds_offset,
                end=end + cds_offset,
                qualifiers=(("note", note), ("label", label)),
            )
        )
    return sorted(features, key=lambda f: (f.start, f.end))


def _merge_violations(
    ordered: Sequence[Violation],
) -> list[tuple[str, str, int, int, int, str]]:
    """Group same-constraint, same-severity violations whose spans touch or overlap.

    Returns ``(constraint, severity, start, end, member_count, first_detail)`` per
    merged group, where ``first_detail`` is the detail of the earliest member (so
    the annotation still names a concrete finding rather than a summary).
    """
    groups: list[tuple[str, str, int, int, int, str]] = []
    for violation in ordered:
        constraint = violation.constraint
        severity = violation.severity.value.upper()
        start = violation.start
        end = violation.end
        detail = violation.detail
        if groups:
            g_con, g_sev, g_start, g_end, g_n, g_detail = groups[-1]
            if g_con == constraint and g_sev == severity and start <= g_end:
                groups[-1] = (
                    g_con,
                    g_sev,
                    g_start,
                    max(g_end, end),
                    g_n + 1,
                    g_detail,
                )
                continue
        groups.append((constraint, severity, start, end, 1, detail))
    return groups


def _comment_lines(result: Result, context: ConstructContext | None) -> list[str]:
    """The COMMENT block: certificate, key metrics, and the provenance stamp."""
    certificate = result.certificate
    metrics = result.metrics
    audit = result.audit
    body: list[str] = [
        "Designed by BT4 (back-translation with constrained optimization).",
        f"Optimality: {certificate.status.value} (solver: {certificate.solver}).",
    ]
    if certificate.relaxed_terms:
        body.append(f"Relaxed rules: {', '.join(certificate.relaxed_terms)}.")
    if certificate.detail:
        body.append(f"Certificate detail: {certificate.detail}")
    cai = audit.get("cai")
    reference_set = audit.get("codon_reference_set")
    if cai is not None:
        line = f"CAI: {float(cai):.4f}"  # type: ignore[arg-type]
        if reference_set:
            line += f" (reference set: {reference_set})"
        body.append(line + ".")
    body.append(f"GC: {metrics.gc * 100:.1f}%. Length: {metrics.length_nt} nt.")
    body.append(
        f"Residual violations: {metrics.hard_violations} hard, "
        f"{metrics.soft_violations} soft -- annotated as misc_feature below."
    )
    if context is not None and not context.is_empty:
        body.append(
            f"Construct context supplied: {len(context.upstream)} nt upstream, "
            f"{len(context.downstream)} nt downstream; the CDS begins at base "
            f"{context.cds_offset + 1}."
        )
    manifest = audit.get("manifest")
    if isinstance(manifest, dict):
        for key in ("bt4_version", "config_hash", "git_commit", "seed"):
            value = manifest.get(key)
            if value not in (None, ""):
                body.append(f"{key}: {value}")
    lines: list[str] = []
    for i, entry in enumerate(body):
        prefix = "COMMENT     " if i == 0 else " " * 12
        wrapped = textwrap.wrap(
            entry, width=_LINE_WIDTH, initial_indent=prefix, subsequent_indent=" " * 12
        )
        lines.extend(wrapped or [prefix.rstrip()])
    return lines


def write_genbank(
    result: Result,
    *,
    context: ConstructContext | None = None,
    locus: str = "BT4",
    definition: str = "",
    organism: str = "synthetic construct",
    whole_construct: bool = True,
) -> str:
    """Render ``result`` as a GenBank record annotated with its residual violations.

    Args:
        result: The optimization result to serialize.
        context: Optional construct context. When supplied (and ``whole_construct``
            is true) the record spans the assembled construct -- flanks included and
            labelled -- with the CDS feature at its real coordinates, so a defect at
            a junction is visible where it actually sits.
        locus: LOCUS name; sanitized to GenBank's fixed-width token rules.
        definition: DEFINITION text. Defaults to a short generated description.
        organism: Value of the ``source`` feature's ``/organism`` qualifier.
        whole_construct: When false, emit only the CDS even if ``context`` is given.

    Returns:
        The GenBank flat file as text, ending in ``//`` and a newline. The output
        is deterministic: no timestamp is written (see the module docstring).
    """
    use_context = (
        context is not None and not context.is_empty and whole_construct
    )
    ctx = context if use_context else None
    sequence = ctx.assemble(result.dna) if ctx is not None else result.dna
    cds_offset = ctx.cds_offset if ctx is not None else 0
    topology = ctx.topology if ctx is not None else LINEAR

    features: list[GenBankFeature] = [
        GenBankFeature(
            key="source",
            start=0,
            end=len(sequence),
            qualifiers=(
                ("organism", organism),
                ("mol_type", "other DNA"),
            ),
        )
    ]
    if ctx is not None:
        if ctx.upstream:
            features.append(
                GenBankFeature(
                    key="misc_feature",
                    start=0,
                    end=len(ctx.upstream),
                    qualifiers=(
                        ("note", "supplied 5' construct context (not designed by BT4)"),
                        ("label", "5' context"),
                    ),
                )
            )
        if ctx.downstream:
            start = cds_offset + len(result.dna)
            features.append(
                GenBankFeature(
                    key="misc_feature",
                    start=start,
                    end=start + len(ctx.downstream),
                    qualifiers=(
                        ("note", "supplied 3' construct context (not designed by BT4)"),
                        ("label", "3' context"),
                    ),
                )
            )
        for span_start, span_end in ctx.masked_spans:
            features.append(
                GenBankFeature(
                    key="misc_feature",
                    start=span_start,
                    end=span_end,
                    qualifiers=(
                        ("note", "masked from the repeat audit (repeat by construction)"),
                        ("label", "masked span"),
                    ),
                )
            )
    features.append(
        GenBankFeature(
            key="CDS",
            start=cds_offset,
            end=cds_offset + len(result.dna),
            qualifiers=(
                ("codon_start", "1"),
                ("product", definition or "BT4 codon-optimized coding sequence"),
                ("translation", result.protein),
                ("note", "designed by BT4"),
            ),
        )
    )
    features.extend(_violation_features(result, cds_offset))

    name = _safe_locus(locus)
    header_topology = CIRCULAR if topology == CIRCULAR else LINEAR
    lines = [
        f"LOCUS       {name:<16} {len(sequence)} bp    DNA     "
        f"{header_topology:<8} SYN",
        f"DEFINITION  {definition or 'BT4 codon-optimized coding sequence.'}",
        f"ACCESSION   {name}",
        f"KEYWORDS    {'.'}",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
    ]
    lines.extend(_comment_lines(result, ctx))
    lines.append("FEATURES             Location/Qualifiers")
    for feature in features:
        lines.extend(_render_feature(feature))
    lines.extend(_render_origin(sequence))
    lines.append("//")
    return "\n".join(lines) + "\n"


_LOCATION = re.compile(r"(\d+)\.\.(\d+)")
_FEATURE_LINE = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
_QUALIFIER_LINE = re.compile(r"^ {21}/(\w+)(?:=(.*))?$")


def parse_genbank(text: str) -> GenBankRecord:
    """Parse a GenBank flat file into a :class:`GenBankRecord`.

    Handles the subset BT4 needs: the LOCUS line's name and topology, DEFINITION,
    a simple feature table (``start..end`` locations, including ``complement(...)``
    which is recorded by span), and the ORIGIN sequence. Join/multi-span locations
    are recorded by their outer bounds.

    Args:
        text: The GenBank file contents.

    Returns:
        The parsed record.

    Raises:
        ValueError: If no LOCUS line or no ORIGIN sequence is present.
    """
    lines = text.splitlines()
    locus = ""
    topology = LINEAR
    definition = ""
    features: list[GenBankFeature] = []
    seq_parts: list[str] = []
    in_features = False
    in_origin = False
    pending: list[tuple[str, str | None]] = []
    current: tuple[str, int, int] | None = None

    def flush() -> None:
        nonlocal current, pending
        if current is not None:
            key, start, end = current
            features.append(
                GenBankFeature(key=key, start=start, end=end, qualifiers=tuple(pending))
            )
        current, pending = None, []

    for raw in lines:
        if raw.startswith("LOCUS"):
            parts = raw.split()
            if len(parts) >= 2:
                locus = parts[1]
            if CIRCULAR in raw.lower():
                topology = CIRCULAR
            continue
        if raw.startswith("DEFINITION"):
            definition = raw[len("DEFINITION") :].strip()
            continue
        if raw.startswith("FEATURES"):
            in_features = True
            continue
        if raw.startswith("ORIGIN"):
            flush()
            in_features, in_origin = False, True
            continue
        if raw.startswith("//"):
            break
        if in_origin:
            seq_parts.append("".join(ch for ch in raw if ch.isalpha()))
            continue
        if in_features:
            match = _FEATURE_LINE.match(raw)
            if match:
                flush()
                key, location = match.group(1), match.group(2)
                bounds = _LOCATION.findall(location)
                if bounds:
                    start = min(int(a) for a, _ in bounds) - 1
                    end = max(int(b) for _, b in bounds)
                else:  # a single-base location such as "42"
                    single = re.findall(r"\d+", location)
                    if not single:
                        continue
                    start, end = int(single[0]) - 1, int(single[0])
                current = (key, start, end)
                continue
            qualifier = _QUALIFIER_LINE.match(raw)
            if qualifier and current is not None:
                name, value = qualifier.group(1), qualifier.group(2)
                pending.append((name, value.strip('"') if value else None))
                continue
            if current is not None and pending and raw.startswith(_QUALIFIER_INDENT):
                # Continuation of the previous qualifier's wrapped value.
                name, value = pending[-1]
                extra = raw.strip().strip('"')
                pending[-1] = (name, f"{value or ''}{extra}")

    flush()
    sequence = "".join(seq_parts).upper()
    if not locus:
        raise ValueError("not a GenBank record: no LOCUS line found")
    if not sequence:
        raise ValueError("GenBank record has no ORIGIN sequence")
    return GenBankRecord(
        locus=locus,
        sequence=sequence,
        topology=topology,
        definition=definition,
        features=tuple(features),
    )


def read_genbank(path: str) -> GenBankRecord:
    """Read and parse a GenBank file from ``path``."""
    with open(path, encoding="utf-8") as handle:
        return parse_genbank(handle.read())


def context_from_genbank(
    record: GenBankRecord | str,
    *,
    insertion_point: int | None = None,
    upstream_nt: int = 0,
    downstream_nt: int = 0,
) -> ConstructContext:
    """Build a :class:`~bt4.domain.context.ConstructContext` from a vector map.

    This is the pairing that makes GenBank reading worth having: the backbone the
    user already has becomes the context the design is optimized *inside*.

    Args:
        record: A parsed record, or GenBank text to parse.
        insertion_point: Base index (0-based) where the CDS will be inserted. When
            omitted, an existing ``CDS`` feature's start is used if present, else
            the sequence is treated as pure upstream context.
        upstream_nt: How many bases before the insertion point to keep (``0`` keeps
            all of them).
        downstream_nt: How many bases after the insertion point to keep (``0`` keeps
            all of them).

    Returns:
        The context, carrying the vector's topology so a circular backbone stays
        circular.
    """
    parsed = parse_genbank(record) if isinstance(record, str) else record
    sequence = parsed.sequence
    if insertion_point is None:
        cds = parsed.feature("CDS")
        insertion_point = cds.start if cds is not None else len(sequence)
    point = max(0, min(insertion_point, len(sequence)))
    upstream = sequence[:point]
    downstream = sequence[point:]
    if upstream_nt > 0:
        upstream = upstream[-upstream_nt:]
    if downstream_nt > 0:
        downstream = downstream[:downstream_nt]
    return ConstructContext(
        upstream=upstream, downstream=downstream, topology=parsed.topology
    )
