"""A measured CDS-variant panel: the file format the expression gate is fed.

:mod:`bt4.biomodels.expression.gate` is model- *and* format-agnostic -- it consumes
in-memory ``(predicted, measured, group)`` triples. That is right for the gate and
wrong for a maintainer, who has to turn a published supplementary table into something
runnable and then be able to prove months later exactly which bytes produced a result.
This module is that bridge: a small tab-separated format, a strict reader, and a
content hash.

**Why strict.** RiboNN silently drops any transcript whose 5'UTR exceeds 1381 nt or
whose CDS+3'UTR exceeds 11937 nt (the caps are applied inside its data module, which
prints a warning and filters the frame). A silently shortened panel is the worst
possible input to an acceptance gate: the gate would report ``n_test`` honestly and
still be answering a question about a different, quietly smaller dataset. So every row
is validated here, **before** any model runs, and a bad row is a refusal naming the row
rather than a filtered one.

**Why a content hash.** A gate result is only meaningful against the exact panel it was
computed on, and thresholds have to be pre-registered *before* the run to stay honest.
:meth:`ExpressionPanel.content_hash` is order-independent and timestamp-free
(invariant #7), so it can be written into a pre-registration file, compared afterwards,
and carried in an attestation.

**The grouping column is the protein, and that is load-bearing.** Synonymous variants
of one protein share length, amino-acid composition, UTR context and assay batch --
they are a dependent cluster, not independent observations. The gate splits and
aggregates by ``group``, so the effective sample size for any cross-group claim is the
number of *proteins*.

Pure standard library: no pandas, no model, nothing lazy to import.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from bt4.domain.sequence import validate_dna

__all__ = [
    "MAX_CDS_UTR3_LEN",
    "MAX_UTR5_LEN",
    "PANEL_COLUMNS",
    "ExpressionPanel",
    "PanelRow",
    "panel_from_rows",
    "read_panel",
]

# RiboNN's own training-time caps, hard-coded in its predict entry point. Rows beyond
# them are dropped by its data module rather than refused, so we refuse here instead.
MAX_UTR5_LEN = 1_381
MAX_CDS_UTR3_LEN = 11_937

_REQUIRED = ("group", "variant_id", "cds", "measured", "utr5", "utr3")
_OPTIONAL = ("readout", "cell_type", "species")
PANEL_COLUMNS = _REQUIRED + _OPTIONAL

_STOP_CODONS = ("TAA", "TGA", "TAG")


@dataclass(frozen=True, slots=True)
class PanelRow:
    """One measured sequence in a CDS-variant panel (immutable, validated).

    Attributes:
        group: The **protein**. Every synonymous variant of one protein shares it; it is
            both the leakage-control unit for the gate's split and the centring unit for
            within-group scoring.
        variant_id: A unique label for this sequence, used only in error messages and
            reports.
        cds: The coding sequence, ACGT, length-3N, ending in a stop codon (RiboNN's
            input contract).
        measured: The measured value, **oriented larger-is-better**, in whatever units
            the source assay used. Log-transform a ratio before putting it here.
        utr5: The 5' UTR actually used in the experiment. Non-empty.
        utr3: The 3' UTR actually used in the experiment. Non-empty.
        readout: What was measured (e.g. ``"mean_ribosome_load"``). Free text, carried
            into the report so a number is never separated from the question it answers.
        cell_type: The cell line or tissue measured. Determines which RiboNN cell-type
            output a comparison should select.
        species: ``"human"`` or ``"mouse"`` -- which RiboNN weight set applies.
    """

    group: str
    variant_id: str
    cds: str
    measured: float
    utr5: str
    utr3: str
    readout: str = ""
    cell_type: str = ""
    species: str = ""

    @property
    def utr_context(self) -> tuple[str, str]:
        """The ``(utr5, utr3)`` pair, the key a predictor must be built for."""
        return (self.utr5, self.utr3)


@dataclass(frozen=True, slots=True)
class ExpressionPanel:
    """A validated set of :class:`PanelRow` plus the identity of the bytes it came from.

    Attributes:
        rows: The panel's rows, in file order.
        source: Where it was read from, for the record (``""`` when built in memory).
    """

    rows: tuple[PanelRow, ...]
    source: str = ""
    _by_context: dict[tuple[str, str], tuple[PanelRow, ...]] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not self.rows:
            raise ValueError("an expression panel needs at least one row")
        seen: set[str] = set()
        for row in self.rows:
            if row.variant_id in seen:
                raise ValueError(f"duplicate variant_id in panel: {row.variant_id!r}")
            seen.add(row.variant_id)

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[PanelRow]:
        return iter(self.rows)

    @property
    def groups(self) -> tuple[str, ...]:
        """Distinct group ids, sorted -- the panel's *effective* sample size.

        The gate's split assigns whole groups, and a cross-group claim is supported by
        this many independent observations, never by ``len(panel)``.
        """
        return tuple(sorted({row.group for row in self.rows}))

    def group_sizes(self) -> dict[str, int]:
        """Rows per group, so an unbalanced panel is visible before it is run."""
        sizes: dict[str, int] = {}
        for row in self.rows:
            sizes[row.group] = sizes.get(row.group, 0) + 1
        return dict(sorted(sizes.items()))

    def contexts(self) -> dict[tuple[str, str], tuple[PanelRow, ...]]:
        """Group rows by their ``(utr5, utr3)`` pair.

        A predictor carries its UTR context on the *model*, not per call, so a panel
        spanning several transcripts with different UTRs cannot be scored in one
        invocation. Callers build one predictor per context and score each bucket --
        this is the split that makes that correct rather than accidental.
        """
        buckets: dict[tuple[str, str], list[PanelRow]] = {}
        for row in self.rows:
            buckets.setdefault(row.utr_context, []).append(row)
        return {context: tuple(rows) for context, rows in buckets.items()}

    def samples(self) -> list[tuple[str, float, str]]:
        """Return ``(cds, measured, group)`` triples, the gate's own input shape."""
        return [(row.cds, row.measured, row.group) for row in self.rows]

    def content_hash(self) -> str:
        """Return a stable SHA-256 over the panel's *content* (not its formatting).

        Rows are canonicalized and sorted, so re-ordering a file, re-quoting it or
        changing its column order does not change the hash, while changing any value
        does. No wall-clock and no RNG (invariant #7), so this is safe to pre-register
        before a gate run and compare against afterwards.
        """
        payload = "\n".join(
            sorted(
                "\t".join(
                    (
                        row.group,
                        row.variant_id,
                        row.cds,
                        repr(row.measured),
                        row.utr5,
                        row.utr3,
                        row.readout,
                        row.cell_type,
                        row.species,
                    )
                )
                for row in self.rows
            )
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def describe(self) -> dict[str, object]:
        """Return a JSON-ready summary, including the sizing facts a gate depends on.

        The gate's own arithmetic makes these load-bearing: a 90% conformal interval
        needs at least 9 calibration rows for a finite half-width, within-group scoring
        needs groups with 2+ members, and a grouped split needs 2+ groups.
        """
        sizes = self.group_sizes()
        rankable = sum(1 for size in sizes.values() if size >= 2)
        return {
            "source": self.source,
            "content_hash": self.content_hash(),
            "n_rows": len(self.rows),
            "n_groups": len(sizes),
            "n_groups_with_2_or_more": rankable,
            "group_sizes": sizes,
            "n_utr_contexts": len(self.contexts()),
            "readouts": sorted({row.readout for row in self.rows if row.readout}),
            "cell_types": sorted({row.cell_type for row in self.rows if row.cell_type}),
            "species": sorted({row.species for row in self.rows if row.species}),
        }


def _validate_flank(value: str, name: str, variant_id: str, cap: int | None) -> str:
    """Return an upper-cased ACGT flank, refusing anything RiboNN would reject."""
    flank = "".join(value.split()).upper()
    if not flank:
        raise ValueError(
            f"panel row {variant_id!r}: {name} is empty. RiboNN reads an all-empty UTR "
            "column as NaN and cannot preprocess it, and the UTRs carry most of its "
            "signal -- an empty-UTR score would not be meaningful."
        )
    bad = sorted({ch for ch in flank if ch not in "ACGT"})
    if bad:
        raise ValueError(f"panel row {variant_id!r}: {name} has non-ACGT {bad}")
    if cap is not None and len(flank) > cap:
        raise ValueError(
            f"panel row {variant_id!r}: {name} is {len(flank)} nt, over RiboNN's "
            f"{cap} nt cap. RiboNN would silently DROP this row, leaving the gate to "
            "answer a question about a quietly smaller panel; fix or remove the row."
        )
    return flank


def _row_from_mapping(raw: dict[str, str], line_number: int) -> PanelRow:
    """Build and fully validate one :class:`PanelRow` from a parsed TSV record."""
    missing = [column for column in _REQUIRED if not (raw.get(column) or "").strip()]
    if missing:
        raise ValueError(f"panel line {line_number}: missing value(s) for {missing}")

    variant_id = raw["variant_id"].strip()
    measured_text = raw["measured"].strip()
    try:
        measured = float(measured_text)
    except ValueError as exc:
        raise ValueError(
            f"panel row {variant_id!r}: measured={measured_text!r} is not a number"
        ) from exc
    if not math.isfinite(measured):
        raise ValueError(
            f"panel row {variant_id!r}: measured={measured_text!r} is not finite"
        )

    try:
        cds = validate_dna(raw["cds"])
    except ValueError as exc:
        raise ValueError(f"panel row {variant_id!r}: {exc}") from exc
    if len(cds) % 3 != 0:
        raise ValueError(
            f"panel row {variant_id!r}: CDS length {len(cds)} is not a multiple of 3 "
            "(RiboNN asserts length-3N and a bad row crashes the whole batch)"
        )
    if cds[-3:] not in _STOP_CODONS:
        raise ValueError(
            f"panel row {variant_id!r}: CDS ends {cds[-3:]!r}, not a stop codon "
            f"({', '.join(_STOP_CODONS)}); RiboNN's input contract requires one"
        )

    utr5 = _validate_flank(raw["utr5"], "utr5", variant_id, MAX_UTR5_LEN)
    utr3 = _validate_flank(raw["utr3"], "utr3", variant_id, None)
    if len(cds) + len(utr3) > MAX_CDS_UTR3_LEN:
        raise ValueError(
            f"panel row {variant_id!r}: CDS+3'UTR is {len(cds) + len(utr3)} nt, over "
            f"RiboNN's {MAX_CDS_UTR3_LEN} nt cap. RiboNN would silently DROP this row; "
            "fix or remove it rather than letting the gate score a smaller panel."
        )

    species = (raw.get("species") or "").strip().lower()
    if species and species not in ("human", "mouse"):
        raise ValueError(
            f"panel row {variant_id!r}: species={species!r}; RiboNN ships human and "
            "mouse weights only"
        )

    return PanelRow(
        group=raw["group"].strip(),
        variant_id=variant_id,
        cds=cds,
        measured=measured,
        utr5=utr5,
        utr3=utr3,
        readout=(raw.get("readout") or "").strip(),
        cell_type=(raw.get("cell_type") or "").strip(),
        species=species,
    )


def _parse(handle: Iterable[str], source: str) -> ExpressionPanel:
    """Parse an open tab-separated panel into a validated :class:`ExpressionPanel`."""
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
        raise ValueError(
            f"{source or 'panel'} has unrecognised column(s) {unknown}; known columns "
            f"are {list(PANEL_COLUMNS)}. A typo'd column would otherwise be ignored "
            "silently."
        )
    rows = [
        _row_from_mapping({key: value for key, value in record.items() if key}, number)
        # DictReader yields data rows, so line 1 is the header.
        for number, record in enumerate(reader, start=2)
    ]
    return ExpressionPanel(rows=tuple(rows), source=source)


def read_panel(path: str | Path) -> ExpressionPanel:
    """Read and fully validate a tab-separated CDS-variant panel.

    The format is one row per measured sequence, with a header naming at least
    ``group``, ``variant_id``, ``cds``, ``measured``, ``utr5`` and ``utr3``; the optional
    ``readout``, ``cell_type`` and ``species`` columns are carried through so a reported
    number is never separated from the question it answers.

    Every row is validated up front and a violation **raises**, naming the row. In
    particular an over-length row is refused rather than dropped, because RiboNN drops
    such rows itself and a quietly-shortened panel would let the gate report an honest
    ``n_test`` for a dataset nobody chose.

    Args:
        path: The ``.tsv`` file to read.

    Returns:
        The validated :class:`ExpressionPanel`.

    Raises:
        ValueError: On a missing/unknown column, an unparseable or non-finite
            ``measured``, a non-ACGT or non-3N or stop-less CDS, an empty or
            over-length UTR, a duplicate ``variant_id``, or an unknown species.
        OSError: If the file cannot be read.
    """
    text = Path(path).read_text(encoding="utf-8")
    return _parse(io.StringIO(text), source=str(path))


def panel_from_rows(rows: Sequence[PanelRow], source: str = "") -> ExpressionPanel:
    """Build a panel from already-constructed rows (for tests and in-memory use)."""
    return ExpressionPanel(rows=tuple(rows), source=source)
