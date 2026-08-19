"""A measured splice-variant panel: the file format the variant half of the gate is fed.

:class:`~bt4.biomodels.splice.gate.SpliceVariantCase` had no way to be built from data.
This module is that bridge, and it is deliberately shaped by the one public benchmark it
has to read: ``kitzmanlab/splicebench2023`` (Smith & Kitzman, *Genome Biol* 24:294, MIT
licensed), 3,616 variants across five genes, already scored by eight tools.

**Scores live beside the labels, not inside them.** A site panel carries sequence and the
gate scores it; a variant panel carries *measurements*, and the predictions may come from
several different sources -- BT4's own adapters, or a benchmark's pre-computed columns.
So the format is labels plus **named score columns**, and
:meth:`SpliceVariantPanel.cases` selects one. That makes the most useful check possible
without any model installed: run the gate on a published benchmark's own scores and
confirm BT4 reproduces its published figures.

**The region is the whole point.** Splice predictors are measured far weaker on exonic
variants than intronic ones -- median prAUC 0.419 vs 0.773 across eight tools -- and BT4
designs coding sequence, so its entire regime is the weaker half. ``region`` is
therefore required on every row, and the gate strata on it.

**Held-out status is checkable here too, and usually fails.** A variant panel's leakage
unit is the **gene**, but whether a gene was held out is a property of its
*chromosome* -- so a row may declare one. It matters more than it looks: 53% of
splicebench2023's variants (BRCA1, FAS, WT1) sit on chromosomes both SpliceAI and
Pangolin trained on. A panel that cannot establish its held-out status says so rather
than assuming.

Pure standard library. Depends only on :mod:`bt4.domain` and this package.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from bt4.biomodels.splice.gate import SpliceVariantCase
from bt4.biomodels.splice.panel import TRAINING_CHROMOSOMES, normalize_chromosome

__all__ = [
    "REGIONS",
    "VARIANT_PANEL_COLUMNS",
    "SpliceVariantPanel",
    "VariantRow",
    "read_variant_panel",
    "variant_panel_from_rows",
]

REGIONS: tuple[str, ...] = ("exonic", "intronic")
"""The two strata. Named rather than free-form, because the exonic/intronic gap is the
finding BT4 has to be able to reproduce and a typo'd stratum would silently split it."""

_REQUIRED = ("variant_id", "group", "region", "label")
_OPTIONAL = ("chromosome", "note")
VARIANT_PANEL_COLUMNS = _REQUIRED + _OPTIONAL
"""Fixed columns. Everything else in the header is read as a **named score column**."""

_TRUTHY = {"1", "true", "yes", "t"}
_FALSY = {"0", "false", "no", "f"}


@dataclass(frozen=True, slots=True)
class VariantRow:
    """One measured variant (immutable, validated).

    Attributes:
        variant_id: Unique label for this variant.
        group: The leakage-control unit -- the **gene**, or the assay. Variants of one
            gene share its sequence context and assay batch; they are a dependent
            cluster, not independent observations.
        region: ``"exonic"`` or ``"intronic"``.
        label: ``1`` if the assay called it splice-disrupting, else ``0``.
        scores: Named predictions for this variant, e.g.
            ``{"spliceai_masked": 0.31, "pangolin_masked": 0.44}``. A score may be
            missing for a variant a tool did not cover; :meth:`SpliceVariantPanel.cases`
            refuses to score a column with gaps rather than quietly dropping rows.
        chromosome: The gene's chromosome, if known. Only used to establish whether the
            panel is held out -- which is a property of the models, not of this file.
        note: Free text carried into reports.
    """

    variant_id: str
    group: str
    region: str
    label: int
    scores: tuple[tuple[str, float], ...] = ()
    chromosome: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        """Validate the row in isolation.

        Raises:
            ValueError: On a blank id or group, an unknown region, or a non-binary label.
        """
        if not self.variant_id.strip():
            raise ValueError("variant row: variant_id is empty")
        if not self.group.strip():
            raise ValueError(f"variant row {self.variant_id!r}: group is empty")
        if self.region not in REGIONS:
            raise ValueError(
                f"variant row {self.variant_id!r}: region={self.region!r}; "
                f"expected one of {list(REGIONS)}"
            )
        if self.label not in (0, 1):
            raise ValueError(
                f"variant row {self.variant_id!r}: label={self.label!r} is not 0 or 1"
            )

    def score(self, column: str) -> float | None:
        """Return this row's score for ``column``, or ``None`` if it has none."""
        return dict(self.scores).get(column)


@dataclass(frozen=True, slots=True)
class SpliceVariantPanel:
    """A validated set of :class:`VariantRow` plus the provenance its numbers need.

    Attributes:
        rows: The panel's rows, in file order.
        negative_construction: How the negative class was built, verbatim. Required for
            the same reason the site panel requires it: average precision's floor is the
            prevalence, and the prevalence is a construction choice.
        assay: What was measured, and under what criterion. **Load-bearing when a
            benchmark pools several assays' differing definitions under one boolean** --
            a pooled figure is then a composite, not a measurement, and this is where
            that gets recorded.
        source: Where it was read from (``""`` when built in memory).
    """

    rows: tuple[VariantRow, ...]
    negative_construction: str
    assay: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        """Validate the panel as a whole.

        Raises:
            ValueError: If it is empty, repeats a ``variant_id``, omits the negative
                construction, or has fewer than two groups.
        """
        if not self.rows:
            raise ValueError("a variant panel needs at least one row")
        if not self.negative_construction.strip():
            raise ValueError(
                "negative_construction is required: average precision's floor is the "
                "prevalence, and the prevalence is a construction choice, so a panel "
                "that does not record how its negatives were built cannot support a "
                "threshold"
            )
        seen: set[str] = set()
        for row in self.rows:
            if row.variant_id in seen:
                raise ValueError(f"duplicate variant_id in panel: {row.variant_id!r}")
            seen.add(row.variant_id)
        if len(self.groups) < 2:
            raise ValueError(
                f"a variant panel needs at least two leakage-control groups (genes), "
                f"got {list(self.groups)}: one gene cannot support a general claim"
            )

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[VariantRow]:
        return iter(self.rows)

    @property
    def groups(self) -> tuple[str, ...]:
        """Distinct genes, sorted -- the panel's *effective* sample size."""
        return tuple(sorted({row.group for row in self.rows}))

    @property
    def score_columns(self) -> tuple[str, ...]:
        """Every named score column present on any row, sorted."""
        return tuple(sorted({name for row in self.rows for name, _ in row.scores}))

    @property
    def training_overlap(self) -> tuple[str, ...]:
        """Groups whose declared chromosome both models trained on, sorted.

        Usually non-empty, and that is the point: 53% of ``splicebench2023``'s variants
        (BRCA1 on chr17, FAS on chr10, WT1 on chr11) sit on training chromosomes. A run
        over them is a tool-ranking benchmark, not a held-out measurement.
        """
        return tuple(
            sorted(
                {
                    row.group
                    for row in self.rows
                    if normalize_chromosome(row.chromosome) in TRAINING_CHROMOSOMES
                }
            )
        )

    @property
    def groups_without_chromosome(self) -> tuple[str, ...]:
        """Groups that declared no usable chromosome, sorted.

        Reported separately from :attr:`training_overlap`, because "not known to be a
        training chromosome" and "known not to be one" are different statements and only
        the second supports a held-out claim.
        """
        return tuple(
            sorted(
                {
                    row.group
                    for row in self.rows
                    if normalize_chromosome(row.chromosome) is None
                }
            )
        )

    @property
    def held_out(self) -> bool:
        """Whether every group is a chromosome neither model trained on."""
        return not self.training_overlap and not self.groups_without_chromosome

    def cases(self, score_column: str) -> list[SpliceVariantCase]:
        """Return gate cases scored from one named column.

        Args:
            score_column: Which prediction to gate. Use
                :attr:`score_columns` to see what is available.

        Returns:
            One :class:`~bt4.biomodels.splice.gate.SpliceVariantCase` per row.

        Raises:
            KeyError: If no row carries ``score_column``.
            ValueError: If some rows carry it and others do not. Scoring the covered
                subset would silently answer a question about a different, smaller panel
                -- the same refusal the expression panel makes about dropped rows.
        """
        if score_column not in self.score_columns:
            raise KeyError(
                f"no score column {score_column!r}; available: {list(self.score_columns)}"
            )
        missing = [row.variant_id for row in self.rows if row.score(score_column) is None]
        if missing:
            raise ValueError(
                f"score column {score_column!r} is missing for {len(missing)} of "
                f"{len(self.rows)} rows (first: {missing[0]!r}). Scoring only the covered "
                "rows would answer a question about a smaller panel than the one "
                "reported; drop those rows deliberately, or use a complete column"
            )
        return [
            SpliceVariantCase(
                predicted=row.score(score_column) or 0.0,
                label=row.label,
                region=row.region,
                group=row.group,
            )
            for row in self.rows
        ]

    def content_hash(self) -> str:
        """Return a stable SHA-256 over the panel's *content* (not its formatting).

        Order-independent and timestamp-free (invariant #7), so it can be pre-registered
        before a gate run and compared afterwards.
        """
        rows = sorted(
            "\t".join(
                (
                    row.variant_id,
                    row.group,
                    row.region,
                    str(row.label),
                    row.chromosome,
                    row.note,
                    ";".join(f"{n}={v!r}" for n, v in sorted(row.scores)),
                )
            )
            for row in self.rows
        )
        payload = "\n".join([self.negative_construction, self.assay, *rows])
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def describe(self) -> dict[str, object]:
        """Return a JSON-ready summary, including the facts a gate result depends on."""
        by_region = {
            region: [row for row in self.rows if row.region == region] for region in REGIONS
        }
        return {
            "source": self.source,
            "assay": self.assay,
            "negative_construction": self.negative_construction,
            "content_hash": self.content_hash(),
            "n_rows": len(self.rows),
            "n_positive": sum(row.label for row in self.rows),
            "groups": list(self.groups),
            "group_sizes": {
                group: sum(1 for row in self.rows if row.group == group)
                for group in self.groups
            },
            "region_sizes": {region: len(rows) for region, rows in by_region.items()},
            "region_prevalence": {
                region: (sum(r.label for r in rows) / len(rows)) if rows else 0.0
                for region, rows in by_region.items()
            },
            "score_columns": list(self.score_columns),
            "training_chromosome_overlap": list(self.training_overlap),
            "groups_without_chromosome": list(self.groups_without_chromosome),
            "held_out": self.held_out,
        }


def _parse_label(raw: str, variant_id: str) -> int:
    """Parse a label cell, accepting ``1/0`` and ``True/False`` (what benchmarks ship)."""
    text = raw.strip().lower()
    if text in _TRUTHY:
        return 1
    if text in _FALSY:
        return 0
    raise ValueError(
        f"variant row {variant_id!r}: label={raw!r} is not a boolean "
        f"(expected one of {sorted(_TRUTHY | _FALSY)})"
    )


def _row_from_mapping(
    raw: dict[str, str], score_names: Sequence[str], line_number: int
) -> VariantRow:
    """Build and fully validate one :class:`VariantRow` from a parsed TSV record."""
    missing = [column for column in _REQUIRED if not (raw.get(column) or "").strip()]
    if missing:
        raise ValueError(f"variant panel line {line_number}: missing value(s) for {missing}")

    variant_id = raw["variant_id"].strip()
    scores: list[tuple[str, float]] = []
    for name in score_names:
        cell = (raw.get(name) or "").strip()
        if not cell:
            continue  # a tool that did not cover this variant; `cases` refuses gaps
        try:
            value = float(cell)
        except ValueError as exc:
            raise ValueError(
                f"variant row {variant_id!r}: {name}={cell!r} is not a number"
            ) from exc
        if not math.isfinite(value):
            raise ValueError(f"variant row {variant_id!r}: {name}={cell!r} is not finite")
        scores.append((name, value))

    return VariantRow(
        variant_id=variant_id,
        group=raw["group"].strip(),
        region=raw["region"].strip().lower(),
        label=_parse_label(raw["label"], variant_id),
        scores=tuple(scores),
        chromosome=(raw.get("chromosome") or "").strip(),
        note=(raw.get("note") or "").strip(),
    )


def variant_panel_from_rows(
    rows: Iterable[VariantRow],
    *,
    negative_construction: str,
    assay: str = "",
    source: str = "",
) -> SpliceVariantPanel:
    """Build a :class:`SpliceVariantPanel` from validated rows."""
    return SpliceVariantPanel(
        rows=tuple(rows),
        negative_construction=negative_construction,
        assay=assay,
        source=source,
    )


def read_variant_panel(
    path: str | Path,
    *,
    negative_construction: str,
    assay: str = "",
) -> SpliceVariantPanel:
    """Read a tab-separated variant panel from ``path``, validating every row.

    Four required columns -- ``variant_id``, ``group``, ``region``, ``label`` -- plus
    optional ``chromosome`` and ``note``. **Every other header column is read as a named
    score column**, so a benchmark's own pre-computed predictions come through as data
    rather than needing a schema change.

    Args:
        path: The panel file.
        negative_construction: How the negative class was built, verbatim.
        assay: What was measured and under what criterion. Record it when the labels
            pool several assays -- a composite is not a measurement.

    Returns:
        A validated :class:`SpliceVariantPanel`.

    Raises:
        ValueError: On a malformed header or a bad row.
        OSError: If ``path`` cannot be read.
    """
    text = Path(path).read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"{path} is empty (no header row)")
    header = [(name or "").strip() for name in reader.fieldnames]
    absent = [column for column in _REQUIRED if column not in header]
    if absent:
        raise ValueError(
            f"{path} header is missing required column(s) {absent}; "
            f"expected at least {list(_REQUIRED)}"
        )
    score_names = [name for name in header if name and name not in VARIANT_PANEL_COLUMNS]
    if not score_names:
        raise ValueError(
            f"{path} has no score column: every header column beyond "
            f"{list(VARIANT_PANEL_COLUMNS)} is read as a named prediction, and a panel "
            "with none cannot be gated"
        )
    rows = [
        _row_from_mapping({k: v for k, v in record.items() if k}, score_names, line_number)
        for line_number, record in enumerate(reader, start=2)
    ]
    if not rows:
        raise ValueError(f"{path} has a header but no rows")
    return variant_panel_from_rows(
        rows, negative_construction=negative_construction, assay=assay, source=str(path)
    )
