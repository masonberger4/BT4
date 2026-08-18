"""An annotated splice-site panel: the file format the splice gate is fed.

:mod:`bt4.biomodels.splice.gate` is model- *and* format-agnostic -- it consumes
in-memory ``(predicted, label, stratum, group)`` cases. That is right for the gate and
wrong for a maintainer, who has to turn a GENCODE release into something runnable and
then prove months later exactly which bytes produced a result. This module is that
bridge: a small tab-separated format, a strict reader, and a content hash. It is the
concrete target for ``docs/DESIGN_splice_cnn_calibration.md`` step B1.

**The position convention is the whole ballgame.** A splice panel has exactly one
catastrophic failure mode, and it is silent: annotate the sites one base off and every
score in the gate is misaligned, the model looks incompetent, and nothing in the
numbers says why. Worse, the conventions in circulation genuinely differ -- an
annotation-derived table names the *exonic boundary* bases, while BT4's own backends
anchor on the *intronic* dinucleotide. So this module pins one convention and then
**verifies it against the sequence itself**:

* a **donor** position is the ``G`` of the intron-opening ``GT``: the first intronic
  base, so ``sequence[i:i + 2] == "GT"``;
* an **acceptor** position is the ``G`` of the intron-closing ``AG``: the last
  intronic base, so ``sequence[i - 1:i + 1] == "AG"``.

That is not an arbitrary pick. It is the anchor
:class:`~bt4.biomodels.splice.baseline.ConsensusPwmSplicePredictor` already uses, which
makes it the one convention verifiable from this repository rather than assumed about
someone else's.

Because ~99% of human introns are canonical U2 ``GT-AG``, a correctly built panel
matches that motif almost everywhere. A panel built to the *other* convention matches
it almost nowhere -- so :func:`read_splice_panel` **refuses** below
:data:`MIN_MOTIF_CONSISTENCY` and reports, per site kind, the shift that *would* have
worked. An off-by-one becomes a message naming the fix instead of a flattering-nonsense
gate result.

**Why a content hash.** A gate result is only meaningful against the exact panel it was
computed on, and thresholds must be pre-registered *before* the run to stay honest.
:meth:`SplicePanel.content_hash` is order-independent and timestamp-free (invariant
#7), so it can be written into a pre-registration file and compared afterwards.

**The grouping column is the chromosome, and that is load-bearing.** Pangolin and
SpliceAI both trained on human chr 2, 4, 6, 8, 10-22, X and Y, leaving chr 1, 3, 5, 7
and 9 held out. Building a panel from training chromosomes produces flattering
nonsense, so the group travels with every window and :meth:`SplicePanel.describe`
reports which chromosomes are present.

Pure standard library: no model, nothing lazy to import. Depends only on
:mod:`bt4.domain`.
"""

from __future__ import annotations

import csv
import hashlib
import io
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from bt4.domain.sequence import validate_dna

__all__ = [
    "ACCEPTOR_MOTIF",
    "DONOR_MOTIF",
    "MIN_MOTIF_CONSISTENCY",
    "OFFSET_SEARCH",
    "PANEL_COLUMNS",
    "TRAINING_CHROMOSOMES",
    "MotifConsistency",
    "SplicePanel",
    "SpliceWindow",
    "canonical_motif_at",
    "panel_from_windows",
    "read_splice_panel",
]

DONOR_MOTIF = "GT"
"""The intron-opening dinucleotide. A donor position is its ``G``."""

ACCEPTOR_MOTIF = "AG"
"""The intron-closing dinucleotide. An acceptor position is its ``G``."""

MIN_MOTIF_CONSISTENCY = 0.90
"""Fraction of sites that must carry the canonical dinucleotide at the declared
position. Roughly 99% of human introns are canonical U2 ``GT-AG`` (the rest are mostly
``GC-AG``, with a small ``AT-AC`` U12 class), so a correct panel clears this easily and
a mis-anchored one cannot come close."""

OFFSET_SEARCH: tuple[int, ...] = (-3, -2, -1, 0, 1, 2, 3)
"""Shifts probed when diagnosing a failed consistency check. Wide enough to name the
common conventions -- the exonic-boundary convention sits at donor ``+1`` /
acceptor ``-1`` -- without implying that any of them is fitted."""

TRAINING_CHROMOSOMES: frozenset[str] = frozenset(
    {"chr2", "chr4", "chr6", "chr8", *(f"chr{n}" for n in range(10, 23)), "chrX", "chrY"}
)
"""Chromosomes both SpliceAI and Pangolin trained on. A panel drawn from these is not
held out, so :meth:`SplicePanel.describe` reports the overlap rather than leaving a
reader to remember the split."""

_REQUIRED = ("window_id", "group", "sequence")
_OPTIONAL = ("donors", "acceptors", "strand", "note")
PANEL_COLUMNS = _REQUIRED + _OPTIONAL

_STRANDS = ("+", "-")


def canonical_motif_at(sequence: str, position: int, kind: str) -> bool:
    """Return whether ``kind``'s canonical dinucleotide sits at ``position``.

    Implements the module's pinned convention: a donor's ``GT`` starts *at* the
    position, an acceptor's ``AG`` *ends* at it. Passing ``"splice"`` asks whether
    *either* sits there, for a backend whose track cannot distinguish the two.

    A position without room for its motif is not canonical (rather than an error), so
    the offset probe can walk off either end and a caller can sweep a whole sequence.

    Args:
        sequence: The window, ACGT.
        position: A 0-based index; out-of-range simply yields ``False``.
        kind: ``"donor"``, ``"acceptor"`` or ``"splice"``.

    Returns:
        Whether the canonical dinucleotide is present at that anchor.

    Raises:
        ValueError: If ``kind`` is not a known site kind.
    """
    if kind == "donor":
        return position >= 0 and sequence[position : position + 2] == DONOR_MOTIF
    if kind == "acceptor":
        return position >= 1 and sequence[position - 1 : position + 1] == ACCEPTOR_MOTIF
    if kind == "splice":
        return canonical_motif_at(sequence, position, "donor") or canonical_motif_at(
            sequence, position, "acceptor"
        )
    raise ValueError(
        f"unknown site kind {kind!r}; expected 'donor', 'acceptor' or 'splice'"
    )


@dataclass(frozen=True, slots=True)
class SpliceWindow:
    """One annotated genomic window: a sequence plus the sites inside it.

    A window is the file's unit; a *position* is the gate's unit. One window yields
    ``len(sequence)`` scored positions, of which only the annotated ones are positive
    -- which is exactly why the negative construction has to be recorded on the panel
    (PR-AUC's floor is the prevalence, and the prevalence is this choice).

    Attributes:
        window_id: Unique label, used in error messages and reports.
        group: The leakage-control group, a **chromosome** (e.g. ``"chr1"``). Cases
            sharing a group are never split across folds.
        sequence: The window, ACGT only. Stored in the orientation the sites are
            annotated on -- a minus-strand gene is reverse-complemented *before* it
            reaches this format, so every position indexes this string directly and no
            consumer has to re-derive an orientation.
        donors: 0-based donor positions, sorted and unique. Each is the ``G`` of the
            intron-opening ``GT`` (see the module docstring).
        acceptors: 0-based acceptor positions, sorted and unique. Each is the ``G`` of
            the intron-closing ``AG``.
        strand: ``"+"`` or ``"-"``, the strand the window was drawn from. Provenance
            only: it never changes scoring, because ``sequence`` is already oriented.
        note: Free text (e.g. ``"deep intron, no annotated site"``), carried into
            reports so a window's purpose is never separated from its numbers.
    """

    window_id: str
    group: str
    sequence: str
    donors: tuple[int, ...] = ()
    acceptors: tuple[int, ...] = ()
    strand: str = "+"
    note: str = ""

    def __post_init__(self) -> None:
        """Validate the window in isolation, refusing anything unscoreable.

        Raises:
            ValueError: On a blank id or group, a non-ACGT sequence, an unknown strand,
                a duplicated or out-of-range position, a position without room for its
                own dinucleotide, or a position annotated as both kinds.
        """
        if not self.window_id.strip():
            raise ValueError("splice window: window_id is empty")
        if not self.group.strip():
            raise ValueError(f"splice window {self.window_id!r}: group is empty")
        if self.strand not in _STRANDS:
            raise ValueError(
                f"splice window {self.window_id!r}: strand={self.strand!r}; "
                f"expected one of {list(_STRANDS)}"
            )
        n = len(self.sequence)
        for kind, positions in (("donor", self.donors), ("acceptor", self.acceptors)):
            if list(positions) != sorted(set(positions)):
                raise ValueError(
                    f"splice window {self.window_id!r}: {kind} positions must be sorted "
                    f"and unique, got {list(positions)}"
                )
            for position in positions:
                if not 0 <= position < n:
                    raise ValueError(
                        f"splice window {self.window_id!r}: {kind} position {position} "
                        f"is outside the {n} nt window"
                    )
                # A site whose dinucleotide would run off the end cannot be checked
                # against the convention, so it cannot be trusted to follow it.
                if kind == "donor" and position + 2 > n:
                    raise ValueError(
                        f"splice window {self.window_id!r}: donor position {position} "
                        f"has no room for its GT in a {n} nt window"
                    )
                if kind == "acceptor" and position < 1:
                    raise ValueError(
                        f"splice window {self.window_id!r}: acceptor position "
                        f"{position} has no room for its AG (needs the base before it)"
                    )
        both = sorted(set(self.donors) & set(self.acceptors))
        if both:
            raise ValueError(
                f"splice window {self.window_id!r}: position(s) {both} are annotated as "
                "both donor and acceptor"
            )

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def n_sites(self) -> int:
        """Total annotated sites -- the window's contribution to the positive class."""
        return len(self.donors) + len(self.acceptors)

    def sites(self) -> tuple[tuple[int, str], ...]:
        """Return ``(position, kind)`` for every annotated site, sorted by position."""
        merged = [(p, "donor") for p in self.donors] + [(p, "acceptor") for p in self.acceptors]
        return tuple(sorted(merged))

    def labels(self, kind: str) -> tuple[int, ...]:
        """Return the 0/1 label of every position for one site ``kind``.

        Args:
            kind: ``"donor"`` or ``"acceptor"``.

        Returns:
            One label per nucleotide, aligned to :attr:`sequence` -- the shape a
            per-position backend track is compared against.

        Raises:
            ValueError: If ``kind`` is not a known site kind.
        """
        if kind == "donor":
            positions = set(self.donors)
        elif kind == "acceptor":
            positions = set(self.acceptors)
        else:
            raise ValueError(f"unknown site kind {kind!r}; expected 'donor' or 'acceptor'")
        return tuple(1 if i in positions else 0 for i in range(len(self.sequence)))


@dataclass(frozen=True, slots=True)
class MotifConsistency:
    """How well a panel's declared positions match the canonical dinucleotides.

    The diagnostic that turns a silent off-by-one into a message naming the fix.

    Attributes:
        n_donor: Annotated donor sites in the panel.
        n_donor_canonical: How many carry ``GT`` at the declared position.
        n_acceptor: Annotated acceptor sites in the panel.
        n_acceptor_canonical: How many carry ``AG`` ending at the declared position.
        donor_offsets: ``(shift, fraction canonical)`` over :data:`OFFSET_SEARCH`.
        acceptor_offsets: The same for acceptors.
    """

    n_donor: int
    n_donor_canonical: int
    n_acceptor: int
    n_acceptor_canonical: int
    donor_offsets: tuple[tuple[int, float], ...]
    acceptor_offsets: tuple[tuple[int, float], ...]

    @property
    def n_sites(self) -> int:
        """Total annotated sites the check ran over."""
        return self.n_donor + self.n_acceptor

    @property
    def fraction(self) -> float:
        """Fraction of all sites carrying their canonical dinucleotide.

        Returns ``0.0`` for a panel with no sites, following the module convention of
        an honest degenerate value rather than a raise -- the panel itself refuses that
        case earlier, for the stronger reason that it has no positive class.
        """
        if not self.n_sites:
            return 0.0
        return (self.n_donor_canonical + self.n_acceptor_canonical) / self.n_sites

    @staticmethod
    def _best(offsets: Sequence[tuple[int, float]]) -> tuple[int, float]:
        """Return the best ``(shift, fraction)``, preferring ``0`` on a tie."""
        return max(offsets, key=lambda item: (item[1], item[0] == 0), default=(0, 0.0))

    @property
    def best_donor_offset(self) -> int:
        """The shift at which donor positions would best match ``GT``."""
        return self._best(self.donor_offsets)[0]

    @property
    def best_acceptor_offset(self) -> int:
        """The shift at which acceptor positions would best match ``AG``."""
        return self._best(self.acceptor_offsets)[0]

    def diagnosis(self) -> str:
        """Return a human-readable explanation of a failed consistency check.

        Names the shift that would have worked per kind, and calls out the
        exonic-boundary convention explicitly when the evidence points at it, because
        that is overwhelmingly the mistake being made.
        """
        donor_shift, donor_best = self._best(self.donor_offsets)
        acceptor_shift, acceptor_best = self._best(self.acceptor_offsets)
        parts = [
            f"{self.n_donor_canonical}/{self.n_donor} donors carry GT at their declared "
            f"position, {self.n_acceptor_canonical}/{self.n_acceptor} acceptors carry AG"
        ]
        if donor_shift or acceptor_shift:
            parts.append(
                f"shifting donors by {donor_shift:+d} would reach {donor_best:.1%} and "
                f"acceptors by {acceptor_shift:+d} would reach {acceptor_best:.1%}"
            )
        if donor_shift == 1 and acceptor_shift == -1:
            parts.append(
                "that pattern is the exonic-boundary convention: your donors look like "
                "the LAST EXONIC base and your acceptors like the FIRST EXONIC base. "
                "BT4 anchors on the intronic dinucleotide -- move each donor +1 and "
                "each acceptor -1"
            )
        elif not donor_shift and not acceptor_shift:
            parts.append(
                "no shift in the probed range helps, so the positions are probably not "
                "off by a fixed offset -- check the window's orientation (a minus-strand "
                "window must be reverse-complemented before it reaches this format) and "
                "whether the sequence and the coordinates came from the same assembly"
            )
        return "; ".join(parts)


def _consistency(windows: Sequence[SpliceWindow]) -> MotifConsistency:
    """Compute the motif-consistency diagnostic over every window in a panel."""
    counts = {"donor": 0, "acceptor": 0}
    canonical = {"donor": 0, "acceptor": 0}
    shifted: dict[str, dict[int, int]] = {
        "donor": dict.fromkeys(OFFSET_SEARCH, 0),
        "acceptor": dict.fromkeys(OFFSET_SEARCH, 0),
    }
    for window in windows:
        for position, kind in window.sites():
            counts[kind] += 1
            if canonical_motif_at(window.sequence, position, kind):
                canonical[kind] += 1
            for offset in OFFSET_SEARCH:
                if canonical_motif_at(window.sequence, position + offset, kind):
                    shifted[kind][offset] += 1

    def profile(kind: str) -> tuple[tuple[int, float], ...]:
        total = counts[kind]
        return tuple(
            (offset, (shifted[kind][offset] / total) if total else 0.0)
            for offset in OFFSET_SEARCH
        )

    return MotifConsistency(
        n_donor=counts["donor"],
        n_donor_canonical=canonical["donor"],
        n_acceptor=counts["acceptor"],
        n_acceptor_canonical=canonical["acceptor"],
        donor_offsets=profile("donor"),
        acceptor_offsets=profile("acceptor"),
    )


@dataclass(frozen=True, slots=True)
class SplicePanel:
    """A validated set of :class:`SpliceWindow` plus the provenance its numbers need.

    Attributes:
        windows: The panel's windows, in file order.
        negative_construction: How the negative class was built, verbatim -- e.g.
            ``"all other positions in the same gene bodies"``. **Required**, because
            PR-AUC's floor is the prevalence and the prevalence is exactly this choice;
            :func:`~bt4.biomodels.splice.gate.verify_splice_gate` demands it for the
            same reason.
        annotation: The gene model the sites came from, e.g.
            ``"GENCODE v44 / GRCh38"``. Recorded because annotation choice alone
            altered SpliceAI's predictions for >10% of variants in some genes
            (Smith & Kitzman 2023).
        source: Where the panel was read from (``""`` when built in memory).
    """

    windows: tuple[SpliceWindow, ...]
    negative_construction: str
    annotation: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        """Validate the panel as a whole.

        Raises:
            ValueError: If it is empty, repeats a ``window_id``, omits the negative
                construction, or contains no annotated site at all.
        """
        if not self.windows:
            raise ValueError("a splice panel needs at least one window")
        if not self.negative_construction.strip():
            raise ValueError(
                "negative_construction is required: PR-AUC's floor is the prevalence, "
                "and the prevalence is a construction choice, so a panel that does not "
                "record how its negatives were built cannot support a threshold"
            )
        seen: set[str] = set()
        for window in self.windows:
            if window.window_id in seen:
                raise ValueError(f"duplicate window_id in panel: {window.window_id!r}")
            seen.add(window.window_id)
        if not self.n_sites:
            raise ValueError(
                "the panel has no annotated splice site: with no positive class every "
                "stratum is unscoreable and the gate would refuse it"
            )

    def __len__(self) -> int:
        return len(self.windows)

    def __iter__(self) -> Iterator[SpliceWindow]:
        return iter(self.windows)

    @property
    def n_positions(self) -> int:
        """Total scored positions -- the gate's actual sample size, per site kind."""
        return sum(len(window.sequence) for window in self.windows)

    @property
    def n_sites(self) -> int:
        """Total annotated sites across every window."""
        return sum(window.n_sites for window in self.windows)

    @property
    def groups(self) -> tuple[str, ...]:
        """Distinct chromosome ids, sorted -- the panel's leakage-control units."""
        return tuple(sorted({window.group for window in self.windows}))

    @property
    def training_overlap(self) -> tuple[str, ...]:
        """Groups that both published models trained on, sorted.

        Non-empty means part of the panel is **not** held out, and any metric computed
        over it is optimistic. Reported rather than refused, because a deliberate
        train-set comparison is a legitimate thing to run -- as long as nobody can
        mistake it for a held-out result.
        """
        return tuple(sorted(g for g in self.groups if g.strip().lower() in TRAINING_CHROMOSOMES))

    def motif_consistency(self) -> MotifConsistency:
        """Check every declared position against its canonical dinucleotide."""
        return _consistency(self.windows)

    def content_hash(self) -> str:
        """Return a stable SHA-256 over the panel's *content* (not its formatting).

        Windows are canonicalized and sorted, so re-ordering a file or changing its
        column order does not move the hash, while changing any sequence, position or
        provenance field does. No wall-clock and no RNG (invariant #7), so it is safe to
        pre-register before a gate run and compare afterwards.
        """
        rows = sorted(
            "\t".join(
                (
                    window.window_id,
                    window.group,
                    window.strand,
                    window.sequence,
                    ",".join(str(p) for p in window.donors),
                    ",".join(str(p) for p in window.acceptors),
                    window.note,
                )
            )
            for window in self.windows
        )
        payload = "\n".join([self.negative_construction, self.annotation, *rows])
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def describe(self) -> dict[str, object]:
        """Return a JSON-ready summary, including the facts a gate result depends on.

        Reports prevalence per site kind (PR-AUC is uninterpretable without its floor),
        the leakage-control groups, any overlap with the models' training chromosomes,
        and the motif-consistency fraction that says whether the positions are anchored
        the way this format means them.
        """
        consistency = self.motif_consistency()
        positions = self.n_positions
        n_donor = sum(len(w.donors) for w in self.windows)
        n_acceptor = sum(len(w.acceptors) for w in self.windows)
        return {
            "source": self.source,
            "annotation": self.annotation,
            "negative_construction": self.negative_construction,
            "content_hash": self.content_hash(),
            "n_windows": len(self.windows),
            "n_positions": positions,
            "n_donor": n_donor,
            "n_acceptor": n_acceptor,
            "donor_prevalence": n_donor / positions if positions else 0.0,
            "acceptor_prevalence": n_acceptor / positions if positions else 0.0,
            "groups": list(self.groups),
            "group_sizes": {
                group: sum(1 for w in self.windows if w.group == group)
                for group in self.groups
            },
            "training_chromosome_overlap": list(self.training_overlap),
            "motif_consistency": consistency.fraction,
        }


def _parse_positions(raw: str, column: str, window_id: str) -> tuple[int, ...]:
    """Parse a comma-separated position list, refusing anything ambiguous.

    Args:
        raw: The cell's text. Empty (or whitespace) means no sites of this kind, which
            is legitimate -- a deep-intron window is a pure-negative control.
        column: The column name, for error messages.
        window_id: The window's id, for error messages.

    Returns:
        Sorted, de-duplicated 0-based positions.

    Raises:
        ValueError: On a non-integer entry, a negative position, or a duplicate. A
            duplicate is refused rather than collapsed: it means the upstream
            extraction emitted the same site twice, which usually means a transcript
            loop double-counted shared exons, and silently de-duplicating would hide
            that while quietly changing the prevalence.
    """
    text = raw.strip()
    if not text:
        return ()
    positions: list[int] = []
    for token in text.split(","):
        entry = token.strip()
        if not entry:
            continue
        try:
            value = int(entry)
        except ValueError as exc:
            raise ValueError(
                f"splice window {window_id!r}: {column} entry {entry!r} is not an integer"
            ) from exc
        if value < 0:
            raise ValueError(
                f"splice window {window_id!r}: {column} position {value} is negative "
                "(positions are 0-based indices into the window)"
            )
        positions.append(value)
    duplicates = sorted({p for p in positions if positions.count(p) > 1})
    if duplicates:
        raise ValueError(
            f"splice window {window_id!r}: {column} repeats position(s) {duplicates}"
        )
    return tuple(sorted(positions))


def _window_from_mapping(raw: dict[str, str], line_number: int) -> SpliceWindow:
    """Build and fully validate one :class:`SpliceWindow` from a parsed TSV record."""
    missing = [column for column in _REQUIRED if not (raw.get(column) or "").strip()]
    if missing:
        raise ValueError(f"panel line {line_number}: missing value(s) for {missing}")

    window_id = raw["window_id"].strip()
    try:
        sequence = validate_dna(raw["sequence"])
    except ValueError as exc:
        # validate_dna rejects N, which is deliberate here: the models' own N-padding is
        # applied at the window's outer edges by the adapter, and an N *inside* a scored
        # window would be an unscoreable position masquerading as a real negative.
        raise ValueError(f"splice window {window_id!r}: {exc}") from exc

    return SpliceWindow(
        window_id=window_id,
        group=raw["group"].strip(),
        sequence=sequence,
        donors=_parse_positions(raw.get("donors") or "", "donors", window_id),
        acceptors=_parse_positions(raw.get("acceptors") or "", "acceptors", window_id),
        strand=(raw.get("strand") or "+").strip() or "+",
        note=(raw.get("note") or "").strip(),
    )


def panel_from_windows(
    windows: Iterable[SpliceWindow],
    *,
    negative_construction: str,
    annotation: str = "",
    source: str = "",
    min_motif_consistency: float = MIN_MOTIF_CONSISTENCY,
) -> SplicePanel:
    """Build a :class:`SplicePanel` and verify its position convention.

    Args:
        windows: The validated windows.
        negative_construction: How the negative class was built, verbatim.
        annotation: The gene model the sites came from.
        source: Where the windows came from, for the record.
        min_motif_consistency: The fraction of sites that must carry their canonical
            dinucleotide. Lower it only for a deliberately non-canonical panel (a U12
            ``AT-AC`` set), and say so in ``annotation`` -- **never** to quiet a failure,
            which is the off-by-one this check exists to catch.

    Returns:
        The validated panel.

    Raises:
        ValueError: From :class:`SpliceWindow` / :class:`SplicePanel` validation, or
            when motif consistency falls below ``min_motif_consistency`` -- with a
            diagnosis naming the shift that would have worked.
    """
    if not 0.0 <= min_motif_consistency <= 1.0:
        raise ValueError(
            f"min_motif_consistency must be in [0, 1], got {min_motif_consistency}"
        )
    panel = SplicePanel(
        windows=tuple(windows),
        negative_construction=negative_construction,
        annotation=annotation,
        source=source,
    )
    consistency = panel.motif_consistency()
    if consistency.fraction < min_motif_consistency:
        raise ValueError(
            f"{source or 'panel'}: only {consistency.fraction:.1%} of annotated sites "
            f"carry their canonical dinucleotide at the declared position, under the "
            f"{min_motif_consistency:.0%} floor. BT4 anchors a donor on the G of the "
            f"intron-opening GT and an acceptor on the G of the intron-closing AG. "
            f"{consistency.diagnosis()}"
        )
    return panel


def _parse(handle: Iterable[str], source: str, **kwargs: object) -> SplicePanel:
    """Parse an open tab-separated panel into a validated :class:`SplicePanel`."""
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
        # A typo'd column would otherwise be silently ignored, and a panel missing its
        # `acceptors` because someone wrote `acceptor` would score every acceptor as a
        # negative -- a quiet, catastrophic relabelling.
        raise ValueError(
            f"{source or 'panel'} has unrecognised column(s) {unknown}; known columns "
            f"are {list(PANEL_COLUMNS)}"
        )
    windows = [
        _window_from_mapping({k: v for k, v in row.items() if k}, line_number)
        for line_number, row in enumerate(reader, start=2)
    ]
    if not windows:
        raise ValueError(f"{source or 'panel'} has a header but no rows")
    return panel_from_windows(windows, source=source, **kwargs)  # type: ignore[arg-type]


def read_splice_panel(
    path: str | Path,
    *,
    negative_construction: str,
    annotation: str = "",
    min_motif_consistency: float = MIN_MOTIF_CONSISTENCY,
) -> SplicePanel:
    """Read a tab-separated splice panel from ``path``, validating every window.

    The format has three required columns -- ``window_id``, ``group``, ``sequence`` --
    and four optional ones: ``donors``, ``acceptors`` (comma-separated 0-based
    positions), ``strand`` and ``note``. A window with neither ``donors`` nor
    ``acceptors`` is a pure-negative control, which is a useful thing to include.

    Args:
        path: The panel file.
        negative_construction: How the negative class was built, verbatim. Required,
            because a PR-AUC threshold is meaningless without a pinned denominator.
        annotation: The gene model the sites came from, e.g. ``"GENCODE v44 / GRCh38"``.
        min_motif_consistency: See :func:`panel_from_windows`.

    Returns:
        A validated :class:`SplicePanel`.

    Raises:
        ValueError: On a malformed header, a bad row, or a failed convention check.
        OSError: If ``path`` cannot be read.
    """
    text = Path(path).read_text(encoding="utf-8")
    return _parse(
        io.StringIO(text),
        str(path),
        negative_construction=negative_construction,
        annotation=annotation,
        min_motif_consistency=min_motif_consistency,
    )
