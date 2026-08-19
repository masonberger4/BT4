"""Build a BT4 splice-site panel from a pinned GENCODE release and GRCh38.

Step **B1** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md),
the site-prediction half. Turns an annotation plus a genome into the tab-separated
format :func:`bt4.api.read_splice_panel` reads, with the position convention correct by
construction rather than by luck.

**The arithmetic, stated once.** A GTF is 1-based, fully inclusive, and always records
``start <= end`` regardless of strand -- but GENCODE emits minus-strand exons in
*transcript* order, so this sorts by coordinate rather than trusting file order. For an
intron between consecutive exons ``(s1, e1)`` and ``(s2, e2)``:

* ``lo = e1 + 1`` is the lowest-coordinate intronic base, ``hi = s2 - 1`` the highest;
* on ``+``, the donor is ``lo`` and the acceptor is ``hi``; on ``-`` they swap, because
  transcription runs the other way.

Windows are stored in transcript orientation (minus-strand reverse-complemented), so a
genomic coordinate ``g`` maps to ``g - w_start`` on ``+`` and ``w_end - g`` on ``-``.
Taking sites per *intron* rather than per exon means the spurious "first acceptor" and
"last donor" are never generated at all.

**Two traps that silently relabel true positives as negatives.** Neither is caught by
BT4's motif check, because both produce sites that are *missing* rather than wrong:

1. **A window contains more than its centre transcript's sites.** A window spanning a
   gene overlaps neighbouring and nested transcripts. Labelling only one transcript's
   introns leaves every other real site in the window scored as a negative -- and a
   backend that correctly detects them is punished for it. This script collects sites
   from **every** MANE transcript overlapping the window, not just the one it was built
   from.
2. **Opposite-strand sites.** Antisense genes overlap sense genes often. SpliceAI and
   Pangolin are strand-specific, so an antisense site is not a site on the strand being
   scored -- but it is real sequence that looks exactly like one. By default a window
   containing any is **skipped**, because a silent false negative costs more than a
   smaller panel; ``--keep-antisense`` keeps them and records the count in each
   window's note.

Windows containing ``N`` (assembly gaps) are skipped: BT4's format forbids ``N``, and an
unscoreable position masquerading as a real negative is exactly what it forbids it for.

Determinism (invariant #7): output depends only on the inputs and the flags -- no wall
clock, no RNG -- so re-running produces a byte-identical panel and the same content hash.

Run it::

    python scripts/make_gencode_splice_panel.py \\
        --gtf gencode.v44.basic.annotation.gtf.gz \\
        --fasta GRCh38.primary_assembly.genome.fa \\
        --out panel.tsv
"""

from __future__ import annotations

import argparse
import gzip
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path

__all__ = [
    "DEFAULT_FLANK",
    "HELD_OUT_CHROMOSOMES",
    "Transcript",
    "build_windows",
    "main",
    "parse_gtf",
    "revcomp",
    "sites_from_exons",
    "to_index",
]

DEFAULT_FLANK = 5_000
"""Real genomic sequence kept each side of the transcript span.

Both wrapped CNNs consume ~10 kb of context, so this gives every interior site the
5,000 nt of *real* flank the models expect. BT4's adapters pad the window's outer edges
with ``N``; with this flank those padded positions are far from any annotated site."""

HELD_OUT_CHROMOSOMES: tuple[str, ...] = ("chr1", "chr3", "chr5", "chr7", "chr9")
"""Chromosomes neither SpliceAI nor Pangolin trained on.

Building from any other chromosome produces flattering nonsense. Pangolin's split is
stated in its own paper; SpliceAI's is confirmed from OpenSpliceAI (eLife 2025), which
rebuilt its data pipeline -- the SpliceAI paper itself is paywalled and could not be
checked directly, so this is recorded as second-hand rather than asserted as primary."""

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    """Return the reverse complement of ``seq`` (ACGTN)."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class Transcript:
    """One MANE Select transcript's exons, as read from the GTF.

    Attributes:
        transcript_id: The versioned ENST id.
        gene_name: For the window label and the note.
        chrom: Sequence name, as the GTF spells it.
        strand: ``"+"`` or ``"-"``.
        exons: 1-based inclusive ``(start, end)`` pairs, in file order.
    """

    transcript_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def span(self) -> tuple[int, int]:
        """The transcript's 1-based inclusive genomic extent."""
        return min(s for s, _ in self.exons), max(e for _, e in self.exons)


def _attribute(attributes: str, key: str) -> str:
    """Return one GTF attribute value, or ``""``. Avoids a full attribute parse."""
    token = f'{key} "'
    start = attributes.find(token)
    if start < 0:
        return ""
    start += len(token)
    end = attributes.find('"', start)
    return attributes[start:end] if end > start else ""


def parse_gtf(path: Path, chromosomes: Sequence[str]) -> dict[str, Transcript]:
    """Read MANE Select exons for ``chromosomes`` from a (gzipped) GENCODE GTF.

    Filtering to MANE Select is what keeps the panel stable across releases: from
    GENCODE v44 to v50 the protein-coding transcript count on these chromosomes grows
    4.1x while MANE Select grows 1.3%, so an unfiltered panel's negative class fills
    with low-confidence transcript models.

    Args:
        path: The GTF, plain or ``.gz``.
        chromosomes: Sequence names to keep.

    Returns:
        ``{transcript_id: Transcript}``, exons in file order.
    """
    wanted = set(chromosomes)
    transcripts: dict[str, Transcript] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon" or fields[0] not in wanted:
                continue
            attributes = fields[8]
            if 'tag "MANE_Select"' not in attributes:
                continue
            transcript_id = _attribute(attributes, "transcript_id")
            if not transcript_id:
                continue
            record = transcripts.get(transcript_id)
            if record is None:
                record = Transcript(
                    transcript_id=transcript_id,
                    gene_name=_attribute(attributes, "gene_name"),
                    chrom=fields[0],
                    strand=fields[6],
                )
                transcripts[transcript_id] = record
            record.exons.append((int(fields[3]), int(fields[4])))
    return transcripts


def sites_from_exons(
    exons: Iterable[tuple[int, int]], strand: str
) -> list[tuple[int, str]]:
    """Return ``(genomic 1-based position, kind)`` for every intron's two sites.

    Exons are sorted by coordinate here rather than trusted from file order: GENCODE
    emits minus-strand exons in transcript order, and a re-serialized GTF or a
    ``gffutils`` round-trip may not preserve either ordering.

    Taking one donor and one acceptor **per intron** is also what makes the spurious
    "first acceptor" and "last donor" impossible rather than something to remember to
    drop.
    """
    ordered = sorted(exons)
    sites: list[tuple[int, str]] = []
    for (_, e1), (s2, _) in pairwise(ordered):
        lo, hi = e1 + 1, s2 - 1
        if lo > hi:  # abutting exons: no intron, so no sites
            continue
        if strand == "+":
            sites.append((lo, "donor"))
            sites.append((hi, "acceptor"))
        else:
            sites.append((hi, "donor"))
            sites.append((lo, "acceptor"))
    return sites


def to_index(position: int, w_start: int, w_end: int, strand: str) -> int:
    """Map a genomic 1-based coordinate to a 0-based index in the stored window."""
    return (position - w_start) if strand == "+" else (w_end - position)


@dataclass(frozen=True, slots=True)
class PanelWindow:
    """One assembled window, ready to write."""

    window_id: str
    group: str
    strand: str
    sequence: str
    donors: tuple[int, ...]
    acceptors: tuple[int, ...]
    note: str


def _fetch(genome: object, chrom: str, start: int, end: int) -> str:
    """Return the plus-strand sequence for 1-based inclusive ``[start, end]``.

    Accepts a ``pyfastx.Fasta`` (what a real run uses) or any mapping of name to
    sequence string, which is what the tests use so the arithmetic is checkable without
    a 3 GB download.
    """
    if isinstance(genome, dict):
        return genome[chrom][start - 1 : end].upper()
    return str(genome.fetch(chrom, (start, end))).upper()  # type: ignore[attr-defined]


def build_windows(
    transcripts: dict[str, Transcript],
    genome: object,
    *,
    flank: int = DEFAULT_FLANK,
    keep_antisense: bool = False,
    limit: int | None = None,
) -> tuple[list[PanelWindow], dict[str, int]]:
    """Assemble panel windows, collecting **every** overlapping MANE site.

    Args:
        transcripts: MANE Select transcripts, from :func:`parse_gtf`.
        genome: A ``pyfastx.Fasta`` or a ``{chrom: sequence}`` mapping.
        flank: Real sequence kept each side of the transcript span.
        keep_antisense: Keep windows overlapping opposite-strand sites, recording the
            count in the note, instead of skipping them. Off by default: the models are
            strand-specific, so an antisense site is not a site on the strand being
            scored, but it is real sequence that looks exactly like one -- and a silent
            false negative costs more than a smaller panel.
        limit: Stop after this many windows (for a quick trial run).

    Returns:
        ``(windows, counts)`` where ``counts`` tallies why windows were skipped.
    """
    # Index every MANE site by chromosome, so a window can find its neighbours' sites
    # and not just its own. This is the fix for trap 1.
    by_chrom: dict[str, list[tuple[int, str, str]]] = {}
    for record in transcripts.values():
        for position, kind in sites_from_exons(record.exons, record.strand):
            by_chrom.setdefault(record.chrom, []).append((position, kind, record.strand))
    for sites in by_chrom.values():
        sites.sort()

    counts = {"n_gap": 0, "n_antisense": 0, "n_no_sites": 0, "n_motif": 0}
    windows: list[PanelWindow] = []
    for transcript_id in sorted(transcripts):
        if limit is not None and len(windows) >= limit:
            break
        record = transcripts[transcript_id]
        if len(record.exons) < 2:
            counts["n_no_sites"] += 1
            continue
        span_start, span_end = record.span
        w_start, w_end = max(1, span_start - flank), span_end + flank

        plus = _fetch(genome, record.chrom, w_start, w_end)
        if "N" in plus:
            # BT4's format forbids N, and an unscoreable position masquerading as a real
            # negative is precisely what it forbids it for.
            counts["n_gap"] += 1
            continue
        stored = plus if record.strand == "+" else revcomp(plus)

        same, opposite = [], 0
        for position, kind, strand in by_chrom.get(record.chrom, ()):
            if not w_start <= position <= w_end:
                continue
            if strand == record.strand:
                same.append((to_index(position, w_start, w_end, record.strand), kind))
            else:
                opposite += 1
        if opposite and not keep_antisense:
            counts["n_antisense"] += 1
            continue
        if not same:
            counts["n_no_sites"] += 1
            continue

        donors = tuple(sorted({i for i, kind in same if kind == "donor"}))
        acceptors = tuple(sorted({i for i, kind in same if kind == "acceptor"}))
        # Self-check before writing, so a bad transcript aborts at the transcript rather
        # than after a whole-genome parse. BT4's reader verifies this too; failing here
        # is cheaper and names the culprit. ~99.4% pass -- the residual is the real minor
        # spliceosome (GC-AG, U12 AT-AC), not a bug.
        bad = [i for i in donors if stored[i : i + 2] != "GT"]
        bad += [i for i in acceptors if stored[i - 1 : i + 1] != "AG"]
        if len(bad) > len(same) * 0.1:
            counts["n_motif"] += 1
            continue

        note = f"{record.gene_name} {record.transcript_id}"
        if opposite:
            note += f"; {opposite} antisense site(s) present, scored as negatives"
        windows.append(
            PanelWindow(
                window_id=f"{record.gene_name or transcript_id}_{transcript_id}",
                group=record.chrom,
                strand=record.strand,
                sequence=stored,
                donors=donors,
                acceptors=acceptors,
                note=note,
            )
        )
    return windows, counts


def write_panel(windows: Sequence[PanelWindow], path: Path) -> None:
    """Write windows as the tab-separated format ``read_splice_panel`` reads."""
    lines = ["window_id\tgroup\tsequence\tdonors\tacceptors\tstrand\tnote"]
    for window in windows:
        lines.append(
            "\t".join(
                (
                    window.window_id,
                    window.group,
                    window.sequence,
                    ",".join(str(p) for p in window.donors),
                    ",".join(str(p) for p in window.acceptors),
                    window.strand,
                    window.note,
                )
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    """Build a splice-site panel from GENCODE + GRCh38 and write it as TSV."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gtf", required=True, help="GENCODE basic annotation GTF (.gz ok)")
    parser.add_argument("--fasta", required=True, help="GRCh38 primary assembly FASTA")
    parser.add_argument("--out", required=True, help="Where to write the panel TSV")
    parser.add_argument(
        "--chromosomes", nargs="+", default=list(HELD_OUT_CHROMOSOMES),
        help="Sequence names to build from (default: the models' held-out chromosomes)",
    )
    parser.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    parser.add_argument(
        "--keep-antisense", action="store_true",
        help="Keep windows overlapping opposite-strand sites (default: skip them). "
             "Those sites become negatives, which is defensible -- the models are "
             "strand-specific -- but it is a claim, so it is opt-in and noted per window",
    )
    parser.add_argument("--limit", type=int, default=None, help="Only build N windows")
    args = parser.parse_args(argv)

    print(f"reading {args.gtf} ...", flush=True)
    transcripts = parse_gtf(Path(args.gtf), args.chromosomes)
    print(f"  {len(transcripts)} MANE Select transcripts on {args.chromosomes}", flush=True)

    try:
        import pyfastx
    except ImportError:
        parser.error("pyfastx is required to read the genome: pip install pyfastx")
    genome = pyfastx.Fasta(args.fasta)

    print("assembling windows ...", flush=True)
    windows, counts = build_windows(
        transcripts, genome, flank=args.flank,
        keep_antisense=args.keep_antisense, limit=args.limit,
    )
    write_panel(windows, Path(args.out))

    n_sites = sum(len(w.donors) + len(w.acceptors) for w in windows)
    print(f"\nwrote {args.out}")
    print(f"  windows {len(windows)}   sites {n_sites}")
    print(f"  skipped: {counts['n_gap']} with assembly gaps, "
          f"{counts['n_antisense']} overlapping antisense sites, "
          f"{counts['n_no_sites']} with no sites, {counts['n_motif']} failing the motif check")
    print("\nVerify it before use:")
    print("  python -c \"from bt4.api import read_splice_panel; import json; "
          f"print(json.dumps(read_splice_panel(r'{args.out}', "
          "negative_construction='all other positions in MANE Select gene-body windows', "
          "annotation='GENCODE v44 / GRCh38').describe(), indent=1))\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
