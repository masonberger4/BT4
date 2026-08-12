"""Build a provenanced codon-usage table from a user-provided CDS set.

This is an honesty feature: rather than trusting a bundled "representative"
table, a user can count codons in their *own* coding sequences and derive a
:class:`~bt4.biomodels.codon.tables.CodonUsageTable` from authentic data, stamped
with content-hashed provenance.

The pipeline is:

1. :func:`count_codons` - validate each CDS and count every codon it contains.
2. :func:`build_table` - turn raw counts into a validated table using Laplace
   (add-``pseudocount``) smoothing so every one of the 64 codons is strictly
   positive (required by :class:`CodonUsageTable`, and it keeps ``log(w)`` finite
   in CAI even for codons unseen in the input).
3. :func:`write_table` - persist the *raw* counts as a ``<organism>.tsv`` plus a
   sibling provenance JSON whose SHA-256 hashes the exact TSV bytes on disk.

This module depends only on :mod:`bt4.domain`, the sibling ``tables`` module, and
the standard library.
"""

from __future__ import annotations

import datetime
import json
import os
from collections.abc import Iterable, Mapping
from pathlib import Path

from bt4.biomodels.codon.tables import CodonUsageTable, sha256_hex
from bt4.domain.genetic_code import CODON_TABLE
from bt4.domain.sequence import validate_dna

__all__ = ["build_table", "count_codons", "write_table"]

_BUILD_NOTE = (
    "Built by counting codons in a user-provided CDS set; "
    "frequencies are raw occurrence counts."
)
_BUILD_METHOD = "codon occurrence counts from a user-provided CDS set"


def count_codons(cds_sequences: Iterable[str]) -> dict[str, int]:
    """Count codons across a set of coding sequences.

    Each sequence is validated as DNA and required to be in-frame (length a
    multiple of three). Nothing is skipped silently: malformed input raises.

    Args:
        cds_sequences: An iterable of coding DNA strings (case-insensitive).

    Returns:
        A mapping from codon to its total occurrence count. Codons that never
        occur are absent from the mapping.

    Raises:
        ValueError: If a sequence contains non-ACGT characters, is empty, or has
            a length that is not a multiple of three (the offending length is
            reported).
    """
    counts: dict[str, int] = {}
    for cds in cds_sequences:
        seq = validate_dna(cds)
        if len(seq) % 3 != 0:
            raise ValueError(f"CDS length {len(seq)} is not a multiple of three")
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            counts[codon] = counts.get(codon, 0) + 1
    return counts


def build_table(
    cds_sequences: Iterable[str],
    *,
    organism: str,
    pseudocount: float = 1.0,
) -> tuple[CodonUsageTable, dict[str, int]]:
    """Build a codon-usage table from a CDS set with Laplace smoothing.

    Codons are counted (via :func:`count_codons`), then every one of the 64
    codons is assigned a smoothed frequency of ``count + pseudocount``. The
    add-``pseudocount`` (Laplace) smoothing keeps each frequency strictly
    positive - which :class:`CodonUsageTable` requires and which avoids
    ``log(0)`` when the resulting table scores a codon never seen in the input.
    The returned counts are the *raw*, unsmoothed occurrence counts.

    Args:
        cds_sequences: An iterable of coding DNA strings.
        organism: Organism label for the resulting table.
        pseudocount: The smoothing constant added to every codon count. Must be
            positive for the smoothing guarantee to hold (the default, ``1.0``,
            is standard Laplace smoothing).

    Returns:
        A ``(table, counts)`` pair: the validated :class:`CodonUsageTable` built
        from smoothed frequencies, and the raw (unsmoothed) codon counts.

    Raises:
        ValueError: If any CDS is malformed (see :func:`count_codons`), or if the
            resulting frequencies fail table validation (e.g. a non-positive
            ``pseudocount`` leaves an unseen codon at zero).
    """
    counts = count_codons(cds_sequences)
    frequency = {codon: counts.get(codon, 0) + pseudocount for codon in CODON_TABLE}
    table = CodonUsageTable(organism=organism, frequency=frequency)
    return table, counts


def write_table(
    counts: Mapping[str, int],
    *,
    organism: str,
    path: str | os.PathLike[str],
    source: str,
    retrieved: str | None = None,
    cds_count: int | None = None,
    pseudocount: float = 0.0,
    build: str | None = None,
    note: str | None = None,
    extra: Mapping[str, object] | None = None,
) -> str:
    """Write codon counts to a TSV plus a sibling provenance JSON.

    Two files are written into the directory ``path``:

    * ``<organism>.tsv`` - the three-column table (``amino_acid``, ``codon``,
      ``frequency``). With ``pseudocount == 0`` (the default) one row is written
      per codon present in ``counts``, with the raw occurrence count. With
      ``pseudocount > 0`` all 64 codons are written with ``count + pseudocount``
      (Laplace smoothing), so the file always covers every amino acid and can be
      loaded straight back with :func:`~bt4.biomodels.codon.tables.load_table_from_file`
      even from a small CDS set.
    * ``<organism>.provenance.json`` - a content-hashed provenance sidecar whose
      ``sha256`` is the digest of the exact TSV bytes written.

    Args:
        counts: Raw codon counts (typically from :func:`count_codons` or the
            second element of :func:`build_table`).
        organism: Organism label; also the TSV/provenance file stem.
        path: Directory in which to write the two files.
        source: Human-readable origin of the CDS set behind the counts.
        retrieved: ISO date the counts were produced. Defaults to today's date.
        cds_count: Number of coding sequences behind the counts, if known.
        pseudocount: When positive, write all 64 codons Laplace-smoothed by this
            amount; when ``0`` write only observed codons with their raw counts.
        build: Overrides the default ``build`` string, for a caller that applied
            its own documented selection/filtering rules (e.g. one representative
            transcript per gene). Describe what was actually done.
        note: Overrides the default honesty caveat. Use it to state precisely
            what the numbers are and are not.
        extra: Additional provenance keys merged into the sidecar (e.g. the
            source URL, assembly, database release, and the SHA-256 of the
            downloaded source file). They make the table **re-derivable by a
            third party**, which is the point of the stamp (CLAUDE.md §8).
            Reserved keys (``source``/``build``/``cds_count``/``retrieved``/
            ``sha256``/``note``) cannot be overwritten -- pass those through
            their own parameters -- so the loaded provenance can never disagree
            with the sidecar.

    Returns:
        The filesystem path to the written TSV, as a string.

    Raises:
        ValueError: If ``extra`` tries to set a reserved provenance key.
    """
    directory = Path(path)
    tsv_path = directory / f"{organism}.tsv"
    provenance_path = directory / f"{organism}.provenance.json"

    # Validate before writing anything: a rejected call must not leave a TSV on
    # disk with no (or a stale) sidecar beside it.
    reserved = {"source", "build", "cds_count", "retrieved", "sha256", "note"}
    clashes = sorted(set(extra or {}) & reserved)
    if clashes:
        raise ValueError(
            f"extra provenance keys {clashes} are reserved; pass them via "
            "their own parameters so the sidecar cannot disagree with itself"
        )

    lines = ["amino_acid\tcodon\tfrequency"]
    if pseudocount > 0:
        for codon in sorted(CODON_TABLE):
            freq = counts.get(codon, 0) + pseudocount
            lines.append(f"{CODON_TABLE[codon]}\t{codon}\t{freq:g}")
    else:
        for codon in sorted(counts):
            lines.append(f"{CODON_TABLE[codon]}\t{codon}\t{counts[codon]}")
    tsv_bytes = ("\n".join(lines) + "\n").encode("utf-8")
    tsv_path.write_bytes(tsv_bytes)

    resolved_retrieved = (
        retrieved if retrieved is not None else datetime.date.today().isoformat()
    )
    provenance: dict[str, object] = {
        "source": source,
        "build": build if build is not None else _BUILD_METHOD,
        "cds_count": cds_count,
        "retrieved": resolved_retrieved,
        "sha256": sha256_hex(tsv_bytes),
        "note": note if note is not None else _BUILD_NOTE,
    }
    if extra:
        provenance.update(extra)
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return str(tsv_path)
