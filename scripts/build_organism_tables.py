#!/usr/bin/env python3
"""Rebuild BT4's bundled codon-usage tables from release-pinned Ensembl CDS sets.

This is the reproducible provenance trail behind the organism tables BT4 ships
(CLAUDE.md §8): every number in ``src/bt4/biomodels/codon/data/<organism>.tsv``
is a **codon count from a real, publicly downloadable coding-sequence set**, and
this script is how anyone re-derives it. Nothing here is hand-authored, and a
table is never invented for an organism whose CDS set could not be fetched --
the run fails loudly instead (the "never fabricate a table" rule).

Usage::

    python scripts/build_organism_tables.py                  # all organisms
    python scripts/build_organism_tables.py mus_musculus     # just one
    python scripts/build_organism_tables.py --cache-dir /tmp/cds --verify

``--verify`` rebuilds into a temporary directory and diffs against the committed
TSVs **and their provenance sidecars** (every field but ``retrieved``, a
wall-clock stamp) instead of overwriting them, so CI or a reviewer can confirm
the shipped tables really are what this script produces from the pinned sources
-- and that each sidecar's re-derivation trail is the true one.

Why these choices (all of them recorded in each table's provenance sidecar):

* **Release-pinned URLs with a pinned digest.** A ``current_fasta`` link moves;
  a release-pinned one does not. Ensembl also reuses the *same filename* across
  releases, so the archive's expected SHA-256 is pinned in the spec and checked
  on every run -- cache hits included -- and a mismatch aborts rather than
  silently counting a different release's genes.
* **One representative CDS per gene.** Ensembl ships every annotated transcript,
  and gene families differ wildly in isoform count -- counting all of them would
  weight codon usage by how finely a gene happens to be annotated rather than by
  the organism's actual coding content. The longest valid CDS per gene is taken
  (ties broken by transcript ID, so the choice is deterministic).
* **Strict validity filtering.** A counted CDS must be ACGT-only, a multiple of
  three, start with ATG, end in a stop codon, and contain no internal stop.
  Partial/ambiguous transcripts are dropped rather than silently miscounted, and
  the drop counts are reported and stamped.
* **The terminal stop codon is counted.** BT4 chooses the stop codon it appends,
  so an organism's real stop-codon usage is decision-relevant. Every other codon
  in the table is a sense codon, exactly as CAI expects.

This script is a maintainer tool: it reaches the network and writes into the
package data directory. It is not imported by the library, and BT4 never fetches
anything at runtime.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from bt4.biomodels.codon.build import count_codons, write_table  # noqa: E402
from bt4.domain.genetic_code import CODON_TABLE  # noqa: E402

DATA_DIR = REPO_ROOT / "src" / "bt4" / "biomodels" / "codon" / "data"

_STOPS = frozenset({"TAA", "TAG", "TGA"})
_ACGT = frozenset("ACGT")

# Ensembl vertebrates/metazoa and Ensembl Plants releases current at build time.
# Pinned rather than "current" so this script keeps producing the shipped bytes.
ENSEMBL_RELEASE = "116"
ENSEMBL_PLANTS_RELEASE = "63"


@dataclass(frozen=True)
class OrganismSpec:
    """One organism's pinned CDS source.

    Attributes:
        key: BT4's canonical organism key (also the TSV/provenance file stem).
        common_name: Human-readable name for the provenance ``source`` string.
        url: Release-pinned URL of the gzipped CDS FASTA.
        assembly: Genome assembly the CDS set is annotated on.
        database: Which Ensembl division served it.
        release: That division's release number.
        source_sha256: Expected SHA-256 of the downloaded archive. Ensembl reuses
            the SAME filename across releases, so a cache keyed on the filename
            alone could silently serve a different release's bytes. Verified on
            every run, cache hit included; a mismatch aborts.
        extra_citation: An additional citation this specific assembly requires
            (e.g. TAIR10 for Arabidopsis), appended to the license note. Empty
            for species whose Ensembl citation alone suffices -- a sidecar should
            not carry a citation irrelevant to its own data.
    """

    key: str
    common_name: str
    url: str
    assembly: str
    database: str
    release: str
    source_sha256: str
    extra_citation: str = ""


def _ensembl(
    key: str, common_name: str, filename: str, assembly: str, source_sha256: str
) -> OrganismSpec:
    """Build a spec for a main-Ensembl (vertebrates/metazoa) species."""
    return OrganismSpec(
        key=key,
        common_name=common_name,
        url=(
            f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/"
            f"{key}/cds/{filename}"
        ),
        assembly=assembly,
        database="Ensembl",
        release=ENSEMBL_RELEASE,
        source_sha256=source_sha256,
    )


SPECS: tuple[OrganismSpec, ...] = (
    _ensembl(
        "mus_musculus", "Mus musculus (house mouse)",
        "Mus_musculus.GRCm39.cds.all.fa.gz", "GRCm39",
        "8c0a53a7fe973adaf6d884884ec86a329deebfb519d0d4f736bd1c0c67f09964",
    ),
    _ensembl(
        "rattus_norvegicus", "Rattus norvegicus (Norway rat)",
        "Rattus_norvegicus.GRCr8.cds.all.fa.gz", "GRCr8",
        "25a90be0825e525afb60eeb5568a6f5c6464216e4c262260a2547432d40acd0d",
    ),
    _ensembl(
        "danio_rerio", "Danio rerio (zebrafish)",
        "Danio_rerio.GRCz11.cds.all.fa.gz", "GRCz11",
        "b5ff5b220650e0ad426df3082ed0c2587c3f283c41d8f92490297b8ef08c5fad",
    ),
    _ensembl(
        "drosophila_melanogaster", "Drosophila melanogaster (fruit fly)",
        "Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz", "BDGP6.54",
        "112c99442c2ad59c6f012ea803cd6c51dbb96eba491e8d33d810873381fb6112",
    ),
    _ensembl(
        "caenorhabditis_elegans", "Caenorhabditis elegans (nematode)",
        "Caenorhabditis_elegans.WBcel235.cds.all.fa.gz", "WBcel235",
        "13a3ebf9bd0bfa3097a3565c979ba88a1c9fccd05cb343a7c6ef21841289e05e",
    ),
    OrganismSpec(
        key="arabidopsis_thaliana",
        common_name="Arabidopsis thaliana (thale cress)",
        url=(
            f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/"
            f"release-{ENSEMBL_PLANTS_RELEASE}/fasta/arabidopsis_thaliana/cds/"
            "Arabidopsis_thaliana.TAIR10.cds.all.fa.gz"
        ),
        assembly="TAIR10",
        database="Ensembl Plants",
        release=ENSEMBL_PLANTS_RELEASE,
        source_sha256="3580fd49dc079892359c1e5b7ffe76b7c2aa13718daa4daa6f397882ad41795d",
        extra_citation=(
            " The TAIR10 assembly/annotation is from TAIR (Lamesch et al., "
            "Nucleic Acids Res 2012, doi:10.1093/nar/gkr1090)."
        ),
    ),
)


@dataclass
class FilterStats:
    """How many CDS records survived each documented filter."""

    records: int = 0
    kept: int = 0
    dropped_non_acgt: int = 0
    dropped_frame: int = 0
    dropped_no_start: int = 0
    dropped_no_stop: int = 0
    dropped_internal_stop: int = 0
    dropped_isoform: int = 0

    def as_dict(self) -> dict[str, int]:
        """The stats as a provenance-ready mapping."""
        return {
            "records_in_source": self.records,
            "cds_counted": self.kept,
            "dropped_non_acgt": self.dropped_non_acgt,
            "dropped_not_multiple_of_three": self.dropped_frame,
            "dropped_no_atg_start": self.dropped_no_start,
            "dropped_no_terminal_stop": self.dropped_no_stop,
            "dropped_internal_stop": self.dropped_internal_stop,
            "dropped_non_representative_isoform": self.dropped_isoform,
        }


def sha256_file(path: Path) -> str:
    """Return the SHA-256 hex digest of ``path`` (streamed, constant memory)."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    """Fetch ``url`` to ``dest`` (skipping the download if it is already there)."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached  {dest.name}")
        return
    print(f"  fetching {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, partial.open("wb") as out:
        shutil.copyfileobj(response, out)
    partial.replace(dest)


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` from a gzipped FASTA, streaming."""
    header: str | None = None
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def parse_ids(header: str) -> tuple[str, str]:
    """Return ``(transcript_id, gene_id)`` from an Ensembl CDS FASTA header.

    Ensembl headers look like ``ENST0001.1 cds chromosome:...  gene:ENSG0002.3
    gene_biotype:protein_coding ...``. When no ``gene:`` field is present the
    transcript id doubles as the gene id, so the record still participates (as
    its own single-isoform gene) rather than being dropped.
    """
    transcript = header.split(None, 1)[0]
    gene = transcript
    for field in header.split():
        if field.startswith("gene:"):
            gene = field[len("gene:") :]
            break
    return transcript, gene


def valid_cds(seq: str, stats: FilterStats) -> bool:
    """Whether ``seq`` is a complete, unambiguous, in-frame coding sequence.

    Every rejection is counted (never silently skipped), and the counts are
    stamped into the table's provenance so the filtering is auditable.
    """
    if not seq or not _ACGT.issuperset(seq):
        stats.dropped_non_acgt += 1
        return False
    if len(seq) % 3 != 0:
        stats.dropped_frame += 1
        return False
    if not seq.startswith("ATG"):
        stats.dropped_no_start += 1
        return False
    if seq[-3:] not in _STOPS:
        stats.dropped_no_stop += 1
        return False
    if any(seq[i : i + 3] in _STOPS for i in range(0, len(seq) - 3, 3)):
        stats.dropped_internal_stop += 1
        return False
    return True


def representative_cds(path: Path) -> tuple[list[str], FilterStats]:
    """Select one valid CDS per gene (the longest) from a gzipped CDS FASTA.

    Returns:
        ``(sequences, stats)`` -- the chosen coding sequences and the filter
        tally behind them.
    """
    stats = FilterStats()
    best: dict[str, tuple[int, str, str]] = {}  # gene -> (length, transcript, seq)
    for header, raw in iter_fasta(path):
        stats.records += 1
        seq = raw.upper()
        if not valid_cds(seq, stats):
            continue
        transcript, gene = parse_ids(header)
        # Longest CDS wins; an exact tie keeps the lexicographically smaller
        # transcript id, so the pick is stable across runs and machines
        # regardless of the order records happen to appear in the file
        # (invariant #7 reaches the data, too).
        candidate = (len(seq), transcript, seq)
        current = best.get(gene)
        if current is None:
            best[gene] = candidate
            continue
        stats.dropped_isoform += 1
        longer = candidate[0] > current[0]
        tie_break = candidate[0] == current[0] and candidate[1] < current[1]
        if longer or tie_break:
            best[gene] = candidate
    stats.kept = len(best)
    return [entry[2] for entry in best.values()], stats


def build_one(spec: OrganismSpec, cache_dir: Path, out_dir: Path) -> Path:
    """Download, filter, count, and write one organism's table + provenance."""
    print(f"{spec.key}:")
    archive = cache_dir / Path(spec.url).name
    download(spec.url, archive)
    source_sha = sha256_file(archive)
    if source_sha != spec.source_sha256:
        # Ensembl reuses the same filename across releases, so a stale cache
        # entry -- or a re-cut upstream file -- would otherwise be counted
        # silently and produce a table that does not match its own provenance.
        raise SystemExit(
            f"{spec.key}: source archive sha256 {source_sha} != pinned "
            f"{spec.source_sha256} ({archive}). Delete the cached file to "
            "re-download, or update the pin deliberately if upstream re-cut it."
        )

    sequences, stats = representative_cds(archive)
    if not sequences:
        raise SystemExit(f"{spec.key}: no valid CDS survived filtering -- refusing to write")
    counts = count_codons(sequences)

    missing = sorted(set(CODON_TABLE) - set(counts))
    if missing:
        # A genome-scale CDS set observes all 64 codons. If any is missing the
        # input is not what we think it is -- refuse rather than paper over it
        # with smoothing, which would put an invented number in a shipped table.
        raise SystemExit(f"{spec.key}: codons unobserved in the source: {missing}")

    total = sum(counts.values())
    print(
        f"  {stats.records} records -> {stats.kept} genes, "
        f"{total} codons counted (sha256 {source_sha[:12]}...)"
    )

    written = write_table(
        counts,
        organism=spec.key,
        path=out_dir,
        source=(
            f"{spec.database} release {spec.release} -- {spec.common_name} "
            f"CDS set ({spec.assembly})"
        ),
        cds_count=stats.kept,
        build=(
            "genome-wide codon occurrence counts recounted from the pinned "
            f"{spec.database} CDS FASTA by scripts/build_organism_tables.py: one "
            "representative CDS per gene (the longest; ties broken by transcript "
            "id for determinism), keeping only ACGT-only sequences that are a "
            "multiple of three, start with ATG, end in a stop codon, and contain "
            "no internal stop. The terminal stop codon is counted; every other "
            "codon is a sense codon."
        ),
        note=(
            "Real genome-wide codon counts from a public, release-pinned CDS set "
            "-- NOT a hand-curated or published summary table. Frequencies are raw "
            "occurrence counts on a positive scale; only ratios within an amino "
            "acid's synonymous group matter for relative adaptiveness and CAI. "
            "Counts reflect one representative (longest) transcript per gene, so "
            "they are not weighted by expression: CAI built on them is a codon-bias "
            "proxy, not a measured expression prediction. Re-derivable from "
            "source_url + source_sha256 by rerunning "
            "scripts/build_organism_tables.py."
        ),
        extra={
            "organism_common_name": spec.common_name,
            "assembly": spec.assembly,
            "database": spec.database,
            "database_release": spec.release,
            "source_url": spec.url,
            "source_sha256": source_sha,
            "total_codons_counted": total,
            "filters": stats.as_dict(),
            "rebuild_command": (
                f"python scripts/build_organism_tables.py {spec.key}"
            ),
            "license_note": (
                "Ensembl annotation and sequence data are made freely available by "
                "EMBL-EBI (see https://www.ebi.ac.uk/about/terms-of-use). Cite "
                "Ensembl (Harrison et al., Nucleic Acids Res 2024, "
                "doi:10.1093/nar/gkad1049)." + spec.extra_citation
            ),
        },
    )
    return Path(written)


def _verify_against_committed(spec: OrganismSpec, rebuilt_tsv: Path) -> list[str]:
    """Diff a rebuilt table AND its sidecar against the committed ones.

    Checking only the TSV would leave the *provenance* unverified -- and the
    provenance is what carries the re-derivation trail, so a sidecar claiming
    the wrong source URL or digest would pass a TSV-only check while telling a
    third party to rebuild from the wrong data. Every sidecar field is compared
    except ``retrieved``, which is a wall-clock stamp of when the rebuild ran and
    is expected to differ.

    Returns:
        A list of human-readable mismatch descriptions (empty when clean).
    """
    problems: list[str] = []
    committed_tsv = DATA_DIR / rebuilt_tsv.name
    if not committed_tsv.is_file():
        return [f"{spec.key}: no committed table"]
    if committed_tsv.read_bytes() != rebuilt_tsv.read_bytes():
        problems.append(f"{spec.key}: committed TSV differs from rebuild")

    name = f"{spec.key}.provenance.json"
    committed_side = DATA_DIR / name
    rebuilt_side = rebuilt_tsv.parent / name
    if not committed_side.is_file():
        return [*problems, f"{spec.key}: no committed provenance sidecar"]
    committed = json.loads(committed_side.read_text(encoding="utf-8"))
    rebuilt = json.loads(rebuilt_side.read_text(encoding="utf-8"))
    for key in sorted(set(committed) | set(rebuilt) - {"retrieved"}):
        if key == "retrieved":
            continue
        if committed.get(key) != rebuilt.get(key):
            problems.append(f"{spec.key}: provenance field {key!r} differs from rebuild")
    return problems


def main(argv: list[str] | None = None) -> int:
    """Rebuild (or verify) the bundled organism tables."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "organisms", nargs="*",
        help="organism keys to build (default: all pinned organisms)",
    )
    parser.add_argument(
        "--cache-dir", default=None,
        help="where to keep downloaded CDS archives (default: a temp directory)",
    )
    parser.add_argument(
        "--verify", action="store_true",
        help="rebuild into a temp directory and diff against the committed TSVs "
             "instead of overwriting them",
    )
    args = parser.parse_args(argv)

    chosen = SPECS
    if args.organisms:
        by_key = {spec.key: spec for spec in SPECS}
        unknown = sorted(set(args.organisms) - set(by_key))
        if unknown:
            parser.error(f"unknown organism(s): {', '.join(unknown)}")
        chosen = tuple(by_key[key] for key in args.organisms)

    with tempfile.TemporaryDirectory() as tmp:
        cache_dir = Path(args.cache_dir) if args.cache_dir else Path(tmp) / "cache"
        out_dir = Path(tmp) / "out" if args.verify else DATA_DIR
        out_dir.mkdir(parents=True, exist_ok=True)

        mismatched: list[str] = []
        for spec in chosen:
            written = build_one(spec, cache_dir, out_dir)
            if args.verify:
                mismatched.extend(_verify_against_committed(spec, written))

        if mismatched:
            print("\nVERIFY FAILED:", file=sys.stderr)
            for line in mismatched:
                print(f"  {line}", file=sys.stderr)
            return 1

    print("\nverified" if args.verify else "\nwrote tables into " + str(DATA_DIR))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
