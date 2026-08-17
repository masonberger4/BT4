#!/usr/bin/env python3
"""Rebuild a bundled tRNA gene-copy-number table from GtRNAdb (CLAUDE.md §8).

BT4's tAI weights are derived from real tRNA gene copy numbers, and §8 requires
every bundled dataset to be **re-derivable**: a third party must be able to rebuild
the exact shipped bytes from a pinned public source. The codon-usage tables have had
that pipeline (``build_organism_tables.py``); the tRNA tables did not -- they were
counted ad hoc. This script closes that gap.

What it does:

* downloads (or reads) a GtRNAdb ``<genome>-tRNAs.fa`` prediction set,
* verifies the source file's SHA-256 against the pinned value, so a silently
  updated upstream release aborts the build instead of changing the shipped table,
* counts genes per anticodon under documented, uniform filtering rules, and
* writes ``<organism>.trna.tsv`` plus a provenance sidecar recording the source URL,
  its SHA-256, the genome build, the totals, and the per-filter drop tally.

It **refuses rather than fabricates**: an unparsable header, a source-hash mismatch,
or a gene total that disagrees with the expected count aborts the write.

Usage::

    python scripts/build_trna_tables.py --organism escherichia_coli
    python scripts/build_trna_tables.py --organism escherichia_coli --verify
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent.parent / "src" / "bt4" / "biomodels" / "codon" / "data"

# Header form: >..._tRNA-<AA>-<ANTICODON>-<n>-<m> (tRNAscan-SE ID: ...) ...
_HEADER = re.compile(r"tRNA-(\w+)-([ACGTN]{3})-")

# Amino-acid labels that are never counted, and why.
#   SeC  selenocysteine: inserted by a dedicated recoding mechanism at UGA, not a
#        synonymous choice, and excluded from the bundled human table too.
#   Und  undetermined: tRNAscan could not assign an identity.
#   Sup  suppressor tRNAs: read stop codons, not sense codons.
_EXCLUDED_AA = frozenset({"SeC", "Und", "Sup", "Undet"})


@dataclass(frozen=True)
class TrnaSpec:
    """A pinned GtRNAdb source for one organism."""

    organism: str
    url: str
    sha256: str
    genome: str
    super_kingdom: str
    expected_genes: int
    retrieved: str
    note: str


SPECS: dict[str, TrnaSpec] = {
    "escherichia_coli": TrnaSpec(
        organism="escherichia_coli",
        url=(
            "https://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/"
            "eschColi_K_12_MG1655-tRNAs.fa"
        ),
        sha256="bd31c3b159caea4e5526a51824f15d5595da34cc2ef5421a2f1967edb253a8df",
        genome="Escherichia coli K-12 MG1655 (U00096.3); chromosomal, no plasmid",
        super_kingdom="bacteria",
        expected_genes=86,
        retrieved="2026-08-17",
        note=(
            "89 predictions minus 1 SeC (selenocysteine) and 2 undetermined = 86 "
            "tRNA genes, matching the standard count for K-12 MG1655. Initiator "
            "fMet, elongator Met and the lysidine-modified tRNA-Ile2 all carry the "
            "anticodon CAT and are counted together under that key -- Ile2's reading "
            "of AUA is applied by the tAI model's prokaryotic (sking=1) term, which "
            "the dos Reis reference implements as a constant independent of copy "
            "number, so no information is lost by the merge."
        ),
    ),
}


def fetch(spec: TrnaSpec, cache: Path | None) -> bytes:
    """Return the source FASTA bytes, from ``cache`` when present, else the URL."""
    if cache is not None and cache.is_file():
        return cache.read_bytes()
    with urllib.request.urlopen(spec.url, timeout=120) as response:
        payload: bytes = response.read()
    if cache is not None:
        cache.write_bytes(payload)
    return payload


def count_anticodons(fasta: str) -> tuple[dict[str, int], dict[str, int]]:
    """Count tRNA genes per anticodon.

    Returns:
        ``(counts, dropped)`` -- the anticodon counts, and a per-reason tally of
        the predictions that were excluded.

    Raises:
        ValueError: If a header cannot be parsed (never guess an identity).
    """
    counts: collections.Counter[str] = collections.Counter()
    dropped: collections.Counter[str] = collections.Counter()
    for line in fasta.splitlines():
        if not line.startswith(">"):
            continue
        match = _HEADER.search(line)
        if match is None:
            raise ValueError(f"unparsable tRNA header (refusing to guess): {line[:120]!r}")
        amino_acid, anticodon = match.group(1), match.group(2)
        if amino_acid in _EXCLUDED_AA:
            dropped[amino_acid] += 1
            continue
        if "N" in anticodon:
            dropped["ambiguous_anticodon"] += 1
            continue
        counts[anticodon] += 1
    return dict(counts), dict(dropped)


def render_tsv(counts: dict[str, int]) -> str:
    """Render the anticodon counts as the bundled ``anticodon<TAB>count`` TSV."""
    lines = ["anticodon\tcount"]
    lines.extend(f"{ac}\t{counts[ac]}" for ac in sorted(counts))
    return "\n".join(lines) + "\n"


def render_provenance(
    spec: TrnaSpec, counts: dict[str, int], dropped: dict[str, int], tsv: str
) -> str:
    """Render the provenance sidecar as deterministic JSON.

    ``sha256`` is the hash of the rendered table bytes (what the loader checks);
    ``source_sha256`` pins the upstream FASTA the table was counted from.
    """
    payload = {
        "source": f"GtRNAdb (UCSC Lowe Lab), {spec.organism} tRNA gene set",
        "source_url": spec.url,
        "source_sha256": spec.sha256,
        "build": (
            "tRNAscan-SE predictions from the pinned GtRNAdb FASTA; genes counted "
            "per anticodon. Excluded: " + ", ".join(sorted(_EXCLUDED_AA)) + ". "
            "Rebuild with scripts/build_trna_tables.py --organism "
            f"{spec.organism} (--verify diffs against the committed bytes)."
        ),
        "genome": spec.genome,
        "super_kingdom": spec.super_kingdom,
        "retrieved": spec.retrieved,
        "sha256": hashlib.sha256(tsv.encode("utf-8")).hexdigest(),
        "total_genes": sum(counts.values()),
        "dropped": dict(sorted(dropped.items())),
        "note": spec.note
        + " Source database: GtRNAdb (Chan & Lowe, Nucleic Acids Res 2016, "
        "doi:10.1093/nar/gkv1309), predictions by tRNAscan-SE (Chan et al. 2021, "
        "doi:10.1093/nar/gkab688). GtRNAdb states no explicit data license: "
        "citation-gated academic use, NOT a CC/public-domain grant.",
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def build(spec: TrnaSpec, *, verify: bool, cache: Path | None) -> int:
    """Build (or verify) one organism's table. Returns a process exit code."""
    payload = fetch(spec, cache)
    digest = hashlib.sha256(payload).hexdigest()
    if digest != spec.sha256:
        print(
            f"ERROR: source hash mismatch for {spec.organism}.\n"
            f"  expected {spec.sha256}\n  got      {digest}\n"
            "  The upstream release changed. Refusing to silently rebuild a "
            "different table; review the new release and update the pin.",
            file=sys.stderr,
        )
        return 2

    counts, dropped = count_anticodons(payload.decode("utf-8"))
    total = sum(counts.values())
    if total != spec.expected_genes:
        print(
            f"ERROR: counted {total} genes for {spec.organism}, expected "
            f"{spec.expected_genes}. Refusing to write a table that disagrees with "
            "its own documented total.",
            file=sys.stderr,
        )
        return 2

    tsv_path = _DATA_DIR / f"{spec.organism}.trna.tsv"
    prov_path = _DATA_DIR / f"{spec.organism}.trna.provenance.json"
    tsv = render_tsv(counts)
    prov = render_provenance(spec, counts, dropped, tsv)

    if verify:
        ok = True
        for path, expected in ((tsv_path, tsv), (prov_path, prov)):
            actual = path.read_text(encoding="utf-8") if path.is_file() else None
            if actual != expected:
                print(f"MISMATCH: {path.name} differs from a fresh rebuild", file=sys.stderr)
                ok = False
        if ok:
            print(f"verified: {spec.organism} reproduces the committed bytes exactly")
        return 0 if ok else 1

    tsv_path.write_text(tsv, encoding="utf-8")
    prov_path.write_text(prov, encoding="utf-8")
    print(
        f"wrote {tsv_path.name} ({total} genes, {len(counts)} anticodons) "
        f"and {prov_path.name}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--organism", required=True, choices=sorted(SPECS))
    parser.add_argument(
        "--verify",
        action="store_true",
        help="rebuild and diff against the committed bytes instead of writing",
    )
    parser.add_argument(
        "--cache", type=Path, default=None, help="reuse/store the source FASTA here"
    )
    args = parser.parse_args(argv)
    return build(SPECS[args.organism], verify=args.verify, cache=args.cache)


if __name__ == "__main__":
    raise SystemExit(main())
