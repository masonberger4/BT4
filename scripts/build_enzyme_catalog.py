#!/usr/bin/env python3
"""Rebuild BT4's restriction-enzyme catalog from a pinned REBASE release.

BT4's enzyme catalog used to be seventeen hand-typed name/site pairs described
only as "textbook-correct" -- no source, no version, no way to check them. This
script replaces that with the same discipline the codon tables get (CLAUDE.md
§8): every recognition sequence is **derived from REBASE**, the authoritative
restriction-enzyme database, and the shipped TSV is re-derivable from a
version-pinned source file whose own SHA-256 is recorded.

Usage::

    python scripts/build_enzyme_catalog.py            # rebuild the catalog
    python scripts/build_enzyme_catalog.py --verify   # diff against the shipped one

What is selected, and why:

* **Restriction enzymes** (REBASE ``ET`` starting with ``R``) -- in this release
  that is Type II (``R2``), Type II bifunctional/IIG (``RM2``), Type II
  modification-dependent/IIM (``R2*``, e.g. **DpnI**), and one Type III
  (``EcoP15I``). Methyltransferases and homing endonucleases are excluded: they
  do not define a site a synthetic coding sequence must avoid. The
  modification-dependent ones are kept on purpose -- avoiding DpnI's ``GATC`` is
  a mainstream goal precisely because a plasmid from a dam+ strain *is*
  Dam-methylated and *is* cut by it.
* **Commercially available only** (REBASE ``CR`` naming at least one supplier).
  An enzyme nobody sells is not one a user will clone with, and the restriction
  is what keeps the shipped subset small and practical rather than a bulk copy
  of the database.
* **A single, fully-specified IUPAC recognition sequence** of 4-20 bases.
  Entries whose site is unknown (``?``) or given as multiple alternatives are
  skipped rather than guessed at.

Isoschizomers are deliberately kept (e.g. both ``KpnI`` and ``Acc65I`` for
``GGTACC``): users name the enzyme they actually own, and BT4 bans the site
either way.

**On re-derivability, honestly.** REBASE publishes only a *moving* URL for the
current release -- there is no versioned permalink for a specific version (the
obvious candidates 404). So the pin that actually holds is the **digest**:
``REBASE_SHA256`` is checked on every run and a mismatch aborts, which detects
drift but cannot resurrect the old bytes. Re-deriving the shipped catalog
therefore requires REBASE version 608 specifically; a user holding a later
release will get an abort telling them so, not a silently different catalog.

**On REBASE and licensing.** REBASE is Copyright (c) Dr. Richard J. Roberts and
is made freely available to the community for use *with citation*; it is not a
CC or public-domain grant. What BT4 ships is a derived subset -- enzyme names
paired with their recognition sequences, each a published fact carried in its
suppliers' own catalogs -- not a redistribution of the REBASE database
(no cut positions, organisms, references, methylation data, or the other fields
that make up the database). The provenance sidecar records the REBASE version,
URL, and source digest so any user can re-derive and re-verify. This mirrors how
BT4 already treats GtRNAdb's tRNA data: attribute precisely, state the terms
accurately, never overclaim a license.

This is a maintainer tool: it reaches the network and writes into the package
data directory. It is not imported by the library, and BT4 never fetches
anything at runtime.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shutil
import sys
import tempfile
import urllib.request
from collections.abc import Iterator
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from bt4.constraints.iupac import is_iupac, reverse_complement_iupac  # noqa: E402

DATA_DIR = REPO_ROOT / "src" / "bt4" / "constraints" / "data"
TSV_NAME = "rebase_enzymes.tsv"
PROVENANCE_NAME = "rebase_enzymes.provenance.json"

# Version-pinned REBASE download. REBASE publishes a monotonically increasing
# version (the "608" in bairoch.608); the file at this URL is the current one, so
# the digest below is what pins the bytes -- a changed digest means REBASE moved
# and the catalog must be rebuilt and re-reviewed deliberately.
REBASE_URL = "http://rebase.neb.com/rebase/link_bairoch"
REBASE_VERSION = "608"
REBASE_SHA256 = "3c60506f3fd49f5c18afb993d1ed680f792d642e90b291ef4427d143ebcfa40c"

_MIN_SITE = 4
# Upper bound is a parse-sanity guard, not a biological limit: the longest real
# commercially-available site in REBASE 608 is XcmI's 15-base CCANNNNNNNNNTGG, so
# anything past 20 is a malformed field rather than an enzyme. An earlier cap of
# 12 silently dropped SfiI (GGCCNNNNNGGCC) and eleven others.
_MAX_SITE = 20

# REBASE enzyme-type codes BT4 ships: every restriction enzyme (ET starting with
# "R"), which in this release resolves to
#   R2   Type II restriction enzyme                              (547)
#   RM2  Type II bifunctional restriction-modification / IIG      (17)
#   R2*  Type II MODIFICATION-DEPENDENT, subtype IIM              (11)
#   R3   Type III (EcoP15I only)                                   (1)
# Methyltransferases (M*) and homing endonucleases (IE) are excluded: they do not
# define a site a synthetic coding sequence must avoid.
#
# The modification-dependent ones are kept deliberately, and it would be a
# mistake to drop them as "cannot cut synthetic DNA". DpnI is R2*, and avoiding
# its GATC site is a mainstream design goal precisely BECAUSE a plasmid grown in
# a standard dam+ strain is Dam-methylated and is cut by it (DpnI template
# digestion is a routine QuikChange/Gibson step). Eight of the R2* entries also
# share a site byte-for-byte with a plain-R2 enzyme already in the catalog
# (DpnI GATC = Bsp143I, GlaI GCGC = HhaI, KroI GCCGGC = NaeI), so excluding them
# would delete an alias a user may reasonably type, not a constraint.
_SHIPPED_TYPES_PREFIX = "R"


def sha256_file(path: Path) -> str:
    """Return the SHA-256 hex digest of ``path`` (streamed, constant memory)."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    """Fetch ``url`` to ``dest`` unless it is already cached there."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  cached  {dest.name}")
        return
    print(f"  fetching {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, partial.open("wb") as out:
        shutil.copyfileobj(response, out)
    partial.replace(dest)


def iter_records(path: Path) -> Iterator[dict[str, str]]:
    """Yield REBASE Bairoch-format records as ``{tag: value}`` mappings.

    The format is one ``TAG   value`` line per field, records separated by
    ``//``. Repeated tags (a multi-line field) are concatenated. The file is
    latin-1: it carries accented author names in its reference lines.
    """
    current: dict[str, str] = {}
    with path.open(encoding="latin-1") as handle:
        for line in handle:
            if line.startswith("//"):
                if current:
                    yield current
                current = {}
                continue
            if len(line) < 5:
                continue
            tag, value = line[:2], line[5:].rstrip()
            if tag in ("ID", "ET", "RS", "CR", "PT"):
                current[tag] = current.get(tag, "") + value
    if current:
        yield current


def recognition_site(
    record: dict[str, str], stats: dict[str, int] | None = None
) -> str | None:
    """Extract a single fully-specified IUPAC recognition site, or ``None``.

    A REBASE ``RS`` field looks like ``GAATTC, 1;`` -- the site, then the cut
    offset. BT4 bans the *site* and does not model cut positions, so the offset
    is dropped.

    An **asymmetric** enzyme is listed as two entries, one per strand, e.g. BsaI
    ``GGTCTC, 7; GAGACC, -5;``. That is one site seen from both directions, not
    two sites, and BT4 bans every site's reverse complement anyway -- so the
    first strand is taken and the second is *verified* to be its reverse
    complement rather than assumed. Dropping these outright would lose exactly
    the Type IIS enzymes modern cloning depends on (BsaI, BsmBI, BbsI, SapI).

    Anything unknown (``?``), non-IUPAC, outside the sane length range, or
    listing two genuinely different sites yields ``None``, so it is skipped
    rather than guessed at.
    """
    def reject(reason: str) -> None:
        if stats is not None:
            stats[f"rejected_{reason}"] = stats.get(f"rejected_{reason}", 0) + 1

    raw = record.get("RS", "").strip()
    if not raw or raw.startswith("?"):
        reject("site_unknown")
        return None
    alternatives = [part for part in raw.split(";") if part.strip()]
    sites = [part.split(",")[0].strip().upper() for part in alternatives]
    if not sites or not all(is_iupac(site) for site in sites if site):
        reject("site_not_iupac")
        return None
    if len(sites) > 2 or (
        len(sites) == 2 and reverse_complement_iupac(sites[0]) != sites[1]
    ):
        # Not a two-strand listing of one site: genuinely several targets.
        reject("site_multiple")
        return None
    site = sites[0]
    if not site or not _MIN_SITE <= len(site) <= _MAX_SITE:
        # Counted, never silent: an unexplained length cut is how SfiI went
        # missing from a shipped catalog without anyone noticing.
        reject("site_length")
        return None
    return site


def select_enzymes(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    """Select the shipped enzyme subset from a REBASE file.

    Returns:
        ``(catalog, stats)`` -- a ``{name: site}`` mapping and the tally of what
        was considered and dropped at each step (stamped into the provenance, so
        the selection is auditable rather than a black box).
    """
    stats = {
        "records_in_source": 0,
        "restriction_enzymes": 0,
        "commercially_available": 0,
        "with_usable_site": 0,
        "duplicate_ids_merged": 0,
        "rejected_site_unknown": 0,
        "rejected_site_not_iupac": 0,
        "rejected_site_multiple": 0,
        "rejected_site_length": 0,
    }
    catalog: dict[str, str] = {}
    for record in iter_records(path):
        stats["records_in_source"] += 1
        if not record.get("ET", "").strip().startswith(_SHIPPED_TYPES_PREFIX):
            continue
        stats["restriction_enzymes"] += 1
        if record.get("CR", ".").strip() in (".", ""):
            continue
        stats["commercially_available"] += 1
        site = recognition_site(record, stats)
        name = record.get("ID", "").strip()
        if not site or not name:
            continue
        stats["with_usable_site"] += 1
        previous = catalog.get(name)
        if previous is not None:
            # REBASE lists some enzymes (the Type IIB ones that cut both sides,
            # e.g. AjuI/AloI/PsrI) as TWO records under one name, each leading
            # with a different strand. That is one double-stranded site seen
            # twice, not two sites -- and BT4 bans both strands regardless. Keep
            # the first-listed strand (the source is digest-pinned, so file order
            # is reproducible) and *verify* the equivalence rather than assuming
            # it: a name that genuinely disagrees with itself is a data problem to
            # resolve deliberately, not to settle by whichever record comes last.
            if previous != site and reverse_complement_iupac(previous) != site:
                raise SystemExit(
                    f"REBASE lists {name} with conflicting sites {previous!r} and "
                    f"{site!r} that are not reverse complements; resolve "
                    "deliberately rather than taking whichever comes last"
                )
            stats["duplicate_ids_merged"] += 1
            continue
        catalog[name] = site
    stats["shipped"] = len(catalog)
    return catalog, stats


def render_tsv(catalog: dict[str, str]) -> bytes:
    """Render the catalog as sorted ``name<TAB>site`` TSV bytes."""
    lines = ["enzyme\tsite"]
    lines.extend(f"{name}\t{catalog[name]}" for name in sorted(catalog))
    return ("\n".join(lines) + "\n").encode("utf-8")


def build(cache_dir: Path, out_dir: Path) -> tuple[Path, Path]:
    """Download REBASE, select the subset, and write the TSV + provenance."""
    archive = cache_dir / "rebase_link_bairoch.txt"
    download(REBASE_URL, archive)
    digest = sha256_file(archive)
    if digest != REBASE_SHA256:
        raise SystemExit(
            "REBASE source digest changed:\n"
            f"  expected {REBASE_SHA256}\n"
            f"  got      {digest}\n"
            "REBASE has published a new version. Update REBASE_VERSION and "
            "REBASE_SHA256 in this script and rebuild deliberately, so the "
            "catalog never changes underneath a release by accident."
        )

    catalog, stats = select_enzymes(archive)
    if not catalog:
        raise SystemExit("no enzymes selected -- refusing to write an empty catalog")

    tsv_bytes = render_tsv(catalog)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / TSV_NAME
    tsv_path.write_bytes(tsv_bytes)

    provenance = {
        "source": f"REBASE (The Restriction Enzyme Database) version {REBASE_VERSION}",
        "source_url": REBASE_URL,
        "source_sha256": digest,
        "rebase_version": REBASE_VERSION,
        "build": (
            "derived by scripts/build_enzyme_catalog.py from the pinned REBASE "
            "Bairoch-format file: restriction enzymes (REBASE ET starting with "
            "R -- here Type II R2, Type II bifunctional/IIG RM2, Type II "
            "modification-dependent/IIM R2* such as DpnI, and one Type III "
            "EcoP15I) that name at least one commercial supplier (CR) and have a "
            f"single fully-specified IUPAC recognition site of {_MIN_SITE}-{_MAX_SITE} "
            "bases. Methyltransferases and homing endonucleases are excluded. Cut "
            "positions, source organisms, references and methylation "
            "data are not extracted -- only the name and its recognition sequence, "
            "which is all BT4 needs to ban a site. Isoschizomers are kept so a "
            "user can name the enzyme they actually own."
        ),
        "enzyme_count": len(catalog),
        "selection": stats,
        "retrieved": datetime.date.today().isoformat(),
        "sha256": hashlib.sha256(tsv_bytes).hexdigest(),
        "rebuild_command": "python scripts/build_enzyme_catalog.py",
        "note": (
            "Real recognition sequences derived from REBASE, NOT hand-typed "
            "values. Each entry is an enzyme name paired with the IUPAC "
            "recognition sequence BT4 bans (together with its reverse complement, "
            "since sites are double-stranded). BT4 does not model cut positions, "
            "star activity, methylation sensitivity, or buffer conditions: a site "
            "being absent means the recognition sequence does not occur, not that "
            "a digest will behave as you expect."
        ),
        "license_note": (
            "REBASE is Copyright (c) Dr. Richard J. Roberts and is made freely "
            "available to the community for use WITH CITATION; it is not a "
            "CC/public-domain grant. What is bundled here is a derived subset "
            "(enzyme name + recognition sequence for commercially available "
            "Type II enzymes -- each a published fact carried in its suppliers' "
            "catalogs), not a redistribution of the REBASE database. Cite Roberts "
            "RJ, Vincze T, Posfai J, Macelis D. REBASE - a database for DNA "
            "restriction and modification: enzymes, genes and genomes. Nucleic "
            "Acids Res (2023), doi:10.1093/nar/gkac975."
        ),
        "rederivation_note": (
            "REBASE publishes only a MOVING url for its current release, with no "
            "versioned permalink, so source_url alone does not pin anything: the "
            "pin that holds is source_sha256, checked on every build. Re-deriving "
            f"these exact bytes therefore requires REBASE version {REBASE_VERSION} "
            "specifically. A later release makes the build abort with the digest "
            "mismatch rather than silently producing a different catalog."
        ),
    }
    provenance_path = out_dir / PROVENANCE_NAME
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"  {stats['records_in_source']} records -> {stats['restriction_enzymes']} type-II -> "
        f"{stats['commercially_available']} commercial -> {len(catalog)} shipped"
    )
    return tsv_path, provenance_path


def main(argv: list[str] | None = None) -> int:
    """Rebuild (or verify) the bundled REBASE-derived enzyme catalog."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--cache-dir", default=None,
        help="where to keep the downloaded REBASE file (default: a temp directory)",
    )
    parser.add_argument(
        "--verify", action="store_true",
        help="rebuild into a temp directory and diff against the shipped catalog "
             "instead of overwriting it",
    )
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory() as tmp:
        cache_dir = Path(args.cache_dir) if args.cache_dir else Path(tmp) / "cache"
        out_dir = Path(tmp) / "out" if args.verify else DATA_DIR
        print("rebase enzyme catalog:")
        tsv_path, provenance_path = build(cache_dir, out_dir)

        if args.verify:
            problems: list[str] = []
            shipped_tsv = DATA_DIR / TSV_NAME
            if not shipped_tsv.is_file():
                problems.append("no shipped catalog TSV")
            elif shipped_tsv.read_bytes() != tsv_path.read_bytes():
                problems.append("shipped catalog TSV differs from the rebuild")
            shipped_prov = DATA_DIR / PROVENANCE_NAME
            if not shipped_prov.is_file():
                problems.append("no shipped provenance sidecar")
            else:
                # 'retrieved' is a wall-clock stamp; everything else must match.
                old = json.loads(shipped_prov.read_text(encoding="utf-8"))
                new = json.loads(provenance_path.read_text(encoding="utf-8"))
                old.pop("retrieved", None)
                new.pop("retrieved", None)
                if old != new:
                    problems.append("shipped provenance differs from the rebuild")
            if problems:
                print("\nVERIFY FAILED:", file=sys.stderr)
                for line in problems:
                    print(f"  {line}", file=sys.stderr)
                return 1
            print("\nverified")
            return 0

    print(f"\nwrote catalog into {DATA_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
