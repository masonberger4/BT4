#!/usr/bin/env python3
"""Build BT4's *highly-expressed* codon-usage reference tables.

This is the reference set CAI was defined on. Sharp & Li (1987) derive the
relative adaptiveness ``w = f/f_max`` from a set of **very highly expressed
genes**, so ``w = 1`` marks the codon that translation in that organism actually
*prefers*. A genome-wide count answers a different question -- which codon is
most *common* across all genes, most of which are lowly expressed -- and under
weak translational selection that answer is set by mutation and GC bias instead.
The two reference sets disagree sharply exactly where it matters: in *E. coli*
the genome-wide table's top Phe codon is ``TTT`` while highly-expressed genes
prefer ``TTC`` (and eight amino acids differ in total; see the table below).

So this script builds a **second, declared reference set** per organism rather
than "correcting" the genome-wide one. Both ship; both are counted from public
pinned sources; every result records which one it used.

How the reference set is chosen (all of it stamped into each sidecar):

* **Expression comes from PaxDb** (v6.1, CC BY 4.0), the whole-organism
  *integrated* dataset -- a weighted consensus of many published quantitative
  proteomics studies, in ppm. Using measured protein abundance instead of a
  hand-picked gene list is the part that is not 1987: the reference set is
  data, re-derivable by anyone, not a curated opinion.
* **Sequences come from the same release-pinned Ensembl CDS sets** the
  genome-wide tables are counted from, filtered by the same documented rules
  (one representative CDS per gene, ACGT-only, in-frame, ATG start, terminal
  stop, no internal stop, no alternate-haplotype/patch loci). Sharing the filter
  is what makes the two tables comparable: they differ in *which genes*, not in
  how a gene is read.
* **PaxDb protein IDs are joined to Ensembl genes with no third-party mapping
  layer** -- against the pinned release's own peptide FASTA. An identifier that
  resolves to two different genes is dropped as ambiguous, never guessed.
* **Organelle-encoded genes are excluded.** Mitochondria and plastids translate
  with a *different genetic code* and their own tRNA pool, and they are never
  the target of a BT4 design (which is a nuclear transgene), so their codon
  counts are not evidence about nuclear translation. The exclusion is justified
  by that, not by a claim about how many would have ranked highly: most organelle
  CDS never reach the exclusion step at all, because under the standard code they
  read as having internal stops and are dropped as invalid first. Both numbers
  are therefore stamped -- ``genes_excluded_organelle_encoded`` (what this filter
  removed: 1 each for mouse, rat and yeast, 0 elsewhere) and
  ``organelle_records_in_cds_source`` (how many were in the annotation at all) --
  so a zero is never mistaken for "this organism has no organelle-encoded
  proteins".
* **Top ``N = 300`` genes by abundance**, ties broken by gene ID. N is not a
  free knob: it was chosen as the smallest size on a tested grid
  (50/100/200/300/500/1000/2000) at which *every* bundled organism observes all
  64 codons, so no shipped table needs smoothing -- an invented number in a
  reference table is exactly what BT4 refuses to ship. Below 300, yeast alone
  leaves ``CGA``/``CGG`` unobserved; far above it the reference set dilutes back
  toward the genome-wide answer (at N=2000 the yeast and mouse tables agree with
  their genome-wide counterparts at every amino acid).

Amino acids whose most-used codon differs from the genome-wide table at N=300
(measured, not asserted -- ``--report`` reprints this):

    C. elegans 11 · E. coli 8 · zebrafish 7 · yeast 5 · mouse 3 · rat 2 ·
    human 2 · fruit fly 2

**Known bias, measured rather than assumed: the shared ATG-start filter.**
Bacteria initiate translation at ``GTG`` and ``TTG`` as well as ``ATG``, and the
inherited validity filter drops those genes -- 409 of 4,239 *E. coli* CDS records
(9.6%), among them *tufA* and *hupB*, both classic highly-expressed genes. The
effect on the shipped table was measured, not guessed: relaxing the filter to
accept all three bacterial starts changes 16 of the 300 selected genes and moves
**no** amino acid's most-used codon, so the *E. coli* design BT4 delivers is the
same either way. The filter is therefore left as it is -- identical to the
genome-wide builder's, which is what keeps the two tables comparable -- and the
drop tally is stamped in every sidecar. Handling alternative starts properly
(which also means deciding whether an initiator codon is a codon *choice* at all,
since the ribosome uses fMet-tRNA regardless) is a change to both builders and is
queued separately.

**Arabidopsis is deliberately absent.** PaxDb identifies *A. thaliana* proteins
by UniProt accession, which the pinned Ensembl Plants annotation does not carry,
so joining them would need an unpinned external mapping. Rather than ship a
table built on a guess, BT4 ships none for that organism and says so -- the
genome-wide table stays its only (honestly labeled) option.

**What this is not.** A high CAI against a highly-expressed reference set is a
better-founded proxy than a genome-wide one, but it is still a proxy, not a
measured expression prediction: Welch et al. (PLoS ONE 2009,
doi:10.1371/journal.pone.0007002) found an *E. coli* variant built by maximizing
exactly this quantity expressed at a fraction of alternatives. That is why BT4
keeps CAI as one axis of a vector rather than the objective (CLAUDE.md §1,
§10.7).

Usage::

    python scripts/build_highly_expressed_tables.py                 # all
    python scripts/build_highly_expressed_tables.py escherichia_coli
    python scripts/build_highly_expressed_tables.py --verify
    python scripts/build_highly_expressed_tables.py --report        # no writes

It is not imported by the library, and no BT4 *table* is
ever fetched at runtime -- the only runtime network access BT4 has at all is
the opt-in, explicitly-consented ASSP splice cross-check (CLAUDE.md §6).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_organism_tables import (  # noqa: E402
    DATA_DIR,
    SPECS,
    FilterStats,
    OrganismSpec,
    cds_region,
    download,
    iter_fasta,
    on_alt_locus,
    parse_ids,
    sha256_file,
    suspicious_kept_region,
    valid_cds,
)

from bt4.biomodels.codon.build import count_codons, write_table  # noqa: E402
from bt4.biomodels.codon.tables import (  # noqa: E402
    HIGHLY_EXPRESSED,
    REFERENCE_SET_SUFFIX,
    sha256_hex,
)
from bt4.domain.genetic_code import CODON_TABLE, STOP  # noqa: E402

# PaxDb release. Pinned, never "latest" -- the same rule the CDS sources follow.
PAXDB_RELEASE = "6.1"
PAXDB_BASE = f"https://pax-db.org/downloads/{PAXDB_RELEASE}/datasets"

# The reference-set size. See the module docstring for how it was chosen; the
# short version is that it is the smallest tested size at which no bundled
# organism leaves a codon unobserved.
DEFAULT_TOP_N = 300

# Region names that denote an organelle genome across the Ensembl divisions BT4
# pulls from: "MT" (human/mouse/rat), "Mito" (yeast), "mitochondrion_genome"
# (fly), plus the plastid spellings. Organelles use a different genetic code and
# their own tRNA pool, so their codon usage is not evidence about the nuclear
# translation machinery a BT4 design will actually meet.
_ORGANELLE_NAMES = ("mt", "mtdna", "mito", "pt", "pltd")
_ORGANELLE_SUBSTRINGS = ("mitochondri", "chloroplast", "plastid", "apicoplast")


def is_organelle_region(region: str) -> bool:
    """Whether an assembly region name denotes a mitochondrial/plastid genome."""
    name = region.strip().lower()
    if not name:
        return False
    return name in _ORGANELLE_NAMES or any(s in name for s in _ORGANELLE_SUBSTRINGS)


@dataclass(frozen=True)
class AbundanceSpec:
    """One organism's pinned protein-abundance source and ID-join inputs.

    Attributes:
        key: BT4's canonical organism key. Must name an
            :data:`build_organism_tables.SPECS` entry -- the CDS set counted here
            is *that* one, so the two reference sets stay comparable.
        taxon: NCBI taxonomy ID, which is also PaxDb's dataset directory and the
            prefix of every ``string_external_id`` in the file.
        paxdb_sha256: Expected SHA-256 of the PaxDb dataset file.
        pep_url: Release-pinned URL of the peptide FASTA used to resolve PaxDb
            protein IDs to Ensembl gene IDs. It must come from the *same*
            release as the CDS set, or the two would disagree about which
            transcripts exist.
        pep_sha256: Expected SHA-256 of that peptide FASTA.
        id_note: What PaxDb's identifiers are for this organism, in one clause.
            Recorded in the sidecar so a reader can see what was joined to what
            without re-downloading anything.
    """

    key: str
    taxon: int
    paxdb_sha256: str
    pep_url: str
    pep_sha256: str
    id_note: str

    @property
    def paxdb_filename(self) -> str:
        """The PaxDb whole-organism integrated dataset's file name."""
        return f"{self.taxon}-WHOLE_ORGANISM-integrated.txt"

    @property
    def paxdb_url(self) -> str:
        """The pinned URL of that dataset."""
        return f"{PAXDB_BASE}/{self.taxon}/{self.paxdb_filename}"


_ENS = "https://ftp.ensembl.org/pub/release-116/fasta"
_BACTERIA = (
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/bacteria/release-63/fasta/"
    "bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
)

SPECS_HE: tuple[AbundanceSpec, ...] = (
    AbundanceSpec(
        key="homo_sapiens",
        taxon=9606,
        paxdb_sha256="dbc2566f94d85117adb686a5acd99fc9fd1bad0ed7db7e803e332e0320f76a4a",
        pep_url=f"{_ENS}/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
        pep_sha256="9b43da92651b35814597af6a8b18f500b768679a49fa4678224f384917ce7668",
        id_note="Ensembl protein IDs (ENSP…), matched without their version suffix",
    ),
    AbundanceSpec(
        key="saccharomyces_cerevisiae",
        taxon=4932,
        paxdb_sha256="6d1b614fa95d47408e3699b03e05cb4670eaa53d8a490fc8954d127e00261f10",
        pep_url=(
            f"{_ENS}/saccharomyces_cerevisiae/pep/"
            "Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
        ),
        pep_sha256="67ae76c720e52cec167378ec2b8cb5c5929360fc30de31abb37ae1035ef8455c",
        id_note="SGD systematic ORF names (YGR192C…), which Ensembl uses as gene IDs",
    ),
    AbundanceSpec(
        key="escherichia_coli",
        taxon=511145,
        paxdb_sha256="e2e52262945a2052c9905609659367f9e019cc9a79b9497be299b93c24d28877",
        pep_url=(
            f"{_BACTERIA}/pep/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
            ".ASM584v2.pep.all.fa.gz"
        ),
        pep_sha256="90a765e1fdc034b7ba261852f5402b3cfd76f0a3ffdd23185a7bf9592f3eadd8",
        id_note="Blattner b-numbers (b3495…), which Ensembl Bacteria uses as gene IDs",
    ),
    AbundanceSpec(
        key="mus_musculus",
        taxon=10090,
        paxdb_sha256="d1679d0d5381b5530f0168b158d1ecf2fd1e14a03a4ae7795d6db98d0a1ac3c8",
        pep_url=f"{_ENS}/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz",
        pep_sha256="480d4a6eb540b1cc26bb4a80ad8f8a50aba514791c576ed41362127f56809a37",
        id_note="Ensembl protein IDs (ENSMUSP…), matched without their version suffix",
    ),
    AbundanceSpec(
        key="rattus_norvegicus",
        taxon=10116,
        paxdb_sha256="b9307d8f8f0e1472ccaf93349e9281addeb59f593370a6a2a06ba0f87bc8f47c",
        pep_url=f"{_ENS}/rattus_norvegicus/pep/Rattus_norvegicus.GRCr8.pep.all.fa.gz",
        pep_sha256="860aef1226c1ac924cac38fb18327e2c8275a8b78c23e2ec20e60f573cffb228",
        id_note="Ensembl protein IDs (ENSRNOP…), matched without their version suffix",
    ),
    AbundanceSpec(
        key="danio_rerio",
        taxon=7955,
        paxdb_sha256="c6d06b562a5f41525ca36ad6d46c44783672ba284cdd94fbb16eda9e94f8e01a",
        pep_url=f"{_ENS}/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz",
        pep_sha256="554ad9e76101d96db674bc9eab1c07116d7c84ecfcd24519f3d70ae810e82ecc",
        id_note="Ensembl protein IDs (ENSDARP…), matched without their version suffix",
    ),
    AbundanceSpec(
        key="drosophila_melanogaster",
        taxon=7227,
        paxdb_sha256="5eb95d67fe05dd9db9d94b17eb8dbbd7f8fbc43a7dcf6439b7b9dc276d089411",
        pep_url=(
            f"{_ENS}/drosophila_melanogaster/pep/"
            "Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz"
        ),
        pep_sha256="a2175de6335a8af53935ccf2ca223a2335d0c0f70cbb2db48159645367f9a7ee",
        id_note="FlyBase polypeptide IDs (FBpp…), which Ensembl uses as protein IDs",
    ),
    AbundanceSpec(
        key="caenorhabditis_elegans",
        taxon=6239,
        paxdb_sha256="6b3438e6e7f392a14c82d9e21bcdc78f23070f54dc9c261c0d8868c974dc4f18",
        pep_url=(
            f"{_ENS}/caenorhabditis_elegans/pep/"
            "Caenorhabditis_elegans.WBcel235.pep.all.fa.gz"
        ),
        pep_sha256="78fb77ec7908a1a7ae2a9ab56b6194b0e7d25b2e1170779230a351fa9e9c378a",
        id_note="WormBase transcript-style protein IDs (ZK1010.1.1…), matched exactly",
    ),
)


@dataclass
class JoinStats:
    """How PaxDb's rows fared on their way to becoming a reference set."""

    paxdb_rows: int = 0
    unmatched_id: int = 0
    ambiguous_id: int = 0
    no_counted_cds: int = 0
    genes_joined: int = 0
    organelle_excluded: int = 0
    organelle_records: int = 0
    eligible_genes: int = 0
    by_protein_id: int = 0
    by_unversioned_protein_id: int = 0
    by_gene_id: int = 0

    def as_dict(self) -> dict[str, int]:
        """The stats as a provenance-ready mapping.

        Three families, and the key names carry which is which because they do
        not reconcile by subtraction:

        * ``rows_matched_*`` + ``rows_unmatched_*`` partition ``paxdb_rows``
          exactly -- every abundance row is one or the other.
        * ``rows_matched_whose_gene_has_no_counted_cds`` is a **subset** of the
          matched rows, not a fourth part of the partition: the identifier
          resolved, the gene just has no valid representative CDS.
        * ``genes_*`` count genes, after rows collapse onto them (one gene can
          carry several protein rows), so they are a different unit entirely.

        An earlier ``joined_via_*`` / ``dropped_*`` naming implied a single flat
        partition and invited arithmetic that does not close.
        """
        return {
            "paxdb_rows": self.paxdb_rows,
            "rows_matched_via_protein_id": self.by_protein_id,
            "rows_matched_via_unversioned_protein_id": self.by_unversioned_protein_id,
            "rows_matched_via_gene_id": self.by_gene_id,
            "rows_unmatched_identifier_ambiguous": self.ambiguous_id,
            "rows_unmatched_identifier_not_in_annotation": self.unmatched_id,
            "rows_matched_whose_gene_has_no_counted_cds": self.no_counted_cds,
            "genes_joined": self.genes_joined,
            "genes_excluded_organelle_encoded": self.organelle_excluded,
            "genes_eligible_for_reference_set": self.eligible_genes,
            "organelle_records_in_cds_source": self.organelle_records,
        }


@dataclass
class ReferenceSet:
    """The selected genes, in rank order, and how they were selected."""

    genes: list[str] = field(default_factory=list)
    names: list[str] = field(default_factory=list)
    sequences: list[str] = field(default_factory=list)
    stats: JoinStats = field(default_factory=JoinStats)
    filters: FilterStats = field(default_factory=FilterStats)

    @property
    def roster_sha256(self) -> str:
        """Digest of the ranked gene roster.

        Pins the exact reference set -- membership *and* order -- in one field,
        so a third party can prove they reproduced the same 300 genes without
        the sidecar having to carry all 300 IDs.
        """
        return sha256_hex(("\n".join(self.genes) + "\n").encode("utf-8"))


def parse_paxdb(path: Path) -> list[tuple[str, str, float]]:
    """Parse a PaxDb dataset file into ``(gene_name, protein_id, abundance)``.

    PaxDb's ``string_external_id`` column is ``<taxon>.<identifier>``; only the
    identifier half joins to an annotation, so the taxon prefix is stripped here.
    Comment lines (the ``#``-prefixed header block, which is where the dataset's
    provenance lives) and malformed rows are skipped.

    Raises:
        SystemExit: If the file yields no usable rows -- a silently empty
            reference set would produce a table built on nothing.
    """
    rows: list[tuple[str, str, float]] = []
    malformed = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3 or "." not in parts[1]:
            malformed += 1
            continue
        try:
            abundance = float(parts[2])
        except ValueError:
            malformed += 1
            continue
        rows.append((parts[0], parts[1].split(".", 1)[1], abundance))
    if not rows:
        raise SystemExit(f"{path.name}: no usable abundance rows -- refusing to continue")
    if malformed:
        print(f"  note: skipped {malformed} malformed PaxDb row(s)")
    return rows


def iter_pep_headers(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(protein_id, gene_id)`` from a gzipped Ensembl peptide FASTA."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            protein = header.split(None, 1)[0]
            for field_ in header.split():
                if field_.startswith("gene:"):
                    yield protein, field_[len("gene:") :]
                    break


@dataclass(frozen=True)
class IdIndex:
    """Peptide-FASTA identifiers resolved to gene IDs, by key kind.

    Attributes:
        exact: Protein IDs as Ensembl writes them, plus gene IDs mapped to
            themselves.
        unversioned: Protein IDs with their ``.N`` version suffix stripped.
        ambiguous: Keys that resolved to more than one gene and were therefore
            removed from both maps. Kept so a dropped row can be *reported* as
            ambiguous instead of being lumped in with identifiers the annotation
            simply does not contain -- those are different failures, and a
            reference set that quietly conflated them would misstate its own
            coverage.
    """

    exact: dict[str, str]
    unversioned: dict[str, str]
    ambiguous: frozenset[str]


def build_id_index(path: Path) -> IdIndex:
    """Index a peptide FASTA by every identifier PaxDb might use.

    Three key kinds are indexed, because PaxDb uses a different one per organism
    (see each spec's ``id_note``): the protein ID exactly as Ensembl writes it,
    the protein ID without its ``.N`` version suffix, and the gene ID. They are
    kept in separate maps so a join can report *which* kind matched rather than
    leaving the reader to guess.

    A key that resolves to more than one gene is dropped, never guessed at.

    **Ambiguity is judged across the maps, not within each one.** The unversioned
    key is synthesized here rather than read from the source, so for
    WormBase-style IDs that natively contain dots it can collide with a real
    identifier: given proteins ``ZK1010.1.1`` (gene B) and ``ZK1010.1`` (gene A),
    the key ``ZK1010.1`` means gene A in ``exact`` and gene B in ``unversioned``.
    Neither map alone sees a conflict, so a per-map check would resolve that key
    to whichever map is consulted first -- a wrong gene, reported as a clean
    join. Pooling the two maps' genes per key is what makes the collision this
    docstring describes actually get dropped.
    """
    exact: dict[str, set[str]] = {}
    unversioned: dict[str, set[str]] = {}
    for protein, gene in iter_pep_headers(path):
        exact.setdefault(protein, set()).add(gene)
        exact.setdefault(gene, set()).add(gene)
        stem = protein.rsplit(".", 1)[0]
        if stem != protein:
            unversioned.setdefault(stem, set()).add(gene)
    pooled: dict[str, set[str]] = {}
    for source in (exact, unversioned):
        for key, genes in source.items():
            pooled.setdefault(key, set()).update(genes)
    ambiguous = frozenset(key for key, genes in pooled.items() if len(genes) > 1)
    return IdIndex(
        exact={k: next(iter(v)) for k, v in exact.items() if k not in ambiguous},
        unversioned={k: next(iter(v)) for k, v in unversioned.items() if k not in ambiguous},
        ambiguous=ambiguous,
    )


def representative_cds_by_gene(
    path: Path,
) -> tuple[dict[str, tuple[str, str]], FilterStats, int]:
    """Select one valid CDS per gene, keeping each gene's assembly region.

    Identical selection to ``build_organism_tables.representative_cds`` -- the
    longest valid CDS per gene, ties broken by transcript ID -- so the two
    reference sets differ only in *which* genes they count. The region comes
    along because organelle-encoded genes have to be excluded downstream.

    Organelle-encoded **records** are tallied here, before the validity filter,
    and reported separately from the genes excluded at selection time. Without
    that, most organelle CDS never reach the exclusion step at all: they are
    translated with a different genetic code, so under the standard code they
    show internal stops and get dropped as invalid first. The exclusion count
    alone would then read ``0`` for human -- whose mitochondrial proteins are
    among the most abundant in any proteomics consensus -- and be mistaken for
    "this organism's abundance data contained no organelle-encoded genes".

    Returns:
        ``(gene -> (sequence, region), stats, organelle_records)``.
    """
    stats = FilterStats()
    organelle_records = 0
    best: dict[str, tuple[int, str, str, str]] = {}
    for header, raw in iter_fasta(path):
        stats.records += 1
        if is_organelle_region(cds_region(header)):
            organelle_records += 1
        if on_alt_locus(header):
            stats.dropped_alt_locus += 1
            continue
        if suspicious_kept_region(header):
            stats.kept_suspicious_region += 1
        seq = raw.upper()
        if not valid_cds(seq, stats):
            continue
        transcript, gene = parse_ids(header)
        candidate = (len(seq), transcript, seq, cds_region(header))
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
    return (
        {gene: (entry[2], entry[3]) for gene, entry in best.items()},
        stats,
        organelle_records,
    )


def select_reference_set(
    rows: list[tuple[str, str, float]],
    index: IdIndex,
    cds: dict[str, tuple[str, str]],
    filters: FilterStats,
    organelle_records: int,
    top_n: int,
) -> ReferenceSet:
    """Rank the joined genes by abundance and keep the top ``top_n``.

    Several PaxDb rows can point at one gene (one row per protein isoform); the
    gene takes the **highest** abundance among them, since the question being
    asked is how heavily that gene's product is made.

    Ranking is ``(-abundance, gene_id)``, so an abundance tie resolves the same
    way on every machine and the roster digest is reproducible (invariant #7
    reaching the data).
    """
    stats = JoinStats(paxdb_rows=len(rows), organelle_records=organelle_records)
    best: dict[str, tuple[float, str]] = {}
    for gene_name, identifier, abundance in rows:
        # Ambiguity is checked FIRST, before either lookup. Checking it only
        # after both maps miss would let an identifier that is ambiguous overall
        # -- but singular in whichever map answers -- through as a clean join,
        # which is exactly the guess this builder claims never to make.
        if identifier in index.ambiguous:
            stats.ambiguous_id += 1
            continue
        gene = index.exact.get(identifier)
        if gene is not None:
            stats.by_protein_id += 1 if identifier != gene else 0
            stats.by_gene_id += 1 if identifier == gene else 0
        else:
            gene = index.unversioned.get(identifier)
            if gene is None:
                # Reported separately from the ambiguous case above: an
                # identifier the annotation does not contain is a different
                # failure from one it resolves two ways.
                stats.unmatched_id += 1
                continue
            stats.by_unversioned_protein_id += 1
        if gene not in cds:
            stats.no_counted_cds += 1
            continue
        current = best.get(gene)
        if current is None or abundance > current[0]:
            best[gene] = (abundance, gene_name)
    stats.genes_joined = len(best)

    eligible = {gene: v for gene, v in best.items() if not is_organelle_region(cds[gene][1])}
    stats.organelle_excluded = len(best) - len(eligible)
    stats.eligible_genes = len(eligible)

    ranked = sorted(eligible.items(), key=lambda item: (-item[1][0], item[0]))[:top_n]
    return ReferenceSet(
        genes=[gene for gene, _ in ranked],
        names=[value[1] for _, value in ranked],
        sequences=[cds[gene][0] for gene, _ in ranked],
        stats=stats,
        filters=filters,
    )


def _fetch_pinned(url: str, dest: Path, expected_sha: str, label: str) -> str:
    """Download ``url`` and abort unless its digest matches ``expected_sha``."""
    download(url, dest)
    actual = sha256_file(dest)
    if actual != expected_sha:
        raise SystemExit(
            f"{label}: sha256 {actual} != pinned {expected_sha} ({dest}). Delete the "
            "cached file to re-download, or update the pin deliberately if upstream "
            "re-cut it."
        )
    return actual


def gather(spec: AbundanceSpec, cache_dir: Path, top_n: int) -> tuple[ReferenceSet, dict[str, str]]:
    """Download every pinned source for one organism and select its genes."""
    cds_spec = _cds_spec(spec.key)
    paxdb_path = cache_dir / f"paxdb-{PAXDB_RELEASE}" / spec.paxdb_filename
    # The peptide FASTA and the CDS FASTA come from the SAME Ensembl release, so
    # they share that release's cache namespace. Ensembl reuses filenames across
    # releases, so caching the peptide file under the PaxDb release number would
    # make a warm cache serve the old release's file after an ENSEMBL_RELEASE
    # bump -- the digest pin catches it, but namespacing means it cannot arise.
    ensembl_root = (
        cache_dir / f"{cds_spec.database.replace(' ', '_')}-{cds_spec.release}"
    )
    pep_path = ensembl_root / Path(spec.pep_url).name
    cds_path = ensembl_root / Path(cds_spec.url).name

    digests = {
        "paxdb": _fetch_pinned(spec.paxdb_url, paxdb_path, spec.paxdb_sha256, spec.key),
        "pep": _fetch_pinned(spec.pep_url, pep_path, spec.pep_sha256, spec.key),
        "cds": _fetch_pinned(cds_spec.url, cds_path, cds_spec.source_sha256, spec.key),
    }

    rows = parse_paxdb(paxdb_path)
    index = build_id_index(pep_path)
    cds, filters, organelle_records = representative_cds_by_gene(cds_path)
    selected = select_reference_set(rows, index, cds, filters, organelle_records, top_n)
    return selected, digests


def _cds_spec(key: str) -> OrganismSpec:
    """Return the pinned CDS spec for ``key``.

    Raises:
        SystemExit: If no genome-wide spec exists -- the two reference sets must
            be counted from the same sequences, so a highly-expressed table
            without its genome-wide counterpart would not be comparable.
    """
    for spec in SPECS:
        if spec.key == key:
            return spec
    raise SystemExit(f"{key}: no pinned CDS spec in build_organism_tables.SPECS")


def build_one(spec: AbundanceSpec, cache_dir: Path, out_dir: Path, top_n: int) -> Path:
    """Download, join, select, count, and write one highly-expressed table."""
    print(f"{spec.key}:")
    reference, digests = gather(spec, cache_dir, top_n)
    cds_spec = _cds_spec(spec.key)

    if len(reference.genes) < top_n:
        raise SystemExit(
            f"{spec.key}: only {len(reference.genes)} genes joined and survived "
            f"filtering, fewer than the requested top {top_n} -- refusing to write a "
            "reference set smaller than the one this table claims to be."
        )

    counts = count_codons(reference.sequences)
    missing = sorted(set(CODON_TABLE) - set(counts))
    if missing:
        # Smoothing an unobserved codon would put an invented number into a
        # shipped reference table. Refuse; raise top_n instead (see the module
        # docstring for how N was chosen).
        raise SystemExit(
            f"{spec.key}: codons unobserved in the top {top_n} genes: {missing}. "
            "Raise --top-n rather than smoothing."
        )

    total = sum(counts.values())
    stats = reference.stats
    print(
        f"  {stats.paxdb_rows} abundance rows -> {stats.genes_joined} genes "
        f"({stats.organelle_excluded} organelle-encoded excluded) -> top {top_n}, "
        f"{total} codons counted"
    )
    print(f"  most abundant: {', '.join(reference.names[:8])}")

    written = write_table(
        counts,
        organism=f"{spec.key}{REFERENCE_SET_SUFFIX[HIGHLY_EXPRESSED]}",
        path=out_dir,
        source=(
            f"PaxDb {PAXDB_RELEASE} whole-organism integrated protein abundances "
            f"({spec.taxon}) over {cds_spec.database} release {cds_spec.release} "
            f"{cds_spec.common_name} CDS (assembly {cds_spec.assembly}, annotation "
            f"{cds_spec.genebuild})"
        ),
        cds_count=len(reference.genes),
        build=(
            f"codon occurrence counts over the {top_n} most abundant proteins in "
            "PaxDb's whole-organism integrated dataset, by "
            "scripts/build_highly_expressed_tables.py. PaxDb protein identifiers "
            "were resolved to Ensembl genes against the pinned peptide FASTA of the "
            "same release (identifiers resolving to more than one gene were dropped, "
            "not guessed); each gene took the highest abundance among its proteins; "
            "organelle-encoded genes were excluded (different genetic code and tRNA "
            "pool); genes were ranked by abundance with ties broken by gene ID. The "
            "counted sequence is the same representative CDS the genome-wide table "
            "uses -- longest valid CDS per gene, ACGT-only, in-frame, ATG start, "
            "terminal stop, no internal stop, no alternate-haplotype or patch loci -- "
            "so the two tables differ only in which genes they count."
        ),
        note=(
            "A HIGHLY-EXPRESSED reference set: the relative adaptiveness w = f/f_max "
            "derived from these counts is CAI in Sharp & Li's original sense, where "
            "w = 1 marks the codon translation prefers rather than the codon that is "
            "merely most common genome-wide. The abundances are WHOLE-ORGANISM "
            "consensus values integrated across many published proteomics studies -- "
            "they are not tissue-, cell-type-, or condition-specific, and this table "
            "must not be presented as any of those. It remains a codon-bias proxy, "
            "NOT a measured expression prediction: maximizing CAI against a "
            "highly-expressed reference set has been shown to underperform "
            "alternatives (Welch et al., PLoS ONE 2009, "
            "doi:10.1371/journal.pone.0007002), which is why BT4 treats it as one "
            "axis of an objective vector. Re-derivable from the pinned source URLs "
            "and digests by rerunning scripts/build_highly_expressed_tables.py."
        ),
        reference_set=HIGHLY_EXPRESSED,
        extra={
            "organism": spec.key,
            "organism_common_name": cds_spec.common_name,
            "assembly": cds_spec.assembly,
            "genebuild": cds_spec.genebuild,
            "database": cds_spec.database,
            "database_release": cds_spec.release,
            "source_url": spec.paxdb_url,
            "source_sha256": digests["paxdb"],
            "abundance_source": (
                f"PaxDb {PAXDB_RELEASE} ({spec.taxon}-WHOLE_ORGANISM-integrated), "
                "abundance in ppm"
            ),
            "abundance_identifiers": spec.id_note,
            "protein_to_gene_url": spec.pep_url,
            "protein_to_gene_sha256": digests["pep"],
            "cds_url": cds_spec.url,
            "cds_sha256": digests["cds"],
            "top_n": top_n,
            "gene_roster_sha256": reference.roster_sha256,
            "most_abundant_genes": reference.names[:20],
            "total_codons_counted": total,
            "join": stats.as_dict(),
            # Named for what it is. This tally describes the whole CDS source --
            # it is byte-identical to the genome-wide sidecar's, because both
            # tables read the same FASTA through the same filters -- so its
            # "kept" number is every gene in the annotation with a valid
            # representative CDS, NOT the 300 counted here. Under the
            # genome-wide sidecar's key name (`cds_counted`) that same number
            # IS the count, so reusing the name would have shipped one key
            # meaning two different things across the two sidecar families.
            "cds_source_filters": {
                **{
                    key: value
                    for key, value in reference.filters.as_dict().items()
                    if key != "cds_counted"
                },
                "genes_with_a_representative_cds": reference.filters.kept,
            },
            "rebuild_command": (
                f"python scripts/build_highly_expressed_tables.py {spec.key}"
            ),
            "license_note": (
                "Protein abundances are from PaxDb (von Mering Lab, SIB / University "
                "of Zurich), released under CC BY 4.0; cite Huang et al., Mol Cell "
                "Proteomics 2023, doi:10.1016/j.mcpro.2023.100640. Coding sequences "
                "and the protein-to-gene mapping are from Ensembl, made freely "
                "available by EMBL-EBI (https://www.ebi.ac.uk/about/terms-of-use); "
                "cite Harrison et al., Nucleic Acids Res 2024, "
                "doi:10.1093/nar/gkad1049." + cds_spec.extra_citation
            ),
        },
    )
    return Path(written)


def _verify_against_committed(spec: AbundanceSpec, rebuilt_tsv: Path) -> list[str]:
    """Diff a rebuilt table and its sidecar against the committed ones.

    Every sidecar field is compared except ``retrieved``, a wall-clock stamp of
    when the rebuild ran that is expected to differ.
    """
    problems: list[str] = []
    committed_tsv = DATA_DIR / rebuilt_tsv.name
    if not committed_tsv.is_file():
        return [f"{spec.key}: no committed table"]
    if committed_tsv.read_bytes() != rebuilt_tsv.read_bytes():
        problems.append(f"{spec.key}: committed TSV differs from rebuild")

    name = f"{rebuilt_tsv.name[: -len('.tsv')]}.provenance.json"
    committed_side = DATA_DIR / name
    rebuilt_side = rebuilt_tsv.parent / name
    if not committed_side.is_file():
        return [*problems, f"{spec.key}: no committed provenance sidecar"]
    committed = json.loads(committed_side.read_text(encoding="utf-8"))
    rebuilt = json.loads(rebuilt_side.read_text(encoding="utf-8"))
    for key in sorted((set(committed) | set(rebuilt)) - {"retrieved"}):
        if committed.get(key) != rebuilt.get(key):
            problems.append(f"{spec.key}: provenance field {key!r} differs from rebuild")
    return problems


def report(spec: AbundanceSpec, cache_dir: Path, grid: tuple[int, ...]) -> None:
    """Print the evidence behind ``DEFAULT_TOP_N`` for one organism.

    For each candidate size: how many codons go unobserved (the constraint that
    sets the floor) and how many amino acids' most-used codon differs from the
    genome-wide table (the signal that fades as the size grows).

    The **stop codon is counted separately** from the amino-acid tally, so this
    output matches the per-organism counts quoted in the module docstring. A stop
    is not an amino acid and its preferred codon moves independently (it moves in
    human, mouse and zebrafish and in none of the others), so folding it into the
    same number would make the tool disagree with the prose that cites it.
    """
    from bt4.biomodels.codon.tables import GENOME_WIDE, load_table

    print(f"{spec.key}:")
    reference, _ = gather(spec, cache_dir, max(grid))
    genome_wide = load_table(spec.key, reference_set=GENOME_WIDE)
    gw_argmax = _argmax_by_aa({c: int(v) for c, v in genome_wide.frequency.items()})
    for size in grid:
        counts = count_codons(reference.sequences[:size])
        missing = sorted(set(CODON_TABLE) - set(counts))
        argmax = _argmax_by_aa(counts)
        differs = sorted(
            aa for aa in argmax if aa != STOP and gw_argmax.get(aa) != argmax[aa]
        )
        stop = (
            f" stop {gw_argmax.get(STOP)}->{argmax[STOP]}"
            if argmax.get(STOP) != gw_argmax.get(STOP)
            else ""
        )
        print(
            f"  N={size:5d}  unobserved codons={len(missing):2d}"
            f"{'' if not missing else ' (' + ','.join(missing) + ')'}"
            f"  amino acids differing from genome-wide={len(differs):2d} "
            f"{''.join(differs)}{stop}"
        )


def _argmax_by_aa(counts: dict[str, int]) -> dict[str, str]:
    """Return each amino acid's most-used codon (ties broken by codon)."""
    best: dict[str, tuple[int, str]] = {}
    for codon, value in counts.items():
        aa = CODON_TABLE[codon]
        candidate = (value, codon)
        if aa not in best or candidate > best[aa]:
            best[aa] = candidate
    return {aa: value[1] for aa, value in best.items()}


def main(argv: list[str] | None = None) -> int:
    """Rebuild, verify, or report on the highly-expressed reference tables."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "organisms", nargs="*",
        help="organism keys to build (default: every organism with a pinned "
             "abundance source)",
    )
    parser.add_argument(
        "--cache-dir", default=None,
        help="where to keep downloaded archives (default: a temp directory)",
    )
    parser.add_argument(
        "--verify", action="store_true",
        help="rebuild into a temp directory and diff against the committed tables "
             "instead of overwriting them",
    )
    parser.add_argument(
        "--report", action="store_true",
        help="print the reference-set-size evidence (codon coverage and divergence "
             "from the genome-wide table across candidate sizes) and write nothing",
    )
    parser.add_argument(
        "--top-n", type=int, default=DEFAULT_TOP_N,
        help=f"reference-set size (default: {DEFAULT_TOP_N}; the shipped tables use "
             "the default, so --verify will fail with any other value)",
    )
    args = parser.parse_args(argv)

    if args.top_n < 1:
        parser.error("--top-n must be positive")

    chosen = SPECS_HE
    if args.organisms:
        by_key = {spec.key: spec for spec in SPECS_HE}
        unknown = sorted(set(args.organisms) - set(by_key))
        if unknown:
            parser.error(
                f"no pinned abundance source for: {', '.join(unknown)} "
                f"(available: {', '.join(sorted(by_key))})"
            )
        chosen = tuple(by_key[key] for key in args.organisms)

    with tempfile.TemporaryDirectory() as tmp:
        cache_dir = Path(args.cache_dir) if args.cache_dir else Path(tmp) / "cache"

        if args.report:
            for spec in chosen:
                report(spec, cache_dir, (50, 100, 200, 300, 500, 1000, 2000))
            return 0

        out_dir = Path(tmp) / "out" if args.verify else DATA_DIR
        out_dir.mkdir(parents=True, exist_ok=True)

        mismatched: list[str] = []
        for spec in chosen:
            written = build_one(spec, cache_dir, out_dir, args.top_n)
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
