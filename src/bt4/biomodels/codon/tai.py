"""tRNA adaptation index (tAI) tables from real tRNA gene copy numbers.

The tAI of dos Reis et al. (2004) measures how well a coding sequence's codons
match the cellular tRNA pool. Unlike CAI -- a proxy read off codon *usage* -- tAI
is grounded in the actual **tRNA gene copy numbers** (tGCN) of a genome and a
fixed wobble model, so it is a more mechanistic prior on translational
efficiency (CLAUDE.md section 6).

This module ports the reference construction faithfully:

* Per codon *i*, the absolute adaptiveness is
  ``W_i = sum_j (1 - s_ij) * tGCN_j`` over the tRNAs *j* that read it, where the
  ``s`` are the dos Reis 2004 optimised wobble penalties
  (:data:`DOSREIS_2004_S`) and the WC/wobble decoding follows ``get.ws`` from the
  ``tai`` R package (github.com/mariodosreis/tai). Methionine (ATG) takes only
  its Watson-Crick reader; stops and Met are excluded from the result; the
  bacterial lysidine reading of Ile ``ATA`` applies only when ``sking == 1``.
* The relative adaptiveness is ``w_i = W_i / max(W)``; codons with no reading
  tRNA (``W_i == 0``) are assigned the geometric mean of the non-zero ``w`` so
  they do not annihilate the index.

The bundled human table is built from GtRNAdb (UCSC) *Homo sapiens* hg38 tRNA
gene copy numbers by anticodon; its content SHA-256 enters the provenance stamp
(CLAUDE.md section 8, invariant #9). This module depends only on
:mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType

from bt4._accel import reverse_complement
from bt4.domain.genetic_code import CODON_TABLE, STOP, translate

__all__ = [
    "DOSREIS_2004_S",
    "TaiProvenance",
    "TaiTable",
    "available_tai_organisms",
    "build_tai_weights",
    "load_tai_provenance",
    "load_tai_table",
    "load_tai_table_from_file",
]

_DATA_PACKAGE = "bt4.biomodels.codon.data"

# Standard genetic-code codon order (base1 outer, base2 middle, base3 inner),
# each cycling T,C,A,G. This ordering is what makes the block-of-four wobble
# decoding below (NNT, NNC, NNA, NNG) line up, exactly as in the reference
# get.ws implementation.
_BASES = "TCAG"
_CODONS: tuple[str, ...] = tuple(
    a + b + c for a in _BASES for b in _BASES for c in _BASES
)

# dos Reis et al. (2004) optimised wobble selective constraints (s), in the order
# the get.ws decoding consumes them. p = 1 - s is each pairing's contribution.
#   index 0  I:U  (A34->inosine reading U3)     WC-equivalent   s=0
#   index 1  G:C                                 Watson-Crick    s=0
#   index 2  U:A                                 Watson-Crick    s=0
#   index 3  C:G                                 Watson-Crick    s=0
#   index 4  G:U  wobble                                         s=0.41
#   index 5  I:C  wobble                                         s=0.28
#   index 6  I:A  wobble                                         s=0.9999
#   index 7  U:G  wobble                                         s=0.68
#   index 8  L:A  lysidine (bacterial Ile ATA only, sking=1)     s=0.89
DOSREIS_2004_S: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)


def build_tai_weights(
    anticodon_counts: Mapping[str, int],
    *,
    s_values: tuple[float, ...] = DOSREIS_2004_S,
    sking: int = 0,
) -> dict[str, float]:
    """Compute relative-adaptiveness ``w`` per codon from tRNA copy numbers.

    A faithful port of ``get.ws`` (tai R package): builds the absolute
    adaptiveness ``W`` for all 64 codons under the dos Reis wobble model, then
    normalizes to ``w = W / max(W)`` and fills zero-``W`` codons with the
    geometric mean of the non-zero ``w``.

    Args:
        anticodon_counts: Mapping ``anticodon -> tRNA gene copy number`` (DNA
            alphabet, uppercase). The reader of codon ``c`` at the Watson-Crick
            wobble position is the tRNA whose anticodon is ``reverse_complement(c)``.
        s_values: The nine wobble ``s`` penalties in :data:`DOSREIS_2004_S` order.
        sking: Super-kingdom -- ``0`` eukaryota (default), ``1`` prokaryota. Only
            ``1`` enables the bacterial lysidine reading of Ile ``ATA``.

    Returns:
        Mapping ``codon -> w`` for the 60 scoreable codons (the 61 sense codons
        minus Met, which -- with the stops -- carries no synonymous choice and is
        excluded exactly as the reference does). Every ``w`` lies in ``(0, 1]``.

    Raises:
        ValueError: If ``s_values`` does not have nine entries.
    """
    if len(s_values) != 9:
        raise ValueError(f"s_values must have 9 entries, got {len(s_values)}")
    p = [1.0 - s for s in s_values]
    # Codon-indexed tRNA vector: t[i] = copy number of the tRNA WC-reading _CODONS[i].
    t = [float(anticodon_counts.get(reverse_complement(codon), 0)) for codon in _CODONS]

    w_abs = [0.0] * 64
    for i in range(0, 61, 4):  # each block: NNT, NNC, NNA, NNG
        w_abs[i] = p[0] * t[i] + p[4] * t[i + 1]        # NNT: I:U (WC-equiv) + G:U
        w_abs[i + 1] = p[1] * t[i + 1] + p[5] * t[i]    # NNC: G:C + I:C
        w_abs[i + 2] = p[2] * t[i + 2] + p[6] * t[i]    # NNA: U:A + I:A
        w_abs[i + 3] = p[3] * t[i + 3] + p[7] * t[i + 2]  # NNG: C:G + U:G

    i_met = _CODONS.index("ATG")
    w_abs[i_met] = p[3] * t[i_met]  # Met: only its Watson-Crick reader, no Ile wobble
    if sking == 1:
        w_abs[_CODONS.index("ATA")] = p[8]  # bacterial lysidine reading of Ile ATA

    # Keep the 60 scoreable codons: drop stops and Met (no synonymous choice).
    kept = [
        codon
        for codon in _CODONS
        if CODON_TABLE[codon] != STOP and codon != "ATG"
    ]
    w_kept = {codon: w_abs[_CODONS.index(codon)] for codon in kept}
    w_max = max(w_kept.values())
    if w_max <= 0.0:
        raise ValueError("no tRNA reads any codon; cannot build a tAI table")
    weights = {codon: value / w_max for codon, value in w_kept.items()}

    nonzero = [v for v in weights.values() if v > 0.0]
    geo_mean = math.exp(sum(math.log(v) for v in nonzero) / len(nonzero))
    return {codon: (v if v > 0.0 else geo_mean) for codon, v in weights.items()}


@dataclass(frozen=True)
class TaiProvenance:
    """Provenance metadata for a bundled tAI (tRNA copy-number) table.

    Attributes:
        source: Human-readable origin of the tRNA gene counts (e.g. GtRNAdb).
        build: How the counts were assembled (genome build, filter, exclusions).
        genome: Genome assembly the counts are drawn from.
        total_genes: Total tRNA genes behind the counts.
        retrieved: ISO date the values were retrieved/curated.
        sha256: Hex SHA-256 of the raw anticodon-count TSV bytes.
        note: Honesty caveat about how the values should (not) be presented.
    """

    source: str
    build: str
    genome: str
    total_genes: int
    retrieved: str
    sha256: str
    note: str


@dataclass(frozen=True, slots=True)
class TaiTable:
    """A per-organism tAI table: tRNA copy numbers plus derived ``w`` values.

    Attributes:
        organism: Canonical organism key.
        anticodon_counts: Read-only mapping ``anticodon -> tRNA gene copy number``.
        sking: Super-kingdom used to build ``w`` (0 eukaryota, 1 prokaryota).
    """

    organism: str
    anticodon_counts: Mapping[str, int]
    sking: int = 0
    _w: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        counts = dict(self.anticodon_counts)
        weights = build_tai_weights(counts, sking=self.sking)
        object.__setattr__(self, "anticodon_counts", MappingProxyType(counts))
        object.__setattr__(self, "_w", MappingProxyType(weights))

    def relative_adaptiveness(self) -> Mapping[str, float]:
        """Return the read-only ``codon -> w`` mapping (60 scoreable codons)."""
        return self._w

    def weight(self, codon: str) -> float:
        """Return the relative adaptiveness ``w`` of ``codon``.

        Raises:
            ValueError: If ``codon`` is a stop or Met (excluded from tAI) or not
                a known codon.
        """
        try:
            return self._w[codon.upper()]
        except KeyError:
            raise ValueError(f"codon {codon!r} has no tAI weight (stop/Met or unknown)") from None

    def tai(self, dna: str) -> float:
        """Compute the tRNA adaptation index of a coding sequence.

        The geometric mean of ``w`` over the sequence's scoreable codons (all
        codons except Met and stops, mirroring the reference), or ``1.0`` when no
        scoreable codon is present.

        Raises:
            ValueError: If ``dna`` has a bad length or an unknown codon.
        """
        translate(dna)  # validate length + codons via domain semantics
        seq = dna.upper()
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if codon in self._w:
                log_sum += math.log(self._w[codon])
                count += 1
        if count == 0:
            return 1.0
        return math.exp(log_sum / count)


def _parse_trna_tsv(text: str) -> dict[str, int]:
    """Parse an ``anticodon<TAB>count`` TSV into a mapping.

    Raises:
        ValueError: On a missing/incorrect header or a malformed row.
    """
    counts: dict[str, int] = {}
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or lines[0].split("\t") != ["anticodon", "count"]:
        raise ValueError("tRNA TSV must start with an 'anticodon<TAB>count' header")
    for lineno, line in enumerate(lines[1:], start=2):
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError(f"malformed tRNA TSV row at line {lineno}: {line!r}")
        anticodon = parts[0].strip().upper()
        try:
            counts[anticodon] = int(parts[1])
        except ValueError:
            raise ValueError(f"non-integer count at line {lineno}: {parts[1]!r}") from None
    return counts


def available_tai_organisms() -> tuple[str, ...]:
    """Return the canonical organism keys with a bundled ``*.trna.tsv``, sorted."""
    data_dir = files(_DATA_PACKAGE)
    names = [
        entry.name[: -len(".trna.tsv")]
        for entry in data_dir.iterdir()
        if entry.name.endswith(".trna.tsv")
    ]
    return tuple(sorted(names))


def load_tai_table_from_file(path: str | Path, *, sking: int = 0) -> TaiTable:
    """Load a :class:`TaiTable` from an ``anticodon<TAB>count`` TSV on disk."""
    p = Path(path)
    counts = _parse_trna_tsv(p.read_text(encoding="utf-8"))
    name = p.name[: -len(".trna.tsv")] if p.name.endswith(".trna.tsv") else p.stem
    return TaiTable(organism=name, anticodon_counts=counts, sking=sking)


def load_tai_table(organism: str) -> TaiTable:
    """Load a bundled tAI table by organism key.

    Args:
        organism: Canonical key (e.g. ``"homo_sapiens"``).

    Returns:
        The :class:`TaiTable` for that organism.

    Raises:
        ValueError: If no bundled ``<organism>.trna.tsv`` exists.
    """
    key = organism.strip().lower().replace(" ", "_")
    resource = files(_DATA_PACKAGE).joinpath(f"{key}.trna.tsv")
    if not resource.is_file():
        raise ValueError(
            f"no bundled tAI table for {organism!r}; "
            f"available: {', '.join(available_tai_organisms())}"
        )
    counts = _parse_trna_tsv(resource.read_text(encoding="utf-8"))
    return TaiTable(organism=key, anticodon_counts=counts, sking=0)


def load_tai_provenance(organism: str) -> TaiProvenance:
    """Load the provenance sidecar for a bundled tAI table.

    Raises:
        ValueError: If no provenance sidecar matches ``organism``.
    """
    key = organism.strip().lower().replace(" ", "_")
    resource = files(_DATA_PACKAGE).joinpath(f"{key}.trna.provenance.json")
    if not resource.is_file():
        raise ValueError(f"no tAI provenance for organism {organism!r}")
    data = json.loads(resource.read_text(encoding="utf-8"))
    return TaiProvenance(
        source=str(data["source"]),
        build=str(data["build"]),
        genome=str(data["genome"]),
        total_genes=int(data["total_genes"]),
        retrieved=str(data["retrieved"]),
        sha256=str(data["sha256"]),
        note=str(data["note"]),
    )
