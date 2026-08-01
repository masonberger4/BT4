"""Data-driven codon-usage tables and the Codon Adaptation Index (CAI).

A :class:`CodonUsageTable` holds per-codon usage frequencies for one organism
and derives the *relative adaptiveness* ``w(codon) = f(codon) / f_max`` where
``f_max`` is the largest frequency among the synonymous codons of the same
amino acid. The CAI of a coding sequence is the geometric mean of ``w`` over its
codons (excluding non-degenerate residues and stops).

Tables are loaded from packaged TSV data (``data/<organism>.tsv``) with a
sidecar ``<organism>.provenance.json``. Bundled values are *representative
published* figures (Kazusa-style), not an authoritative per-genome recount.

This module depends only on :mod:`bt4.domain` and the standard library.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType

from bt4.domain.genetic_code import AMINO_ACIDS, CODON_TABLE, STOP, translate

__all__ = [
    "ALIASES",
    "CodonUsageTable",
    "TableProvenance",
    "available_organisms",
    "load_provenance",
    "load_table",
    "load_table_from_file",
    "sha256_hex",
]

_DATA_PACKAGE = "bt4.biomodels.codon.data"

# Alias -> canonical organism key (canonical == TSV basename without extension).
ALIASES: dict[str, str] = {
    "human": "homo_sapiens",
    "homo sapiens": "homo_sapiens",
    "hsapiens": "homo_sapiens",
    "h_sapiens": "homo_sapiens",
    "hs": "homo_sapiens",
}


def sha256_hex(data: bytes) -> str:
    """Return the lowercase hex SHA-256 digest of ``data``."""
    return hashlib.sha256(data).hexdigest()


@dataclass(frozen=True)
class TableProvenance:
    """Provenance metadata for a codon-usage table.

    Attributes:
        source: Human-readable origin of the counts (e.g. a database name).
        build: How the shipped table was assembled.
        cds_count: Number of coding sequences behind the counts, or ``None`` when
            the table is a representative published summary rather than a recount.
        retrieved: ISO date the values were retrieved/curated.
        sha256: Hex SHA-256 of the raw TSV bytes the provenance describes.
        note: Honesty caveat about how the values should (not) be presented.
    """

    source: str
    build: str
    cds_count: int | None
    retrieved: str
    sha256: str
    note: str


@dataclass(frozen=True, slots=True)
class CodonUsageTable:
    """Per-organism codon usage with derived relative adaptiveness and CAI.

    Attributes:
        organism: Canonical organism key.
        frequency: Read-only mapping codon -> usage frequency (any positive scale).
    """

    organism: str
    frequency: Mapping[str, float]
    _w: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        freq = dict(self.frequency)
        self._validate(freq)
        w = self._compute_relative_adaptiveness(freq)
        object.__setattr__(self, "frequency", MappingProxyType(freq))
        object.__setattr__(self, "_w", MappingProxyType(w))

    @staticmethod
    def _validate(freq: Mapping[str, float]) -> None:
        """Validate coverage: all 20 amino acids plus at least one stop.

        Raises:
            ValueError: If an amino acid (or stop) is uncovered, or a frequency
                is non-positive or a codon is unknown.
        """
        for codon, value in freq.items():
            if codon not in CODON_TABLE:
                raise ValueError(f"unknown codon {codon!r} in codon table")
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(
                    f"frequency for codon {codon!r} must be finite and positive, got {value!r}"
                )
        covered = {CODON_TABLE[codon] for codon in freq}
        missing = sorted(aa for aa in AMINO_ACIDS if aa not in covered)
        if missing:
            raise ValueError(
                "codon-usage table is missing amino acid(s): " + ", ".join(missing)
            )
        if STOP not in covered:
            raise ValueError("codon-usage table is missing a stop codon")

    @staticmethod
    def _compute_relative_adaptiveness(freq: Mapping[str, float]) -> dict[str, float]:
        by_aa: dict[str, float] = {}
        for codon, value in freq.items():
            aa = CODON_TABLE[codon]
            if value > by_aa.get(aa, 0.0):
                by_aa[aa] = value
        return {codon: value / by_aa[CODON_TABLE[codon]] for codon, value in freq.items()}

    def weight(self, codon: str) -> float:
        """Return the relative adaptiveness ``w`` of ``codon``.

        Args:
            codon: A DNA codon (case-insensitive).

        Returns:
            ``w(codon) = f(codon) / f_max`` within its amino acid.

        Raises:
            ValueError: If ``codon`` is not present in the table.
        """
        key = codon.upper()
        try:
            return self._w[key]
        except KeyError:
            raise ValueError(f"codon {codon!r} not in table for {self.organism!r}") from None

    def relative_adaptiveness(self) -> Mapping[str, float]:
        """Return the read-only codon -> relative adaptiveness ``w`` mapping."""
        return self._w

    def cai(self, dna: str) -> float:
        """Compute the Codon Adaptation Index of a coding sequence.

        The CAI is the geometric mean of ``w`` over the sequence's codons,
        counting only codons whose amino acid has more than one synonymous codon
        (i.e. excluding Met ``M``, Trp ``W``, and stop codons).

        Args:
            dna: Coding DNA whose length is a multiple of three.

        Returns:
            The CAI in ``(0, 1]``, or ``1.0`` when no scoreable codon is present.

        Raises:
            ValueError: If ``dna`` has bad length or contains an unknown codon.
        """
        translate(dna)  # validate length + codons via domain semantics.
        seq = dna.upper()
        log_sum = 0.0
        count = 0
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            aa = CODON_TABLE[codon]
            if aa == STOP or aa in ("M", "W"):
                continue
            log_sum += math.log(self._w[codon])
            count += 1
        if count == 0:
            return 1.0
        return math.exp(log_sum / count)


def _parse_tsv(text: str) -> dict[str, float]:
    """Parse the codon-usage TSV text into a codon -> frequency mapping.

    Raises:
        ValueError: On malformed rows or a missing/incorrect header.
    """
    freq: dict[str, float] = {}
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("empty codon-usage TSV")
    header = lines[0].split("\t")
    if header != ["amino_acid", "codon", "frequency"]:
        raise ValueError(f"unexpected TSV header: {header!r}")
    for lineno, line in enumerate(lines[1:], start=2):
        parts = line.split("\t")
        if len(parts) != 3:
            raise ValueError(f"malformed TSV row at line {lineno}: {line!r}")
        _aa, codon, freq_str = parts
        codon = codon.strip().upper()
        try:
            value = float(freq_str)
        except ValueError:
            raise ValueError(f"non-numeric frequency at line {lineno}: {freq_str!r}") from None
        if codon in freq:
            raise ValueError(f"duplicate codon {codon!r} at line {lineno}")
        freq[codon] = value
    return freq


def _canonical(organism: str) -> str:
    key = organism.strip().lower().replace(" ", "_")
    return ALIASES.get(organism.strip().lower(), ALIASES.get(key, key))


def available_organisms() -> tuple[str, ...]:
    """Return the canonical organism keys with a bundled TSV, sorted."""
    data_dir = files(_DATA_PACKAGE)
    names = [
        entry.name[: -len(".tsv")]
        for entry in data_dir.iterdir()
        # Exclude the tRNA tables (``<organism>.trna.tsv``): those are tAI data,
        # not codon-usage tables, and are surfaced by ``available_tai_organisms``.
        if entry.name.endswith(".tsv") and not entry.name.endswith(".trna.tsv")
    ]
    return tuple(sorted(names))


def load_table_from_file(path: str | Path) -> CodonUsageTable:
    """Load a :class:`CodonUsageTable` from a TSV file on disk.

    Args:
        path: Path to a ``.tsv`` file in the standard three-column format. The
            organism name is taken from the file stem.

    Returns:
        The validated table.

    Raises:
        ValueError: If the file is malformed or fails table validation.
    """
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    freq = _parse_tsv(text)
    return CodonUsageTable(organism=p.stem, frequency=freq)


def load_table(organism: str) -> CodonUsageTable:
    """Load a bundled codon-usage table by organism name (alias-aware).

    Args:
        organism: A canonical key (e.g. ``"homo_sapiens"``) or alias (``"human"``).

    Returns:
        The validated table.

    Raises:
        ValueError: If no bundled table matches ``organism``.
    """
    key = _canonical(organism)
    resource = files(_DATA_PACKAGE).joinpath(f"{key}.tsv")
    if not resource.is_file():
        raise ValueError(
            f"unknown organism {organism!r}; available: {', '.join(available_organisms())}"
        )
    freq = _parse_tsv(resource.read_text(encoding="utf-8"))
    return CodonUsageTable(organism=key, frequency=freq)


def load_provenance(organism: str) -> TableProvenance:
    """Load the provenance sidecar for a bundled organism table (alias-aware).

    Args:
        organism: A canonical key or alias.

    Returns:
        The parsed :class:`TableProvenance`.

    Raises:
        ValueError: If no provenance sidecar matches ``organism``.
    """
    key = _canonical(organism)
    resource = files(_DATA_PACKAGE).joinpath(f"{key}.provenance.json")
    if not resource.is_file():
        raise ValueError(f"no provenance for organism {organism!r}")
    data = json.loads(resource.read_text(encoding="utf-8"))
    return TableProvenance(
        source=str(data["source"]),
        build=str(data["build"]),
        cds_count=None if data.get("cds_count") is None else int(data["cds_count"]),
        retrieved=str(data["retrieved"]),
        sha256=str(data["sha256"]),
        note=str(data["note"]),
    )
