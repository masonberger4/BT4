"""Data-driven codon-usage tables and the Codon Adaptation Index (CAI).

A :class:`CodonUsageTable` holds per-codon usage frequencies for one organism
and derives the *relative adaptiveness* ``w(codon) = f(codon) / f_max`` where
``f_max`` is the largest frequency among the synonymous codons of the same
amino acid. The CAI of a coding sequence is the geometric mean of ``w`` over its
codons (excluding non-degenerate residues and stops).

Every table is loaded from packaged TSV data with a sidecar
``.provenance.json``, and every one is a **real count from a public,
release-pinned source** -- never a hand-authored or "representative" summary.
Read the sidecar rather than assuming; it is the record of what was counted.

Only ratios within an amino acid's synonymous group are ever used
(``w = f/f_max``), so the absolute scale never matters. Consumers that need a
scale-free reading must normalize the same way; see
:class:`~bt4.objectives.minmax.MinMaxTerm`, which does.

**Reference sets: which genes ``w`` was counted over.** ``w`` only means
something relative to a reference set, so BT4 ships two and every table says
which one it is (:attr:`CodonUsageTable.reference_set`, and the manifest of
every run):

``highly_expressed`` (``data/<organism>.highly_expressed.tsv``)
    Codon counts over the 300 most abundant proteins in PaxDb's whole-organism
    integrated dataset, built by ``scripts/build_highly_expressed_tables.py``.
    This is CAI in **Sharp & Li's (1987) original sense** -- ``w = 1`` marks the
    codon translation *prefers* -- with the reference set derived from measured
    protein abundance instead of a hand-picked gene list. It is the **default**
    wherever it is bundled (eight of the nine organisms).

``genome_wide`` (``data/<organism>.tsv``)
    Codon counts over one representative CDS per gene across the whole
    release-pinned Ensembl CDS set, built by ``scripts/build_organism_tables.py``.
    This marks the codon that is most *common*, which under weak translational
    selection is set by mutation and GC bias rather than by translation. Bundled
    for all nine organisms, and the default only for *A. thaliana*, whose PaxDb
    identifiers cannot be joined to the pinned annotation without an unpinned
    external mapping (so BT4 ships no highly-expressed table for it rather than
    one built on a guess).

The two disagree exactly where translational selection is strong. In *E. coli*
eight amino acids change their most-used codon between the two tables --
``TTT``→``TTC`` (Phe), ``CGC``→``CGT`` (Arg), ``GGC``→``GGT`` (Gly),
``ATT``→``ATC`` (Ile) and four more, all of them recovering the classic *E.
coli* optimal codons -- against eleven in *C. elegans*, five in yeast and two in
human (``AGA``→``CGC`` for Arg and ``AGC``→``TCC`` for Ser, plus the preferred
stop codon moving ``TGA``→``TAA``). Neither table is a *measured expression
prediction*: a highly-expressed reference set makes CAI a better-founded proxy,
not a validated one (CLAUDE.md §10.7).

Requesting a reference set an organism does not have raises rather than quietly
substituting the other one; :func:`available_reference_sets` reports what is
actually bundled.

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
    "CUSTOM_REFERENCE_SET",
    "GENOME_WIDE",
    "HIGHLY_EXPRESSED",
    "REFERENCE_SETS",
    "REFERENCE_SET_SUFFIX",
    "CodonUsageTable",
    "TableProvenance",
    "available_organisms",
    "available_reference_sets",
    "default_reference_set",
    "load_provenance",
    "load_table",
    "load_table_from_file",
    "sha256_hex",
]

_DATA_PACKAGE = "bt4.biomodels.codon.data"

#: Counts over the most abundant proteins -- CAI's original reference set.
HIGHLY_EXPRESSED = "highly_expressed"

#: Counts over every gene in the annotation -- codon *commonness*.
GENOME_WIDE = "genome_wide"

#: Known reference sets, **in preference order**: the first one an organism
#: actually has is its default. Highly-expressed leads because that is the
#: reference set CAI is defined on; genome-wide is the fallback for an organism
#: whose abundance data cannot be joined to the pinned annotation.
REFERENCE_SETS: tuple[str, ...] = (HIGHLY_EXPRESSED, GENOME_WIDE)

#: Reference set -> the bundled file's stem suffix. Genome-wide owns the bare
#: ``<organism>`` stem for backwards compatibility with tables already shipped.
REFERENCE_SET_SUFFIX: Mapping[str, str] = MappingProxyType(
    {HIGHLY_EXPRESSED: ".highly_expressed", GENOME_WIDE: ""}
)

#: What a table built outside the bundled data reports as its reference set:
#: a user's own CDS set, whose composition BT4 cannot characterize.
CUSTOM_REFERENCE_SET = "custom"

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
        reference_set: Which gene set the counts were taken over --
            :data:`HIGHLY_EXPRESSED` or :data:`GENOME_WIDE`. This is what makes
            ``w = f/f_max`` interpretable, so it is carried alongside the digest
            rather than left to the reader to infer from the file name.
    """

    source: str
    build: str
    cds_count: int | None
    retrieved: str
    sha256: str
    note: str
    reference_set: str


@dataclass(frozen=True, slots=True)
class CodonUsageTable:
    """Per-organism codon usage with derived relative adaptiveness and CAI.

    Attributes:
        organism: Canonical organism key.
        frequency: Read-only mapping codon -> usage frequency (any positive scale).
        reference_set: Which gene set the frequencies were counted over --
            :data:`HIGHLY_EXPRESSED`, :data:`GENOME_WIDE`, or
            :data:`CUSTOM_REFERENCE_SET` for a table built from a caller's own
            CDS. ``w = f/f_max`` means "the codon this reference set uses most",
            so a table that did not carry its reference set would be a number
            without a question; consumers report it rather than assume it.
    """

    organism: str
    frequency: Mapping[str, float]
    reference_set: str = CUSTOM_REFERENCE_SET
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
            try:
                w = self._w[codon]
            except KeyError:
                # A codon valid in the genetic code but absent from a sparse table
                # (validation only checks amino-acid coverage). Raise the documented
                # ValueError -- matching weight() -- not a bare KeyError.
                raise ValueError(
                    f"codon {codon!r} has no weight in the table for {self.organism!r}"
                ) from None
            log_sum += math.log(w)
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
    """Return the canonical organism keys with a bundled codon table, sorted.

    A bundled *organism* table is exactly ``<organism>.tsv`` -- one dot, no
    interior qualifier. Every other TSV in the data directory belongs to some
    other axis of the same organism (``<organism>.trna.tsv`` is tAI data, listed
    by ``available_tai_organisms``; ``<organism>.highly_expressed.tsv`` is a
    reference set, listed by :func:`available_reference_sets`), and an organism
    is listed once regardless of how many of those it has.

    The rule is stated as "a stem with no dot" rather than as a list of suffixes
    to exclude, because a suffix blocklist fails *open*: the next data file added
    beside these would silently appear as an organism named
    ``homo_sapiens.something``.
    """
    data_dir = files(_DATA_PACKAGE)
    names = [
        entry.name[: -len(".tsv")]
        for entry in data_dir.iterdir()
        if entry.name.endswith(".tsv") and entry.name.count(".") == 1
    ]
    return tuple(sorted(names))


def _stem(organism: str, reference_set: str) -> str:
    """Return the bundled file stem for one organism/reference-set pair."""
    return f"{organism}{REFERENCE_SET_SUFFIX[reference_set]}"


def _has(organism: str, reference_set: str) -> bool:
    resource = files(_DATA_PACKAGE).joinpath(f"{_stem(organism, reference_set)}.tsv")
    return bool(resource.is_file())


def available_reference_sets(organism: str) -> tuple[str, ...]:
    """Return the reference sets bundled for ``organism``, in preference order.

    Args:
        organism: A canonical key (e.g. ``"homo_sapiens"``) or alias (``"human"``).

    Returns:
        A subset of :data:`REFERENCE_SETS` ordered as that tuple is, so
        ``[0]`` is the organism's default.

    Raises:
        ValueError: If no bundled table matches ``organism`` at all.
    """
    key = _canonical(organism)
    found = tuple(name for name in REFERENCE_SETS if _has(key, name))
    if not found:
        raise ValueError(
            f"unknown organism {organism!r}; available: {', '.join(available_organisms())}"
        )
    return found


def default_reference_set(organism: str) -> str:
    """Return the reference set :func:`load_table` uses for ``organism``.

    Highly-expressed wherever it is bundled -- that is the reference set CAI is
    defined on -- and genome-wide otherwise.

    Raises:
        ValueError: If no bundled table matches ``organism``.
    """
    return available_reference_sets(organism)[0]


def _resolve_reference_set(organism: str, reference_set: str | None) -> str:
    """Validate an explicit reference-set request, or pick the default.

    A request for a reference set the organism does not have is an error, never
    a quiet substitution of the other one: the two answer different questions,
    so silently swapping them would make a run's CAI mean something other than
    what the caller asked for while still looking like a success.
    """
    if reference_set is None:
        return default_reference_set(organism)
    if reference_set not in REFERENCE_SET_SUFFIX:
        raise ValueError(
            f"unknown reference set {reference_set!r}; known: {', '.join(REFERENCE_SETS)}"
        )
    have = available_reference_sets(organism)
    if reference_set not in have:
        raise ValueError(
            f"no {reference_set!r} codon table is bundled for {organism!r}; "
            f"available for it: {', '.join(have)}"
        )
    return reference_set


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
    # A caller's own file: BT4 has no idea which genes it was counted over, and
    # must not imply either bundled reference set.
    return CodonUsageTable(
        organism=p.stem, frequency=freq, reference_set=CUSTOM_REFERENCE_SET
    )


def load_table(organism: str, *, reference_set: str | None = None) -> CodonUsageTable:
    """Load a bundled codon-usage table by organism name (alias-aware).

    Args:
        organism: A canonical key (e.g. ``"homo_sapiens"``) or alias (``"human"``).
        reference_set: :data:`HIGHLY_EXPRESSED` or :data:`GENOME_WIDE`. Defaults
            to the organism's :func:`default_reference_set`.

    Returns:
        The validated table, carrying the reference set it was loaded from.

    Raises:
        ValueError: If no bundled table matches ``organism``, or the requested
            reference set is unknown or not bundled for it.
    """
    key = _canonical(organism)
    resolved = _resolve_reference_set(key, reference_set)
    resource = files(_DATA_PACKAGE).joinpath(f"{_stem(key, resolved)}.tsv")
    freq = _parse_tsv(resource.read_text(encoding="utf-8"))
    return CodonUsageTable(organism=key, frequency=freq, reference_set=resolved)


def load_provenance(organism: str, *, reference_set: str | None = None) -> TableProvenance:
    """Load the provenance sidecar for a bundled organism table (alias-aware).

    Args:
        organism: A canonical key or alias.
        reference_set: Which reference set's sidecar to read. Defaults to the
            organism's :func:`default_reference_set`, matching
            :func:`load_table`.

    Returns:
        The parsed :class:`TableProvenance`.

    Raises:
        ValueError: If no provenance sidecar matches the organism/reference set.
    """
    key = _canonical(organism)
    resolved = _resolve_reference_set(key, reference_set)
    resource = files(_DATA_PACKAGE).joinpath(f"{_stem(key, resolved)}.provenance.json")
    if not resource.is_file():
        raise ValueError(f"no provenance for organism {organism!r} ({resolved})")
    data = json.loads(resource.read_text(encoding="utf-8"))
    stamped = str(data.get("reference_set", resolved))
    if stamped != resolved:
        # The sidecar names a different reference set than the file it sits
        # beside. Something is mis-filed, and every downstream honesty claim
        # ("this run used the highly-expressed table") rests on this label.
        raise ValueError(
            f"provenance for {organism!r} claims reference set {stamped!r} but sits "
            f"beside the {resolved!r} table"
        )
    return TableProvenance(
        source=str(data["source"]),
        build=str(data["build"]),
        cds_count=None if data.get("cds_count") is None else int(data["cds_count"]),
        retrieved=str(data["retrieved"]),
        sha256=str(data["sha256"]),
        note=str(data["note"]),
        reference_set=resolved,
    )
