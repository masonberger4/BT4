"""Codon-usage tables and the Codon Adaptation Index (CAI).

Re-exports the public API from :mod:`bt4.biomodels.codon.tables`.
"""

from __future__ import annotations

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.pairs import CodonPairTable, build_codon_pair_table
from bt4.biomodels.codon.tables import (
    ALIASES,
    CodonUsageTable,
    TableProvenance,
    available_organisms,
    load_provenance,
    load_table,
    load_table_from_file,
    sha256_hex,
)

__all__ = [
    "ALIASES",
    "CodonPairTable",
    "CodonUsageTable",
    "TableProvenance",
    "available_organisms",
    "build_codon_pair_table",
    "build_table",
    "count_codons",
    "load_provenance",
    "load_table",
    "load_table_from_file",
    "sha256_hex",
    "write_table",
]
