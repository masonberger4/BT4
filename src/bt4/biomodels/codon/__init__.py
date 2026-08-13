"""Codon-usage tables and the Codon Adaptation Index (CAI).

Re-exports the public API from :mod:`bt4.biomodels.codon.tables`.
"""

from __future__ import annotations

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.pairs import CodonPairTable, build_codon_pair_table
from bt4.biomodels.codon.tables import (
    ALIASES,
    CUSTOM_REFERENCE_SET,
    GENOME_WIDE,
    HIGHLY_EXPRESSED,
    REFERENCE_SET_SUFFIX,
    REFERENCE_SETS,
    CodonUsageTable,
    TableProvenance,
    available_organisms,
    available_reference_sets,
    default_reference_set,
    load_provenance,
    load_table,
    load_table_from_file,
    sha256_hex,
)

__all__ = [
    "ALIASES",
    "CUSTOM_REFERENCE_SET",
    "GENOME_WIDE",
    "HIGHLY_EXPRESSED",
    "REFERENCE_SETS",
    "REFERENCE_SET_SUFFIX",
    "CodonPairTable",
    "CodonUsageTable",
    "TableProvenance",
    "available_organisms",
    "available_reference_sets",
    "build_codon_pair_table",
    "build_table",
    "count_codons",
    "default_reference_set",
    "load_provenance",
    "load_table",
    "load_table_from_file",
    "sha256_hex",
    "write_table",
]
