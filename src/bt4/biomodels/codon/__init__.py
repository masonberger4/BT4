"""Codon-usage tables and the Codon Adaptation Index (CAI).

Re-exports the public API from :mod:`bt4.biomodels.codon.tables`.
"""

from __future__ import annotations

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
    "CodonUsageTable",
    "TableProvenance",
    "available_organisms",
    "load_provenance",
    "load_table",
    "load_table_from_file",
    "sha256_hex",
]
