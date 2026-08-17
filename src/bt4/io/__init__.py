"""BT4 io layer: serialization of results to FASTA and versioned JSON.

Pure stdlib, depends only on ``bt4.domain``. This is the only place that turns a
:class:`~bt4.domain.Result` into external text formats; every serializer is
deterministic so its output can enter the provenance stamp.
"""

from __future__ import annotations

from .fasta import parse_fasta, read_fasta
from .genbank import (
    GenBankFeature,
    GenBankRecord,
    context_from_genbank,
    parse_genbank,
    read_genbank,
    write_genbank,
)
from .serialize import result_to_dict, result_to_json, to_fasta

__all__ = [
    "GenBankFeature",
    "GenBankRecord",
    "context_from_genbank",
    "parse_fasta",
    "parse_genbank",
    "read_fasta",
    "read_genbank",
    "result_to_dict",
    "result_to_json",
    "to_fasta",
    "write_genbank",
]
