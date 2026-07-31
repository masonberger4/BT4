"""Serialization for BT4 results: FASTA text and a self-describing JSON schema.

Pure stdlib (``json`` + ``textwrap``), depends only on ``bt4.domain``. Every
serializer here is deterministic: identical inputs produce byte-identical output,
so a run's JSON can enter the provenance stamp and be diffed across runs.

The JSON schema is versioned (``schema_version``) and self-describing — it names
every field it carries, mirrors the recomputed metrics, and preserves the
optimality certificate and violation audit so a serialized result stays honest
about how optimal it is and what (if anything) was relaxed.
"""

from __future__ import annotations

import json
import textwrap
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from bt4.domain import Result

__all__ = ["result_to_dict", "result_to_json", "to_fasta"]

_FASTA_WIDTH = 60


def to_fasta(dna: str, *, header: str = "bt4") -> str:
    """Render a DNA sequence as a single FASTA record.

    Args:
        dna: The nucleotide sequence to serialize (may be empty).
        header: The record identifier placed after ``>`` on the header line.

    Returns:
        A FASTA record: a ``>{header}`` line followed by ``dna`` wrapped at 60
        columns per line, terminated by a trailing newline. An empty ``dna``
        yields just the header line plus a newline.
    """
    lines = [f">{header}"]
    if dna:
        lines.extend(textwrap.wrap(dna, _FASTA_WIDTH))
    return "\n".join(lines) + "\n"


def result_to_dict(result: Result) -> dict[str, object]:
    """Convert a :class:`~bt4.domain.Result` into a JSON-ready dictionary.

    The layout is self-describing and versioned via ``schema_version``. Objective
    terms are emitted in sorted order and metrics are mirrored verbatim, so the
    dictionary is deterministic and independently checkable.

    Args:
        result: The optimization result to serialize.

    Returns:
        A JSON-serializable ``dict`` (assuming ``result.audit`` values are
        themselves JSON-serializable).
    """
    metrics = result.metrics
    objective = {term: metrics.objective.get(term) for term in sorted(metrics.objective.terms())}
    certificate = result.certificate
    return {
        "schema_version": "1",
        "protein": result.protein,
        "dna": result.dna,
        "length_nt": metrics.length_nt,
        "gc": metrics.gc,
        "objective": objective,
        "certificate": {
            "status": certificate.status.value,
            "solver": certificate.solver,
            "gap": certificate.gap,
            "relaxed_terms": list(certificate.relaxed_terms),
            "detail": certificate.detail,
        },
        "violations": [
            {
                "constraint": v.constraint,
                "severity": v.severity.value,
                "start": v.start,
                "end": v.end,
                "detail": v.detail,
            }
            for v in result.violations
        ],
        "metrics": {
            "gc": metrics.gc,
            "length_nt": metrics.length_nt,
            "hard_violations": metrics.hard_violations,
            "soft_violations": metrics.soft_violations,
        },
        "audit": dict(result.audit),
    }


def result_to_json(result: Result, *, indent: int = 2) -> str:
    """Serialize a :class:`~bt4.domain.Result` to a deterministic JSON string.

    Args:
        result: The optimization result to serialize.
        indent: Number of spaces to indent each nesting level.

    Returns:
        A pretty-printed JSON string with keys sorted, so the same result always
        serializes to the identical string.
    """
    return json.dumps(result_to_dict(result), indent=indent, sort_keys=True)
