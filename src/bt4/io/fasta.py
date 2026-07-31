"""Minimal FASTA parsing for the BT4 IO layer.

Parses FASTA text (or a file) into ``(header, sequence)`` records with no
external dependencies. Sequences are concatenated across their (possibly wrapped)
lines, upper-cased, and stripped of all whitespace; headers are the ``>`` line
with the leading marker removed and trimmed. Malformed input - sequence data
before any header, or a record with no sequence - raises ``ValueError`` rather
than being silently dropped.

This module depends only on the standard library.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["parse_fasta", "read_fasta"]


def _record(header: str, chunks: list[str]) -> tuple[str, str]:
    """Assemble one ``(header, sequence)`` record, rejecting an empty sequence.

    Args:
        header: The record header (already stripped of its ``>`` marker).
        chunks: The whitespace-free, upper-cased sequence line fragments.

    Returns:
        The ``(header, sequence)`` pair with ``sequence`` the concatenation of
        ``chunks``.

    Raises:
        ValueError: If the concatenated sequence is empty.
    """
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(f"FASTA record {header!r} has an empty sequence")
    return header, sequence


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``(header, sequence)`` records.

    Each record starts at a ``>`` header line; its sequence is the concatenation
    of the following non-header lines, upper-cased with all whitespace removed.
    Blank lines are ignored. Headers are the ``>`` line without the marker,
    trimmed.

    Args:
        text: The full FASTA document as a single string.

    Returns:
        A list of ``(header, sequence)`` pairs in file order.

    Raises:
        ValueError: If sequence data appears before any header, or if any record
            has an empty sequence.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header is not None:
                records.append(_record(header, chunks))
            header = stripped[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("FASTA sequence data appeared before any header")
            chunks.append("".join(raw.split()).upper())
    if header is not None:
        records.append(_record(header, chunks))
    return records


def read_fasta(path: str | os.PathLike[str]) -> list[tuple[str, str]]:
    """Read a FASTA file from disk and parse it.

    Args:
        path: Filesystem path to a FASTA file (UTF-8 encoded).

    Returns:
        A list of ``(header, sequence)`` pairs in file order.

    Raises:
        ValueError: If the file's contents are malformed (see :func:`parse_fasta`).
        OSError: If the file cannot be read.
    """
    return parse_fasta(Path(path).read_text(encoding="utf-8"))
