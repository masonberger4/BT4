"""Shared CSV reading helpers for panel formats that carry whole sequences.

BT4's panel formats put a **sequence** in a single tab-separated field, and Python's
:mod:`csv` refuses any field longer than 131,072 characters.

For a **splice** panel that cap is fatal: a window is a gene span plus 5,000 nt of
flank each side, which routinely exceeds it, so :func:`bt4.api.read_splice_panel`
could not read a single panel ``scripts/make_gencode_splice_panel.py`` produced from
real genomic sequence.

For an **expression** panel it is only a matter of which error you get -- a valid row's
CDS+3'UTR is capped at 11,937 nt by RiboNN's input width, far below csv's limit -- but
an over-long row should be refused by the panel's own check, which names that limit,
not by a stdlib error about a different one.

The limit is **process-global module state**, so raising it permanently from library
code would change parsing behaviour for anything else in the caller's process. The
context manager below raises it only for the duration of one parse and restores the
caller's value afterward, including on an exception.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from contextlib import contextmanager

__all__ = ["MAX_FIELD_SIZE", "relaxed_field_size"]

MAX_FIELD_SIZE = 2**31 - 1
"""The largest field a panel reader will accept, in characters.

:func:`csv.field_size_limit` takes a C ``long``, and the C standard guarantees
``LONG_MAX >= 2147483647`` -- so this exact value is representable everywhere,
including on Windows where ``long`` is 32-bit even in a 64-bit build. Passing
``sys.maxsize`` would raise :class:`OverflowError` there.

Two gigabytes is far above any real panel field -- the longest human gene span is
about 2.3 Mb -- so this bounds the format without constraining it, while still
refusing a runaway field rather than reading forever.
"""


@contextmanager
def relaxed_field_size() -> Iterator[None]:
    """Raise :mod:`csv`'s field-size cap to :data:`MAX_FIELD_SIZE` for one parse.

    Restores the caller's previous limit on the way out, including on an exception,
    so a panel read never leaves the process's ``csv`` configuration changed.
    """
    previous = csv.field_size_limit(MAX_FIELD_SIZE)
    try:
        yield
    finally:
        csv.field_size_limit(previous)
