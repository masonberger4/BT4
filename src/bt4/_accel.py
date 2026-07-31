"""Optional Rust acceleration with a pure-Python fallback.

The hot-loop primitives are implemented in the ``bt4_native`` Rust extension
(built with maturin/PyO3, see ``rust/bt4_core``). When that extension is not
built, this module provides byte-for-byte identical pure-Python implementations,
so ``import bt4`` and the full test suite work with or without the compiler.

``ACCELERATED`` reports which path is active; a property test pins the two
implementations to identical outputs.
"""

from __future__ import annotations

from collections.abc import Callable

from .domain.sequence import validate_dna

__all__ = ["ACCELERATED", "gc_count", "max_homopolymer_run", "reverse_complement"]

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _py_reverse_complement(seq: str) -> str:
    return validate_dna(seq).translate(_COMPLEMENT)[::-1]


def _py_gc_count(seq: str) -> int:
    return sum(1 for ch in validate_dna(seq) if ch in ("G", "C"))


def _py_max_homopolymer_run(seq: str) -> int:
    s = validate_dna(seq) if seq else ""
    if not s:
        return 0
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


# Public functions are declared with explicit signatures so that binding either
# the Rust extension (typed ``Any`` once the stubless module is ignored) or the
# pure-Python fallbacks stays clean under ``mypy --strict``.
reverse_complement: Callable[[str], str]
gc_count: Callable[[str], int]
max_homopolymer_run: Callable[[str], int]

try:  # pragma: no cover - exercised only when the extension is built
    import bt4_native as _native  # type: ignore[import-not-found]

    reverse_complement = _native.reverse_complement
    gc_count = _native.gc_count
    max_homopolymer_run = _native.max_homopolymer_run
    ACCELERATED = True
except ImportError:
    reverse_complement = _py_reverse_complement
    gc_count = _py_gc_count
    max_homopolymer_run = _py_max_homopolymer_run
    ACCELERATED = False
