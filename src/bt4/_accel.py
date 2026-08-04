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

__all__ = [
    "ACCELERATED",
    "gc_count",
    "longest_repeat",
    "max_gc_run",
    "max_homopolymer_run",
    "reverse_complement",
    "trellis_solve",
]

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


def _py_max_gc_run(seq: str) -> int:
    """Length of the longest run of consecutive ``{G, C}`` bases (mixed allowed).

    ``GCGC`` counts as a run of four -- the "GC length" semantics, distinct from
    a single-base homopolymer. Returns ``0`` for the empty sequence. Byte-for-byte
    identical to the Rust ``bt4_native.max_gc_run``.
    """
    s = validate_dna(seq) if seq else ""
    best = run = 0
    for ch in s:
        if ch in ("G", "C"):
            run += 1
            if run > best:
                best = run
        else:
            run = 0
    return best


def _py_longest_repeat(seq: str) -> int:
    """Length of the longest reverse-complement-aware repeat in ``seq``.

    The largest ``L`` such that some length-``L`` substring occurs at two distinct
    start positions (a *direct* repeat) or some length-``L`` substring's reverse
    complement occurs anywhere in ``seq`` (an *inverted* repeat; a *palindrome*
    when a substring equals its own reverse complement). Returns ``0`` when no
    substring of length ``>= 1`` repeats in that sense.

    Computed as the maximum of two longest-common-substring dynamic programs:
    ``seq`` against itself off the main diagonal (direct repeats, overlaps
    allowed) and ``seq`` against its own reverse complement (inverted repeats).
    Because "offending" is monotone in length, this matches
    :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint`'s notion exactly:
    ``_py_longest_repeat(seq) > max_length`` iff that constraint's ``validate``
    yields a hard violation. Byte-for-byte identical to the Rust
    ``bt4_native.longest_repeat``.
    """
    s = validate_dna(seq) if seq else ""
    n = len(s)
    if n == 0:
        return 0
    rc = reverse_complement(s)
    best = 0
    prev = [0] * n
    cur = [0] * n
    # Direct repeats: longest common substring of ``s`` with itself, excluding the
    # ``i == j`` main diagonal so the two copies sit at distinct positions.
    for i in range(n):
        si = s[i]
        for j in range(n):
            if i != j and si == s[j]:
                v = (prev[j - 1] if j > 0 else 0) + 1
                cur[j] = v
                if v > best:
                    best = v
            else:
                cur[j] = 0
        prev, cur = cur, prev
    # Inverted / palindromic repeats: longest common substring of ``s`` and ``rc``.
    for j in range(n):
        prev[j] = 0
    for i in range(n):
        si = s[i]
        for j in range(n):
            if si == rc[j]:
                v = (prev[j - 1] if j > 0 else 0) + 1
                cur[j] = v
                if v > best:
                    best = v
            else:
                cur[j] = 0
        prev, cur = cur, prev
    return best


def _py_trellis_solve(
    codons: list[str],
    layer_from: list[list[int]],
    layer_to: list[list[int]],
    layer_codon: list[list[int]],
    layer_delta: list[list[float]],
    beam: int | None = None,
) -> tuple[str, float, bool] | None:
    """Pure-Python twin of the Rust ``bt4_native.trellis_solve`` layer DP.

    Runs the exact codon-trellis layer DP over precomputed, position-independent
    transition tables (see the Rust docstring). Context id ``0`` is the empty
    start context. Each layer's four parallel lists give the *allowed*
    transitions ``from -> to`` placing ``codons[codon_id]`` with pre-summed scalar
    ``delta``. The best ``(score, dna)`` per context is kept -- highest scalar,
    ties toward the lexicographically smaller DNA -- exactly as
    :func:`bt4.optimize.exact_dp.solve_exact` does, so the two are byte-for-byte
    identical. Returns ``None`` when a layer has no reachable context (infeasible),
    else ``(best_dna, best_scalar, pruned)`` where ``pruned`` reports whether any
    layer was beam-truncated.
    """
    n_layers = len(layer_from)
    if not (len(layer_to) == len(layer_codon) == len(layer_delta) == n_layers):
        raise ValueError(
            "layer_from/layer_to/layer_codon/layer_delta must have equal length"
        )
    cur: dict[int, tuple[float, str]] = {0: (0.0, "")}
    pruned = False
    for li in range(n_layers):
        froms = layer_from[li]
        tos = layer_to[li]
        cods = layer_codon[li]
        dels = layer_delta[li]
        if not (len(tos) == len(cods) == len(dels) == len(froms)):
            raise ValueError("per-layer transition lists must have equal length")
        nxt: dict[int, tuple[float, str]] = {}
        for t in range(len(froms)):
            state = cur.get(froms[t])
            if state is None:
                continue
            score, dna = state
            ns = score + dels[t]
            nd = dna + codons[cods[t]]
            to = tos[t]
            current = nxt.get(to)
            if current is None or ns > current[0] or (ns == current[0] and nd < current[1]):
                nxt[to] = (ns, nd)
        if not nxt:
            return None
        if beam is not None and len(nxt) > beam:
            kept = sorted(nxt.items(), key=lambda kv: (-kv[1][0], kv[1][1]))[:beam]
            nxt = dict(kept)
            pruned = True
        cur = nxt
    best_score, best_dna = min(cur.values(), key=lambda sv: (-sv[0], sv[1]))
    return best_dna, best_score, pruned


# Public functions are declared with explicit signatures so that binding either
# the Rust extension (typed ``Any`` once the stubless module is ignored) or the
# pure-Python fallbacks stays clean under ``mypy --strict``.
reverse_complement: Callable[[str], str]
gc_count: Callable[[str], int]
max_homopolymer_run: Callable[[str], int]
max_gc_run: Callable[[str], int]
longest_repeat: Callable[[str], int]
trellis_solve: Callable[
    [
        list[str],
        list[list[int]],
        list[list[int]],
        list[list[int]],
        list[list[float]],
        int | None,
    ],
    tuple[str, float, bool] | None,
]

try:  # pragma: no cover - exercised only when the extension is built
    import bt4_native as _native  # type: ignore[import-not-found]

    reverse_complement = _native.reverse_complement
    gc_count = _native.gc_count
    max_homopolymer_run = _native.max_homopolymer_run
    max_gc_run = _native.max_gc_run
    longest_repeat = _native.longest_repeat
    trellis_solve = _native.trellis_solve
    ACCELERATED = True
except ImportError:
    reverse_complement = _py_reverse_complement
    gc_count = _py_gc_count
    max_homopolymer_run = _py_max_homopolymer_run
    max_gc_run = _py_max_gc_run
    longest_repeat = _py_longest_repeat
    trellis_solve = _py_trellis_solve
    ACCELERATED = False
