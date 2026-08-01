"""Small, dependency-light benchmark: naive back-translation vs BT4.

Runs mostly through the stable :mod:`bt4.api` surface. For each protein in a
panel it builds a naive baseline (the first synonymous codon per residue) and
compares it against the BT4-optimized sequence on cheap, locally recomputed
metrics: GC fraction, longest homopolymer run, CpG dinucleotide count, longest
tandem-repeat span, and the mean sliding-window %MinMax profile, plus the BT4
result's CAI and optimality certificate.

Every number here is recomputed from the delivered DNA of the two sequences BT4
itself produced (the naive straw-man and the optimized answer) — this is a
naive-vs-BT4 report on locally computed metrics only. It deliberately does NOT
compare against published tools (GeneOptimizer / IDT / Twist); that external
comparison is a documented Phase 2 roadmap item (CLAUDE.md §7, §9) and no
competitor numbers are fabricated here.

This is a report, not a test: it never asserts a winner, it just tabulates the
numbers so a human can see the trade BT4 made. Run it directly to print a table::

    python scripts/benchmark.py
    python scripts/benchmark.py --json

Being a standalone script (not part of the ``bt4`` package import graph), it may
reach past :mod:`bt4.api` into the ``biomodels``/``objectives`` subpackages for
the pieces the public API does not surface (codon frequencies and the %MinMax
profile helper).
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence
from typing import cast

from bt4 import api
from bt4._accel import max_homopolymer_run
from bt4.biomodels.codon.tables import load_table
from bt4.domain.genetic_code import STOP, synonymous_codons
from bt4.domain.sequence import gc_fraction
from bt4.objectives import min_max_profile

__all__ = ["DEFAULT_PANEL", "benchmark", "main", "naive_backtranslate"]

# A small, built-in panel spanning short/medium/rare-residue/hydrophobic/long.
DEFAULT_PANEL: dict[str, str] = {
    "short": "MAAL",
    "medium": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIE",
    "rare_residues": "MWCWMHWQWYWMW",
    "hydrophobic": "AVLIFMAVLIFMAVLIF",
    "long": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEKQANGTWPADEFHLMNCYVGR",
}

# Sliding-window size (in codons) for the %MinMax profile mean (Clarke & Clark
# classically use 17-18); matches the objectives layer's reporting default.
_MINMAX_WINDOW = 18

# Longest tandem-repeat span is scanned over unit lengths 1..6 -- mononucleotide
# through hexanucleotide periods, the classic microsatellite/STR range.
_MAX_TANDEM_PERIOD = 6

# Column order for the readable table and the header labels shown per column.
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("name", "name"),
    ("len_nt", "len_nt"),
    ("naive_gc", "naive_gc"),
    ("bt4_gc", "bt4_gc"),
    ("naive_cpg", "naive_cpg"),
    ("bt4_cpg", "bt4_cpg"),
    ("naive_maxhomo", "naive_hp"),
    ("bt4_maxhomo", "bt4_hp"),
    ("naive_maxtandem", "naive_tand"),
    ("bt4_maxtandem", "bt4_tand"),
    ("naive_minmax_mean", "naive_mm"),
    ("bt4_minmax_mean", "bt4_mm"),
    ("bt4_cai", "bt4_cai"),
    ("bt4_certificate", "certificate"),
)


def naive_backtranslate(protein: str) -> str:
    """Back-translate ``protein`` by taking the first synonymous codon per residue.

    This is the deliberately unoptimized straw-man baseline: it ignores codon
    usage entirely and simply picks the lexicographically first codon for each
    amino acid, then appends the first stop codon.

    Args:
        protein: A stop-free single-letter amino-acid string.

    Returns:
        The naive coding DNA, including a trailing stop codon.

    Raises:
        KeyError: If ``protein`` contains a non-amino-acid character.
    """
    codons = [synonymous_codons(residue)[0] for residue in protein]
    codons.append(synonymous_codons(STOP)[0])
    return "".join(codons)


def _cpg_count(dna: str) -> int:
    """Return the number of ``CG`` (CpG) dinucleotides in ``dna``.

    A simple, position-overlapping count of the substring ``"CG"``; deterministic
    and computed directly from the delivered sequence.

    Args:
        dna: The coding sequence to scan.
    """
    return dna.upper().count("CG")


def _longest_tandem_span(dna: str, max_period: int = _MAX_TANDEM_PERIOD) -> int:
    """Return the length of the longest tandem repeat in ``dna``.

    A *tandem repeat of period* ``p`` is a substring whose bases satisfy
    ``seq[k] == seq[k + p]`` throughout (its period is ``p``); it is counted only
    when it spans at least two full copies (span ``>= 2 * p``). This scans every
    period ``p`` in ``1..max_period`` -- mononucleotide (``p == 1``, i.e. a
    homopolymer run) through hexanucleotide by default -- and returns the maximum
    span found, or ``0`` when no length-``>= 2 * p`` tandem repeat exists for any
    period. For a run of ``run`` consecutive period-``p`` matches the covered span
    is ``run + p`` bases (``run / p + 1`` copies).

    Args:
        dna: The coding sequence to scan.
        max_period: Largest repeat-unit length to consider (``>= 1``).

    Returns:
        The longest tandem-repeat span in bases, or ``0`` if none.
    """
    seq = dna.upper()
    n = len(seq)
    best = 0
    for period in range(1, max_period + 1):
        run = 0
        for k in range(n - period):
            if seq[k] == seq[k + period]:
                run += 1
                span = run + period
                if span >= 2 * period and span > best:
                    best = span
            else:
                run = 0
    return best


def _minmax_mean(dna: str, frequencies: Mapping[str, float]) -> float | None:
    """Return the mean sliding-window %MinMax of ``dna`` (``None`` if too short).

    Averages :func:`bt4.objectives.min_max_profile` over the sequence with a
    window of :data:`_MINMAX_WINDOW` codons. When the sequence is shorter than one
    window the profile is empty and there is no meaningful mean, so ``None`` is
    returned (rendered as ``null`` in JSON, keeping ``--json`` valid).

    Args:
        dna: Coding DNA whose length is a multiple of three.
        frequencies: Mapping ``codon -> frequency`` covering every codon of
            ``dna``.

    Returns:
        The mean %MinMax rounded to six places, or ``None`` when the sequence is
        shorter than one window.
    """
    profile = min_max_profile(dna, frequencies, window=_MINMAX_WINDOW)
    if not profile:
        return None
    return round(sum(profile) / len(profile), 6)


def benchmark(
    proteins: Mapping[str, str], config: api.OptimizeConfig | None = None
) -> list[dict[str, object]]:
    """Compare the naive baseline against BT4 over ``proteins``.

    Every metric is recomputed locally from the delivered DNA of the naive
    straw-man and the BT4 answer; nothing is compared against external tools.

    Args:
        proteins: Mapping of panel name to stop-free protein string.
        config: Optimization configuration; defaults to :class:`bt4.api.OptimizeConfig`.
            Its ``organism`` also selects the codon-frequency table used for the
            %MinMax metric.

    Returns:
        One row dict per protein with keys ``name``, ``len_nt``, ``naive_gc``,
        ``bt4_gc``, ``naive_cpg``, ``bt4_cpg``, ``naive_maxhomo``,
        ``bt4_maxhomo``, ``naive_maxtandem``, ``bt4_maxtandem``,
        ``naive_minmax_mean``, ``bt4_minmax_mean``, ``bt4_cai``, and
        ``bt4_certificate``. Floats are rounded to six places; a ``*_minmax_mean``
        is ``None`` when the sequence is shorter than one %MinMax window.
    """
    cfg = config if config is not None else api.OptimizeConfig()
    frequencies = load_table(cfg.organism).frequency
    rows: list[dict[str, object]] = []
    for name, protein in proteins.items():
        naive_dna = naive_backtranslate(protein)
        result = api.optimize(protein, cfg)
        rows.append(
            {
                "name": name,
                "len_nt": result.metrics.length_nt,
                "naive_gc": round(gc_fraction(naive_dna), 6),
                "bt4_gc": round(result.metrics.gc, 6),
                "naive_cpg": _cpg_count(naive_dna),
                "bt4_cpg": _cpg_count(result.dna),
                "naive_maxhomo": max_homopolymer_run(naive_dna),
                "bt4_maxhomo": max_homopolymer_run(result.dna),
                "naive_maxtandem": _longest_tandem_span(naive_dna),
                "bt4_maxtandem": _longest_tandem_span(result.dna),
                "naive_minmax_mean": _minmax_mean(naive_dna, frequencies),
                "bt4_minmax_mean": _minmax_mean(result.dna, frequencies),
                "bt4_cai": round(float(cast("float", result.audit["cai"])), 6),
                "bt4_certificate": result.certificate.status.value,
            }
        )
    return rows


def _format_table(rows: Sequence[Mapping[str, object]]) -> str:
    """Render benchmark rows as a fixed-width, human-readable table.

    Args:
        rows: The row dicts produced by :func:`benchmark`.

    Returns:
        A newline-joined table string (header plus one line per row).
    """
    headers = [label for _, label in _COLUMNS]
    cells: list[list[str]] = [
        [_render_cell(row[key]) for key, _ in _COLUMNS] for row in rows
    ]
    widths = [
        max(len(headers[i]), *(len(line[i]) for line in cells)) if cells else len(headers[i])
        for i in range(len(_COLUMNS))
    ]
    lines = ["  ".join(h.ljust(widths[i]) for i, h in enumerate(headers))]
    lines.append("  ".join("-" * widths[i] for i in range(len(_COLUMNS))))
    lines.extend(
        "  ".join(line[i].ljust(widths[i]) for i in range(len(_COLUMNS))) for line in cells
    )
    return "\n".join(lines)


def _render_cell(value: object) -> str:
    """Format a single cell: floats to three decimals, everything else as text."""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the benchmark over the built-in panel and print a table (or JSON).

    Args:
        argv: Optional argument vector (defaults to ``sys.argv[1:]``).

    Returns:
        Process exit code (``0`` on success).
    """
    parser = argparse.ArgumentParser(description="Benchmark naive back-translation vs BT4.")
    parser.add_argument(
        "--organism",
        default="homo_sapiens",
        help="Codon-usage table key (default: homo_sapiens).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the rows as JSON instead of a formatted table.",
    )
    args = parser.parse_args(argv)

    config = api.OptimizeConfig(organism=args.organism)
    rows = benchmark(DEFAULT_PANEL, config)

    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        print(_format_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
