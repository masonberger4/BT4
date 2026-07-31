"""Small, dependency-light benchmark: naive back-translation vs BT4.

Runs entirely through the stable :mod:`bt4.api` surface. For each protein in a
panel it builds a naive baseline (the first synonymous codon per residue) and
compares it against the BT4-optimized sequence on cheap, recomputed metrics
(GC fraction, longest homopolymer run) plus the BT4 result's CAI and optimality
certificate.

This is a report, not a test: it never asserts a winner, it just tabulates the
numbers so a human can see the trade BT4 made. Run it directly to print a table::

    python scripts/benchmark.py
    python scripts/benchmark.py --json
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence

from bt4 import api
from bt4._accel import max_homopolymer_run
from bt4.domain.genetic_code import STOP, synonymous_codons
from bt4.domain.sequence import gc_fraction

__all__ = ["DEFAULT_PANEL", "benchmark", "main", "naive_backtranslate"]

# A small, built-in panel spanning short/medium/rare-residue/hydrophobic/long.
DEFAULT_PANEL: dict[str, str] = {
    "short": "MAAL",
    "medium": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIE",
    "rare_residues": "MWCWMHWQWYWMW",
    "hydrophobic": "AVLIFMAVLIFMAVLIF",
    "long": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEKQANGTWPADEFHLMNCYVGR",
}

# Column order for the readable table and the header labels shown per column.
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("name", "name"),
    ("len_nt", "len_nt"),
    ("naive_gc", "naive_gc"),
    ("bt4_gc", "bt4_gc"),
    ("naive_maxhomo", "naive_hp"),
    ("bt4_maxhomo", "bt4_hp"),
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


def benchmark(
    proteins: Mapping[str, str], config: api.OptimizeConfig | None = None
) -> list[dict[str, object]]:
    """Compare the naive baseline against BT4 over ``proteins``.

    Args:
        proteins: Mapping of panel name to stop-free protein string.
        config: Optimization configuration; defaults to :class:`bt4.api.OptimizeConfig`.

    Returns:
        One row dict per protein with keys ``name``, ``len_nt``, ``naive_gc``,
        ``bt4_gc``, ``naive_maxhomo``, ``bt4_maxhomo``, ``bt4_cai``, and
        ``bt4_certificate``. Floats are rounded to six places.
    """
    rows: list[dict[str, object]] = []
    for name, protein in proteins.items():
        naive_dna = naive_backtranslate(protein)
        result = api.optimize(protein, config)
        rows.append(
            {
                "name": name,
                "len_nt": result.metrics.length_nt,
                "naive_gc": round(gc_fraction(naive_dna), 6),
                "bt4_gc": round(result.metrics.gc, 6),
                "naive_maxhomo": max_homopolymer_run(naive_dna),
                "bt4_maxhomo": max_homopolymer_run(result.dna),
                "bt4_cai": round(float(result.audit["cai"]), 6),
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
