"""Cited, license-clean comparison of BT4 vs real commercial codon optimizers.

This standalone report loads a fixed panel of ten published KRas4B coding
sequences (native human CDS plus one variant per commercial/academic optimizer)
from Ranaghan et al. (2021), *BMC Biology* 19:36, DOI 10.1186/s12915-021-00968-8,
reused here under CC BY 4.0. Full attribution, the tool-name mapping (GeneArt =
Thermo GeneOptimizer, DNA2.0 = ATUM, etc.), and the DNA2.0-truncation caveat are
in the sidecar ``scripts/data/kras_ranaghan2021.LICENSE.md``.

Honest framing (read this before reading the numbers):

* **Every value in the table is recomputed by BT4's own functions directly from
  each delivered nucleotide sequence** - CAI and tAI from the bundled
  ``homo_sapiens`` tables, GC% and the longest homopolymer run from
  ``bt4._accel``, and the CpG count by direct dinucleotide counting. No number is
  taken from the source paper or from any tool's own report; nothing is
  fabricated. The tool sequences are real, published output; only the *metrics*
  are BT4's recomputation of them.
* **BT4 is not claimed to be "better" than any tool.** BT4 optimizes its *own*
  objective (a CAI-weighted, GC-aware trellis solve), which is a different target
  than most of these tools pursue and than what actually predicts expression
  (the Ranaghan paper's own finding is that CAI is a poor yield predictor). This
  table only shows *where BT4 lands* on these recomputed axes relative to the
  panel - it is a placement, not a scoreboard.
* **The DNA2.0 (ATUM) sequence encodes a different, C-terminally truncated
  protein** (169 aa vs 188 aa). Its per-codon and compositional metrics are still
  well-defined, but it is not a synonymous variant of the reference, so its row
  is flagged as a length/protein mismatch and must not be read as directly
  comparable to the full-length rows.

Being a standalone script (not part of the ``bt4`` package import graph), it runs
through the stable :mod:`bt4.api` / :mod:`bt4.io` surfaces and reaches into the
``biomodels`` subpackage only for the codon and tRNA tables the public API does
not surface. Run it directly to print a table::

    python scripts/compare_tools.py
    python scripts/compare_tools.py --json
"""

from __future__ import annotations

import argparse
import json
import pathlib
from collections.abc import Mapping, Sequence

from bt4 import api
from bt4._accel import gc_count, max_homopolymer_run
from bt4.biomodels.codon.tables import CodonUsageTable, load_table
from bt4.biomodels.codon.tai import TaiTable, load_tai_table
from bt4.domain.genetic_code import STOP, translate

__all__ = ["DEFAULT_PANEL_PATH", "TOOL_LABELS", "compare", "load_panel", "main"]

# The staged, CC BY 4.0 panel (see the sidecar .LICENSE.md next to it).
DEFAULT_PANEL_PATH = pathlib.Path(__file__).resolve().parent / "data" / "kras_ranaghan2021.fasta"

# FASTA header the native reference CDS is expected under.
_NATIVE_HEADER = "Native"

# Row name used for BT4's own optimized sequence.
_BT4_ROW = "BT4"

# FASTA header -> human-readable optimizer/vendor label. Mirrors the mapping in
# the sidecar license file; GeneArt is Thermo's GeneOptimizer, DNA2.0 is ATUM.
TOOL_LABELS: dict[str, str] = {
    "Native": "native human KRAS CDS (reference)",
    "GeneArt": "Thermo GeneArt (GeneOptimizer)",
    "GeneWiz": "GENEWIZ (Azenta)",
    "DNA2.0": "ATUM (formerly DNA2.0)",
    "IDT": "Integrated DNA Technologies",
    "Genscript": "GenScript (OptimumGene)",
    "Twist": "Twist Bioscience",
    "JCAT": "Java Codon Adaptation Tool",
    "OPTIMIZER": "OPTIMIZER",
    "COOL": "Codon Optimization OnLine",
    _BT4_ROW: "BT4 (this tool)",
}

# Column order for the readable table: (row-dict key, header label).
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("name", "name"),
    ("optimizer", "optimizer"),
    ("len_nt", "len_nt"),
    ("aa", "aa"),
    ("matches_native", "match"),
    ("cai", "cai"),
    ("tai", "tai"),
    ("gc_pct", "gc_pct"),
    ("cpg", "cpg"),
    ("max_homo", "maxhomo"),
    ("certificate", "certificate"),
    ("note", "note"),
)

# Honest header printed above the table (never present a recomputed placement as
# a claim of superiority; flag the truncated row).
_BANNER = (
    "BT4 vs published codon-optimization tools - KRas4B panel\n"
    "Ranaghan et al. (2021) BMC Biology 19:36, DOI 10.1186/s12915-021-00968-8 (CC BY 4.0)\n"
    "\n"
    "All metrics are RECOMPUTED by BT4 from each delivered sequence (CAI/tAI from the\n"
    "bundled homo_sapiens tables, GC%/homopolymer from bt4._accel, CpG by direct count).\n"
    "BT4 optimizes its OWN objective and is NOT claimed 'better' than any tool - this\n"
    "only shows where BT4 lands on these axes. DNA2.0 is a 169-aa C-terminal truncation\n"
    "(a different protein) and is flagged; its row is not directly comparable."
)


def load_panel(path: str | pathlib.Path = DEFAULT_PANEL_PATH) -> list[tuple[str, str]]:
    """Load the KRas panel FASTA into ``(header, sequence)`` records.

    Args:
        path: Path to the panel FASTA; defaults to the staged Ranaghan panel.

    Returns:
        The ``(header, sequence)`` records in file order (via :func:`bt4.api.read_fasta`).

    Raises:
        ValueError: If the FASTA is malformed (see :func:`bt4.io.parse_fasta`).
        OSError: If the file cannot be read.
    """
    return api.read_fasta(path)


def _protein(dna: str) -> str:
    """Return the encoded protein of ``dna`` with any trailing stop removed.

    Args:
        dna: A coding DNA sequence whose length is a multiple of three.

    Returns:
        The single-letter amino-acid string, excluding a trailing stop codon.

    Raises:
        ValueError: If ``dna`` has a bad length or an unknown codon.
    """
    aa = translate(dna)
    return aa[:-1] if aa.endswith(STOP) else aa


def _recompute_metrics(
    dna: str, table: CodonUsageTable, tai_table: TaiTable
) -> dict[str, object]:
    """Recompute every reported metric for ``dna`` from the sequence itself.

    All values are derived here by BT4's own functions; none is read from any
    external source. This is the ``reported == computed`` guarantee applied to
    the whole panel.

    Args:
        dna: The coding sequence to score.
        table: Codon-usage table supplying the CAI weights.
        tai_table: tAI table supplying the tRNA-adaptation weights.

    Returns:
        Mapping with ``cai``, ``tai`` (both rounded to six places), ``gc_pct``
        (rounded to three), ``cpg`` (``CG`` dinucleotide count), and ``max_homo``
        (longest homopolymer run).
    """
    length = len(dna)
    return {
        "cai": round(table.cai(dna), 6),
        "tai": round(tai_table.tai(dna), 6),
        "gc_pct": round(100.0 * gc_count(dna) / length, 3) if length else 0.0,
        "cpg": dna.upper().count("CG"),
        "max_homo": max_homopolymer_run(dna),
    }


def _mismatch_note(aa_len: int, native_aa_len: int) -> str:
    """Return a caveat string when a sequence's protein differs from the native.

    Args:
        aa_len: Amino-acid length of this sequence's protein.
        native_aa_len: Amino-acid length of the native reference protein.

    Returns:
        A flag describing the mismatch, or the empty string when lengths agree
        (a same-length but non-identical protein is reported as a residue
        mismatch, since the panel's optimizers are meant to be synonymous).
    """
    if aa_len != native_aa_len:
        return f"length mismatch: {aa_len} aa vs {native_aa_len} aa (encodes a different protein)"
    return "residue mismatch vs native (not a synonymous variant)"


def compare(
    records: Sequence[tuple[str, str]], config: api.OptimizeConfig | None = None
) -> list[dict[str, object]]:
    """Tabulate recomputed metrics for every panel sequence plus BT4's output.

    Each panel record is translated (flagging any whose protein differs from the
    native reference) and scored by :func:`_recompute_metrics`. BT4 is then run
    on the native protein and scored the same way, so the BT4 row is like-for-like
    with the tools. Nothing is claimed about which is "best" - see the module
    docstring.

    Args:
        records: The panel ``(header, sequence)`` records; one must be the native
            reference (header ``"Native"``), else the first record is used.
        config: BT4 optimization config; defaults to :class:`bt4.api.OptimizeConfig`.
            Its ``organism`` selects the codon and tRNA tables used to recompute
            *every* row's CAI and tAI (the panel is human, so the default
            ``homo_sapiens`` is the right table).

    Returns:
        One row dict per panel sequence in file order, followed by the BT4 row.
        Keys: ``name``, ``optimizer``, ``len_nt``, ``aa``, ``matches_native``,
        ``cai``, ``tai``, ``gc_pct``, ``cpg``, ``max_homo``, ``certificate``
        (``None`` for the published tools, whose optimality BT4 cannot vouch for,
        and the certificate status for the BT4 row), and ``note`` (a mismatch
        flag, empty when the protein matches the native reference).

    Raises:
        ValueError: If ``records`` is empty, or a sequence has a bad length or
            unknown codon, or the organism has no bundled table.
    """
    if not records:
        raise ValueError("panel has no records")
    cfg = config if config is not None else api.OptimizeConfig()
    table = load_table(cfg.organism)
    tai_table = load_tai_table(cfg.organism)

    native = next((seq for header, seq in records if header == _NATIVE_HEADER), records[0][1])
    native_protein = _protein(native)
    native_aa = len(native_protein)

    rows: list[dict[str, object]] = []
    for header, dna in records:
        protein = _protein(dna)
        matches = protein == native_protein
        row: dict[str, object] = {
            "name": header,
            "optimizer": TOOL_LABELS.get(header, header),
            "len_nt": len(dna),
            "aa": len(protein),
            "matches_native": matches,
            "certificate": None,
            "note": "" if matches else _mismatch_note(len(protein), native_aa),
        }
        row.update(_recompute_metrics(dna, table, tai_table))
        rows.append(row)

    result = api.optimize(native_protein, cfg)
    bt4_row: dict[str, object] = {
        "name": _BT4_ROW,
        "optimizer": TOOL_LABELS[_BT4_ROW],
        "len_nt": result.metrics.length_nt,
        "aa": len(_protein(result.dna)),
        "matches_native": _protein(result.dna) == native_protein,
        "certificate": result.certificate.status.value,
        "note": "",
    }
    bt4_row.update(_recompute_metrics(result.dna, table, tai_table))
    rows.append(bt4_row)
    return rows


def _render_cell(value: object) -> str:
    """Format one table cell: floats to three decimals, everything else as text."""
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _format_table(rows: Sequence[Mapping[str, object]]) -> str:
    """Render comparison rows as a fixed-width, human-readable table.

    Args:
        rows: The row dicts produced by :func:`compare`.

    Returns:
        A newline-joined table string (header, rule, one line per row).
    """
    headers = [label for _, label in _COLUMNS]
    cells: list[list[str]] = [[_render_cell(row[key]) for key, _ in _COLUMNS] for row in rows]
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


def main(argv: Sequence[str] | None = None) -> int:
    """Run the comparison over the staged panel and print a table (or JSON).

    Args:
        argv: Optional argument vector (defaults to ``sys.argv[1:]``).

    Returns:
        Process exit code (``0`` on success).
    """
    parser = argparse.ArgumentParser(
        description="Compare BT4 vs published codon-optimization tools on recomputed metrics."
    )
    parser.add_argument(
        "--panel",
        default=str(DEFAULT_PANEL_PATH),
        help="Path to the panel FASTA (default: the staged Ranaghan 2021 KRas panel).",
    )
    parser.add_argument(
        "--organism",
        default="homo_sapiens",
        help="Codon/tRNA table key for BT4 and for recomputing every row (default: homo_sapiens).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the rows as JSON instead of a formatted table.",
    )
    args = parser.parse_args(argv)

    records = load_panel(args.panel)
    config = api.OptimizeConfig(organism=args.organism)
    rows = compare(records, config)

    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        print(_BANNER)
        print()
        print(_format_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
