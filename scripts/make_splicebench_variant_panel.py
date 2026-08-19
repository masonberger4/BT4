"""Convert ``kitzmanlab/splicebench2023`` into a BT4 variant panel.

Step **B1** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md),
the variant-effect half. Smith & Kitzman (*Genome Biol* 24:294, 2023) published 3,616
variants across five genes, measured by MPSA and **already scored by eight tools** --
under an MIT licence. That combination makes it the one panel that can check BT4's gate
machinery *before* any model is installed: run the gate on the benchmark's own
pre-computed columns and see whether it reproduces the published figures.

**The data is not in the repository.** The GitHub repo is four notebooks and a LICENSE;
everything lives in one Zenodo archive whose top directory must be renamed::

    curl -L -o splicebench_data.tar.gz \\
      "https://zenodo.org/records/8351879/files/splicebench_data.tar.gz?download=1"
    # md5 e628ca38209064be73d28d5bddf1ae80  (334,223,475 bytes)
    tar -xzf splicebench_data.tar.gz
    mv for_zenodo data

Labels and scores alone are 11.5 MB:
``tar -xzf splicebench_data.tar.gz for_zenodo/scored_data``.

**The columns this reads**, verified against the archive rather than the paper:

===================  =========================================================
``sdv_fc2``          the label, ``"True"``/``"False"`` (intermediates already
                     dropped upstream)
``exon``             the exonic/intronic stratifier, ``"True"``/``"False"``.
                     This, and nothing else, is what the published 0.419/0.773
                     split uses
``DS_maxm``          SpliceAI, **masked**; ``DS_max`` is unmasked
``pang_max_abs``     Pangolin, **masked**; ``pang_max_nomask_abs`` unmasked
``varlist``          unique key
``gene_name``        absent in the MLH1 file -- the gene comes from the filename
===================  =========================================================

**Over half of this benchmark is not held out, and the panel says so.** The genes sit on
chr17 (BRCA1), chr10 (FAS), chr11 (WT1), chr3 (POU1F1, MST1R, MLH1) -- so **2,077 of
3,616 variants are on chromosomes both SpliceAI and Pangolin trained on**, including
BRCA1, which is otherwise the closest public thing to BT4's synonymous-CDS regime. Each
row carries its chromosome, so :attr:`SpliceVariantPanel.held_out` reports it rather than
leaving a reader to remember. ``--held-out-only`` keeps just the chr3 genes.

Note also that ``mlh1_final_scored.txt`` (296 variants) is a **separate** manually
curated clinical set: 3,616 + 296 = 3,912, the paper's benchmark total. Both numbers are
true about different things, so this script excludes MLH1 unless asked.

Run it::

    python scripts/make_splicebench_variant_panel.py --data data/scored_data --out variants.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Sequence
from pathlib import Path

__all__ = ["GENE_CHROMOSOME", "SCORE_COLUMNS", "SOURCE_FILES", "convert", "main"]

SOURCE_FILES: dict[str, str] = {
    "brca1_findlay_scored.txt": "BRCA1",
    "fas_ex6_snvs_scored.txt": "FAS",
    "pou1f1_snvs_scored.txt": "POU1F1",
    "ron_ex11_scored.txt": "MST1R",
    "wt1_ex9_scored.txt": "WT1",
}
"""The five files whose row counts sum to exactly 3,616. ``mlh1_final_scored.txt`` is a
separate clinical set and is opt-in."""

MLH1_FILE = "mlh1_final_scored.txt"

GENE_CHROMOSOME: dict[str, str] = {
    "BRCA1": "17",
    "FAS": "10",
    "POU1F1": "3",
    "MST1R": "3",
    "WT1": "11",
    "MLH1": "3",
}
"""Each gene's chromosome, so held-out status is checkable.

Only the chr3 genes are held out; chr10, chr11 and chr17 are in both models' training
set. This is the single most consequential fact about this benchmark for BT4's purposes,
which is why it is a table here rather than a sentence in a doc."""

SCORE_COLUMNS: dict[str, str] = {
    "DS_maxm": "spliceai_masked",
    "DS_max": "spliceai_unmasked",
    "pang_max_abs": "pangolin_masked",
    "pang_max_nomask_abs": "pangolin_unmasked",
}
"""Upstream column -> the name BT4's panel will carry.

Masked and unmasked are kept as separate columns rather than one being chosen here: they
answer different questions (masked suppresses scores at annotated sites) and the gate
scores one named column at a time, so the choice belongs at gate time and on the record."""

_LABEL = "sdv_fc2"
_REGION = "exon"
_KEY = "varlist"


def _boolean(cell: str) -> bool | None:
    """Parse the archive's ``"True"``/``"False"`` strings; ``None`` if unparseable."""
    text = cell.strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    return None


def convert(
    directory: Path, *, include_mlh1: bool = False, held_out_only: bool = False
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """Read the scored-data TSVs and return BT4 panel rows plus a per-file tally.

    Args:
        directory: The archive's ``scored_data`` directory.
        include_mlh1: Also read the separate 296-variant clinical set.
        held_out_only: Keep only genes on chromosomes neither model trained on.

    Returns:
        ``(rows, counts)``. ``counts`` maps each gene to how many rows it contributed,
        plus ``"skipped_unparseable"`` for rows whose label or region could not be read.

    Raises:
        FileNotFoundError: If an expected file is absent, naming it -- a partial panel
            silently answers a question about a different dataset.
    """
    wanted = dict(SOURCE_FILES)
    if include_mlh1:
        wanted[MLH1_FILE] = "MLH1"

    rows: list[dict[str, str]] = []
    counts: dict[str, int] = {"skipped_unparseable": 0}
    for filename, gene in sorted(wanted.items()):
        chromosome = GENE_CHROMOSOME[gene]
        if held_out_only and chromosome != "3":
            continue
        path = directory / filename
        if not path.is_file():
            raise FileNotFoundError(
                f"missing {path}. Extract the Zenodo archive and rename its top "
                "directory: tar -xzf splicebench_data.tar.gz && mv for_zenodo data"
            )
        kept = 0
        with path.open("r", encoding="utf-8", newline="") as handle:
            for record in csv.DictReader(handle, delimiter="\t"):
                label = _boolean(record.get(_LABEL, ""))
                exonic = _boolean(record.get(_REGION, ""))
                if label is None or exonic is None:
                    counts["skipped_unparseable"] += 1
                    continue
                row = {
                    "variant_id": (record.get(_KEY) or f"{gene}_{kept}").strip(),
                    "group": gene,
                    "region": "exonic" if exonic else "intronic",
                    "label": "1" if label else "0",
                    "chromosome": chromosome,
                    "note": filename,
                }
                for upstream, name in SCORE_COLUMNS.items():
                    cell = (record.get(upstream) or "").strip()
                    row[name] = cell if cell.lower() not in ("", "na", "nan") else ""
                rows.append(row)
                kept += 1
        counts[gene] = kept
    return rows, counts


def write_panel(rows: Sequence[dict[str, str]], path: Path) -> None:
    """Write rows as the tab-separated format ``read_variant_panel`` reads."""
    columns = ["variant_id", "group", "region", "label", "chromosome", "note"]
    columns += list(SCORE_COLUMNS.values())
    lines = ["\t".join(columns)]
    lines += ["\t".join(row.get(column, "") for column in columns) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    """Convert the splicebench2023 scored data into a BT4 variant panel."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--data", required=True, help="the archive's scored_data directory")
    parser.add_argument("--out", required=True, help="where to write the panel TSV")
    parser.add_argument(
        "--include-mlh1", action="store_true",
        help="also read the separate 296-variant clinical set (not part of the 3,616)",
    )
    parser.add_argument(
        "--held-out-only", action="store_true",
        help="keep only the chr3 genes. The other three (BRCA1, FAS, WT1) are on "
             "chromosomes both models trained on, so a run including them is a "
             "tool-ranking benchmark rather than a held-out measurement",
    )
    args = parser.parse_args(argv)

    rows, counts = convert(
        Path(args.data), include_mlh1=args.include_mlh1, held_out_only=args.held_out_only
    )
    write_panel(rows, Path(args.out))

    print(f"wrote {args.out}  ({len(rows)} variants)")
    for gene in sorted(k for k in counts if k != "skipped_unparseable"):
        chromosome = GENE_CHROMOSOME[gene]
        status = "held out" if chromosome == "3" else "TRAINING chromosome"
        print(f"  {gene:8} chr{chromosome:<3} {counts[gene]:>5}   {status}")
    if counts["skipped_unparseable"]:
        print(f"  skipped {counts['skipped_unparseable']} row(s) with an unreadable "
              f"{_LABEL}/{_REGION}")
    print("\nVerify it before use:")
    print("  python -c \"from bt4.api import read_variant_panel; import json; "
          f"print(json.dumps(read_variant_panel(r'{args.out}', "
          "negative_construction='assayed variants called non-disruptive', "
          "assay='MPSA sdv_fc2 composite').describe(), indent=1))\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
