"""``bt4`` command-line interface -- a thin, print-only shell over :mod:`bt4.api`.

Subcommands:

* ``bt4 optimize PROTEIN`` -- back-translate and optimize a protein.
* ``bt4 validate DNA`` -- audit a coding sequence against the constraints.
* ``bt4 organisms`` -- list bundled codon-usage tables.
* ``bt4 --version`` -- print the single-sourced BT4 version.

Only this module prints; everything else returns data.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from bt4 import __version__, api
from bt4.optimize import InfeasibleError

__all__ = ["main"]


def _build_config(args: argparse.Namespace) -> api.OptimizeConfig:
    motifs = tuple(m.strip().upper() for m in args.forbid) if args.forbid else ()
    enzymes = tuple(e.strip() for e in args.enzyme) if args.enzyme else ()
    return api.OptimizeConfig(
        organism=args.organism,
        gc_target=args.gc_target,
        cai_weight=args.cai_weight,
        gc_weight=args.gc_weight,
        max_homopolymer=None if args.max_homopolymer <= 0 else args.max_homopolymer,
        forbidden_motifs=motifs,
        restriction_enzymes=enzymes,
        ramp_weight=args.ramp_weight,
        ramp_codons=args.ramp_codons,
        cpg_weight=args.cpg_weight,
        cpg_mode=args.cpg_mode,
        minmax_weight=args.minmax_weight,
        minmax_direction=args.minmax_direction,
        tandem_unit=args.tandem_unit,
        tandem_copies=args.tandem_copies,
        inverted_stem=args.inverted_stem,
        inverted_loop=args.inverted_loop,
        avoid_internal_start=args.avoid_internal_start,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        beam=None if args.beam <= 0 else args.beam,
        seed=args.seed,
    )


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--organism", default="homo_sapiens", help="codon-usage table (alias ok)")
    parser.add_argument("--gc-target", type=float, default=0.55, dest="gc_target")
    parser.add_argument("--cai-weight", type=float, default=1.0, dest="cai_weight")
    parser.add_argument("--gc-weight", type=float, default=0.0, dest="gc_weight")
    parser.add_argument("--max-homopolymer", type=int, default=6, dest="max_homopolymer")
    parser.add_argument("--forbid", action="append", metavar="MOTIF", help="forbidden motif")
    parser.add_argument(
        "--enzyme", action="append", metavar="NAME", help="forbid a restriction site (repeatable)"
    )
    parser.add_argument("--ramp-weight", type=float, default=0.0, dest="ramp_weight",
                        help="5' translation-ramp weight (0 = off)")
    parser.add_argument("--ramp-codons", type=int, default=35, dest="ramp_codons")
    parser.add_argument("--cpg-weight", type=float, default=0.0, dest="cpg_weight",
                        help="CpG-dinucleotide term weight (0 = off)")
    parser.add_argument("--cpg-mode", choices=("deplete", "elevate"), default="deplete",
                        dest="cpg_mode")
    parser.add_argument("--minmax-weight", type=float, default=0.0, dest="minmax_weight",
                        help="%%MinMax codon-commonness term weight (0 = off)")
    parser.add_argument("--minmax-direction", choices=("max", "min"), default="max",
                        dest="minmax_direction", help="favour common (max) or rare (min) codons")
    parser.add_argument("--tandem-unit", type=int, default=None, dest="tandem_unit",
                        help="ban tandem repeats of this unit length (off unless set)")
    parser.add_argument("--tandem-copies", type=int, default=3, dest="tandem_copies",
                        help="copies that constitute a banned tandem repeat (default 3)")
    parser.add_argument("--inverted-stem", type=int, default=None, dest="inverted_stem",
                        help="ban hairpins with this stem/arm length (off unless set)")
    parser.add_argument("--inverted-loop", type=int, default=0, dest="inverted_loop",
                        help="max hairpin loop length between arms (default 0)")
    parser.add_argument("--avoid-internal-start", action="store_true",
                        dest="avoid_internal_start",
                        help="forbid internal ATG in a strong Kozak context")
    parser.add_argument("--gc-min", type=int, default=None, dest="gc_min",
                        help="min total GC count (CP-SAT, or Lagrangian with local/pairwise terms)")
    parser.add_argument("--gc-max", type=int, default=None, dest="gc_max",
                        help="max total GC count (CP-SAT, or Lagrangian with local/pairwise terms)")
    parser.add_argument("--beam", type=int, default=0, help="beam width (0 = exact DP)")
    parser.add_argument("--seed", type=int, default=0)


def _cmd_optimize(args: argparse.Namespace) -> int:
    config = _build_config(args)
    result = api.optimize(args.protein, config)
    if args.fasta:
        sys.stdout.write(api.to_fasta(result.dna, header=args.header))
        return 0
    if args.json:
        sys.stdout.write(api.result_to_json(result) + "\n")
        return 0
    cai = float(result.audit.get("cai", 0.0))  # type: ignore[arg-type]
    cert = result.certificate
    print(f"protein   {result.protein}")
    print(f"dna       {result.dna}")
    print(f"length    {result.metrics.length_nt} nt")
    print(f"CAI       {cai:.4f}")
    print(f"GC        {result.metrics.gc * 100:.1f}%")
    print(f"optimality {cert.status.value} ({cert.solver})")
    hard, soft = result.metrics.hard_violations, result.metrics.soft_violations
    print(f"violations {hard} hard / {soft} soft")
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    config = _build_config(args)
    report = api.validate(args.dna, config)
    print(f"dna        {report.dna}")
    print(f"length     {report.metrics.length_nt} nt")
    print(f"GC         {report.metrics.gc * 100:.1f}%")
    print(f"feasible   {report.is_feasible}")
    for v in report.violations:
        print(f"  {v.severity.value:4} {v.constraint} [{v.start}:{v.end}] {v.detail}")
    return 0


def _cmd_organisms(_args: argparse.Namespace) -> int:
    for name in api.available_organisms():
        print(name)
    return 0


def _cmd_enzymes(_args: argparse.Namespace) -> int:
    for name in api.available_enzymes():
        print(name)
    return 0


def _cmd_build_table(args: argparse.Namespace) -> int:
    records = api.read_fasta(args.cds)
    sequences = [seq for _header, seq in records]
    if not sequences:
        print("error: no sequences in the CDS FASTA", file=sys.stderr)
        return 2
    _table, counts = api.build_table(sequences, organism=args.organism)
    path = api.write_table(
        counts,
        organism=args.organism,
        path=args.out,
        source=args.source,
        cds_count=len(sequences),
        pseudocount=args.pseudocount,
    )
    print(f"wrote {path}")
    print(f"organism  {args.organism}")
    print(f"CDS count {len(sequences)}")
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bt4", description="BT4 back-translation optimizer")
    parser.add_argument("--version", action="version", version=f"bt4 {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_opt = sub.add_parser("optimize", help="optimize a protein into coding DNA")
    p_opt.add_argument("protein", help="stop-free amino-acid string")
    p_opt.add_argument("--fasta", action="store_true", help="emit FASTA")
    p_opt.add_argument("--json", action="store_true", help="emit JSON (with manifest)")
    p_opt.add_argument("--header", default="bt4", help="FASTA header")
    _add_common(p_opt)
    p_opt.set_defaults(func=_cmd_optimize)

    p_val = sub.add_parser("validate", help="audit a coding sequence")
    p_val.add_argument("dna", help="ACGT coding sequence")
    _add_common(p_val)
    p_val.set_defaults(func=_cmd_validate)

    p_org = sub.add_parser("organisms", help="list bundled codon-usage tables")
    p_org.set_defaults(func=_cmd_organisms)

    p_enz = sub.add_parser("enzymes", help="list known restriction enzymes")
    p_enz.set_defaults(func=_cmd_enzymes)

    p_bt = sub.add_parser("build-table", help="build a codon table from a CDS FASTA")
    p_bt.add_argument("cds", help="path to a FASTA of coding sequences")
    p_bt.add_argument("--organism", required=True, help="name for the built table")
    p_bt.add_argument("--out", default=".", help="output directory (default: .)")
    p_bt.add_argument("--source", default="user-provided CDS set", help="provenance source label")
    p_bt.add_argument(
        "--pseudocount",
        type=float,
        default=1.0,
        help="Laplace smoothing added to every codon so the table always loads (default: 1.0)",
    )
    p_bt.set_defaults(func=_cmd_build_table)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse ``argv`` and dispatch a subcommand; return a process exit code."""
    args = _parser().parse_args(argv)
    try:
        exit_code: int = args.func(args)
    except (ValueError, InfeasibleError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
