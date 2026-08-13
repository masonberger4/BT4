"""``bt4`` command-line interface -- a thin, print-only shell over :mod:`bt4.api`.

Subcommands:

* ``bt4 optimize PROTEIN`` -- back-translate and optimize a protein.
* ``bt4 library PROTEIN --n N`` -- sample a library from the codon distribution
  (stochastic; SAMPLED certificate, not an optimum).
* ``bt4 validate DNA`` -- audit a coding sequence against the constraints.
* ``bt4 organisms`` -- list bundled codon-usage tables and their reference sets.
* ``bt4 enzymes`` -- list known restriction enzymes.
* ``bt4 build-table CDS.fasta`` -- build a codon table from a CDS FASTA.
* ``bt4 --version`` -- print the single-sourced BT4 version.

Only this module prints; everything else returns data.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence

from bt4 import __version__, api
from bt4.api import InfeasibleError

__all__ = ["main"]

_SPLICE_BACKENDS: tuple[str, ...] = ("assp", "pwm", "pangolin", "spliceai")
"""Backends selectable by ``--check-splice`` / ``--splice-backend``.

``assp`` is the opt-in **online** cross-check (needs the ``bt4[assp]`` extra);
``pwm`` is the offline PWM baseline; ``pangolin`` / ``spliceai`` are the wrapped
CNNs (used only when the user's own install and weights are present).
"""


def _print_splice_crosscheck(dna: str, backend: str) -> None:
    """Run the opt-in splice cross-check on ``dna`` and print it to stderr.

    Written to **stderr** (never stdout) so stdout stays a clean FASTA / JSON
    artifact and the cross-check's numbers -- network-derived for ASSP -- are kept
    out of any captured, reproducible-from-manifest output. The pass is advisory and
    can never fail the run (an unavailable backend degrades gracefully).
    """
    cc = api.splice_crosscheck(dna, backend=backend)
    tags = ["network-derived"] if cc.network_derived else ["local"]
    tags.append("calibrated" if cc.calibrated else "UNCALIBRATED")
    lines = [
        f"--- splice cross-check [{cc.backend}] ---",
        "  " + ", ".join(tags) + "; advisory only, NOT part of the run manifest",
    ]
    if not cc.available:
        lines.append(f"  unavailable: {cc.reason}")
        lines.append("  (an opt-in splice cross-check outage never fails the run)")
    else:
        lines.append(f"  pooled risk {cc.pooled_risk:.3f} (top-{cc.top_k} log-odds; uncalibrated)")
        lines.append(f"  sites       {len(cc.sites)} predicted")
        for site in cc.sites:
            cls = f"  {site.site_class}" if site.site_class else ""
            lines.append(
                f"    {site.kind:8} pos {site.position:>5}  score {site.score:.3f}{cls}"
            )
    print("\n".join(lines), file=sys.stderr)


def _resolve_dinuc_budget(
    args: argparse.Namespace,
) -> tuple[str | None, int | None, int | None]:
    """Map the ``--cpg-*`` / ``--upa-*`` convenience flags to a dinucleotide budget.

    ``--cpg-min``/``--cpg-max`` target the CpG (``CG``) 2-mer and ``--upa-*`` the
    UpA (``TA``) 2-mer. Only one dinucleotide family may be budgeted per run (the
    engine tracks a single count budget at a time), so setting bounds for both is
    a clear error.

    Args:
        args: The parsed CLI namespace (the flags are optional; ``getattr`` keeps
            this usable from subcommands that don't declare them, e.g. validate).

    Returns:
        ``(dinuc_budget, dinuc_min, dinuc_max)`` -- all ``None`` when no flag set.

    Raises:
        ValueError: If both a ``--cpg-*`` and a ``--upa-*`` bound are given.
    """
    cpg_min, cpg_max = getattr(args, "cpg_min", None), getattr(args, "cpg_max", None)
    upa_min, upa_max = getattr(args, "upa_min", None), getattr(args, "upa_max", None)
    cpg_set = cpg_min is not None or cpg_max is not None
    upa_set = upa_min is not None or upa_max is not None
    if cpg_set and upa_set:
        raise ValueError(
            "only one dinucleotide budget at a time: use either --cpg-min/--cpg-max "
            "or --upa-min/--upa-max, not both"
        )
    if cpg_set:
        return "CG", cpg_min, cpg_max
    if upa_set:
        return "TA", upa_min, upa_max
    return None, None, None


def _build_config(args: argparse.Namespace) -> api.OptimizeConfig:
    motifs = tuple(m.strip().upper() for m in args.forbid) if args.forbid else ()
    enzymes = tuple(e.strip() for e in args.enzyme) if args.enzyme else ()
    presets = tuple(p.strip() for p in args.forbid_preset) if args.forbid_preset else ()
    dinuc_budget, dinuc_min, dinuc_max = _resolve_dinuc_budget(args)
    cpb_cds: tuple[str, ...] = ()
    if args.cpb_cds:
        # Read the reference CDS FASTA now (CLI layer); the engine builds the
        # codon-pair table from it. No default table is bundled (§8 honesty).
        cpb_cds = tuple(seq for _header, seq in api.read_fasta(args.cpb_cds))
    return api.OptimizeConfig(
        organism=args.organism,
        reference_set=args.reference_set,
        gc_target=args.gc_target,
        cai_weight=args.cai_weight,
        tai_weight=args.tai_weight,
        gc_weight=args.gc_weight,
        cpb_weight=args.cpb_weight,
        cpb_reference_cds=cpb_cds,
        max_homopolymer=None if args.max_homopolymer <= 0 else args.max_homopolymer,
        max_gc_run=None if args.max_gc_run <= 0 else args.max_gc_run,
        max_repeat_length=None if args.max_repeat_length <= 0 else args.max_repeat_length,
        forbidden_motifs=motifs,
        forbidden_presets=presets,
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
        avoid_splice_sites=args.avoid_splice_sites,
        avoid_internal_start=args.avoid_internal_start,
        avoid_uorf=args.avoid_uorf,
        uorf_region_nt=args.uorf_region_nt,
        refine=args.refine,
        refine_iterations=args.refine_iterations,
        folding_weight=args.folding_weight,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        dinuc_budget=dinuc_budget,
        dinuc_min=dinuc_min,
        dinuc_max=dinuc_max,
        beam=None if args.beam <= 0 else args.beam,
        seed=args.seed,
    )


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--organism", default="homo_sapiens", help="codon-usage table (alias ok)")
    parser.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
        help="which gene set the CAI weights come from: highly_expressed (the "
             "reference set CAI is defined on; the default wherever bundled) or "
             "genome_wide (codon commonness). Omit to use the organism's default; "
             "see `bt4 organisms`",
    )
    parser.add_argument("--gc-target", type=float, default=0.55, dest="gc_target")
    parser.add_argument("--cai-weight", type=float, default=1.0, dest="cai_weight")
    parser.add_argument("--tai-weight", type=float, default=0.0, dest="tai_weight",
                        help="tRNA-adaptation-index weight (0 = off; human tRNA data only)")
    parser.add_argument("--gc-weight", type=float, default=0.0, dest="gc_weight")
    parser.add_argument("--cpb-weight", type=float, default=0.0, dest="cpb_weight",
                        help="codon-pair-bias weight (0 = off; needs --cpb-cds). "
                        "Negative deoptimizes pairs (attenuated-vaccine design)")
    parser.add_argument("--cpb-cds", default=None, dest="cpb_cds", metavar="FASTA",
                        help="reference CDS FASTA to build the codon-pair table from "
                        "(required when --cpb-weight is set; no default is bundled)")
    parser.add_argument("--max-homopolymer", type=int, default=6, dest="max_homopolymer",
                        help="longest allowed single-base run (<=0 = off)")
    parser.add_argument("--max-gc-run", type=int, default=0, dest="max_gc_run",
                        help="longest allowed run of consecutive G/C bases (<=0 = off)")
    parser.add_argument("--max-repeat-length", type=int, default=0, dest="max_repeat_length",
                        help="longest allowed repeated substring anywhere (<=0 = off; "
                        "non-local, refinement-enforced -> HEURISTIC result)")
    parser.add_argument("--forbid", action="append", metavar="MOTIF", help="forbidden motif")
    parser.add_argument(
        "--forbid-preset", action="append", metavar="KEY", dest="forbid_preset",
        help="forbid a named preset's motifs (repeatable; see 'bt4 presets')"
    )
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
    parser.add_argument("--avoid-splice-sites", action="store_true",
                        dest="avoid_splice_sites",
                        help="forbid strong splice-consensus donor/acceptor motifs "
                        "(sense strand; structural heuristic, not a CNN)")
    parser.add_argument("--avoid-internal-start", action="store_true",
                        dest="avoid_internal_start",
                        help="forbid internal ATG in a strong Kozak context")
    parser.add_argument("--avoid-uorf", action="store_true", dest="avoid_uorf",
                        help="suppress out-of-frame internal ATG..in-frame-stop uORFs "
                        "(non-local, refinement-enforced -> HEURISTIC; structural, "
                        "not a calibrated expression claim)")
    parser.add_argument("--uorf-region-nt", type=int, default=100, dest="uorf_region_nt",
                        help="5' scan window (nt) for uORF-opening ATGs (default 100)")
    parser.add_argument("--refine", action="store_true",
                        help="run a 5'-folding-aware SA refinement pass (HEURISTIC result)")
    parser.add_argument("--refine-iterations", type=int, default=2000, dest="refine_iterations",
                        help="SA proposals when --refine is set (default 2000)")
    parser.add_argument("--folding-weight", type=float, default=1.0, dest="folding_weight",
                        help="weight on the 5' folding score during --refine")
    parser.add_argument("--gc-min", type=int, default=None, dest="gc_min",
                        help="min total GC count (CP-SAT, or Lagrangian with local/pairwise terms)")
    parser.add_argument("--gc-max", type=int, default=None, dest="gc_max",
                        help="max total GC count (CP-SAT, or Lagrangian with local/pairwise terms)")
    parser.add_argument("--beam", type=int, default=0, help="beam width (0 = exact DP)")
    parser.add_argument("--seed", type=int, default=0)


def _cmd_optimize(args: argparse.Namespace) -> int:
    config = _build_config(args)
    result = api.optimize(args.protein, config)
    # Opt-in, out-of-loop splice cross-check on the DELIVERED sequence (never in the
    # optimizer). Printed to stderr so it runs in every output mode without touching
    # the stdout FASTA/JSON artifact (its ASSP numbers stay out of the manifest).
    if args.check_splice:
        _print_splice_crosscheck(result.dna, args.check_splice)
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
    # Name the reference set on the same line as the number: a CAI of 1.0 against
    # highly-expressed counts and one against genome-wide counts are different
    # claims, and the caller may not have passed --reference-set at all.
    reference_set = result.audit.get("codon_reference_set", "")
    print(f"CAI       {cai:.4f}" + (f"  (vs {reference_set})" if reference_set else ""))
    if "tai" in result.audit:
        print(f"tAI       {float(result.audit['tai']):.4f}")  # type: ignore[arg-type]
    print(f"GC        {result.metrics.gc * 100:.1f}%")
    for key in ("cg_count", "ta_count"):
        if key in result.audit:
            label = key[:2].upper()  # "CG" (CpG) or "TA" (UpA)
            print(f"{label} count  {result.audit[key]}")
    print(f"optimality {cert.status.value} ({cert.solver})")
    hard, soft = result.metrics.hard_violations, result.metrics.soft_violations
    print(f"violations {hard} hard / {soft} soft")
    if "max_repeat_enforced" in result.audit:
        enforced = result.audit["max_repeat_enforced"]
        limit = result.audit.get("max_repeat_length")
        residual = result.audit.get("max_repeat_residual")
        print(f"max-repeat {enforced} (limit {limit}, residual {residual})")
        if enforced == "partial":
            print("  NOTE: refinement could not remove every long repeat (the protein may "
                  "force some); residual repeats are reported in the violations above.")
    if "uorf_enforced" in result.audit:
        enforced = result.audit["uorf_enforced"]
        region = result.audit.get("uorf_region_nt")
        residual = result.audit.get("uorf_residual")
        print(f"uORF       {enforced} (5' window {region} nt, residual {residual})")
        if enforced == "partial":
            print("  NOTE: refinement could not remove every out-of-frame uORF (the protein "
                  "may force some); residual uORFs are reported in the violations above. "
                  "This is a structural flag, not a calibrated expression prediction.")
    if "folding_model" in result.audit:
        model = result.audit.get("folding_model")
        dg = float(result.audit.get("folding_dg", 0.0))  # type: ignore[arg-type]
        calibrated = bool(result.audit.get("folding_calibrated", False))
        units = "kcal/mol" if calibrated else "arbitrary units (UNCALIBRATED proxy)"
        print(f"folding    5' dG {dg:.3f} {units} [{model}]")
        if not calibrated:
            print("  NOTE: folding is a labeled baseline proxy, not real thermodynamics; "
                  "install bt4[fold] (ViennaRNA) for calibrated dG.")
    return 0


def _cmd_library(args: argparse.Namespace) -> int:
    config = _build_config(args)
    # seed is threaded via _build_config -> config.seed (from --seed); pass None so
    # run_library honors it. This is a stochastic sampler, not an optimizer.
    lib = api.library(args.protein, config, n=args.n, temperature=args.temperature)
    if args.json:
        payload = {
            "protein": args.protein,
            "n": len(lib.results),
            "temperature": args.temperature,
            "seed": args.seed,
            "certificate": "sampled",
            "distinct": lib.distinct,
            "mean_pairwise_hamming": lib.mean_pairwise_hamming,
            "manifest": lib.manifest.to_dict(),
            "sequences": [api.result_to_dict(r) for r in lib.results],
        }
        sys.stdout.write(json.dumps(payload) + "\n")
        return 0
    # Default: one FASTA record per sampled sequence (clean stdout, so a pipe to a
    # FASTA reader still works). Honest framing goes to stderr, never stdout.
    for i, r in enumerate(lib.results, start=1):
        sys.stdout.write(api.to_fasta(r.dna, header=f"{args.header}_{i}"))
    print(
        f"note: {len(lib.results)} sequence(s) SAMPLED from the codon distribution "
        f"(temperature {args.temperature}); distinct {lib.distinct}, mean pairwise "
        f"Hamming {lib.mean_pairwise_hamming:.3f}. These are sampled, NOT optimized, "
        "and are not an expression prediction.",
        file=sys.stderr,
    )
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
    if args.splice_backend:
        _print_splice_crosscheck(args.dna, args.splice_backend)
    return 0


def _cmd_tracks(args: argparse.Namespace) -> int:
    result = api.tracks(
        args.dna,
        args.organism,
        reference_set=args.reference_set,
        nt_window=args.nt_window,
        codon_window=args.codon_window,
    )
    if args.json:
        payload = {
            "dna": result.dna,
            "tracks": [
                {
                    "name": t.name,
                    "window": t.window,
                    "window_unit": t.window_unit,
                    "unit": t.unit,
                    "values": list(t.values),
                }
                for t in result.tracks
            ],
        }
        sys.stdout.write(json.dumps(payload) + "\n")
        return 0
    print(f"dna     {len(result.dna)} nt")
    print(f"{'track':14} {'window':>8} {'n':>5} {'min':>8} {'max':>8} {'mean':>8}")
    for row in api.summarize(result.tracks):
        win = f"{row['window']} {row['window_unit']}"
        mn = "-" if row["min"] is None else f"{float(row['min']):.3f}"  # type: ignore[arg-type]
        mx = "-" if row["max"] is None else f"{float(row['max']):.3f}"  # type: ignore[arg-type]
        mean = "-" if row["mean"] is None else f"{float(row['mean']):.3f}"  # type: ignore[arg-type]
        print(f"{row['name']!s:14} {win:>8} {row['n']:>5} {mn:>8} {mx:>8} {mean:>8}")
    return 0


def _cmd_organisms(_args: argparse.Namespace) -> int:
    print(f"{'organism':26} {'default':17} reference sets")
    for name in api.available_organisms():
        sets = api.available_reference_sets(name)
        print(f"{name:26} {sets[0]:17} {', '.join(sets)}")
    return 0


def _cmd_enzymes(_args: argparse.Namespace) -> int:
    for name in api.available_enzymes():
        print(name)
    return 0


def _cmd_presets(_args: argparse.Namespace) -> int:
    for preset in api.available_forbidden_presets():
        motifs = ", ".join(preset.motifs)
        print(f"{preset.key:26} {preset.label}")
        print(f"{'':26} {preset.description}")
        print(f"{'':26} motifs: {motifs}")
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
    p_opt.add_argument("--cpg-min", type=int, default=None, dest="cpg_min",
                       help="min total CpG (CG) count over the CDS (exact budget DP)")
    p_opt.add_argument("--cpg-max", type=int, default=None, dest="cpg_max",
                       help="max total CpG (CG) count -- CpG depletion for stealth")
    p_opt.add_argument("--upa-min", type=int, default=None, dest="upa_min",
                       help="min total UpA (TA) count over the CDS (exact budget DP)")
    p_opt.add_argument("--upa-max", type=int, default=None, dest="upa_max",
                       help="max total UpA (TA) count over the CDS")
    p_opt.add_argument("--check-splice", choices=_SPLICE_BACKENDS, default=None,
                       dest="check_splice",
                       help="opt-in, out-of-loop splice cross-check of the delivered "
                       "sequence (advisory, printed to stderr; 'assp' is an online "
                       "service and needs the bt4[assp] extra)")
    _add_common(p_opt)
    p_opt.set_defaults(func=_cmd_optimize)

    p_lib = sub.add_parser(
        "library",
        help="sample a library of N sequences from the codon distribution (SAMPLED, "
        "not optimized)",
    )
    p_lib.add_argument("protein", help="stop-free amino-acid string")
    p_lib.add_argument("--n", type=int, required=True, dest="n",
                       help="number of sequences to sample (>= 1)")
    p_lib.add_argument("--temperature", type=float, default=1.0, dest="temperature",
                       help="sampling temperature (>0; ->0 argmax, 1 natural, large uniform)")
    p_lib.add_argument("--json", action="store_true",
                       help="emit JSON (sequences + manifest + diversity)")
    p_lib.add_argument("--header", default="bt4_lib", help="FASTA header prefix")
    _add_common(p_lib)
    p_lib.set_defaults(func=_cmd_library)

    p_val = sub.add_parser("validate", help="audit a coding sequence")
    p_val.add_argument("dna", help="ACGT coding sequence")
    p_val.add_argument("--splice-backend", choices=_SPLICE_BACKENDS, default=None,
                       dest="splice_backend",
                       help="opt-in, out-of-loop splice cross-check of the sequence "
                       "(advisory, printed to stderr; 'assp' is an online service and "
                       "needs the bt4[assp] extra)")
    _add_common(p_val)
    p_val.set_defaults(func=_cmd_validate)

    p_org = sub.add_parser(
        "organisms", help="list bundled codon-usage tables and their reference sets"
    )
    p_org.set_defaults(func=_cmd_organisms)

    p_enz = sub.add_parser("enzymes", help="list known restriction enzymes")
    p_enz.set_defaults(func=_cmd_enzymes)

    p_pre = sub.add_parser("presets", help="list named forbidden-sequence presets")
    p_pre.set_defaults(func=_cmd_presets)

    p_trk = sub.add_parser("tracks", help="per-site composition tracks (GC/CpG/%MinMax)")
    p_trk.add_argument("dna", help="ACGT coding sequence")
    p_trk.add_argument("--organism", default="homo_sapiens", help="table for the %MinMax track")
    p_trk.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
        help="that organism's reference set for the %%MinMax frequencies "
             "(default: the organism's own default)",
    )
    p_trk.add_argument("--nt-window", type=int, default=50, dest="nt_window",
                       help="nucleotide window for the GC/CpG tracks")
    p_trk.add_argument("--codon-window", type=int, default=18, dest="codon_window",
                       help="codon window for the %%MinMax track")
    p_trk.add_argument("--json", action="store_true", help="emit full per-window arrays as JSON")
    p_trk.set_defaults(func=_cmd_tracks)

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
