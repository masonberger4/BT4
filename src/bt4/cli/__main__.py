"""``bt4`` command-line interface -- a thin, print-only shell over :mod:`bt4.api`.

Subcommands:

* ``bt4 optimize PROTEIN`` -- back-translate and optimize a protein.
* ``bt4 library PROTEIN --n N`` -- sample a library from the codon distribution
  (stochastic; SAMPLED certificate, not an optimum).
* ``bt4 validate DNA`` -- audit a coding sequence against the constraints.
* ``bt4 organisms`` -- list bundled codon-usage tables and their reference sets.
* ``bt4 enzymes`` -- list known restriction enzymes.
* ``bt4 build-table CDS.fasta`` -- build a codon table from a CDS FASTA.
* ``bt4 expression-gate PANEL.tsv`` -- run the expression acceptance gate on a
  measured CDS-variant panel, against the mandatory baselines. Reports; never
  promotes.
* ``bt4 splice-gate PANEL.tsv`` -- run the splice gate on an annotated splice-site
  panel, against four permanent baselines.
* ``bt4 variant-gate PANEL.tsv`` -- run the splice gate on a measured variant panel.
  Needs no model installed: a benchmark's own pre-computed scores are data.
* ``bt4 --version`` -- print the single-sourced BT4 version.

Only this module prints; everything else returns data.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Sequence

from bt4 import __version__, api
from bt4.api import InfeasibleError

__all__ = ["main"]

_SPLICE_BACKENDS: tuple[str, ...] = ("assp", "pwm", "pangolin", "spliceai")
_OFFLINE_SPLICE_BACKENDS: tuple[str, ...] = ("pwm", "pangolin", "spliceai")
"""Backends selectable by ``--check-splice`` / ``--splice-backend``.

``assp`` is the opt-in **online** cross-check (needs the ``bt4[assp]`` extra);
``pwm`` is the offline PWM baseline; ``pangolin`` / ``spliceai`` are the wrapped
CNNs (used only when the user's own install and weights are present).
"""


def _enable_attested_splice(args: argparse.Namespace) -> None:
    """Opt this process into honoring bundled splice attestations, if asked.

    The flag sets the same environment variable the library consults, so one
    switch governs every path (CLI, api, Studio) rather than each threading its
    own parameter through.
    """
    if getattr(args, "use_attested_splice", False):
        os.environ[api.USE_ATTESTED_SPLICE_ENV_VAR] = "1"


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


# Every OptimizeConfig field an application preset may set, mapped to the CLI
# option that overrides it and the argparse dest it lands in. A preset must never
# be able to set something the user cannot then override from the command line, so
# a test asserts this table covers every field the bundled presets touch.
_PRESET_FIELD_TO_FLAG: dict[str, tuple[str, str]] = {
    "gc_window_nt": ("--gc-window", "gc_window"),
    "gc_window_min": ("--gc-window-min", "gc_window_min"),
    "gc_window_max": ("--gc-window-max", "gc_window_max"),
    "max_homopolymer": ("--max-homopolymer", "max_homopolymer"),
    "max_gc_run": ("--max-gc-run", "max_gc_run"),
    "max_repeat_length": ("--max-repeat-length", "max_repeat_length"),
    "cpg_weight": ("--cpg-weight", "cpg_weight"),
    "cpg_mode": ("--cpg-mode", "cpg_mode"),
    "avoid_splice_sites": ("--avoid-splice-sites", "avoid_splice_sites"),
    "avoid_polya": ("--avoid-polya", "avoid_polya"),
    "avoid_uorf": ("--avoid-uorf", "avoid_uorf"),
    "inverted_stem": ("--inverted-stem", "inverted_stem"),
    "inverted_loop": ("--inverted-loop", "inverted_loop"),
    "refine": ("--refine", "refine"),
}


def _apply_preset_to_args(args: argparse.Namespace, argv: Sequence[str] | None) -> None:
    """Fold ``args.preset``'s values into ``args`` without clobbering explicit flags.

    A preset supplies starting values; a flag the user actually typed always wins,
    so a preset is a starting point rather than a cage. (The values cannot simply
    be pre-seeded into the argparse namespace: a subparser parses into a fresh
    namespace and would overwrite them with its own defaults.)

    Raises:
        KeyError: If ``args.preset`` is not a known preset.
    """
    preset = api.resolve_preset(args.preset)
    tokens = {
        token.split("=", 1)[0]
        for token in (sys.argv[1:] if argv is None else argv)
    }
    for field, value in preset.overrides.items():
        flag, dest = _PRESET_FIELD_TO_FLAG[field]
        if flag not in tokens:  # user did not name it -> the preset supplies it
            setattr(args, dest, value)


def _read_flank(value: str | None) -> str:
    """Read a flank given either literally or as a path to a FASTA file."""
    if not value:
        return ""
    candidate = value.strip()
    if set(candidate.upper()) <= set("ACGTN"):
        return candidate
    # Not a bare sequence -> treat it as a FASTA path (a clear error if it is not).
    return "".join(seq for _header, seq in api.read_fasta(candidate))


def _build_context(args: argparse.Namespace) -> api.ConstructContext | None:
    """Build the construct context from --utr5/--utr3, or None when neither is set."""
    upstream = _read_flank(getattr(args, "utr5", None))
    downstream = _read_flank(getattr(args, "utr3", None))
    if not upstream and not downstream:
        return None
    return api.ConstructContext(upstream=upstream, downstream=downstream)


def _build_config(args: argparse.Namespace) -> api.OptimizeConfig:
    motifs = tuple(m.strip().upper() for m in args.forbid) if args.forbid else ()
    enzymes = tuple(e.strip() for e in args.enzyme) if args.enzyme else ()
    extra_sites = (
        tuple(s.strip().upper() for s in args.enzyme_site) if args.enzyme_site else ()
    )
    presets = tuple(p.strip() for p in args.forbid_preset) if args.forbid_preset else ()
    dinuc_budget, dinuc_min, dinuc_max = _resolve_dinuc_budget(args)
    context = _build_context(args)
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
        gc_window_nt=None if args.gc_window <= 0 else args.gc_window,
        gc_window_min=args.gc_window_min,
        gc_window_max=args.gc_window_max,
        max_repeat_length=None if args.max_repeat_length <= 0 else args.max_repeat_length,
        forbidden_motifs=motifs,
        forbidden_presets=presets,
        restriction_enzymes=enzymes,
        restriction_extra_sites=extra_sites,
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
        avoid_polya=args.avoid_polya,
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
        application_preset=getattr(args, "preset", None) or "",
        context=context,
        context_provenance=getattr(args, "context_provenance", "omit"),
    )


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--preset", default=None, metavar="KEY",
        help="start from a named application preset (see 'bt4 presets'). NO preset "
        "is applied by default -- BT4 stays regime-agnostic. A preset only supplies "
        "starting values: any flag you pass explicitly still wins"
    )
    parser.add_argument(
        "--utr5", default=None, metavar="SEQ|FASTA",
        help="known sequence immediately 5' of the CDS (5'UTR / vector backbone), "
        "as literal ACGT or a FASTA path. With it, every LOCAL rule is also checked "
        "ACROSS the junction and uORF pairing can see a leader ATG. Use N for "
        "unknown bases: each flank is truncated at the N nearest the CDS"
    )
    parser.add_argument(
        "--utr3", default=None, metavar="SEQ|FASTA",
        help="known sequence immediately 3' of the CDS (same accepted forms)"
    )
    parser.add_argument(
        "--context-provenance", default="omit", dest="context_provenance",
        choices=("omit", "hash"),
        help="what the run manifest records about --utr5/--utr3: 'omit' (default) "
        "stores only their lengths; 'hash' stores a content hash, which makes the "
        "run fully reproducible from its stamp but fingerprints your backbone. The "
        "context is never transmitted anywhere either way"
    )
    parser.add_argument("--organism", default="homo_sapiens", help="codon-usage table (alias ok)")
    parser.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
        help="which gene set the CAI weights come from: highly_expressed (the "
             "reference set CAI is defined on; the default wherever bundled) or "
             "genome_wide (codon commonness). Omit to use the organism's default; "
             "see `bt4 organisms`",
    )
    parser.add_argument("--gc-target", type=float, default=0.55, dest="gc_target",
                        help="target GC fraction in [0,1] for the GC-proximity objective. "
                        "This is a SOFT objective: on `bt4 optimize` it only steers the "
                        "sequence when --gc-weight > 0 (default 0 = off, so --gc-target "
                        "alone has no effect); the frontier always sweeps it. For a HARD "
                        "GC bound use --gc-min/--gc-max")
    parser.add_argument("--cai-weight", type=float, default=1.0, dest="cai_weight")
    parser.add_argument("--gc-weight", type=float, default=0.0, dest="gc_weight",
                        help="weight on the GC-proximity objective in a single solve "
                        "(0 = off; must be > 0 for --gc-target to steer `bt4 optimize`). "
                        "For a hard GC-count window use --gc-min/--gc-max instead")
    parser.add_argument("--tai-weight", type=float, default=0.0, dest="tai_weight",
                        help="tRNA-adaptation-index weight (0 = off; human tRNA data only)")
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
    parser.add_argument("--gc-window", type=int, default=0, dest="gc_window",
                        help="sliding-window length in nt for the windowed GC bound "
                        "(<=0 = off). The rule synthesis vendors specify: bound GC in "
                        "EVERY window, e.g. --gc-window 50 --gc-window-min 0.25 "
                        "--gc-window-max 0.65. LOCAL, so it stays exact in the trellis")
    parser.add_argument("--gc-window-min", type=float, default=0.0, dest="gc_window_min",
                        help="minimum GC fraction in any --gc-window window (0-1)")
    parser.add_argument("--gc-window-max", type=float, default=1.0, dest="gc_window_max",
                        help="maximum GC fraction in any --gc-window window (0-1)")
    parser.add_argument("--max-repeat-length", type=int, default=0, dest="max_repeat_length",
                        help="longest allowed repeated substring anywhere (<=0 = off; "
                        "non-local, refinement-enforced -> HEURISTIC result)")
    parser.add_argument("--forbid", action="append", metavar="MOTIF", help="forbidden motif")
    parser.add_argument(
        "--forbid-preset", action="append", metavar="KEY", dest="forbid_preset",
        help="forbid a named preset's motifs (repeatable; see 'bt4 presets')"
    )
    parser.add_argument(
        "--enzyme-site", action="append", metavar="SITE", dest="enzyme_site",
        help="forbid a recognition SITE directly, as IUPAC (repeatable; e.g. GANTC). "
        "Use this for an enzyme the bundled catalog lacks -- --forbid takes only "
        "literal ACGT, so a degenerate site can only be given here"
    )
    parser.add_argument(
        "--enzyme", action="append", metavar="NAME", help="forbid a restriction site (repeatable)"
    )
    parser.add_argument("--ramp-weight", type=float, default=0.0, dest="ramp_weight",
                        help="5' shaping-prior weight (0 = off): favours less-adapted "
                        "codons early. A PRIOR, not a mechanism -- the 5' benefit is "
                        "driven by reduced RNA structure, not codon rarity (Goodman "
                        "2013), so use --refine/--folding-weight for that lever")
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
    parser.add_argument("--avoid-polya", action="store_true", dest="avoid_polya",
                        help="forbid FUNCTIONAL poly(A) signals: an AATAAA/ATTAAA "
                        "hexamer paired with a downstream U/GU-rich element (the "
                        "bipartite signal the cleavage machinery recognises). More "
                        "permissive than --forbid-preset poly_a_signal, which bans "
                        "every bare hexamer. Refinement-enforced -> HEURISTIC result")
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
    if args.genbank:
        # The annotated map: residual violations the optimizer could not remove
        # ride along as misc_feature spans, so they are visible in SnapGene /
        # Benchling rather than only in the JSON audit.
        sys.stdout.write(
            api.write_genbank(
                result, context=config.context, locus=args.header,
            )
        )
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
        splice=args.splice_track,
    )
    if args.json:
        payload = {
            "dna": result.dna,
            "organism": result.organism,
            "codon_reference_set": result.reference_set,
            "splice_model": result.splice_model,
            "splice_calibrated": result.splice_calibrated,
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
    if result.reference_set:
        # %MinMax is reference-set-dependent, so the track table is not
        # self-describing without this line.
        print(f"tables  {result.organism} / {result.reference_set}")
    if result.splice_model:
        state = "calibrated" if result.splice_calibrated else "UNCALIBRATED (advisory)"
        print(f"splice  {result.splice_model} [{state}]")
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
    print("Application presets (--preset KEY) -- starting points, none applied by default:")
    print()
    for app in api.available_presets():
        settings = ", ".join(f"{k}={v!r}" for k, v in sorted(app.overrides.items()))
        print(f"{app.key:26} {app.label}  [{app.regime}]")
        print(f"{'':26} {app.description}")
        print(f"{'':26} why: {app.rationale}")
        print(f"{'':26} sets: {settings}")
        print()
    print("Forbidden-sequence presets (--forbid-preset KEY):")
    print()
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


def _cmd_expression_gate(args: argparse.Namespace) -> int:
    """Run the expression acceptance gate over a measured panel and print the verdict."""
    panel = api.read_panel(args.panel)
    comparison = api.expression_gate(
        panel,
        args.backend,
        settings=api.GateSettings(
            within_group=args.within_group,
            recalibrate=args.recalibrate,
            target_coverage=args.target_coverage,
            coverage_tolerance=args.coverage_tolerance,
            min_spearman=args.min_spearman,
            calibration_fraction=args.calibration_fraction,
            bootstrap_resamples=args.bootstrap_resamples,
            seed=args.seed,
        ),
        species=args.species,
        cell_types=tuple(args.cell_types or ()),
        num_workers=args.num_workers,
        organism=args.organism,
        reference_set=args.reference_set,
    )

    if not args.within_group:
        print(
            "warning: pooled mode credits between-protein skill, which is NOT the "
            "regime BT4 deploys in. Pass --within-group for the strict bar.",
            file=sys.stderr,
        )

    summary = panel.describe()
    flag = "calibrated" if comparison.backend_calibrated else "UNCALIBRATED"
    print(f"panel:    {summary['n_rows']} rows / {summary['n_groups']} groups")
    print(f"          sha256 {comparison.panel_hash[:16]}...")
    print(f"backend:  {comparison.backend}  [{flag}]")
    print(f"mode:     {'within-protein' if args.within_group else 'POOLED'}")
    print()
    print(f"{'':<14}{'spearman':>10}{'CI low':>9}{'coverage':>10}{'width/IQR':>11}")
    rows = [("HEAD", comparison.head), *((f"  {n}", r) for n, r in comparison.baselines)]
    for label, report in rows:
        print(
            f"{label:<14}{report.spearman:>10.3f}{report.spearman_ci_low:>9.3f}"
            f"{report.empirical_coverage:>10.3f}{report.width_over_iqr:>11.3f}"
        )
    print()
    print(f"gate passed (thresholds) : {comparison.head.passed}")
    print(
        f"beats every baseline     : {comparison.beats_every_baseline} "
        f"(best: {comparison.best_baseline} at {comparison.best_baseline_spearman:.3f})"
    )
    print(f"interval is informative  : {comparison.interval_is_informative}")
    print(f"PROMOTABLE on this panel : {comparison.promotable}")
    print()
    print(
        "This command flips nothing: 'promotable' means the pre-registered conditions "
        "held on this panel. Promotion is a separate, recorded step, and min_spearman "
        "is a pre-commitment rather than a community standard (none exists)."
    )
    return 0


def _resolve_anchor_offsets(args: argparse.Namespace) -> int | dict[str, int]:
    """Combine the anchor flags into one scalar or per-kind mapping.

    Per-kind wins over the scalar, and ``--cnn-anchors`` is the measured default for the
    wrapped CNNs. A scalar remains valid because a kind-separated panel is a legitimate
    way to run them.
    """
    base = 0 if args.anchor_offset is None else args.anchor_offset
    if args.cnn_anchors:
        offsets = dict(api.CNN_ANCHOR_OFFSETS)
    elif args.donor_offset is None and args.acceptor_offset is None:
        return base
    else:
        offsets = dict.fromkeys(("donor", "acceptor"), base)
    if args.donor_offset is not None:
        offsets["donor"] = args.donor_offset
    if args.acceptor_offset is not None:
        offsets["acceptor"] = args.acceptor_offset
    return offsets


def _splice_gate_progress(index: int, total: int, window_id: str, length: int) -> None:
    """Report which window is being scored, to stderr.

    A wrapped CNN reads ~10 kb of context per position, so a real panel runs for tens of
    minutes with nothing to distinguish it from a hang. This goes to **stderr** so the
    report on stdout stays pipeable, and it names the window's length because that -- not
    the count -- is what the remaining wait is proportional to.
    """
    print(
        f"  scoring {index}/{total}  {window_id}  ({length:,} nt)",
        file=sys.stderr,
        flush=True,
    )


def _cmd_splice_gate(args: argparse.Namespace) -> int:
    """Run the splice acceptance gate over an annotated panel and print the verdict."""
    panel = api.read_splice_panel(
        args.panel,
        negative_construction=args.negative_construction,
        annotation=args.annotation,
        min_motif_consistency=args.min_motif_consistency,
    )
    comparison = api.splice_panel_gate(
        panel,
        args.backend,
        settings=api.SpliceGateSettings(
            threshold=args.threshold,
            min_pr_auc=args.min_pr_auc,
            min_pr_auc_skill=args.min_pr_auc_skill,
            max_ece=args.max_ece,
            n_bins=args.bins,
            seed=args.seed,
        ),
        anchor_offset=_resolve_anchor_offsets(args),
        progress=None if args.quiet else _splice_gate_progress,
    )

    summary = panel.describe()
    flag = "fidelity-attested" if comparison.backend_calibrated else "UNCALIBRATED"
    print(f"panel:    {summary['n_windows']} windows / {summary['n_positions']} positions")
    print(f"          {panel.n_sites} annotated sites, groups {list(panel.groups)}")
    print(f"          motif consistency {summary['motif_consistency']:.1%}")
    print(f"          sha256 {comparison.panel_hash[:16]}...")
    print(f"backend:  {comparison.backend}  [{flag}]")
    for note in comparison.notes:
        print(f"note:     {note}")
    print()

    header = f"{'':<14}{'stratum':<10}{'AP':>8}{'skill':>8}{'ROC':>8}{'top-k':>8}{'ECE':>8}"
    print(header)
    rows = [("HEAD", comparison.head), *((f"  {n}", r) for n, r in comparison.baselines)]
    for label, report in rows:
        for stratum in report.strata:
            print(
                f"{label:<14}{stratum.name:<10}{stratum.pr_auc:>8.3f}"
                f"{stratum.pr_auc_skill:>8.3f}{stratum.roc_auc:>8.3f}"
                f"{stratum.top_k_accuracy:>8.3f}{stratum.ece:>8.3f}"
            )
    print()
    for name, baseline, skill in comparison.best_baseline:
        print(f"best baseline ({name}): {baseline} at skill {skill:.3f}")
    print()
    print(f"gate passed (thresholds) : {comparison.head.passed}")
    for reason in comparison.head.reasons:
        print(f"  - {reason}")
    print(f"beats every baseline     : {comparison.beats_every_baseline}")
    print(f"panel is held out        : {comparison.held_out}")
    print(f"a bar was declared       : {comparison.thresholds_declared}")
    print(f"PROMOTABLE on this panel : {comparison.promotable}")
    print()
    print(
        "This command flips nothing, and it answers a DIFFERENT question from the "
        "fidelity attestation: that one proves BT4's wrapper reproduces the published "
        "model, this one asks whether the numbers mean what they claim. A backend needs "
        "both, and the thresholds here are pre-commitments set at gate time.",
        file=sys.stderr,
    )
    return 0


def _cmd_variant_gate(args: argparse.Namespace) -> int:
    """Run the splice gate over a measured variant panel and print the verdict."""
    panel = api.read_variant_panel(
        args.panel,
        negative_construction=args.negative_construction,
        assay=args.assay,
    )
    summary = panel.describe()

    if args.list_scores or not args.score:
        print(f"panel:    {summary['n_rows']} variants / {len(panel.groups)} genes")
        print(f"          {summary['n_positive']} positive, groups {list(panel.groups)}")
        print(f"          sha256 {panel.content_hash()[:16]}...")
        print(f"\nscore columns: {list(panel.score_columns)}")
        print("\nPick one with --score. Masked and unmasked answer different questions "
              "(masked\nsuppresses scores at annotated sites), so the choice belongs on "
              "the record.")
        return 0

    report = api.verify_splice_gate(
        panel.cases(args.score),
        negative_construction=panel.negative_construction,
        panel_note=panel.assay,
        min_pr_auc=args.min_pr_auc,
        min_pr_auc_skill=args.min_pr_auc_skill,
        max_ece=args.max_ece,
    )

    held = "held out" if panel.held_out else "NOT HELD OUT"
    print(f"panel:    {summary['n_rows']} variants / {len(panel.groups)} genes  [{held}]")
    print(f"          sha256 {panel.content_hash()[:16]}...")
    print(f"score:    {args.score}")
    if panel.training_overlap:
        print(f"\nWARNING: {list(panel.training_overlap)} are on chromosomes both models "
              "trained on.\n         These metrics are optimistic and cannot support "
              "promotion. Rebuild with\n         --held-out-only for a held-out run.")
    print()
    print(
        f"{'stratum':<12}{'n':>7}{'prev':>8}{'AP':>9}{'skill':>9}{'ROC':>9}{'ECE':>9}"
    )
    for stratum in report.strata:
        print(f"{stratum.name:<12}{stratum.n_cases:>7}{stratum.prevalence:>8.3f}"
              f"{stratum.pr_auc:>9.3f}{stratum.pr_auc_skill:>9.3f}"
              f"{stratum.roc_auc:>9.3f}{stratum.ece:>9.3f}")
    print()
    print("'prev' is average precision's FLOOR. Compare strata on 'skill' (AP rescaled")
    print("so 0 is no-skill at any prevalence), never on raw AP.")
    print()

    # The published anchor this panel exists to reproduce.
    by_name = {s.name: s.pr_auc for s in report.strata}
    n_rows, n_groups = len(panel.rows), len(panel.groups)
    if {"exonic", "intronic"} <= set(by_name):
        print("Smith & Kitzman 2023 report, as a MEDIAN ACROSS TOOLS pooled over ALL SIX")
        print("of their datasets (MLH1 included):")
        print(f"  exonic  0.419   (this run: {by_name['exonic']:.3f})")
        print(f"  intronic 0.773  (this run: {by_name['intronic']:.3f})")
        print()
        print(f"This panel has {n_rows} variants over {n_groups} gene(s).")
        if n_groups < 6:
            print("--include-mlh1 matches the published composition (six datasets).")
        print()
        print("NEITHER the levels NOR the gap are directly comparable to those numbers.")
        print("A single tool is not a median over tools: a strong one sits above the")
        print("median on both strata, and by more on the harder one, which COMPRESSES")
        print("the gap. Measured here, matching the composition moved both figures")
        print("FURTHER from the published pair, not closer.")
        print()
        print("What is comparable is the ORDERING -- and it must be read on SKILL, not on")
        print("AP, because the two strata rarely share a prevalence. (The archive carries")
        print("all eight tools' scores, so an exact reproduction is possible; it needs the")
        print("other six mapped in.)")
        skill = {s.name: s.pr_auc_skill for s in report.strata}
        prevalence = {s.name: s.prevalence for s in report.strata}
        ap_inverted = by_name["exonic"] >= by_name["intronic"]
        skill_ordered = skill["exonic"] < skill["intronic"]
        print()
        if ap_inverted and skill_ordered:
            # The documented failure mode, caught rather than left to mislead: raw AP
            # can invert purely because one stratum has more positives.
            print("NOTE: raw AP is INVERTED here (exonic above intronic), but that is a")
            exonic_prev = f"{prevalence['exonic']:.1%}"
            intronic_prev = f"{prevalence['intronic']:.1%}"
            print(f"      prevalence artifact, not a finding -- exonic {exonic_prev} vs "
                  f"intronic {intronic_prev}")
            print("      positives, so exonic AP starts from a much higher floor. On skill")
            print(f"      the expected ordering holds ({skill['exonic']:.3f} < "
                  f"{skill['intronic']:.3f}). The panel is fine.")
        elif not skill_ordered:
            print("WARNING: exonic skill is NOT below intronic. That is the one result that")
            print("         should prompt suspicion of the panel build before the model.")
    print()
    print(f"gate passed (thresholds) : {report.passed}")
    for reason in report.reasons:
        print(f"  - {reason}")
    print(
        "\nThis command flips nothing. It measures a backend's scores against a measured "
        "assay;\nit does not promote anything, and the exonic stratum is the half BT4 "
        "actually operates in.",
        file=sys.stderr,
    )
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bt4", description="BT4 back-translation optimizer")
    parser.add_argument("--version", action="version", version=f"bt4 {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_opt = sub.add_parser("optimize", help="optimize a protein into coding DNA")
    p_opt.add_argument("protein", help="stop-free amino-acid string")
    p_opt.add_argument("--fasta", action="store_true", help="emit FASTA")
    p_opt.add_argument("--json", action="store_true", help="emit JSON (with manifest)")
    p_opt.add_argument("--genbank", action="store_true",
                       help="emit an annotated GenBank record: residual violations "
                       "become misc_feature spans, so a defect the optimizer could "
                       "not remove is visible on the map you open. Includes the "
                       "construct context (--utr5/--utr3) when given")
    p_opt.add_argument("--header", default="bt4", help="FASTA header")
    p_opt.add_argument("--cpg-min", type=int, default=None, dest="cpg_min",
                       help="min total CpG (CG) count over the CDS (exact budget DP)")
    p_opt.add_argument("--cpg-max", type=int, default=None, dest="cpg_max",
                       help="max total CpG (CG) count -- CpG depletion for stealth")
    p_opt.add_argument("--upa-min", type=int, default=None, dest="upa_min",
                       help="min total UpA (TA) count over the CDS (exact budget DP)")
    p_opt.add_argument("--upa-max", type=int, default=None, dest="upa_max",
                       help="max total UpA (TA) count over the CDS")
    p_opt.add_argument("--use-attested-splice", action="store_true",
                       dest="use_attested_splice",
                       help="Honor the bundled splice fidelity attestation, promoting a "
                            "wrapped CNN to calibrated=True. Off by default: an attestation "
                            "records that BT4's wrapper is faithful to the published model, "
                            "not that its scores are calibrated probabilities for coding "
                            "sequence (see docs/DESIGN_splice_cnn_calibration.md).")
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
    p_val.add_argument("--use-attested-splice", action="store_true",
                       dest="use_attested_splice",
                       help="Honor the bundled splice fidelity attestation (see "
                            "'optimize --use-attested-splice').")
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

    # argparse %-formats every help string as it RENDERS help, so a literal
    # percent sign must be escaped as "%%": a bare "%M" raises ValueError from
    # inside the help formatter, which took out both `bt4 --help` (the subparser
    # entry below is listed there) and `bt4 tracks --help`.
    p_trk = sub.add_parser("tracks", help="per-site composition tracks (GC/CpG/%%MinMax)")
    p_trk.add_argument("dna", help="ACGT coding sequence")
    p_trk.add_argument(
        "--organism", default="homo_sapiens", help="table for the %%MinMax track"
    )
    p_trk.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
        help="that organism's reference set for the %%MinMax frequencies "
             "(default: the organism's own default)",
    )
    p_trk.add_argument("--splice-track", action="store_true", dest="splice_track",
                       help="add a per-nucleotide cryptic-splice-site track. Opt-in "
                       "because it runs a model: with the default PWM baseline the "
                       "values are UNCALIBRATED pseudo-scores showing where "
                       "consensus-like positions sit, not probabilities")
    p_trk.add_argument("--nt-window", type=int, default=50, dest="nt_window",
                       help="nucleotide window for the GC/CpG tracks")
    p_trk.add_argument("--codon-window", type=int, default=18, dest="codon_window",
                       help="codon window for the %%MinMax track")
    p_trk.add_argument("--json", action="store_true", help="emit full per-window arrays as JSON")
    p_trk.set_defaults(func=_cmd_tracks)

    p_gate = sub.add_parser(
        "expression-gate",
        help="run the expression acceptance gate on a measured CDS-variant panel",
    )
    p_gate.add_argument("panel", help="panel TSV (group/variant_id/cds/measured/utr5/utr3)")
    p_gate.add_argument("--backend", default="ribonn", help="expression backend name")
    p_gate.add_argument("--species", default="human", choices=("human", "mouse"))
    p_gate.add_argument(
        "--cell-type", action="append", dest="cell_types", default=None,
        help="restrict the head to this cell type (repeatable); match it to the panel, "
             "since averaging every tissue against a single-cell-line measurement is a "
             "scope error",
    )
    p_gate.add_argument(
        "--within-group", action="store_true", dest="within_group",
        help="score inside each protein -- the strict bar, and the regime BT4 actually "
             "deploys in. Without it, a head that merely recognises which gene it is "
             "looking at scores well",
    )
    p_gate.add_argument(
        "--recalibrate", action="store_true",
        help="fit measured ~ a*pred + b on the calibration fold before residuals; "
             "required whenever the head's units differ from the assay's",
    )
    p_gate.add_argument("--min-spearman", type=float, default=0.30, dest="min_spearman")
    p_gate.add_argument(
        "--target-coverage", type=float, default=0.90, dest="target_coverage"
    )
    p_gate.add_argument(
        "--coverage-tolerance", type=float, default=0.05, dest="coverage_tolerance"
    )
    p_gate.add_argument(
        "--calibration-fraction", type=float, default=0.50, dest="calibration_fraction"
    )
    p_gate.add_argument(
        "--bootstrap-resamples", type=int, default=1000, dest="bootstrap_resamples"
    )
    p_gate.add_argument("--num-workers", type=int, default=0, dest="num_workers")
    p_gate.add_argument("--organism", default="homo_sapiens")
    p_gate.add_argument(
        "--reference-set", default=None, dest="reference_set",
        choices=list(api.REFERENCE_SETS),
    )
    p_gate.add_argument("--seed", type=int, default=0)
    p_gate.set_defaults(func=_cmd_expression_gate)

    p_sgate = sub.add_parser(
        "splice-gate",
        help="run the splice acceptance gate on an annotated splice-site panel",
    )
    p_sgate.add_argument(
        "panel", help="panel TSV (window_id/group/sequence/donors/acceptors)"
    )
    p_sgate.add_argument(
        "--negative-construction", required=True, dest="negative_construction",
        help="how the negative class was built, verbatim (e.g. 'all other positions in "
             "the same gene bodies'). REQUIRED: average precision's floor is the "
             "prevalence, which is a construction choice, so a threshold without a "
             "pinned denominator is passable by sampling fewer negatives",
    )
    p_sgate.add_argument(
        "--annotation", default="",
        help="the gene model the sites came from, e.g. 'GENCODE v44 / GRCh38'. "
             "Annotation choice alone moved SpliceAI's predictions for >10%% of "
             "variants in some genes (Smith & Kitzman 2023)",
    )
    p_sgate.add_argument(
        "--backend", default="pwm", choices=_OFFLINE_SPLICE_BACKENDS,
        help="ASSP is deliberately absent: it is network-derived and excluded from the "
             "reproducible-from-manifest guarantee, so it cannot support a gate result",
    )
    p_sgate.add_argument(
        "--anchor-offset", type=int, default=None, dest="anchor_offset",
        help="one offset for every site kind. Fine for a kind-separated panel; NOT "
             "sufficient for a mixed one scored by a real backend -- see --cnn-anchors",
    )
    p_sgate.add_argument(
        "--donor-offset", type=int, default=None, dest="donor_offset",
        help="anchor offset for donors only (overrides --anchor-offset)",
    )
    p_sgate.add_argument(
        "--acceptor-offset", type=int, default=None, dest="acceptor_offset",
        help="anchor offset for acceptors only (overrides --anchor-offset)",
    )
    p_sgate.add_argument(
        "--cnn-anchors", action="store_true", dest="cnn_anchors",
        help="use the measured SpliceAI/Pangolin anchors: donor -1, acceptor +1. Both "
             "models score a site on the EXONIC boundary base while BT4's panel anchors "
             "on the intronic dinucleotide, so the two kinds need OPPOSITE offsets and "
             "no single value is correct for a mixed panel",
    )
    p_sgate.add_argument(
        "--min-motif-consistency", type=float, default=api.MIN_SPLICE_MOTIF_CONSISTENCY,
        dest="min_motif_consistency",
        help="fraction of sites that must carry their canonical dinucleotide. Lower it "
             "only for a deliberately non-canonical (U12 AT-AC) panel, never to quiet "
             "the off-by-one this check exists to catch",
    )
    p_sgate.add_argument(
        "--threshold", type=float, default=api.DEFAULT_SITE_PROBABILITY,
        help="the operating point the MCC is computed at (BT4's own, by default)",
    )
    p_sgate.add_argument("--min-pr-auc", type=float, default=0.0, dest="min_pr_auc")
    p_sgate.add_argument(
        "--min-pr-auc-skill", type=float, default=0.0, dest="min_pr_auc_skill",
        help="per-stratum floor on average precision rescaled so no-skill is 0 at any "
             "prevalence. NOT comparable across panels -- it still moves with how the "
             "negatives were sampled, which is why --negative-construction is required",
    )
    p_sgate.add_argument(
        "--max-ece", type=float, default=1.0, dest="max_ece",
        help="per-stratum calibration ceiling. Reported, but NOT a bar on its own: a "
             "base-rate predictor scores ECE 0.0 at splice prevalence, so declaring "
             "only this leaves the run un-promotable. Pair it with --min-pr-auc-skill",
    )
    p_sgate.add_argument("--bins", type=int, default=10, help="reliability bins")
    p_sgate.add_argument("--seed", type=int, default=0)
    p_sgate.add_argument(
        "--quiet", action="store_true",
        help="suppress the per-window scoring progress on stderr. A CNN-backed run takes "
             "tens of minutes and reports nothing without it; the report on stdout is "
             "unaffected either way",
    )
    p_sgate.set_defaults(func=_cmd_splice_gate)

    p_vgate = sub.add_parser(
        "variant-gate",
        help="run the splice gate on a measured variant panel (no model needed)",
    )
    p_vgate.add_argument("panel", help="variant panel TSV (from make_splicebench_variant_panel.py)")
    p_vgate.add_argument(
        "--score", default=None,
        help="which score column to gate. Omit to list what the panel carries",
    )
    p_vgate.add_argument(
        "--list-scores", action="store_true", dest="list_scores",
        help="list the panel's score columns and exit",
    )
    p_vgate.add_argument(
        "--negative-construction", required=True, dest="negative_construction",
        help="how the negative class was built, verbatim. REQUIRED: average precision's "
             "floor is the prevalence, which is a construction choice",
    )
    p_vgate.add_argument(
        "--assay", default="",
        help="what was measured and under what criterion. Record it when the labels pool "
             "several assays -- a composite is not a measurement",
    )
    p_vgate.add_argument("--min-pr-auc", type=float, default=0.0, dest="min_pr_auc")
    p_vgate.add_argument(
        "--min-pr-auc-skill", type=float, default=0.0, dest="min_pr_auc_skill"
    )
    p_vgate.add_argument("--max-ece", type=float, default=1.0, dest="max_ece")
    p_vgate.set_defaults(func=_cmd_variant_gate)

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
    parser = _parser()
    args = parser.parse_args(argv)
    if getattr(args, "preset", None):
        try:
            _apply_preset_to_args(args, argv)
        except KeyError as exc:
            print(f"error: {exc.args[0] if exc.args else exc}", file=sys.stderr)
            return 2
    try:
        exit_code: int = args.func(args)
    except (ValueError, InfeasibleError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
