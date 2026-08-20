"""Measure a splice backend's **detection floor** inside designed coding sequence.

The question this answers, and the one it does not
-------------------------------------------------
On a label-free designed-CDS panel, "the model is blind here" and "the model is
correctly silent because there is no strong site here" produce **identical output**.
Nothing in BT4 could separate them, and the difference decides whether a splice
backend is usable at all in the regime BT4 designs for.

This separates them the only way available without an assay: plant a site whose
position you chose, and see whether the model finds it. That bounds **inertness from
below**. It is a negative control, and it is *not* evidence of correct silence --
detecting a site BT4 planted says nothing about sites nobody put there. Inferring
"the CDS is clean" from "the planted site was found" is affirming the consequent, and
this script deliberately reports nothing that would support it.

Design
------
* **Substitution, never insertion.** A planted 9-mer replaces 9 bases, so CDS length
  and reading frame are preserved and the probe stays a coding sequence. A probe that
  introduces an in-frame stop is **skipped and reported**, never scored silently.
* **A graded ladder, not one motif.** Full consensus, weakened, weaker, a ``GT``->``CT``
  ablation, and a composition-matched scramble. The ablation is the load-bearing
  control: it keeps 7 of the 9 bases and destroys only the invariant dinucleotide, so a
  model that still responds is reacting to *an edit* rather than to a splice signal.
* **The scramble is constrained, and the constraint is stated.** It must be a
  permutation of the full consensus that introduces no in-frame stop and does not
  recreate ``GT`` at the junction. A first pass used an unconstrained shuffle that hit a
  stop codon at every plant site and was skipped everywhere -- a control that never ran.
* **Real flanking context by default.** The adapters pad with 5,000 literal ``N``, which
  is measured to deflate scores inside the CDS; ``--flank-fasta`` supplies real sequence
  instead. ``--n-padded`` reproduces the shipped path for comparison.
* **Anchor-aware.** A donor's score sits one base 5' of the intronic ``G``
  (``CNN_ANCHOR_OFFSETS``), so the reported offset is measured against that prediction
  rather than the plant coordinate. An offset that is not 0 is a finding, not noise.

Licensing
---------
Prints **derived scalars only** -- peak, offset, and whether a threshold was cleared.
Per-position score arrays are the licensed model's output (Pangolin GPL-3.0, SpliceAI
CC BY-NC 4.0) and are never written out. Run it where the weights live.

Determinism (invariant #7): output depends only on the inputs, the ladder constants,
and ``--seed``; no timestamp is emitted.

Run it::

    python scripts/probe_splice_detection_floor.py --panel designed.tsv \\
        --flank-fasta chr1_window.fa --groups KRas4B Beclin1
"""

from __future__ import annotations

import argparse
import itertools
import statistics
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT / "src") not in sys.path:  # pragma: no cover - script convenience
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from bt4.biomodels.splice.base import score_in_context  # noqa: E402
from bt4.domain.genetic_code import CODON_TABLE  # noqa: E402

__all__ = ["LADDER", "Rung", "build_scramble", "main", "plant", "probe_host"]

FULL_CONSENSUS = "CAGGTAAGT"
"""The canonical 5' splice site as a 9-mer: exon ``MAG`` | intron ``GTRAGT``.

The junction falls between index 2 and 3, so the invariant ``GT`` sits at 3..4. Every
other rung is defined by how it differs from this one.
"""

JUNCTION = 3
"""Index within the 9-mer where the intron begins (the ``G`` of the invariant ``GT``)."""


@dataclass(frozen=True, slots=True)
class Rung:
    """One step of the ladder.

    Attributes:
        label: Human-readable name, printed verbatim.
        motif: The 9-mer planted, or ``None`` for the unmodified host baseline.
        note: What this rung controls for -- printed so a reader of the output knows
            why the rung exists without reading this file.
    """

    label: str
    motif: str | None
    note: str


LADDER: tuple[Rung, ...] = (
    Rung("L0 host (unmodified)", None, "baseline: what the host scores untouched"),
    Rung("L1 full consensus", FULL_CONSENSUS, "strongest possible donor"),
    Rung("L2 weakened", "CAGGTATGT", "R->T at +5; degraded but keeps the invariant GT"),
    Rung("L3 weaker", "AAAGTATCT", "exonic and distal positions degraded"),
    Rung(
        "L4 GT->CT ablation",
        "CAGCTAAGT",
        "7 of 9 bases identical to L1, invariant GT destroyed -- the decisive control",
    ),
)


def build_scramble(consensus: str = FULL_CONSENSUS) -> str:
    """Return a composition-matched permutation that is a *usable* control.

    A bare shuffle is not good enough: it can introduce an in-frame stop (making the
    probe non-coding, so it is skipped and the control silently never runs) or
    recreate ``GT`` at the junction (making it a weak positive rather than a negative).

    Args:
        consensus: The 9-mer whose base composition the scramble must match.

    Returns:
        A permutation of ``consensus`` with no in-frame stop, no ``GT`` at the
        junction, and not equal to ``consensus``. Chosen as the median of the sorted
        candidates so the result is deterministic.

    Raises:
        ValueError: If no permutation satisfies the constraints.
    """
    candidates = sorted(
        {
            "".join(p)
            for p in itertools.permutations(consensus)
            if not _has_stop("".join(p))
            and "".join(p)[JUNCTION : JUNCTION + 2] != "GT"
            and "".join(p) != consensus
        }
    )
    if not candidates:
        raise ValueError(f"no usable composition-matched scramble of {consensus!r}")
    return candidates[len(candidates) // 2]


def _has_stop(seq: str) -> bool:
    """Whether ``seq`` contains an in-frame stop codon, reading from index 0."""
    return any(CODON_TABLE.get(seq[i : i + 3]) == "*" for i in range(0, len(seq) - 2, 3))


def plant(host: str, motif: str, at: int) -> str:
    """Substitute ``motif`` into ``host`` at ``at``, preserving length and frame.

    Args:
        host: The coding sequence to plant into.
        motif: The 9-mer to substitute.
        at: Start index, which the caller must keep a multiple of 3 so the motif's own
            codons stay in the host's frame.

    Returns:
        The probe sequence, the same length as ``host``.
    """
    probe = host[:at] + motif + host[at + len(motif) :]
    assert len(probe) == len(host), "planting must not change CDS length"
    return probe


def probe_host(
    predictor: object,
    host: str,
    positions: Sequence[int],
    upstream: str,
    downstream: str,
    *,
    threshold: float,
    window: int = 6,
) -> list[tuple[int, Rung, float | None, int | None]]:
    """Score every ladder rung at every position in one host.

    Returns:
        ``(position, rung, peak, offset)`` per probe. ``peak`` and ``offset`` are
        ``None`` when the probe was skipped for introducing an in-frame stop -- a
        skipped rung is reported, never dropped.
    """
    from bt4.biomodels.splice.base import SpliceResult  # noqa: F401  (typing only)

    rows: list[tuple[int, Rung, float | None, int | None]] = []
    scramble = Rung(
        "L5 scrambled",
        build_scramble(),
        "composition-matched to L1, order destroyed -- controls base content",
    )
    for at in positions:
        # A donor's score anchors one base 5' of the intronic G (CNN_ANCHOR_OFFSETS).
        expected = at + JUNCTION - 1
        for rung in (*LADDER, scramble):
            probe = host if rung.motif is None else plant(host, rung.motif, at)
            if rung.motif is not None and _has_stop(probe[at : at + len(rung.motif)]):
                rows.append((at, rung, None, None))
                continue
            result = score_in_context(predictor, probe, upstream, downstream)  # type: ignore[arg-type]
            lo = max(0, expected - window)
            hi = min(len(result.donor), expected + window + 1)
            local = list(result.donor[lo:hi])
            peak = max(local) if local else 0.0
            rows.append((at, rung, peak, lo + local.index(peak) - expected))
    return rows


def _read_flank(path: Path, flank: int) -> tuple[str, str]:
    """Return ``(upstream, downstream)`` of ``flank`` nt each from a FASTA."""
    seq = "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    ).upper()
    if len(seq) < 2 * flank:
        raise ValueError(f"{path} holds {len(seq)} nt, need {2 * flank}")
    return seq[:flank], seq[flank : 2 * flank]


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ladder and print the summary."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--panel", required=True, help="designed-CDS panel TSV")
    parser.add_argument("--provenance", required=True, help="the panel's provenance line")
    parser.add_argument("--backend", default="pangolin", choices=("pangolin", "spliceai"))
    parser.add_argument("--groups", nargs="*", help="proteins to probe (default: all)")
    parser.add_argument(
        "--flank-fasta",
        help="real flanking sequence. Omit to reproduce the shipped N-padded path, "
        "which is measured to deflate scores inside the CDS",
    )
    parser.add_argument("--flank", type=int, default=5000, help="nt of flank per side")
    parser.add_argument("--threshold", type=float, default=0.5, help="display cutoff only")
    args = parser.parse_args(argv)

    from bt4.api import read_designed_cds_panel
    from bt4.pipeline.splice_crosscheck import resolve_splice_backend

    panel = read_designed_cds_panel(args.panel, provenance=args.provenance)
    predictor = resolve_splice_backend(args.backend)
    up, dn = _read_flank(Path(args.flank_fasta), args.flank) if args.flank_fasta else ("", "")
    if not up:
        print("NOTE: no --flank-fasta, so the adapter pads with literal N. That path is")
        print("      measured to deflate scores inside the CDS; read peaks as a lower bound.\n")

    groups = tuple(args.groups) if args.groups else panel.groups
    print(f"{'host':10s} {'pos':>6s} {'rung':22s} {'peak':>7s} {'off':>4s} {'>thr':>5s}")
    print("-" * 60)
    pooled: dict[str, list[float]] = {}
    for group in groups:
        member = next(m for m in panel.group_members(group) if m.role == "designed")
        positions = [len(member.cds) * f // 3 * 3 for f in (1, 2, 3)]
        positions = [p // 4 for p in positions]
        for at, rung, peak, offset in probe_host(
            predictor, member.cds, positions, up, dn, threshold=args.threshold
        ):
            if peak is None:
                print(f"{group:10s} {at:6d} {rung.label:22s}  skipped: in-frame stop")
                continue
            pooled.setdefault(rung.label, []).append(peak)
            mark = "YES" if peak > args.threshold else "no"
            print(f"{group:10s} {at:6d} {rung.label:22s} {peak:7.4f} {offset:+4d} {mark:>5s}")
        print()

    print("=" * 60)
    print(f"Pooled over hosts x positions (threshold {args.threshold:g}):")
    for rung in (*LADDER, Rung("L5 scrambled", "", "")):
        values = pooled.get(rung.label)
        if not values:
            continue
        cleared = sum(v > args.threshold for v in values)
        print(
            f"  {rung.label:22s} n={len(values)} median {statistics.median(values):.4f} "
            f"cleared {cleared}/{len(values)}"
        )
    print(
        "\nThis bounds INERTNESS from below and nothing else. Detecting a site planted\n"
        "here is not evidence that an unmodified designed CDS is free of cryptic sites --\n"
        "that inference is affirming the consequent, and no label-free probe supplies it."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
