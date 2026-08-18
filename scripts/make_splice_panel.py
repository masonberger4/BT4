"""Build the deterministic sequence panel for the splice integration-fidelity gate.

Step **A4a** of [`docs/DESIGN_splice_cnn_calibration.md`](../docs/DESIGN_splice_cnn_calibration.md).
The gate (:func:`~bt4.biomodels.splice.verify_pangolin_fidelity` /
:func:`~bt4.biomodels.splice.verify_spliceai_fidelity`) needs two halves that must
come from *different* sources:

1. **the sequences** -- this script, which may freely use ``bt4``; DNA is not a
   licensed artifact and the panel is committed-friendly and reproducible;
2. **the expected per-position scores** -- ``scripts/capture_pangolin_panel.py``,
   which must **never** import ``bt4``, because capturing expectations with the
   very adapter under test would make the gate pass trivially.

Keeping the halves in separate files is what makes that separation auditable
rather than a promise in a docstring.

**Why the panel's composition matters.** The gate is a numerical equivalence
check, so in principle any sequences would do. In practice a panel on which the
model scores ~0 everywhere is *weak*: a wrong output channel, a wrong fold set, or
a transposed one-hot would still "match" within tolerance because both sides are
near zero. A useful panel therefore has to make the model produce a **wide spread**
of scores. This script builds one accordingly, and
``scripts/run_splice_fidelity_gate.py`` reports the realized spread so a weak panel
is visible rather than silently flattering.

Panel composition (every sequence pure ``ACGT`` -- ``validate_dna`` rejects ``N``):

* ``cds`` -- BT4-designed coding sequences, the regime BT4 actually ships;
* ``cds_variant`` -- sampled synonymous variants of those same proteins, the
  *exact* regime the splice audit runs on (designed, not natural);
* ``donor`` / ``acceptor`` -- windows carrying a canonical U2 consensus
  (``MAG|GTRAGT`` / polypyrimidine + ``YAG|G``) centred in flanking sequence, so
  the model has something it should score **high**. These are synthetic probes,
  **not** claimed to be real annotated splice sites -- they exist to move the
  score range, and no biological claim is made about them;
* ``edge`` -- a very short sequence and GC-extreme sequences, which guard the
  padding / context-crop arithmetic at the boundaries.

Determinism (invariant #7): sequence generation is seeded and uses no wall clock,
so re-running produces a byte-identical panel and the same ``content_hash``.

Run it::

    python scripts/make_splice_panel.py --out panel_sequences.json
    python scripts/make_splice_panel.py --summary        # inspect without writing
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections.abc import Sequence
from dataclasses import asdict, dataclass

__all__ = ["PanelSequence", "build_panel", "main", "panel_content_hash"]

# Proteins used to generate the designed-CDS members. Short, real-ish, and
# deliberately degenerate (many synonymous choices) so sampled variants differ.
_PROTEINS: tuple[tuple[str, str], ...] = (
    ("insulin_b", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"),
    ("gfp_frag", "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTL"),
    (
        "degenerate",
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    ),
)

# Canonical U2 consensus probes. These are *textbook consensus* strings, the same
# motifs BT4's PWM baseline and SpliceSiteMotifConstraint encode -- included to
# push the model's score range up, never as a claim that they are real sites.
_DONOR_CONSENSUS = "CAGGTAAGT"      # MAG | GT RAGT, invariant GT at +1,+2
_ACCEPTOR_CONSENSUS = "TTTTTCTTTTTCCAGG"  # polypyrimidine tract, then YAG | G

_BASES = "ACGT"


@dataclass(frozen=True, slots=True)
class PanelSequence:
    """One sequence in the fidelity panel (no scores -- those are captured later).

    Attributes:
        id: Stable identifier, unique within the panel.
        kind: ``cds`` | ``cds_variant`` | ``donor`` | ``acceptor`` | ``edge``.
        sequence: The pure-``ACGT`` sequence to score.
        note: Human-readable provenance of this member.
    """

    id: str
    kind: str
    sequence: str
    note: str


def _rng(seed: int) -> random.Random:
    """Return a seeded RNG (no global state, so callers cannot perturb each other)."""
    return random.Random(seed)


def _random_dna(rng: random.Random, length: int, gc: float = 0.5) -> str:
    """Return random ACGT of ``length`` with an approximate GC fraction ``gc``."""
    out = []
    for _ in range(length):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def _centred_probe(rng: random.Random, motif: str, flank: int, seed_offset: int) -> str:
    """Embed ``motif`` in the middle of random flanking sequence.

    Centring matters: both wrapped CNNs crop a 10,000 nt context and BT4 pads the
    remainder with ``N``, so a motif at the very edge of a short sequence sits in
    the padding-artifact regime that OpenSpliceAI documents. Putting it in the
    middle keeps the probe's score driven by real sequence on both sides.
    """
    left = _random_dna(_rng(seed_offset), flank)
    right = _random_dna(_rng(seed_offset + 1), flank)
    return f"{left}{motif}{right}"


def _designed_sequences(seed: int) -> list[PanelSequence]:
    """Generate BT4-designed CDSs and sampled synonymous variants of them.

    Imports :mod:`bt4.api` lazily and *only here*: this is the one half of the
    panel that is allowed to touch BT4 (see the module docstring).
    """
    from bt4 import api

    members: list[PanelSequence] = []
    for name, protein in _PROTEINS:
        result = api.optimize(protein)
        members.append(
            PanelSequence(
                id=f"cds_{name}",
                kind="cds",
                sequence=result.dna,
                note=f"api.optimize({name}), {len(result.dna)} nt, default organism",
            )
        )
        lib = api.library(protein, n=2, seed=seed, temperature=1.5)
        for i, member in enumerate(lib.results):
            members.append(
                PanelSequence(
                    id=f"var_{name}_{i}",
                    kind="cds_variant",
                    sequence=member.dna,
                    note=f"api.library({name}, seed={seed}, T=1.5) member {i} -- SAMPLED",
                )
            )
    return members


def build_panel(seed: int = 0, *, include_designed: bool = True) -> list[PanelSequence]:
    """Build the full deterministic panel.

    Args:
        seed: Master seed. Same seed -> byte-identical panel (invariant #7).
        include_designed: When ``False``, skip the members that require ``bt4``
            (useful for generating probe-only panels in an environment where BT4
            is not installed).

    Returns:
        The panel, ordered deterministically.
    """
    rng = _rng(seed)
    members: list[PanelSequence] = []

    if include_designed:
        members.extend(_designed_sequences(seed))

    # Consensus probes -- the members that should make the model light up.
    for i in range(3):
        members.append(
            PanelSequence(
                id=f"donor_probe_{i}",
                kind="donor",
                sequence=_centred_probe(rng, _DONOR_CONSENSUS, 300, seed + 100 + i * 2),
                note="synthetic canonical donor consensus MAG|GTRAGT centred in random flanks",
            )
        )
        members.append(
            PanelSequence(
                id=f"acceptor_probe_{i}",
                kind="acceptor",
                sequence=_centred_probe(rng, _ACCEPTOR_CONSENSUS, 300, seed + 200 + i * 2),
                note=(
                    "synthetic canonical acceptor consensus (polyY + YAG|G) "
                    "centred in random flanks"
                ),
            )
        )

    # Edge cases -- these guard the padding / crop arithmetic, not the biology.
    members.append(
        PanelSequence(
            id="edge_short",
            kind="edge",
            sequence=_random_dna(_rng(seed + 300), 99),
            note="99 nt -- far shorter than the 10,000 nt context crop",
        )
    )
    members.append(
        PanelSequence(
            id="edge_gc_high",
            kind="edge",
            sequence=_random_dna(_rng(seed + 301), 600, gc=0.85),
            note="600 nt at ~85% GC",
        )
    )
    members.append(
        PanelSequence(
            id="edge_gc_low",
            kind="edge",
            sequence=_random_dna(_rng(seed + 302), 600, gc=0.15),
            note="600 nt at ~15% GC",
        )
    )
    return members


def panel_content_hash(members: Sequence[PanelSequence]) -> str:
    """Return a deterministic SHA-256 over the panel's sequences.

    Lets a captured score file be bound to the exact panel it was captured from,
    so a mismatched pairing is detectable rather than silently compared.
    """
    canonical = json.dumps([asdict(m) for m in members], sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def main(argv: Sequence[str] | None = None) -> int:
    """Build the panel and write it as JSON (or summarize it)."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--out", default=None, help="Write the panel JSON here.")
    parser.add_argument("--seed", type=int, default=0, help="Master seed (default: 0).")
    parser.add_argument(
        "--no-designed",
        action="store_true",
        help="Skip members that require bt4 (probe/edge members only).",
    )
    parser.add_argument("--summary", action="store_true", help="Print a summary table.")
    args = parser.parse_args(argv)

    members = build_panel(args.seed, include_designed=not args.no_designed)
    payload = {
        "schema_version": 1,
        "seed": args.seed,
        "content_hash": panel_content_hash(members),
        "sequences": [asdict(m) for m in members],
    }

    if args.summary or not args.out:
        print(f"{'id':22} {'kind':12} {'nt':>6}  note")
        for m in members:
            print(f"{m.id:22} {m.kind:12} {len(m.sequence):>6}  {m.note}")
        by_kind: dict[str, int] = {}
        for m in members:
            by_kind[m.kind] = by_kind.get(m.kind, 0) + 1
        kinds = ", ".join(f"{k}={v}" for k, v in sorted(by_kind.items()))
        print(f"\n{len(members)} sequences: {kinds}")
        print(f"content_hash: {payload['content_hash']}")

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
