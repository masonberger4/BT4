"""Build the designed synonymous CDS panel — BT4's own regime.

Converts a grouped CDS FASTA into the format
:func:`bt4.api.read_designed_cds_panel` reads. The bundled default is
``scripts/data/ranaghan2021_tab4.fasta``: **93 sequences, three proteins, each with the
native CDS plus 30 designs** (three anonymized commercial algorithms x ten repeat runs),
from Ranaghan et al. 2021, CC BY 4.0, already cited elsewhere in this repository.

**Why that source rather than sequences generated here.** The panel exists to ask how
splice models behave on designed coding sequence. Generating the designs with BT4 would
make the answer partly a fact about BT4, and BT4 is the tool that would consume it -- the
same circularity the splice gate refuses when it runs BT4's own PWM as a permanent
baseline. Ranaghan's designs come from independent commercial optimizers, are published,
and are citable. ``--include-bt4`` adds BT4's own design for the same protein as one more
member when the comparison is wanted explicitly; it is off by default, and the format
records who designed each sequence in the ``note`` rather than in a field the analysis
can branch on.

**What this panel cannot do.** It carries no splice labels, because designed coding
sequence has no splice ground truth: nothing here has been assayed, none of it appears in
any annotation, and a motif is not a site. See
``bt4.biomodels.splice.designed_panel`` for what is measurable without labels --
response to synonymous change, cross-backend agreement, and Δ against the native.

Determinism (invariant #7): output depends only on the input FASTA and the flags. With
``--include-bt4`` it additionally depends on BT4's optimizer, which is itself
deterministic from its config and seed.

Run it::

    python scripts/make_designed_cds_panel.py --out designed.tsv
    python scripts/make_designed_cds_panel.py --out designed.tsv --include-bt4
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT / "src") not in sys.path:  # pragma: no cover - script convenience
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from bt4.biomodels.splice.designed_panel import (  # noqa: E402
    DesignedMember,
    panel_from_members,
    write_designed_cds_panel,
)

__all__ = ["DEFAULT_FASTA", "PROVENANCE", "build_members", "main", "read_fasta"]

DEFAULT_FASTA = _REPO_ROOT / "scripts" / "data" / "ranaghan2021_tab4.fasta"
"""Ranaghan et al. 2021 Table 4: 3 proteins x (1 native + 3 algorithms x 10 runs)."""

PROVENANCE = (
    "Ranaghan et al. 2021 (CC BY 4.0) Table 4 -- three proteins, each the native human "
    "CDS plus 30 designs from three ANONYMIZED commercial optimizers over ten repeat "
    "runs each. Tools are anonymized upstream and are not re-identified here; the "
    "repeat runs make these a determinism axis as well as a design axis."
)
"""Recorded verbatim into the panel. Required by the format."""


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return ``(header, sequence)`` pairs, in file order.

    A local reader rather than a dependency: the panel builder must run wherever BT4
    runs, and this file format is two rules.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records:
        raise ValueError(f"{path} contains no FASTA records")
    return records


def _parse_header(header: str) -> tuple[str, str, str]:
    """Return ``(group, member_id, role)`` from a Ranaghan-style header.

    The format is pipe-separated, ``Protein|Source|run|...``. Anything whose source
    field is ``Native`` (case-insensitively) is the group's reference; everything else
    is a design.
    """
    fields = [part.strip() for part in header.split("|") if part.strip()]
    if not fields:
        raise ValueError(f"unparseable FASTA header: {header!r}")
    group = fields[0]
    source = fields[1] if len(fields) > 1 else "unknown"
    run = fields[2] if len(fields) > 2 and fields[2].lower().startswith("run") else ""
    role = "native" if source.lower() == "native" else "designed"
    member_id = "_".join(part for part in (group, source, run) if part)
    return group, member_id, role


def build_members(
    records: Sequence[tuple[str, str]], *, include_bt4: bool = False
) -> list[DesignedMember]:
    """Turn FASTA records into panel members, optionally adding BT4's own designs.

    Args:
        records: ``(header, sequence)`` pairs.
        include_bt4: Add one BT4-designed member per group. Imported lazily and only
            here, so the panel builds without a working optimizer when it is off.

    Returns:
        The members, in file order, with any BT4 designs appended per group.
    """
    members: list[DesignedMember] = []
    for header, sequence in records:
        group, member_id, role = _parse_header(header)
        members.append(
            DesignedMember(
                group=group,
                member_id=member_id,
                role=role,
                cds=sequence.upper(),
                note=header,
            )
        )

    if not include_bt4:
        return members

    from bt4 import api
    from bt4.domain.genetic_code import translate

    for group in sorted({m.group for m in members}):
        native = next(m for m in members if m.group == group and m.role == "native")
        protein = translate(native.cds)
        # BT4 emits the stop itself, so hand it the protein without one.
        result = api.optimize(protein.rstrip("*"))
        members.append(
            DesignedMember(
                group=group,
                member_id=f"{group}_BT4",
                role="designed",
                cds=result.dna,
                note=f"bt4.api.optimize, {len(result.dna)} nt, default organism",
            )
        )
    return members


def main(argv: Sequence[str] | None = None) -> int:
    """Build the panel and write it."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--fasta", default=str(DEFAULT_FASTA),
        help="grouped CDS FASTA, headers 'Protein|Source|run' (default: Ranaghan 2021)",
    )
    parser.add_argument("--out", required=True, help="where to write the panel TSV")
    parser.add_argument(
        "--include-bt4", action="store_true",
        help="add BT4's own design per protein. OFF by default: the panel asks how "
             "splice models behave on designed CDS, and generating the designs with the "
             "tool that would consume the answer makes it partly a fact about BT4",
    )
    args = parser.parse_args(argv)

    records = read_fasta(Path(args.fasta))
    print(f"read {len(records)} records from {args.fasta}")

    members = build_members(records, include_bt4=args.include_bt4)
    provenance = PROVENANCE
    if args.include_bt4:
        provenance += " Plus one bt4.api.optimize design per protein, added locally."
    panel = panel_from_members(members, source=args.out, provenance=provenance)
    write_designed_cds_panel(members, Path(args.out))

    print(f"\nwrote {args.out}")
    print(f"  members {len(panel)}   groups {list(panel.groups)}")
    for group in panel.groups:
        native = panel.native(group)
        print(
            f"    {group:12s} {len(panel.designed(group)):3d} designs, "
            f"native {len(native.cds):5d} nt, protein {len(native.protein)} aa"
        )
    print(f"  content_hash {panel.content_hash()}")
    print(
        "\nThis panel carries NO splice labels, and cannot: designed coding sequence has "
        "\nno splice ground truth. It supports label-free measurements only -- response "
        "\nto synonymous change, cross-backend agreement, and Δ against the native."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
