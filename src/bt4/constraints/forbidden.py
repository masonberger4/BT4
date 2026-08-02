"""Named forbidden-sequence presets (CLAUDE.md §6, §6.6).

BT3's ``almost-there`` branch shipped a flat ``DEFAULT_FORBIDDEN_LIST`` of motifs
its optimizer banned. BT4 keeps that capability but makes it *legible*: instead of
one opaque list, it offers a small catalog of **named, documented presets**, each
a group of motifs with a plain-language description. BT4 Studio renders one
checkbox per preset (with the description as its hover tooltip); the CLI takes
their keys by name. Every preset resolves to literal DNA motifs that are enforced
*exactly* by :class:`~bt4.constraints.rules.ForbiddenMotifConstraint` (with
reverse-complement banning), so nothing here is an unvalidated claim -- only a
convenient, attributable grouping of motifs a user may choose to forbid.

Users who need something outside the catalog pass their own motifs directly
(``forbidden_motifs`` in the config / the CLI ``--forbid`` flag / the Studio
"custom motifs" field); the presets are additive with those.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "FORBIDDEN_PRESETS",
    "ForbiddenPreset",
    "available_forbidden_presets",
    "resolve_forbidden_motifs",
]


@dataclass(frozen=True, slots=True)
class ForbiddenPreset:
    """A named group of forbidden DNA motifs.

    Attributes:
        key: Stable identifier used by the CLI and config (e.g. ``"poly_a_signal"``).
        label: Short human-readable name for the UI checkbox.
        description: One-line explanation shown as the control's hover tooltip.
        motifs: The literal DNA motifs this preset forbids (ACGT, upper-case).
    """

    key: str
    label: str
    description: str
    motifs: tuple[str, ...] = field(default=())


FORBIDDEN_PRESETS: tuple[ForbiddenPreset, ...] = (
    ForbiddenPreset(
        key="poly_a_signal",
        label="Poly-A signals",
        description=(
            "Canonical (AATAAA) and common variant (ATTAAA) polyadenylation "
            "signals -- avoid premature 3'-end cleavage/polyadenylation of the mRNA."
        ),
        motifs=("AATAAA", "ATTAAA"),
    ),
    ForbiddenPreset(
        key="tata_box",
        label="TATA box",
        description=(
            "Core TATA-box promoter element (TATAAA) -- avoid a cryptic internal "
            "promoter inside the coding sequence."
        ),
        motifs=("TATAAA",),
    ),
    ForbiddenPreset(
        key="telomere_repeat",
        label="Telomere repeat",
        description=(
            "Human telomeric repeat unit (TTAGGG) -- a common unwanted motif in "
            "synthetic constructs."
        ),
        motifs=("TTAGGG",),
    ),
    ForbiddenPreset(
        key="bt3_synthesis_artifacts",
        label="BT3 synthesis artifacts",
        description=(
            "Cloning/synthesis-artifact motifs inherited from BT3's default "
            "forbidden-sequence set."
        ),
        motifs=(
            "GCTGGTGG",
            "GTTGTAAC",
            "TTATCCACA",
            "GCCGTCTGAA",
            "AAGTGCGGT",
            "ACAAGCGGTC",
        ),
    ),
)
"""The ordered catalog of forbidden-sequence presets (stable display order)."""

_BY_KEY: dict[str, ForbiddenPreset] = {p.key: p for p in FORBIDDEN_PRESETS}


def available_forbidden_presets() -> tuple[ForbiddenPreset, ...]:
    """Return the forbidden-sequence presets in their stable catalog order."""
    return FORBIDDEN_PRESETS


def resolve_forbidden_motifs(keys: tuple[str, ...]) -> tuple[str, ...]:
    """Return the de-duplicated motifs for the named presets, in catalog order.

    Args:
        keys: Preset keys to resolve (e.g. ``("poly_a_signal", "tata_box")``).

    Returns:
        The union of the selected presets' motifs, de-duplicated, ordered by the
        catalog then by motif order within each preset.

    Raises:
        KeyError: If any key is not a known preset.
    """
    wanted = set(keys)
    unknown = wanted - _BY_KEY.keys()
    if unknown:
        known = ", ".join(sorted(_BY_KEY))
        raise KeyError(f"unknown forbidden preset(s): {sorted(unknown)}; known: {known}")
    seen: set[str] = set()
    out: list[str] = []
    for preset in FORBIDDEN_PRESETS:
        if preset.key in wanted:
            for motif in preset.motifs:
                if motif not in seen:
                    seen.add(motif)
                    out.append(motif)
    return tuple(out)
