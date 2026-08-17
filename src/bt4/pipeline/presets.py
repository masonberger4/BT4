"""Application presets: named starting points for a design run (CLAUDE.md §6).

A first-time user faces ~30 independent knobs and no guidance about which values
suit *their* construct. These presets are the guidance -- but they are deliberately
**data, not engine behaviour**: each is a small mapping of
:class:`~bt4.pipeline.optimize.OptimizeConfig` field names to values, applied by
:func:`apply_preset` and freely overridable afterwards. The solver has no idea a
preset was used.

Three honesty rules shape this module.

* **No preset is applied by default** (a maintainer decision). BT4 stays
  regime-agnostic: the CpG direction, the length ceilings and the structure
  handling that suit an AAV transgene are wrong for an IVT mRNA, and guessing
  which the user meant would be exactly the kind of unstated assumption §10.6
  forbids. A run uses a preset only when it is asked for by name.
* **A preset is a *convention*, not a validated prediction.** The values below are
  drawn from published synthesis-vendor complexity guidance and from standard
  vector-design practice. They are a defensible place to *start*; none of them is
  a calibrated claim that a given sequence will express better, and nothing here
  changes what BT4 reports (§10.6, §10.7).
* **Presets never enable a rule that cannot be satisfied.** In particular none of
  them sets ``avoid_internal_start``: it is genuinely infeasible on most real
  proteins (an internal Met's Kozak context can be synonymously forced), so
  bundling it would hand almost every user a ``RELAXED`` certificate they did not
  ask for. It stays an explicit opt-in, where the relaxation path
  (:mod:`bt4.domain.relax`) makes it honest.

Note which certificate a preset implies: presets that set ``max_repeat_length`` or
``avoid_uorf`` engage a **non-local** rule, so a run whose exact-DP seed violates it
is refined and reports ``HEURISTIC``; the purely LOCAL rules (windowed GC, GC run,
homopolymer, splice motifs) stay exact in the trellis and keep ``PROVEN_OPTIMAL``.

This module depends only on its sibling :mod:`bt4.pipeline.optimize` and the
standard library.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType

from bt4.pipeline.optimize import OptimizeConfig

__all__ = [
    "APPLICATION_PRESETS",
    "ApplicationPreset",
    "apply_preset",
    "available_presets",
    "resolve_preset",
]


@dataclass(frozen=True, slots=True)
class ApplicationPreset:
    """One named starting point for a design run.

    Attributes:
        key: Stable identifier used on the CLI and in the manifest.
        label: Short human-readable name.
        regime: The construct class this preset targets (``"vector"``,
            ``"plasmid"``, ``"mrna"``, or ``"synthesis"``). Recorded so a user can
            see at a glance that regimes differ; BT4 itself stays regime-agnostic.
        description: What the preset is for, in one sentence.
        rationale: Why these particular values -- the provenance of the numbers.
        overrides: :class:`~bt4.pipeline.optimize.OptimizeConfig` field names mapped
            to the values this preset sets. Every key must be a real field.
    """

    key: str
    label: str
    regime: str
    description: str
    rationale: str
    overrides: Mapping[str, object]


# Vendor-style synthesis manufacturability, shared by most presets below. GC is
# bounded per 50-nt window (a *local* GC extreme is what fails synthesis, which a
# whole-sequence GC number cannot see), plus caps on homopolymers, GC runs, and
# dispersed repeats. These mirror the complexity screens IDT and Twist publish;
# they are conventions, not measured thresholds for any one vendor's process.
_SYNTHESIS: dict[str, object] = {
    "gc_window_nt": 50,
    "gc_window_min": 0.25,
    "gc_window_max": 0.65,
    "max_homopolymer": 6,
    "max_gc_run": 6,
    "max_repeat_length": 20,
}


APPLICATION_PRESETS: tuple[ApplicationPreset, ...] = (
    ApplicationPreset(
        key="synthesis",
        label="Gene synthesis (vendor-style)",
        regime="synthesis",
        description=(
            "Manufacturability only: keeps a sequence inside the complexity limits "
            "commercial gene-synthesis vendors screen for."
        ),
        rationale=(
            "Windowed GC (25-65% per 50 nt) rather than a whole-sequence GC number, "
            "because local GC extremes are what break oligo synthesis and assembly; "
            "plus homopolymer, GC-run and dispersed-repeat caps. No expression claim."
        ),
        overrides=dict(_SYNTHESIS),
    ),
    ApplicationPreset(
        key="aav_transgene",
        label="AAV transgene",
        regime="vector",
        description=(
            "AAV-packaged transgene: synthesis limits plus CpG depletion and "
            "avoidance of strong cryptic splice motifs and uORFs."
        ),
        rationale=(
            "CpG depletion is standard practice for AAV transgenes (unmethylated CpG "
            "is sensed by TLR9 and is associated with transgene silencing and "
            "inflammatory response). Strong splice-consensus motifs and out-of-frame "
            "uORFs are avoided structurally. NOTE: BT4 optimizes the CDS only -- it "
            "does not check the ~4.7 kb AAV packaging limit, which depends on the "
            "whole construct BT4 cannot see."
        ),
        overrides={
            **_SYNTHESIS,
            "cpg_weight": 1.0,
            "cpg_mode": "deplete",
            "avoid_splice_sites": True,
            "avoid_polya": True,
            "avoid_uorf": True,
        },
    ),
    ApplicationPreset(
        key="lvv_transgene",
        label="Lentiviral (LVV) transgene",
        regime="vector",
        description=(
            "Lentivirally-delivered transgene: synthesis limits plus the strongest "
            "available cryptic-splice avoidance and tight repeat limits."
        ),
        rationale=(
            "Cryptic splicing inside a lentiviral genome is a documented clinical "
            "failure mode (aberrant splicing producing truncated transcripts), so the "
            "structural splice-motif rule is on and repeats are capped tighter than "
            "the synthesis default (long repeats also promote recombination during "
            "vector production). The real splice audit is the wrapped SpliceAI / "
            "Pangolin pass, which is out-of-loop and currently UNCALIBRATED."
        ),
        overrides={
            **_SYNTHESIS,
            "max_repeat_length": 16,
            "avoid_splice_sites": True,
            "avoid_polya": True,
            "avoid_uorf": True,
        },
    ),
    ApplicationPreset(
        key="plasmid_production",
        label="Plasmid production",
        regime="plasmid",
        description=(
            "Plasmid propagated in E. coli: synthesis limits plus tighter repeat and "
            "homopolymer caps for replication stability."
        ),
        rationale=(
            "Long direct and inverted repeats promote recombination and plasmid "
            "instability during propagation, so the dispersed-repeat cap is tightened "
            "and homopolymers are held shorter. The cap is reverse-complement aware, "
            "so it already covers inverted and palindromic repeats -- the separate "
            "LOCAL hairpin rule (inverted_stem) is deliberately NOT set here: a stem "
            "long enough to matter needs a trellis context wide enough to make the "
            "exact DP intractable (measured: minutes per run), and it would be "
            "redundant with the cap. Set --organism escherichia_coli as well if the "
            "protein is expressed in E. coli rather than merely cloned."
        ),
        overrides={
            **_SYNTHESIS,
            "max_repeat_length": 16,
            "max_homopolymer": 5,
        },
    ),
    ApplicationPreset(
        key="mrna_ivt",
        label="IVT mRNA",
        regime="mrna",
        description=(
            "In-vitro-transcribed mRNA: synthesis limits plus 5' structure "
            "refinement and uORF avoidance."
        ),
        rationale=(
            "An IVT transcript never passes through the spliceosome, so the cryptic-"
            "splice rules are deliberately NOT enabled here; the relevant levers are "
            "start-proximal structure (refinement) and uORFs that divert scanning "
            "ribosomes. CpG direction is left unset on purpose: depletion suits a "
            "stealth therapeutic while elevation suits an immunostimulatory vaccine, "
            "and BT4 will not guess which you are building."
        ),
        overrides={
            **_SYNTHESIS,
            "avoid_uorf": True,
            "refine": True,
        },
    ),
)


_BY_KEY: Mapping[str, ApplicationPreset] = MappingProxyType(
    {preset.key: preset for preset in APPLICATION_PRESETS}
)


def available_presets() -> tuple[ApplicationPreset, ...]:
    """Return every bundled application preset, in a stable order."""
    return APPLICATION_PRESETS


def resolve_preset(key: str) -> ApplicationPreset:
    """Return the preset named ``key``.

    Args:
        key: A preset key (see :func:`available_presets`).

    Raises:
        KeyError: If ``key`` is unknown -- the message names the known keys rather
            than silently falling back to a default (there is no default preset).
    """
    preset = _BY_KEY.get(key)
    if preset is None:
        known = ", ".join(sorted(_BY_KEY))
        raise KeyError(f"unknown preset {key!r}; known presets: {known}")
    return preset


def apply_preset(key: str, config: OptimizeConfig | None = None) -> OptimizeConfig:
    """Return ``config`` with preset ``key``'s fields applied.

    The preset supplies a starting point; any field the caller sets *after* this
    call wins, so a preset is never a cage. The preset key is recorded on the
    returned config's ``application_preset`` field so it reaches the run manifest
    (invariant #9): two runs that differ only by preset must not stamp the same
    provenance.

    Args:
        key: A preset key (see :func:`available_presets`).
        config: Base configuration; a default :class:`OptimizeConfig` when omitted.

    Returns:
        A new :class:`OptimizeConfig` -- the input is never mutated.

    Raises:
        KeyError: If ``key`` is unknown.
        TypeError: If the preset names a field the config does not have (guards a
            typo in the preset table from silently doing nothing).
    """
    preset = resolve_preset(key)
    base = config if config is not None else OptimizeConfig()
    fields = {f for f in OptimizeConfig.__dataclass_fields__}
    unknown = sorted(set(preset.overrides) - fields)
    if unknown:
        raise TypeError(
            f"preset {key!r} sets unknown OptimizeConfig field(s): {', '.join(unknown)}"
        )
    return replace(base, **preset.overrides, application_preset=key)  # type: ignore[arg-type]
