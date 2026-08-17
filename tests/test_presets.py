"""Tests for the application presets (:mod:`bt4.pipeline.presets`).

Presets are *data*, and the honesty rules that make them safe are testable: no
preset is applied by default, every preset names real config fields, a preset key
reaches the manifest, and -- the load-bearing one -- every bundled preset actually
produces a feasible run on a realistic protein panel rather than handing the user
an infeasible configuration.
"""

from __future__ import annotations

import pytest

from bt4 import api
from bt4.domain.genetic_code import translate
from bt4.pipeline.optimize import OptimizeConfig
from bt4.pipeline.presets import (
    APPLICATION_PRESETS,
    apply_preset,
    available_presets,
    resolve_preset,
)

# A small but realistic panel: a fluorescent-protein N-terminus, a repeat-heavy
# tagged construct (linkers + His tag + FLAG), and a hydrophobic stretch.
_PANEL = {
    "gfp_like": "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF",
    "linker_rich": "M" + "EAAAK" * 4 + "GGGGSGGGGS" + "HHHHHH" + "DYKDDDDK",
    "hydrophobic": "MLLLVVAAIIGGFFWWYYMLLVVAAIIGGFFWWYY",
}


def test_no_preset_is_applied_by_default() -> None:
    """BT4 stays regime-agnostic: a default config carries no preset."""
    assert OptimizeConfig().application_preset == ""


def test_every_preset_sets_only_real_config_fields() -> None:
    fields = set(OptimizeConfig.__dataclass_fields__)
    for preset in APPLICATION_PRESETS:
        unknown = set(preset.overrides) - fields
        assert not unknown, f"{preset.key} sets unknown field(s): {sorted(unknown)}"


def test_presets_have_distinct_keys_and_describe_themselves() -> None:
    keys = [p.key for p in APPLICATION_PRESETS]
    assert len(keys) == len(set(keys))
    for preset in APPLICATION_PRESETS:
        assert preset.label and preset.description and preset.rationale
        assert preset.regime in {"vector", "plasmid", "mrna", "synthesis"}


def test_unknown_preset_raises_and_names_the_known_keys() -> None:
    # It must refuse, never quietly fall back -- there is no default preset.
    with pytest.raises(KeyError) as exc_info:
        resolve_preset("no_such_preset")
    message = str(exc_info.value)
    assert "synthesis" in message


def test_apply_preset_records_the_key_for_provenance() -> None:
    config = apply_preset("aav_transgene")
    assert config.application_preset == "aav_transgene"
    assert config.cpg_mode == "deplete"
    assert config.avoid_splice_sites is True


def test_apply_preset_does_not_mutate_the_caller_config() -> None:
    base = OptimizeConfig(organism="escherichia_coli")
    applied = apply_preset("synthesis", base)
    assert base.gc_window_nt is None
    assert applied.gc_window_nt == 50
    # Fields the preset does not set are carried through untouched.
    assert applied.organism == "escherichia_coli"


def test_caller_overrides_beat_the_preset() -> None:
    from dataclasses import replace

    config = replace(apply_preset("synthesis"), max_homopolymer=4)
    assert config.max_homopolymer == 4  # a preset is a starting point, not a cage


def test_preset_key_changes_the_manifest() -> None:
    # Invariant #9: two runs differing only by preset must not stamp the same
    # provenance, even if they happened to deliver the same sequence.
    plain = api.optimize("MAALKHETQWY", OptimizeConfig(max_homopolymer=6))
    preset = api.optimize(
        "MAALKHETQWY", OptimizeConfig(max_homopolymer=6, application_preset="synthesis")
    )
    assert plain.audit["manifest"] != preset.audit["manifest"]


def test_no_preset_enables_the_infeasible_internal_start_rule() -> None:
    """`avoid_internal_start` fails on most real proteins, so no preset ships it.

    Bundling it would hand nearly every user a RELAXED certificate they never
    asked for; it stays an explicit opt-in.
    """
    for preset in APPLICATION_PRESETS:
        assert preset.overrides.get("avoid_internal_start") is not True


@pytest.mark.parametrize("key", [p.key for p in APPLICATION_PRESETS])
@pytest.mark.parametrize("protein_name", sorted(_PANEL))
def test_every_preset_is_feasible_on_the_panel(key: str, protein_name: str) -> None:
    """A shipped preset must produce a real result, not an infeasible config.

    The refinement budget is trimmed here: the question is whether the preset is
    *satisfiable and honestly reported*, not how close annealing gets, and the
    audit assertions below hold at any iteration count.
    """
    from dataclasses import replace

    protein = _PANEL[protein_name]
    result = api.optimize(protein, replace(apply_preset(key), refine_iterations=250))
    # Invariant #1 holds regardless of preset.
    assert translate(result.dna) == protein + "*"
    # Any rule the preset could not fully enforce is disclosed, never hidden.
    for audit_key, value in result.audit.items():
        if audit_key.endswith("_enforced"):
            assert value in {"clean", "partial"}
            residual = result.audit[audit_key.replace("_enforced", "_residual")]
            assert (value == "clean") == (residual == 0)


def test_available_presets_is_stable_and_public() -> None:
    assert available_presets() == APPLICATION_PRESETS
    assert api.available_presets() == APPLICATION_PRESETS


def test_no_preset_creates_an_intractable_trellis_context() -> None:
    """A preset must not silently make a run take minutes.

    The exact DP keys its state on the literal trailing context, so its state space
    grows roughly exponentially in the widest ``context_len`` among the LOCAL
    constraints. This bites in practice: an ``inverted_stem=12, inverted_loop=4``
    hairpin rule (context 27) was measured at 115-339 s per run, versus 1.6-8.4 s
    for every other preset. Presets therefore keep the trellis context small and
    push wide rules (windowed GC, dispersed repeats) to the refinement layer, where
    cost is linear and the certificate says so.

    MUTATION THAT MUST FAIL THIS: put ``inverted_stem=12, inverted_loop=4`` back in
    the plasmid preset.
    """
    from bt4.pipeline.optimize import _build_constraints

    for preset in APPLICATION_PRESETS:
        config = apply_preset(preset.key)
        widest = max((c.context_len() for c in _build_constraints(config)), default=0)
        assert widest <= 16, (
            f"preset {preset.key} builds a trellis context of {widest} nt, wide "
            "enough to make the exact DP intractable"
        )


def test_every_preset_field_has_a_cli_flag_that_overrides_it() -> None:
    """A preset must never set something the CLI cannot then override.

    MUTATION THAT MUST FAIL THIS: add a field to a preset's ``overrides`` without
    adding it to ``_PRESET_FIELD_TO_FLAG``; the preset would then silently pin a
    value the user has no way to change from the command line.
    """
    from bt4.cli.__main__ import _PRESET_FIELD_TO_FLAG

    used = {field for preset in APPLICATION_PRESETS for field in preset.overrides}
    missing = sorted(used - set(_PRESET_FIELD_TO_FLAG))
    assert not missing, f"preset field(s) with no overriding CLI flag: {missing}"


def test_cli_preset_supplies_values_but_explicit_flags_win() -> None:
    from bt4.cli.__main__ import _apply_preset_to_args, _parser

    parser = _parser()
    argv = ["optimize", "MAALKHETQWY", "--preset", "synthesis", "--max-homopolymer", "4"]
    args = parser.parse_args(argv)
    _apply_preset_to_args(args, argv)
    assert args.gc_window == 50  # supplied by the preset
    assert args.max_homopolymer == 4  # explicit flag beats the preset's 6

    argv2 = ["optimize", "MAALKHETQWY", "--preset", "synthesis"]
    args2 = parser.parse_args(argv2)
    _apply_preset_to_args(args2, argv2)
    assert args2.max_homopolymer == 6  # not named -> preset supplies it

    with pytest.raises(KeyError):
        bad = parser.parse_args(["optimize", "MA", "--preset", "nope"])
        _apply_preset_to_args(bad, ["optimize", "MA", "--preset", "nope"])


def test_cli_without_preset_is_unchanged() -> None:
    """No preset named -> nothing is applied (there is no default preset)."""
    from bt4.cli.__main__ import _build_config, _parser

    args = _parser().parse_args(["optimize", "MAALKHETQWY"])
    config = _build_config(args)
    assert config.application_preset == ""
    assert config.gc_window_nt is None
    assert config.max_repeat_length is None
