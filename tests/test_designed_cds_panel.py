"""Tests for the designed synonymous CDS panel — BT4's own regime.

Every splice measurement before this one was *recall on natural sites in natural
genes*. BT4 emits neither: it emits a synonymous re-encoding of one protein's coding
sequence. This panel is that regime, and its most important property is what it
**refuses** to claim.

There is no splice ground truth for designed coding sequence — nothing here has been
assayed, none of it is annotated, and a motif is not a site. So the format carries no
label column, the reader refuses one, and nothing built on it reports `passed`.

What is tested here is that the panel's *defining* property is enforced rather than
assumed: within a group every member must encode the same protein. A panel whose
members are not synonymous is not a weaker synonymous panel, it is a different
experiment, and every number drawn from it would be about the wrong thing.
"""

from __future__ import annotations

import pathlib

import pytest

from bt4.biomodels.splice.designed_panel import (
    PANEL_COLUMNS,
    DesignedMember,
    panel_from_members,
    read_designed_cds_panel,
    write_designed_cds_panel,
)

# Two synonymous encodings of Met-Lys-stop, and one that is NOT synonymous.
_NATIVE = "ATGAAATAA"
_SYNONYMOUS = "ATGAAGTAA"  # AAA -> AAG, both Lys
_DIFFERENT = "ATGGGGTAA"  # Lys -> Gly


def _members() -> list[DesignedMember]:
    return [
        DesignedMember("P", "P_native", "native", _NATIVE, "the real CDS"),
        DesignedMember("P", "P_v1", "designed", _SYNONYMOUS, "vendor A"),
    ]


def test_a_valid_panel_round_trips(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "designed.tsv"
    write_designed_cds_panel(_members(), path)
    panel = read_designed_cds_panel(path, provenance="test fixture")
    assert len(panel) == 2
    assert panel.groups == ("P",)
    assert panel.native("P").member_id == "P_native"
    assert [m.member_id for m in panel.designed("P")] == ["P_v1"]


def test_the_panel_declares_that_it_has_no_labels() -> None:
    """The honest field, and the one a reader most needs."""
    panel = panel_from_members(_members(), provenance="test fixture")
    assert panel.describe()["carries_splice_labels"] is False
    assert "label" not in PANEL_COLUMNS


def test_a_label_column_is_refused_by_name(tmp_path: pathlib.Path) -> None:
    """Someone adding one believed this panel has ground truth. It does not."""
    path = tmp_path / "labelled.tsv"
    path.write_text(
        "group\tmember_id\trole\tcds\tlabel\n"
        f"P\tP_native\tnative\t{_NATIVE}\t1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="NO label column"):
        read_designed_cds_panel(path, provenance="test fixture")


def test_a_non_synonymous_group_is_refused() -> None:
    """The defining property, verified rather than trusted."""
    members = [
        DesignedMember("P", "P_native", "native", _NATIVE),
        DesignedMember("P", "P_bad", "designed", _DIFFERENT),
    ]
    with pytest.raises(ValueError, match="not synonymous"):
        panel_from_members(members, provenance="test fixture")


def test_the_refusal_names_the_offending_member() -> None:
    """A 93-member panel is not debuggable from 'a group is wrong'."""
    members = [
        DesignedMember("P", "P_native", "native", _NATIVE),
        DesignedMember("P", "P_ok", "designed", _SYNONYMOUS),
        DesignedMember("P", "P_bad", "designed", _DIFFERENT),
    ]
    with pytest.raises(ValueError, match="P_bad"):
        panel_from_members(members, provenance="test fixture")


@pytest.mark.parametrize("count", [0, 2])
def test_a_group_needs_exactly_one_native(count: int) -> None:
    """It is the reference every Δ in the group is measured against."""
    members = [DesignedMember("P", "P_v1", "designed", _SYNONYMOUS)]
    members += [
        DesignedMember("P", f"P_native{i}", "native", _NATIVE) for i in range(count)
    ]
    with pytest.raises(ValueError, match="native"):
        panel_from_members(members, provenance="test fixture")


def test_a_sequence_that_is_not_whole_codons_is_refused() -> None:
    members = [DesignedMember("P", "P_native", "native", "ATGAAAT")]
    with pytest.raises(ValueError, match="whole number of codons"):
        panel_from_members(members, provenance="test fixture")


def test_provenance_is_required() -> None:
    """A designed-CDS panel is uninterpretable without knowing who designed it.

    Unlike a label, this one *is* knowable, so it is demanded rather than omitted.
    """
    with pytest.raises(ValueError, match="provenance is required"):
        panel_from_members(_members(), provenance="   ")


def test_content_hash_is_order_independent_and_moves_on_a_real_change() -> None:
    """Invariant #7: the same panel hashes the same however it is ordered."""
    forward = panel_from_members(_members(), provenance="x")
    backward = panel_from_members(list(reversed(_members())), provenance="x")
    assert forward.content_hash() == backward.content_hash()

    changed = panel_from_members(
        [
            DesignedMember("P", "P_native", "native", _NATIVE),
            DesignedMember("P", "P_v1", "designed", "ATGAAGTAG"),  # different stop
        ],
        provenance="x",
    )
    assert changed.content_hash() != forward.content_hash()


def test_the_bundled_ranaghan_panel_is_genuinely_synonymous() -> None:
    """The real panel, built from the committed CC BY 4.0 FASTA.

    Guards the fixture that matters: 3 proteins, each the native CDS plus 30 designs
    from three anonymized commercial optimizers over ten repeat runs.
    """
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(
        "bt4_make_designed_cds_panel",
        pathlib.Path(__file__).parent.parent / "scripts" / "make_designed_cds_panel.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["bt4_make_designed_cds_panel"] = module
    spec.loader.exec_module(module)

    records = module.read_fasta(module.DEFAULT_FASTA)
    assert len(records) == 93
    # panel_from_members enforces synonymy, so building it IS the assertion.
    panel = panel_from_members(
        module.build_members(records), provenance=module.PROVENANCE
    )
    assert panel.groups == ("Beclin1", "KRas4B", "PDE3A")
    for group in panel.groups:
        assert len(panel.designed(group)) == 30
        proteins = {m.protein for m in panel.group_members(group)}
        assert len(proteins) == 1, group


def test_the_probe_reports_progress_and_is_silent_by_default() -> None:
    """A silent multi-minute run is indistinguishable from a hang.

    Fixed once for `splice-gate` and then reintroduced here: the probe accepted a
    progress callback and the CLI did not pass one. Per *group* rather than per
    sequence, because the probe hands a whole group to `backend_agreement`, which
    scores its members internally — the group is the finest boundary this layer can
    honestly report.
    """
    from bt4.pipeline.splice_gate import probe_designed_cds

    panel = panel_from_members(_members(), provenance="test fixture")
    seen: list[tuple[int, int, str, int]] = []
    probe_designed_cds(panel, ("pwm",), progress=lambda *call: seen.append(call))
    assert [entry[2] for entry in seen] == list(panel.groups)
    assert all(entry[1] == len(panel.groups) for entry in seen)
    assert [entry[0] for entry in seen] == list(range(1, len(panel.groups) + 1))

    # The API default stays print-free and callback-free (section 3: only `cli` prints).
    probe_designed_cds(panel, ("pwm",))


def test_the_probe_reports_the_backends_own_names_not_the_aliases() -> None:
    """`pwm` is a registry alias; `consensus-pwm-baseline` is what ran.

    `backend_agreement` keys its results by the predictor's own name, so looking them
    up by alias raised KeyError — found by running the probe. The resolved names also
    carry configuration an alias loses, such as Pangolin's tissue set.
    """
    from bt4.pipeline.splice_gate import probe_designed_cds

    panel = panel_from_members(_members(), provenance="test fixture")
    probe = probe_designed_cds(panel, ("pwm",))[0]
    assert probe.backends == ("consensus-pwm-baseline",)
    assert set(probe.delta_spread) == {"consensus-pwm-baseline"}
    assert set(probe.delta_range) == {"consensus-pwm-baseline"}


def test_the_probe_has_no_pass_or_promotable_field() -> None:
    """Structural: this panel has no labels, so nothing built on it may render a verdict."""
    from bt4.pipeline.splice_gate import DesignedCdsProbe

    fields = set(DesignedCdsProbe.__dataclass_fields__)
    assert not fields & {"passed", "promotable", "thresholds_declared", "reasons"}
