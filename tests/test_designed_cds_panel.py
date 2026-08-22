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


def test_zero_and_nearly_zero_stay_distinguishable_in_the_report() -> None:
    """`%.4f` printed an exactly-zero Δ spread and a 1e-7 one identically.

    On the first real run Pangolin's spread displayed as `0.0000` for two of three
    proteins, and the two readings that display permits are opposite: exactly zero
    means the backend gave the native and all 30 designs the same pooled risk, while
    1e-7 means it separates them and the signal is merely tiny. Rounding that away
    decides the reader's conclusion for them.

    **Correction.** This docstring used to finish "…and cannot rank them at all", and
    that inference was wrong — it is what the exactly-zero reading was then used to
    conclude. An exactly-zero *pooled-risk* spread does not imply the backend cannot
    rank: measured, Pangolin's scores varied more than twofold across those designs and
    BT4's pooling hinge floored all of it, because no position reached the 0.5
    background. Keeping the digits distinguishable is still right; the conclusion drawn
    from them needed `DesignedCdsProbe.degenerate` to be sound.
    """
    from bt4.cli.__main__ import _signal

    assert _signal(0.0) == "0"
    assert _signal(0.0, signed=True) == "+0"
    assert _signal(1e-7) == "1.00e-07"
    assert _signal(-3.2e-6, signed=True) == "-3.20e-06"
    # Ordinary magnitudes keep the readable fixed-point form.
    assert _signal(0.8349) == "0.8349"
    assert _signal(-1.0885, signed=True) == "-1.0885"


# --------------------------------------------------------------------------
# A zero spread must say WHICH zero it is


def test_the_probe_carries_a_background_free_response_beside_the_risk() -> None:
    """The risk pooling is degenerate in this panel's regime; the response is not.

    Measured against the hash-verified Pangolin weights, only 6 of the 93 sequences
    carried any position above the 0.5 background — all six designs of a single protein —
    so for the other two proteins every risk Δ was exactly zero. Without a second,
    background-free statistic the probe's central question — does a synonymous change move
    the score at all — is unanswerable for those groups by construction.
    """
    from bt4.pipeline.splice_gate import probe_designed_cds

    panel = panel_from_members(_members(), provenance="test fixture")
    probe = probe_designed_cds(panel, ("pwm",))[0]
    for field in ("response_spread", "response_range", "sub_background", "max_score"):
        assert set(getattr(probe, field)) == {"consensus-pwm-baseline"}
    assert isinstance(probe.response_sign_agreement, float)


def test_degenerate_distinguishes_a_floored_zero_from_a_measured_one() -> None:
    """The distinction the probe's whole conclusion rests on."""
    from bt4.pipeline.splice_gate import DesignedCdsProbe

    def _probe(sub_background: int) -> DesignedCdsProbe:
        return DesignedCdsProbe(
            group="G",
            n_designs=30,
            backends=("b",),
            delta_spread={"b": 0.0},
            delta_range={"b": (0.0, 0.0)},
            rank_correlations={},
            sign_agreement=1.0,
            response_spread={"b": 3.9},
            response_range={"b": (-2.1, 1.7)},
            response_rank_correlations={},
            response_sign_agreement=1.0,
            sub_background={"b": sub_background},
            max_score={"b": 0.323},
        )

    # Every one of the 31 sequences floored => the zero is by construction.
    assert _probe(31).degenerate("b") is True
    # Even one sequence clearing background means the zero is a real measurement.
    assert _probe(30).degenerate("b") is False
    assert _probe(0).degenerate("b") is False


def test_the_cli_never_calls_a_floored_zero_a_failure_to_rank() -> None:
    """Structural: the sentence that turned BT4's silence into the model's.

    `designed-probe` printed "This backend cannot rank these candidates at all" on any
    exactly-zero spread. On the panel it was written for, that spread is produced by
    BT4's own pooling hinge for every backend whose scores stay under 0.5 — so the line
    reported a property of `pool_log_odds` as a property of the CNN.
    """
    from pathlib import Path as _Path

    source = (
        _Path(__file__).resolve().parent.parent / "src" / "bt4" / "cli" / "__main__.py"
    ).read_text(encoding="utf-8")
    assert "cannot rank these candidates at all" not in source
    assert "ZERO BY CONSTRUCTION" in source
