"""Tests for the annotated splice-panel format.

A splice panel has exactly one catastrophic failure mode, and it is silent: annotate
the sites one base off and every score in the gate is misaligned, the model looks
incompetent, and nothing in the numbers says why. The conventions genuinely differ --
an annotation-derived table names the *exonic boundary* bases, BT4's backends anchor on
the *intronic* dinucleotide -- so the format pins one and verifies it against the
sequence itself.

These tests pin:

* the convention is the one BT4's own PWM baseline uses, so it is verifiable from this
  repository rather than assumed about someone else's;
* a panel built to the other convention is **refused**, with a message naming the exact
  shift that would have worked;
* the negative construction is mandatory, because PR-AUC's floor is the prevalence;
* overlap with the models' training chromosomes is reported, not silently scored;
* the content hash is order-independent and moves on any real change (invariant #7).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bt4.biomodels.splice import ConsensusPwmSplicePredictor
from bt4.biomodels.splice.panel import (
    MIN_MOTIF_CONSISTENCY,
    SplicePanel,
    SpliceWindow,
    canonical_motif_at,
    panel_from_windows,
    read_splice_panel,
)

_NEG = "all other positions in the same gene bodies"

# One intron: exon ... CAG | GTAAGT ... polypyrimidine ... CAG | G ... exon
# Donor  = the G of GT (index 9). Acceptor = the G of AG (index 39).
_SEQ = "CCCCCCCAG" "GTAAGT" "CCCCCCCCCCCCCTTTTTTTTTCAG" "GCCCCCCCCC"
_DONOR = 9
_ACCEPTOR = 39


def _window(window_id: str = "w0", group: str = "chr1") -> SpliceWindow:
    return SpliceWindow(window_id, group, _SEQ, donors=(_DONOR,), acceptors=(_ACCEPTOR,))


def _panel(*windows: SpliceWindow, **kwargs: object) -> SplicePanel:
    return panel_from_windows(
        windows or (_window(),), negative_construction=_NEG, **kwargs
    )  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# The convention, and that it is BT4's own


def test_the_fixture_is_anchored_the_way_the_format_says() -> None:
    """Sanity: the shared sequence really does carry GT/AG at the declared bases."""
    assert _SEQ[_DONOR : _DONOR + 2] == "GT"
    assert _SEQ[_ACCEPTOR - 1 : _ACCEPTOR + 1] == "AG"
    assert canonical_motif_at(_SEQ, _DONOR, "donor")
    assert canonical_motif_at(_SEQ, _ACCEPTOR, "acceptor")


def test_the_convention_matches_bt4s_own_pwm_baseline() -> None:
    """The anchor is verifiable from this repo, not assumed about someone else's.

    This is *why* this convention was chosen over the annotation-derived one: the
    shipped baseline peaks on exactly these positions, so a panel and the backend BT4
    already returns from ``default()`` agree by construction.
    """
    scores = ConsensusPwmSplicePredictor().score_sequence(_SEQ)
    assert max(range(len(_SEQ)), key=lambda i: scores.donor[i]) == _DONOR
    assert max(range(len(_SEQ)), key=lambda i: scores.acceptor[i]) == _ACCEPTOR


def test_splice_kind_asks_whether_either_motif_is_present() -> None:
    """A combined-track backend cannot distinguish the two, so neither does the check."""
    assert canonical_motif_at(_SEQ, _DONOR, "splice")
    assert canonical_motif_at(_SEQ, _ACCEPTOR, "splice")
    assert not canonical_motif_at(_SEQ, 0, "splice")


def test_an_unknown_site_kind_is_refused() -> None:
    with pytest.raises(ValueError, match="unknown site kind"):
        canonical_motif_at(_SEQ, 0, "branchpoint")


def test_a_position_without_room_for_its_motif_is_not_canonical() -> None:
    """The probe walks off both ends, so an offset sweep never raises."""
    assert not canonical_motif_at(_SEQ, -1, "donor")
    assert not canonical_motif_at(_SEQ, len(_SEQ), "donor")
    assert not canonical_motif_at(_SEQ, 0, "acceptor")


# --------------------------------------------------------------------------
# The off-by-one this format exists to catch


def test_the_exonic_boundary_convention_is_refused_by_name() -> None:
    """The single most likely mistake, turned into a message naming the fix.

    Donors one base early and acceptors one base late is exactly what you get from an
    annotation table's exon boundaries. Without this check the gate would report a
    competent model as hopeless and nothing would say why.
    """
    windows = [
        SpliceWindow(f"w{i}", "chr1", _SEQ, donors=(_DONOR - 1,), acceptors=(_ACCEPTOR + 1,))
        for i in range(10)
    ]
    with pytest.raises(ValueError) as excinfo:
        panel_from_windows(windows, negative_construction=_NEG)
    message = str(excinfo.value)
    assert "exonic-boundary convention" in message
    assert "move each donor +1 and each acceptor -1" in message


def test_the_diagnosis_reports_the_shift_that_would_have_worked() -> None:
    """Even for a shift with no named convention, the number is in the message."""
    windows = [
        SpliceWindow(f"w{i}", "chr1", _SEQ, donors=(_DONOR + 2,), acceptors=(_ACCEPTOR + 2,))
        for i in range(10)
    ]
    with pytest.raises(ValueError, match=r"shifting donors by -2"):
        panel_from_windows(windows, negative_construction=_NEG)


def test_a_correctly_anchored_panel_scores_full_consistency() -> None:
    consistency = _panel().motif_consistency()
    assert consistency.fraction == 1.0
    assert consistency.best_donor_offset == 0
    assert consistency.best_acceptor_offset == 0


def test_the_floor_tolerates_the_real_non_canonical_minority() -> None:
    """~1% of human introns are GC-AG or AT-AC; the floor must not refuse them."""
    windows = [_window(f"w{i}") for i in range(19)]
    odd = _SEQ[:_DONOR] + "GC" + _SEQ[_DONOR + 2 :]
    windows.append(SpliceWindow("w19", "chr1", odd, donors=(_DONOR,), acceptors=(_ACCEPTOR,)))
    panel = panel_from_windows(windows, negative_construction=_NEG)
    assert MIN_MOTIF_CONSISTENCY <= panel.motif_consistency().fraction < 1.0


def test_a_scrambled_panel_says_the_positions_are_not_off_by_a_fixed_offset() -> None:
    """No shift helps, so the message points at orientation/assembly instead."""
    seq = "A" * 60
    windows = [SpliceWindow(f"w{i}", "chr1", seq, donors=(10,), acceptors=(40,)) for i in range(6)]
    with pytest.raises(ValueError, match="no shift in the probed range helps"):
        panel_from_windows(windows, negative_construction=_NEG)


# --------------------------------------------------------------------------
# Window validation


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"window_id": " "}, "window_id is empty"),
        ({"group": ""}, "group is empty"),
        ({"strand": "?"}, "strand="),
        ({"donors": (500,)}, "outside the"),
        ({"donors": (5, 3)}, "sorted and unique"),
        ({"acceptors": (0,)}, "no room for its AG"),
    ],
)
def test_a_malformed_window_is_refused(kwargs: dict[str, object], match: str) -> None:
    base: dict[str, object] = {
        "window_id": "w0",
        "group": "chr1",
        "sequence": _SEQ,
        "donors": (_DONOR,),
        "acceptors": (_ACCEPTOR,),
    }
    with pytest.raises(ValueError, match=match):
        SpliceWindow(**{**base, **kwargs})  # type: ignore[arg-type]


def test_a_donor_at_the_very_end_has_no_room_for_its_gt() -> None:
    with pytest.raises(ValueError, match="no room for its GT"):
        SpliceWindow("w0", "chr1", "ACGTACG", donors=(6,))


def test_a_position_cannot_be_both_kinds() -> None:
    with pytest.raises(ValueError, match="both donor and acceptor"):
        SpliceWindow("w0", "chr1", _SEQ, donors=(_DONOR,), acceptors=(_DONOR,))


def test_a_pure_negative_window_is_legitimate() -> None:
    """A deep-intron control carries no sites and must not be refused."""
    window = SpliceWindow("deep", "chr3", "A" * 100, note="deep intron, no annotated site")
    panel = panel_from_windows([_window(), window], negative_construction=_NEG)
    assert panel.n_sites == 2
    assert window.n_sites == 0


def test_labels_align_to_the_sequence() -> None:
    window = _window()
    labels = window.labels("donor")
    assert len(labels) == len(_SEQ)
    assert labels[_DONOR] == 1
    assert sum(labels) == 1
    with pytest.raises(ValueError, match="unknown site kind"):
        window.labels("splice")


def test_sites_are_returned_in_position_order() -> None:
    assert _window().sites() == ((_DONOR, "donor"), (_ACCEPTOR, "acceptor"))


# --------------------------------------------------------------------------
# Panel-level provenance


def test_negative_construction_is_mandatory() -> None:
    """PR-AUC's floor is the prevalence, and the prevalence is this choice."""
    with pytest.raises(ValueError, match="negative_construction is required"):
        SplicePanel(windows=(_window(),), negative_construction="  ")


def test_a_panel_with_no_sites_is_refused() -> None:
    """With no positive class every stratum is unscoreable."""
    with pytest.raises(ValueError, match="no annotated splice site"):
        panel_from_windows(
            [SpliceWindow("w0", "chr1", "A" * 50)], negative_construction=_NEG
        )


def test_duplicate_window_ids_are_refused() -> None:
    with pytest.raises(ValueError, match="duplicate window_id"):
        panel_from_windows([_window("w"), _window("w")], negative_construction=_NEG)


def test_training_chromosome_overlap_is_reported() -> None:
    """chr2 is a training chromosome for both models; chr1 and chr3 are held out."""
    panel = panel_from_windows(
        [_window("a", "chr1"), _window("b", "chr2"), _window("c", "chr3")],
        negative_construction=_NEG,
    )
    assert panel.training_overlap == ("chr2",)
    held_out = panel_from_windows(
        [_window("a", "chr1"), _window("c", "chr3")], negative_construction=_NEG
    )
    assert held_out.training_overlap == ()


def test_describe_reports_prevalence_beside_the_counts() -> None:
    """A PR-AUC without its floor is not interpretable."""
    summary = _panel().describe()
    assert summary["n_positions"] == len(_SEQ)
    assert summary["donor_prevalence"] == pytest.approx(1 / len(_SEQ))
    assert summary["motif_consistency"] == 1.0
    assert summary["negative_construction"] == _NEG


# --------------------------------------------------------------------------
# The content hash (invariant #7)


def test_content_hash_is_order_independent() -> None:
    a = panel_from_windows([_window("a"), _window("b")], negative_construction=_NEG)
    b = panel_from_windows([_window("b"), _window("a")], negative_construction=_NEG)
    assert a.content_hash() == b.content_hash()


def test_content_hash_moves_on_any_real_change() -> None:
    base = _panel().content_hash()
    moved = panel_from_windows(
        [SpliceWindow("w0", "chr3", _SEQ, donors=(_DONOR,), acceptors=(_ACCEPTOR,))],
        negative_construction=_NEG,
    )
    relabelled = panel_from_windows(
        [_window()], negative_construction="a different denominator"
    )
    annotated = panel_from_windows(
        [_window()], negative_construction=_NEG, annotation="GENCODE v44"
    )
    hashes = {base, moved.content_hash(), relabelled.content_hash(), annotated.content_hash()}
    assert len(hashes) == 4


def test_content_hash_is_stable_across_calls() -> None:
    panel = _panel()
    assert panel.content_hash() == panel.content_hash()


# --------------------------------------------------------------------------
# The reader


def _write(
    tmp_path: Path,
    rows: str,
    header: str = "window_id\tgroup\tsequence\tdonors\tacceptors",
) -> Path:
    path = tmp_path / "panel.tsv"
    path.write_text(f"{header}\n{rows}", encoding="utf-8")
    return path


def test_reader_round_trips_a_well_formed_panel(tmp_path: Path) -> None:
    path = _write(tmp_path, f"w0\tchr1\t{_SEQ}\t{_DONOR}\t{_ACCEPTOR}\n")
    panel = read_splice_panel(path, negative_construction=_NEG, annotation="GENCODE v44")
    assert len(panel) == 1
    assert panel.windows[0].donors == (_DONOR,)
    assert panel.annotation == "GENCODE v44"
    assert panel.source == str(path)


def test_reader_refuses_a_typod_column(tmp_path: Path) -> None:
    """A silently-ignored ``acceptor`` column would relabel every acceptor a negative."""
    path = _write(
        tmp_path,
        f"w0\tchr1\t{_SEQ}\t{_DONOR}\t{_ACCEPTOR}\n",
        header="window_id\tgroup\tsequence\tdonors\tacceptor",
    )
    with pytest.raises(ValueError, match="unrecognised column"):
        read_splice_panel(path, negative_construction=_NEG)


def test_reader_refuses_an_n_inside_a_window(tmp_path: Path) -> None:
    """An N would be an unscoreable position masquerading as a real negative."""
    path = _write(tmp_path, f"w0\tchr1\t{'N' + _SEQ[1:]}\t{_DONOR}\t{_ACCEPTOR}\n")
    with pytest.raises(ValueError, match="w0"):
        read_splice_panel(path, negative_construction=_NEG)


def test_reader_refuses_a_repeated_position(tmp_path: Path) -> None:
    """Usually a transcript loop double-counting shared exons -- it changes prevalence."""
    path = _write(tmp_path, f"w0\tchr1\t{_SEQ}\t{_DONOR},{_DONOR}\t{_ACCEPTOR}\n")
    with pytest.raises(ValueError, match="repeats position"):
        read_splice_panel(path, negative_construction=_NEG)


def test_reader_refuses_a_non_integer_position(tmp_path: Path) -> None:
    path = _write(tmp_path, f"w0\tchr1\t{_SEQ}\tnine\t{_ACCEPTOR}\n")
    with pytest.raises(ValueError, match="is not an integer"):
        read_splice_panel(path, negative_construction=_NEG)


def test_reader_accepts_an_empty_position_list(tmp_path: Path) -> None:
    path = _write(
        tmp_path, f"w0\tchr1\t{_SEQ}\t{_DONOR}\t{_ACCEPTOR}\nw1\tchr3\t{'A' * 40}\t\t\n"
    )
    panel = read_splice_panel(path, negative_construction=_NEG)
    assert panel.windows[1].n_sites == 0


def test_reader_refuses_a_header_only_file(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no rows"):
        read_splice_panel(_write(tmp_path, ""), negative_construction=_NEG)


def test_reader_refuses_a_missing_required_column(tmp_path: Path) -> None:
    path = _write(tmp_path, f"w0\t{_SEQ}\n", header="window_id\tsequence")
    with pytest.raises(ValueError, match="missing required column"):
        read_splice_panel(path, negative_construction=_NEG)


def test_lowering_the_consistency_floor_is_possible_but_deliberate(tmp_path: Path) -> None:
    """A U12 AT-AC panel is legitimate; it just has to be declared, never quietly."""
    seq = "CCCCCCCAG" "ATATCC" "CCCCCCCCCCCCCTTTTTTTTTCAC" "GCCCCCCCCC"
    path = _write(tmp_path, f"w0\tchr1\t{seq}\t{_DONOR}\t{_ACCEPTOR}\n")
    with pytest.raises(ValueError, match="canonical dinucleotide"):
        read_splice_panel(path, negative_construction=_NEG)
    panel = read_splice_panel(
        path, negative_construction=_NEG, annotation="U12 AT-AC set", min_motif_consistency=0.0
    )
    assert panel.motif_consistency().fraction == 0.0


# --------------------------------------------------------------------------
# Chromosome naming: the panel builder's most likely trip-up


@pytest.mark.parametrize(
    ("group", "normalized"),
    [
        ("chr2", "2"), ("CHR2", "2"), ("2", "2"), ("chrX", "X"), ("x", "X"),
        ("chr1", "1"), (" chr9 ", "9"),
        ("NC_000002.12", None), ("scaffold_7", None), ("", None), ("chrZ", None),
    ],
)
def test_chromosome_names_normalize_across_conventions(
    group: str, normalized: str | None
) -> None:
    """GENCODE writes ``chr2``, Ensembl writes ``2``, and both are the same chromosome."""
    from bt4.biomodels.splice.panel import normalize_chromosome

    assert normalize_chromosome(group) == normalized


def test_training_overlap_matches_either_spelling() -> None:
    """A training panel must be caught whichever convention named it.

    Matching only the ``chr``-prefixed spelling meant an Ensembl-named panel drawn
    entirely from the models' training chromosomes reported no overlap at all.
    """
    for group in ("chr2", "2", "CHR2", "chrX", "X"):
        panel = panel_from_windows([_window("w", group)], negative_construction=_NEG)
        assert panel.training_overlap == (group,), group
        assert panel.unclassified_groups == ()


def test_held_out_chromosomes_match_in_either_spelling() -> None:
    for group in ("chr1", "1", "chr9", "9"):
        panel = panel_from_windows([_window("w", group)], negative_construction=_NEG)
        assert panel.training_overlap == ()
        assert panel.unclassified_groups == ()


def test_an_unclassifiable_group_is_reported_as_unknown_not_clean() -> None:
    """Absence from ``training_overlap`` must not be read as evidence of anything."""
    panel = panel_from_windows(
        [_window("w", "NC_000002.12")], negative_construction=_NEG
    )
    assert panel.training_overlap == ()
    assert panel.unclassified_groups == ("NC_000002.12",)
    assert panel.describe()["unclassified_groups"] == ["NC_000002.12"]


# --------------------------------------------------------------------------
# Sites the backend structurally cannot score


def test_a_site_at_a_window_edge_is_reported_as_a_forced_miss() -> None:
    """The labels are right and the case is still unwinnable.

    A donor at position 0 carries a real ``GT``, so the motif check passes and the panel
    is accepted -- but no backend has flanking sequence there, so the PWM returns exactly
    ``0.0``. That is a ``label=1`` case the model cannot get right, depressing every
    metric through no fault of its own. Reported separately precisely because nothing is
    wrong with the panel's labels; the cure is a wider window, not a dropped site.
    """
    from bt4.biomodels.splice.panel import DEFAULT_EDGE_MARGIN

    # A donor at position 0, and an identical one comfortably in the interior.
    filler = "".join("ACTC"[i % 4] for i in range(60))
    edge_seq = "GT" + filler + "CAG" + "GTAAGT" + filler
    interior = edge_seq.index("CAG" + "GTAAGT") + 3
    window = SpliceWindow("edgy", "chr1", edge_seq, donors=(0, interior))
    panel = panel_from_windows([window], negative_construction=_NEG)

    assert panel.motif_consistency().fraction == 1.0  # the labels are correct
    assert panel.edge_sites() == (("edgy", 0, "donor"),)
    assert panel.describe()["n_edge_sites"] == 1

    scores = ConsensusPwmSplicePredictor().score_sequence(edge_seq)
    assert scores.donor[0] == 0.0  # structurally unscoreable
    assert scores.donor[interior] > 0.9  # the same motif, scored properly
    assert DEFAULT_EDGE_MARGIN > 0


def test_the_edge_margin_is_caller_settable_for_a_cnn_panel() -> None:
    """A 10 kb-context CNN wants orders of magnitude more margin than the PWM."""
    panel = _panel()
    assert panel.edge_sites(margin=0) == ()
    # With a CNN-sized margin every site in a short window is near an edge.
    assert len(panel.edge_sites(margin=len(_SEQ))) == panel.n_sites
