"""Tests for the measured CDS-variant panel format (``expression/panel.py``).

The panel is the bridge between a published supplementary table and the acceptance
gate, so its job is to **refuse** rather than to cope. The load-bearing property: a row
RiboNN would silently drop must be a hard error here, because a quietly-shortened panel
lets the gate report an honest ``n_test`` for a dataset nobody chose.
"""

from __future__ import annotations

import pathlib

import pytest

from bt4.biomodels.expression import (
    MAX_CDS_UTR3_LEN,
    MAX_UTR5_LEN,
    PANEL_COLUMNS,
    ExpressionPanel,
    PanelRow,
    panel_from_rows,
    read_panel,
)

_HEADER = "group\tvariant_id\tcds\tmeasured\tutr5\tutr3\n"


def _row(**overrides: object) -> PanelRow:
    base = dict(
        group="P1",
        variant_id="v1",
        cds="ATGAAATAA",
        measured=1.0,
        utr5="GCCACC",
        utr3="GCTAAT",
    )
    base.update(overrides)
    return PanelRow(**base)  # type: ignore[arg-type]


def _write(tmp_path: pathlib.Path, body: str, header: str = _HEADER) -> pathlib.Path:
    path = tmp_path / "panel.tsv"
    path.write_text(header + body, encoding="utf-8")
    return path


# --- reading ------------------------------------------------------------------


def test_reads_a_minimal_panel(tmp_path: pathlib.Path) -> None:
    path = _write(
        tmp_path,
        "P1\tv1\tATGAAATAA\t1.5\tGCCACC\tGCTAAT\n"
        "P1\tv2\tATGAAGTAA\t2.5\tGCCACC\tGCTAAT\n"
        "P2\tv3\tATGGATTGA\t-0.5\tGCCACC\tGCTAAT\n",
    )
    panel = read_panel(path)
    assert len(panel) == 3
    assert panel.groups == ("P1", "P2")
    assert panel.group_sizes() == {"P1": 2, "P2": 1}
    assert panel.samples()[0] == ("ATGAAATAA", 1.5, "P1")
    assert panel.source == str(path)


def test_optional_columns_are_carried_through(tmp_path: pathlib.Path) -> None:
    # A measurement without its readout and cell type is a number with no question
    # attached, so these travel with it into every report.
    header = _HEADER.rstrip("\n") + "\treadout\tcell_type\tspecies\n"
    path = _write(
        tmp_path,
        "P1\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTAAT\tmean_ribosome_load\tHEK293T\thuman\n"
        "P1\tv2\tATGAAGTAA\t2.0\tGCCACC\tGCTAAT\tmean_ribosome_load\tHEK293T\thuman\n",
        header=header,
    )
    panel = read_panel(path)
    assert panel.describe()["readouts"] == ["mean_ribosome_load"]
    assert panel.describe()["cell_types"] == ["HEK293T"]
    assert panel.rows[0].species == "human"


def test_case_and_whitespace_are_normalised(tmp_path: pathlib.Path) -> None:
    path = _write(
        tmp_path,
        "P1\tv1\tatgaaataa\t1.0\t gccacc \tgctaat\n"
        "P1\tv2\tATGAAGTAA\t2.0\tGCCACC\tGCTAAT\n",
    )
    panel = read_panel(path)
    assert panel.rows[0].cds == "ATGAAATAA"
    assert panel.rows[0].utr5 == "GCCACC"


# --- refusals: the whole point ------------------------------------------------


def test_an_over_length_utr5_is_refused_not_dropped(tmp_path: pathlib.Path) -> None:
    # RiboNN filters such rows inside its data module. If we passed them through, the
    # gate would score a smaller panel than the one that was pre-registered.
    long_utr5 = "A" * (MAX_UTR5_LEN + 1)
    path = _write(tmp_path, f"P1\tv1\tATGAAATAA\t1.0\t{long_utr5}\tGCTAAT\n")
    with pytest.raises(ValueError, match="over RiboNN's 1381 nt cap"):
        read_panel(path)
    with pytest.raises(ValueError, match="silently DROP"):
        read_panel(path)


def test_an_over_length_cds_plus_utr3_is_refused(tmp_path: pathlib.Path) -> None:
    cds = "ATG" + "AAA" * ((MAX_CDS_UTR3_LEN // 3) + 1) + "TAA"
    path = _write(tmp_path, f"P1\tv1\t{cds}\t1.0\tGCCACC\tGCTAAT\n")
    with pytest.raises(ValueError, match="over RiboNN's 11937 nt cap"):
        read_panel(path)


@pytest.mark.parametrize(
    ("body", "match"),
    [
        ("P1\tv1\tATGAAATA\t1.0\tGCCACC\tGCTAAT\n", "not a multiple of 3"),
        ("P1\tv1\tATGAAAAAA\t1.0\tGCCACC\tGCTAAT\n", "not a stop codon"),
        ("P1\tv1\tATGXXXTAA\t1.0\tGCCACC\tGCTAAT\n", "non-ACGT"),
        ("P1\tv1\tATGAAATAA\tabc\tGCCACC\tGCTAAT\n", "is not a number"),
        ("P1\tv1\tATGAAATAA\tnan\tGCCACC\tGCTAAT\n", "is not finite"),
        ("P1\tv1\tATGAAATAA\tinf\tGCCACC\tGCTAAT\n", "is not finite"),
        ("P1\tv1\tATGAAATAA\t1.0\t\tGCTAAT\n", "missing value"),
        ("P1\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTNNT\n", "non-ACGT"),
        ("\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTAAT\n", "missing value"),
    ],
)
def test_bad_rows_raise_with_the_row_named(
    tmp_path: pathlib.Path, body: str, match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        read_panel(_write(tmp_path, body))


def test_a_typod_column_is_refused_rather_than_ignored(tmp_path: pathlib.Path) -> None:
    # Silently ignoring an unknown column is how a mislabelled `measurement` column
    # ends up unused while the gate runs happily on nothing.
    header = _HEADER.rstrip("\n") + "\tmeasurment\n"
    path = _write(tmp_path, "P1\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTAAT\tx\n", header=header)
    with pytest.raises(ValueError, match="unrecognised column"):
        read_panel(path)


def test_a_missing_required_column_names_what_is_expected(
    tmp_path: pathlib.Path,
) -> None:
    path = _write(tmp_path, "P1\tv1\tATGAAATAA\t1.0\n", header="group\tvariant_id\tcds\tmeasured\n")
    with pytest.raises(ValueError, match="missing required column"):
        read_panel(path)


def test_an_empty_file_is_refused(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="no header row"):
        read_panel(path)


def test_duplicate_variant_ids_are_refused() -> None:
    with pytest.raises(ValueError, match="duplicate variant_id"):
        panel_from_rows([_row(variant_id="v1"), _row(variant_id="v1", measured=2.0)])


def test_an_empty_panel_is_refused() -> None:
    with pytest.raises(ValueError, match="at least one row"):
        ExpressionPanel(rows=())


def test_an_unknown_species_is_refused(tmp_path: pathlib.Path) -> None:
    header = _HEADER.rstrip("\n") + "\tspecies\n"
    path = _write(
        tmp_path, "P1\tv1\tATGAAATAA\t1.0\tGCCACC\tGCTAAT\tzebrafish\n", header=header
    )
    with pytest.raises(ValueError, match="human and mouse weights only"):
        read_panel(path)


# --- content hash -------------------------------------------------------------


def test_content_hash_ignores_row_order_but_not_values() -> None:
    a = _row(variant_id="v1", measured=1.0)
    b = _row(variant_id="v2", cds="ATGAAGTAA", measured=2.0)

    forward = panel_from_rows([a, b]).content_hash()
    reversed_ = panel_from_rows([b, a]).content_hash()
    assert forward == reversed_  # re-ordering a file is not a different panel

    changed = panel_from_rows([a, _row(variant_id="v2", cds="ATGAAGTAA", measured=2.5)])
    assert changed.content_hash() != forward  # a changed measurement is


def test_content_hash_is_stable_across_calls_and_source() -> None:
    rows = [_row(variant_id="v1"), _row(variant_id="v2", measured=2.0)]
    first = panel_from_rows(rows, source="a.tsv")
    second = panel_from_rows(rows, source="b.tsv")
    # Timestamp-free and path-free, so it can be pre-registered before a gate run and
    # compared afterwards (invariant #7).
    assert first.content_hash() == first.content_hash() == second.content_hash()


# --- the sizing facts a gate depends on ---------------------------------------


def test_describe_surfaces_the_gate_sizing_facts() -> None:
    rows = [
        _row(group="P1", variant_id="a1"),
        _row(group="P1", variant_id="a2", measured=2.0),
        _row(group="P2", variant_id="b1", measured=3.0),
    ]
    described = panel_from_rows(rows).describe()
    assert described["n_rows"] == 3
    assert described["n_groups"] == 2
    # Within-group scoring needs 2+ members, so a panel of singletons is visibly unfit
    # before a single model runs.
    assert described["n_groups_with_2_or_more"] == 1
    assert described["n_utr_contexts"] == 1


def test_contexts_bucket_rows_by_their_utr_pair() -> None:
    # A predictor carries its UTR context on the model, so a multi-transcript panel
    # cannot be scored in one invocation; this is the split that makes that explicit.
    rows = [
        _row(variant_id="a", utr5="GCCACC", utr3="GCTAAT"),
        _row(variant_id="b", utr5="GCCACC", utr3="GCTAAT"),
        _row(variant_id="c", utr5="AAACCC", utr3="TTTGGG"),
    ]
    contexts = panel_from_rows(rows).contexts()
    assert len(contexts) == 2
    assert len(contexts[("GCCACC", "GCTAAT")]) == 2
    assert len(contexts[("AAACCC", "TTTGGG")]) == 1


def test_panel_feeds_the_gate_directly() -> None:
    from bt4.biomodels.expression import verify_expression_gate
    from bt4.biomodels.expression.gate import ExpressionEvalCase

    rows = [
        _row(group=f"P{i // 3}", variant_id=f"v{i}", measured=float(i))
        for i in range(12)
    ]
    panel = panel_from_rows(rows)
    cases = [
        ExpressionEvalCase(predicted=measured, measured=measured, group=group)
        for _cds, measured, group in panel.samples()
    ]
    report = verify_expression_gate(cases, bootstrap_resamples=0)
    assert report.n_calibration + report.n_test == 12
    assert report.n_groups == 4


def test_panel_columns_are_the_documented_set() -> None:
    assert PANEL_COLUMNS == (
        "group", "variant_id", "cds", "measured", "utr5", "utr3",
        "readout", "cell_type", "species",
    )


def test_an_over_cap_cds_is_refused_by_the_panel_not_by_the_stdlib(
    tmp_path: pathlib.Path,
) -> None:
    """A row longer than csv's 131,072-char field cap must still reach the real check.

    A *valid* row can never approach that cap -- `MAX_CDS_UTR3_LEN` holds CDS+3'UTR to
    11,937 nt, RiboNN's input width -- so unlike the splice panel, this is not about
    reading a legitimately huge field. It is about which error an over-long row gets:
    the panel's own message, naming RiboNN's limit and what to do about it, rather than
    a bare `_csv.Error` citing a field limit the format does not have.
    """
    import csv

    cds = "ATG" + "AAA" * 50_000 + "TAA"  # 150,006 nt, a well-formed but far-too-long ORF
    assert len(cds) > csv.field_size_limit(), "fixture must exceed the csv cap"

    path = _write(tmp_path, f"P1\tv1\t{cds}\t1.5\tGCCACC\tGCTAAT\n")
    with pytest.raises(ValueError, match="over RiboNN's"):
        read_panel(path)
    assert csv.field_size_limit() == 131_072  # the caller's limit is left alone
