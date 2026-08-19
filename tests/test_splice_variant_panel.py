"""Tests for the measured splice-variant panel.

``SpliceVariantCase`` shipped with no way to build one from data, which made the
variant half of the gate unreachable. This is that bridge, and it is shaped by the one
public benchmark it has to read: ``kitzmanlab/splicebench2023``.

These tests pin the properties that keep a variant gate honest:

* **scores are named columns**, so a benchmark's own pre-computed predictions come
  through as data -- which makes the most useful zero-model check possible: run the gate
  on a published benchmark's numbers and see whether BT4 reproduces its figures;
* a score column with **gaps is refused**, not silently scored on its covered subset;
* the **region is required and enumerated**, because the exonic/intronic gap is the
  finding BT4 has to reproduce and a typo'd stratum would split it silently;
* **held-out status is checkable and usually fails** -- over half of splicebench2023 is
  on chromosomes both models trained on;
* one gene cannot support a general claim, so two groups are required.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bt4.biomodels.splice.gate import VARIANT_EFFECT, verify_splice_gate
from bt4.biomodels.splice.variant_panel import (
    REGIONS,
    SpliceVariantPanel,
    VariantRow,
    read_variant_panel,
    variant_panel_from_rows,
)

_NEG = "all assayed variants the assay called non-disruptive"
_ASSAY = "six-assay sdv_fc2 composite"


def _rows() -> list[VariantRow]:
    """A small panel spanning both regions, two genes, and two scorers."""
    out = []
    for gene, chrom in (("BRCA1", "17"), ("POU1F1", "3")):
        for i in range(10):
            disruptive = i % 5 == 0
            region = "exonic" if i % 2 == 0 else "intronic"
            out.append(
                VariantRow(
                    variant_id=f"{gene}_{i}",
                    group=gene,
                    region=region,
                    label=int(disruptive),
                    scores=(
                        ("spliceai_masked", 0.9 if disruptive else 0.1),
                        ("pangolin_masked", 0.8 if disruptive else 0.2),
                    ),
                    chromosome=chrom,
                )
            )
    return out


def _panel(**kwargs: object) -> SpliceVariantPanel:
    return variant_panel_from_rows(
        _rows(), negative_construction=_NEG, assay=_ASSAY, **kwargs
    )  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# Scores are data, not schema


def test_score_columns_are_discovered_from_the_file() -> None:
    """A benchmark's own predictions come through without a schema change."""
    panel = _panel()
    assert panel.score_columns == ("pangolin_masked", "spliceai_masked")


def test_cases_are_built_from_one_named_column() -> None:
    panel = _panel()
    cases = panel.cases("spliceai_masked")
    assert len(cases) == len(panel)
    assert {c.region for c in cases} == set(REGIONS)
    assert {c.group for c in cases} == {"BRCA1", "POU1F1"}


def test_the_gate_recognises_them_as_the_variant_task() -> None:
    """The bridge actually reaches the gate, which is the point of the module."""
    report = verify_splice_gate(
        _panel().cases("pangolin_masked"), negative_construction=_NEG
    )
    assert report.task == VARIANT_EFFECT
    assert {s.name for s in report.strata} == set(REGIONS)


def test_an_unknown_score_column_names_what_is_available() -> None:
    with pytest.raises(KeyError, match="spliceai_masked"):
        _panel().cases("spliceai_unmasked")


def test_a_partially_covered_score_column_is_refused() -> None:
    """Scoring the covered subset answers a question about a smaller panel.

    A tool that did not cover every variant is a real situation; quietly dropping the
    rows it missed while reporting the panel's name is the dishonest response to it.
    """
    rows = _rows()
    rows[3] = VariantRow(
        variant_id=rows[3].variant_id,
        group=rows[3].group,
        region=rows[3].region,
        label=rows[3].label,
        scores=(("pangolin_masked", 0.5),),  # no spliceai score
        chromosome=rows[3].chromosome,
    )
    panel = variant_panel_from_rows(rows, negative_construction=_NEG)
    assert len(panel.cases("pangolin_masked")) == len(rows)
    with pytest.raises(ValueError, match="missing for 1 of 20 rows"):
        panel.cases("spliceai_masked")


# --------------------------------------------------------------------------
# Row validation


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"variant_id": " "}, "variant_id is empty"),
        ({"group": ""}, "group is empty"),
        ({"region": "deep_intronic"}, "region="),
        ({"label": 2}, "is not 0 or 1"),
    ],
)
def test_a_malformed_row_is_refused(kwargs: dict[str, object], match: str) -> None:
    base: dict[str, object] = {
        "variant_id": "v", "group": "G", "region": "exonic", "label": 1,
    }
    with pytest.raises(ValueError, match=match):
        VariantRow(**{**base, **kwargs})  # type: ignore[arg-type]


def test_region_is_enumerated_not_free_text() -> None:
    """The exonic/intronic gap is the finding; a typo'd stratum would split it."""
    assert REGIONS == ("exonic", "intronic")
    with pytest.raises(ValueError, match="region="):
        VariantRow("v", "G", "Exonic", 1)


# --------------------------------------------------------------------------
# Panel-level provenance


def test_negative_construction_is_mandatory() -> None:
    with pytest.raises(ValueError, match="negative_construction is required"):
        SpliceVariantPanel(rows=tuple(_rows()), negative_construction="  ")


def test_one_gene_cannot_support_a_general_claim() -> None:
    single = [r for r in _rows() if r.group == "BRCA1"]
    with pytest.raises(ValueError, match="two leakage-control groups"):
        variant_panel_from_rows(single, negative_construction=_NEG)


def test_duplicate_variant_ids_are_refused() -> None:
    rows = _rows()
    with pytest.raises(ValueError, match="duplicate variant_id"):
        variant_panel_from_rows(rows + rows[:1], negative_construction=_NEG)


def test_training_chromosome_overlap_is_reported_per_gene() -> None:
    """The measured reality of the recommended benchmark: over half is not held out.

    BRCA1 is on chr17 -- a chromosome both models trained on -- and it is the largest
    gene set in splicebench2023 and the closest public thing to BT4's regime.
    """
    panel = _panel()
    assert panel.training_overlap == ("BRCA1",)
    assert panel.held_out is False


def test_a_chr3_only_panel_is_held_out() -> None:
    """The usable subset: POU1F1, MST1R and MLH1 are all on chr3."""
    rows = [r for r in _rows() if r.group == "POU1F1"]
    rows += [
        VariantRow(f"MST1R_{i}", "MST1R", "exonic" if i % 2 else "intronic",
                   int(i % 5 == 0), (("spliceai_masked", 0.5),), "3")
        for i in range(10)
    ]
    panel = variant_panel_from_rows(rows, negative_construction=_NEG)
    assert panel.training_overlap == ()
    assert panel.groups_without_chromosome == ()
    assert panel.held_out is True


def test_a_gene_without_a_chromosome_is_unknown_not_clean() -> None:
    """Absence from the training list must not be read as evidence of anything."""
    rows = [
        VariantRow(r.variant_id, r.group, r.region, r.label, r.scores, "")
        for r in _rows()
    ]
    panel = variant_panel_from_rows(rows, negative_construction=_NEG)
    assert panel.training_overlap == ()
    assert set(panel.groups_without_chromosome) == {"BRCA1", "POU1F1"}
    assert panel.held_out is False


def test_the_assay_definition_travels_with_the_panel() -> None:
    """A composite over several assays' criteria is not a measurement."""
    assert _panel().describe()["assay"] == _ASSAY


def test_describe_reports_prevalence_per_region() -> None:
    summary = _panel().describe()
    assert summary["region_sizes"] == {"exonic": 10, "intronic": 10}
    assert summary["n_positive"] == 4
    assert summary["held_out"] is False


# --------------------------------------------------------------------------
# Content hash (invariant #7)


def test_content_hash_is_order_independent_and_moves_on_change() -> None:
    rows = _rows()
    a = variant_panel_from_rows(rows, negative_construction=_NEG, assay=_ASSAY)
    b = variant_panel_from_rows(list(reversed(rows)), negative_construction=_NEG, assay=_ASSAY)
    assert a.content_hash() == b.content_hash()

    relabelled = variant_panel_from_rows(rows, negative_construction="something else", assay=_ASSAY)
    reassayed = variant_panel_from_rows(rows, negative_construction=_NEG, assay="one assay only")
    assert len({a.content_hash(), relabelled.content_hash(), reassayed.content_hash()}) == 3


# --------------------------------------------------------------------------
# The reader


def _write(tmp_path: Path, body: str, header: str) -> Path:
    path = tmp_path / "variants.tsv"
    path.write_text(f"{header}\n{body}", encoding="utf-8")
    return path


_HEADER = "variant_id\tgroup\tregion\tlabel\tchromosome\tspliceai_masked\tpangolin_masked"


def test_reader_round_trips_and_accepts_boolean_labels(tmp_path: Path) -> None:
    """Benchmarks ship ``True``/``False``; refusing them would be pedantry."""
    body = (
        "v1\tBRCA1\texonic\tTrue\t17\t0.9\t0.8\n"
        "v2\tPOU1F1\tintronic\tFalse\t3\t0.1\t0.2\n"
    )
    panel = read_variant_panel(
        _write(tmp_path, body, _HEADER), negative_construction=_NEG, assay=_ASSAY
    )
    assert [r.label for r in panel] == [1, 0]
    assert panel.score_columns == ("pangolin_masked", "spliceai_masked")
    assert panel.training_overlap == ("BRCA1",)


def test_reader_refuses_a_panel_with_no_score_column(tmp_path: Path) -> None:
    """Every column beyond the fixed set is a prediction; none means nothing to gate."""
    header = "variant_id\tgroup\tregion\tlabel"
    with pytest.raises(ValueError, match="no score column"):
        read_variant_panel(
            _write(tmp_path, "v1\tA\texonic\t1\n", header), negative_construction=_NEG
        )


def test_reader_refuses_a_non_numeric_score(tmp_path: Path) -> None:
    body = "v1\tBRCA1\texonic\tTrue\t17\tNA\t0.8\n"
    with pytest.raises(ValueError, match="is not a number"):
        read_variant_panel(_write(tmp_path, body, _HEADER), negative_construction=_NEG)


def test_reader_refuses_a_non_boolean_label(tmp_path: Path) -> None:
    body = "v1\tBRCA1\texonic\tmaybe\t17\t0.9\t0.8\n"
    with pytest.raises(ValueError, match="is not a boolean"):
        read_variant_panel(_write(tmp_path, body, _HEADER), negative_construction=_NEG)


def test_reader_treats_an_empty_score_cell_as_uncovered(tmp_path: Path) -> None:
    """Not an error at read time -- ``cases`` is where the gap is refused."""
    body = (
        "v1\tBRCA1\texonic\tTrue\t17\t\t0.8\n"
        "v2\tPOU1F1\tintronic\tFalse\t3\t0.1\t0.2\n"
    )
    panel = read_variant_panel(_write(tmp_path, body, _HEADER), negative_construction=_NEG)
    assert len(panel.cases("pangolin_masked")) == 2
    with pytest.raises(ValueError, match="missing for 1 of 2 rows"):
        panel.cases("spliceai_masked")


# --------------------------------------------------------------------------
# The CLI surface


def _write_cli_panel(tmp_path: Path) -> str:
    header = "variant_id\tgroup\tregion\tlabel\tchromosome\tspliceai_masked"
    lines = [header]
    for gene, chrom in (("BRCA1", "17"), ("POU1F1", "3")):
        for i in range(20):
            sdv = i % 4 == 0
            region = "exonic" if (i // 4) % 2 == 0 else "intronic"
            lines.append(
                f"{gene}_{i}\t{gene}\t{region}\t{'True' if sdv else 'False'}\t{chrom}\t"
                f"{0.9 if sdv else 0.1}"
            )
    path = tmp_path / "variants.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def test_cli_lists_score_columns_when_none_is_chosen(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Masked and unmasked answer different questions, so the pick is deliberate."""
    from bt4.cli.__main__ import main

    assert main(["variant-gate", _write_cli_panel(tmp_path), "--negative-construction", _NEG]) == 0
    out = capsys.readouterr().out
    assert "spliceai_masked" in out
    assert "Pick one with --score" in out


def test_cli_warns_that_a_training_chromosome_panel_is_optimistic(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """BRCA1 is on chr17. Over half of splicebench2023 is in this position."""
    from bt4.cli.__main__ import main

    assert (
        main(
            [
                "variant-gate",
                _write_cli_panel(tmp_path),
                "--negative-construction",
                _NEG,
                "--score",
                "spliceai_masked",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "NOT HELD OUT" in out
    assert "BRCA1" in out
    assert "cannot support promotion" in out


def test_cli_prints_the_published_anchor_beside_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The exonic/intronic gap is what this panel exists to let BT4 reproduce."""
    from bt4.cli.__main__ import main

    main(
        [
            "variant-gate",
            _write_cli_panel(tmp_path),
            "--negative-construction",
            _NEG,
            "--score",
            "spliceai_masked",
        ]
    )
    out = capsys.readouterr().out
    assert "0.419" in out and "0.773" in out
    assert "suspect the panel" in out


def test_cli_requires_the_negative_construction(tmp_path: Path) -> None:
    from bt4.cli.__main__ import main

    with pytest.raises(SystemExit):
        main(["variant-gate", _write_cli_panel(tmp_path)])
