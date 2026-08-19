"""Tests for the splicebench2023 -> BT4 variant-panel converter.

Smith & Kitzman's benchmark is the one panel that can check BT4's gate machinery before
any model is installed: 3,616 MPSA-measured variants, already scored by eight tools,
MIT-licensed. The converter's job is to get the column mapping right and to carry the
one fact about the benchmark that most affects how its numbers may be read -- that over
half of it is on chromosomes both models trained on.

The tests run against synthetic TSVs carrying the archive's real column names, so the
mapping is checkable without the 334 MB download.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name: str, filename: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def converter() -> ModuleType:
    return _load("bt4_make_splicebench_variant_panel", "make_splicebench_variant_panel.py")


def _write_archive(tmp_path: Path, converter: ModuleType, *, rows: int = 4) -> Path:
    """Write a synthetic `scored_data` directory using the archive's real column names."""
    directory = tmp_path / "scored_data"
    directory.mkdir()
    header = "\t".join(
        ["varlist", "gene_name", "sdv_fc2", "exon", "DS_maxm", "DS_max",
         "pang_max_abs", "pang_max_nomask_abs"]
    )
    for filename, gene in {**converter.SOURCE_FILES, converter.MLH1_FILE: "MLH1"}.items():
        lines = [header]
        for i in range(rows):
            disruptive = "True" if i % 2 == 0 else "False"
            exonic = "True" if i < rows // 2 else "False"
            lines.append(
                "\t".join(
                    [f"{gene.lower()}:{1000 + i}:A:T", gene, disruptive, exonic,
                     "0.91", "0.93", "0.44", "0.46"]
                )
            )
        (directory / filename).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return directory


# --------------------------------------------------------------------------
# The column mapping


def test_the_label_and_stratifier_columns_are_the_published_ones(
    converter: ModuleType,
) -> None:
    """``sdv_fc2`` and ``exon`` are what the paper's 0.419/0.773 split uses."""
    source = (_SCRIPTS / "make_splicebench_variant_panel.py").read_text(encoding="utf-8")
    assert '_LABEL = "sdv_fc2"' in source
    assert '_REGION = "exon"' in source


def test_masked_and_unmasked_scores_are_kept_apart(converter: ModuleType) -> None:
    """They answer different questions, so the choice belongs at gate time.

    Masked suppresses scores at annotated sites; collapsing the two here would make the
    gate's number depend on an invisible decision taken during conversion.
    """
    assert converter.SCORE_COLUMNS == {
        "DS_maxm": "spliceai_masked",
        "DS_max": "spliceai_unmasked",
        "pang_max_abs": "pangolin_masked",
        "pang_max_nomask_abs": "pangolin_unmasked",
    }


def test_conversion_produces_a_readable_bt4_panel(
    converter: ModuleType, tmp_path: Path
) -> None:
    """End-to-end: the converter's output is what BT4's reader accepts."""
    from bt4.api import read_variant_panel

    directory = _write_archive(tmp_path, converter)
    rows, counts = converter.convert(directory)
    out = tmp_path / "variants.tsv"
    converter.write_panel(rows, out)

    panel = read_variant_panel(
        out,
        negative_construction="assayed variants the assay called non-disruptive",
        assay="MPSA sdv_fc2 composite",
    )
    assert len(panel) == 4 * len(converter.SOURCE_FILES)
    assert set(panel.groups) == set(converter.SOURCE_FILES.values())
    assert panel.score_columns == (
        "pangolin_masked", "pangolin_unmasked", "spliceai_masked", "spliceai_unmasked"
    )
    assert {c.region for c in panel.cases("spliceai_masked")} == {"exonic", "intronic"}
    assert counts["skipped_unparseable"] == 0


# --------------------------------------------------------------------------
# The fact that most affects how the numbers may be read


def test_every_genes_chromosome_is_recorded(converter: ModuleType) -> None:
    """Held-out status is a property of the chromosome, so each gene declares one."""
    assert converter.GENE_CHROMOSOME == {
        "BRCA1": "17", "FAS": "10", "POU1F1": "3",
        "MST1R": "3", "WT1": "11", "MLH1": "3",
    }


def test_the_default_panel_is_not_held_out(converter: ModuleType, tmp_path: Path) -> None:
    """Measured: BRCA1, FAS and WT1 are on chromosomes both models trained on.

    That is 2,077 of the real benchmark's 3,616 variants -- including BRCA1, otherwise
    the closest public thing to BT4's synonymous-CDS regime.
    """
    from bt4.api import read_variant_panel

    directory = _write_archive(tmp_path, converter)
    rows, _ = converter.convert(directory)
    out = tmp_path / "all.tsv"
    converter.write_panel(rows, out)
    panel = read_variant_panel(out, negative_construction="x")

    assert set(panel.training_overlap) == {"BRCA1", "FAS", "WT1"}
    assert panel.held_out is False


def test_held_out_only_keeps_the_chr3_genes(converter: ModuleType, tmp_path: Path) -> None:
    """The usable subset for a claim that depends on being held out."""
    from bt4.api import read_variant_panel

    directory = _write_archive(tmp_path, converter)
    rows, _ = converter.convert(directory, held_out_only=True)
    out = tmp_path / "heldout.tsv"
    converter.write_panel(rows, out)
    panel = read_variant_panel(out, negative_construction="x")

    assert set(panel.groups) == {"POU1F1", "MST1R"}
    assert panel.training_overlap == ()
    assert panel.held_out is True


def test_mlh1_is_excluded_unless_asked(converter: ModuleType, tmp_path: Path) -> None:
    """3,616 and 3,912 are both true about different things; conflating them is the bug."""
    directory = _write_archive(tmp_path, converter)
    default, _ = converter.convert(directory)
    with_mlh1, _ = converter.convert(directory, include_mlh1=True)
    assert {r["group"] for r in default} == set(converter.SOURCE_FILES.values())
    assert "MLH1" in {r["group"] for r in with_mlh1}
    assert len(with_mlh1) > len(default)


# --------------------------------------------------------------------------
# Refusals


def test_a_missing_file_names_the_extraction_step(
    converter: ModuleType, tmp_path: Path
) -> None:
    """A partial panel silently answers a question about a different dataset.

    The most likely cause is the archive's ``for_zenodo`` top directory, which the
    upstream notebooks themselves require renaming to ``data``.
    """
    directory = _write_archive(tmp_path, converter)
    (directory / "wt1_ex9_scored.txt").unlink()
    with pytest.raises(FileNotFoundError, match="mv for_zenodo data"):
        converter.convert(directory)


def test_an_unreadable_label_is_counted_not_guessed(
    converter: ModuleType, tmp_path: Path
) -> None:
    """A row whose label cannot be read is dropped and tallied, never defaulted."""
    directory = _write_archive(tmp_path, converter)
    path = directory / "fas_ex6_snvs_scored.txt"
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = lines[1].split("\t")
    fields[2] = "intermediate"
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    _, counts = converter.convert(directory)
    assert counts["skipped_unparseable"] == 1
    assert counts["FAS"] == 3


# --------------------------------------------------------------------------
# Second review round


def test_a_ragged_row_is_counted_not_crashed_on(
    converter: ModuleType, tmp_path: Path
) -> None:
    """``csv.DictReader`` fills a short row's missing fields with ``None``.

    A 972-column published supplement can easily have one, and ``_boolean(None)`` used to
    raise ``AttributeError`` mid-conversion instead of the row being counted unparseable.
    """
    directory = _write_archive(tmp_path, converter)
    path = directory / "fas_ex6_snvs_scored.txt"
    lines = path.read_text(encoding="utf-8").splitlines()
    lines.append("truncated_row\tFAS")  # far fewer fields than the header
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    _, counts = converter.convert(directory)
    assert counts["skipped_unparseable"] == 1
    assert counts["FAS"] == 4


def test_a_whitespace_only_key_falls_back_instead_of_emitting_an_empty_id(
    converter: ModuleType, tmp_path: Path
) -> None:
    """Truthiness was tested before stripping, so '   ' defeated the fallback.

    The row was written with an empty ``variant_id``, which makes the **whole** converted
    panel unreadable rather than just that row.
    """
    from bt4.api import read_variant_panel

    directory = _write_archive(tmp_path, converter)
    path = directory / "wt1_ex9_scored.txt"
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = lines[1].split("\t")
    fields[0] = "   "
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    rows, _ = converter.convert(directory)
    assert all(row["variant_id"].strip() for row in rows)
    out = tmp_path / "v.tsv"
    converter.write_panel(rows, out)
    assert len(read_variant_panel(out, negative_construction="x")) == len(rows)


@pytest.mark.parametrize("null", ["NA", " na ", "NaN", "None", "-", ".", ""])
def test_every_null_spelling_reads_as_missing_not_as_a_number(
    converter: ModuleType, tmp_path: Path, null: str
) -> None:
    """A missing score that parsed as a number would corrupt the gate silently."""
    directory = _write_archive(tmp_path, converter)
    path = directory / "fas_ex6_snvs_scored.txt"
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = lines[1].split("\t")
    fields[4] = null  # DS_maxm
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    rows, _ = converter.convert(directory)
    blanked = [r for r in rows if r["group"] == "FAS" and r["spliceai_masked"] == ""]
    assert len(blanked) == 1, null
