"""Tests for FASTA parsing and building codon tables from a CDS set."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bt4.biomodels.codon.build import build_table, count_codons, write_table
from bt4.biomodels.codon.tables import load_table_from_file, sha256_hex
from bt4.domain.genetic_code import CODON_TABLE
from bt4.io.fasta import parse_fasta, read_fasta

# --- parse_fasta / read_fasta ------------------------------------------------


def test_parse_fasta_multi_record_multiline() -> None:
    text = (
        ">seq1 first record\n"
        "ATG AAA\n"  # internal whitespace removed
        "cccggg\n"  # lower-case upper-cased
        "\n"  # blank line ignored
        ">  seq2\n"  # header trimmed
        "TTTTTT\n"
    )
    records = parse_fasta(text)
    assert records == [
        ("seq1 first record", "ATGAAACCCGGG"),
        ("seq2", "TTTTTT"),
    ]
    # The first record's sequence round-trips its (normalized) source lines.
    assert records[0][1] == "ATGAAACCCGGG"


def test_parse_fasta_sequence_before_header_raises() -> None:
    with pytest.raises(ValueError):
        parse_fasta("ATGAAA\n>seq1\nATGAAA\n")


def test_parse_fasta_empty_record_raises() -> None:
    # seq1 has no sequence lines before the next header.
    with pytest.raises(ValueError):
        parse_fasta(">seq1\n>seq2\nATGAAA\n")


def test_read_fasta_from_file(tmp_path: Path) -> None:
    p = tmp_path / "in.fasta"
    p.write_text(">a\nATGAAA\n>b\nGGGCCC\n", encoding="utf-8")
    assert read_fasta(p) == [("a", "ATGAAA"), ("b", "GGGCCC")]


# --- count_codons ------------------------------------------------------------


def test_count_codons_exact() -> None:
    # CDS1: ATG AAA ATG -> ATG x2, AAA x1
    # CDS2: AAA GGG     -> AAA x1, GGG x1
    counts = count_codons(["ATGAAAATG", "AAAGGG"])
    assert counts == {"ATG": 2, "AAA": 2, "GGG": 1}


def test_count_codons_non_multiple_of_three_raises() -> None:
    with pytest.raises(ValueError, match="5"):
        count_codons(["ATGAA"])  # length 5


# --- build_table -------------------------------------------------------------


def test_build_table_weights_and_full_coverage() -> None:
    # Leucine: CTG x3, CTT x1 -> CTG dominates. Alanine: GCT x1, GCC x1.
    cds = ["CTGCTGCTGCTT", "GCTGCC"]
    table, counts = build_table(cds, organism="toy")

    # Raw counts are returned unsmoothed.
    assert counts["CTG"] == 3
    assert counts["CTT"] == 1
    assert counts["GCT"] == 1
    assert counts["GCC"] == 1

    # The most-frequent codon within an amino acid has weight 1.0.
    assert table.weight("CTG") == pytest.approx(1.0)
    # Smoothed: CTT=(1+1)=2 over CTG=(3+1)=4 -> 0.5.
    assert table.weight("CTT") == pytest.approx(0.5)

    # Laplace smoothing covers all 64 codons.
    assert set(table.frequency) == set(CODON_TABLE)

    # CAI is a genuine geometric mean in (0, 1].
    cai = table.cai("CTTCTGGCC")
    assert 0.0 < cai <= 1.0


def test_build_table_returns_codon_usage_table_and_counts() -> None:
    table, counts = build_table(["ATGAAA"], organism="tiny")
    assert table.organism == "tiny"
    assert counts == {"ATG": 1, "AAA": 1}


# --- write_table -------------------------------------------------------------


def test_write_table_roundtrips_via_load_table_from_file(tmp_path: Path) -> None:
    # Full coverage so the written TSV loads; a known Alanine skew for spot-check.
    counts: dict[str, int] = dict.fromkeys(CODON_TABLE, 1)
    counts["GCC"] = 5  # Alanine group maximum
    counts["GCT"] = 2  # another Alanine codon

    tsv_path_str = write_table(
        counts,
        organism="toy",
        path=tmp_path,
        source="unit-test CDS set",
        cds_count=3,
    )
    tsv_path = Path(tsv_path_str)
    assert tsv_path.name == "toy.tsv"
    assert tsv_path.parent == tmp_path

    # The raw count is written verbatim as the frequency.
    assert "A\tGCC\t5\n" in tsv_path.read_text(encoding="utf-8")

    # An equivalent table loads back with the expected weights.
    table = load_table_from_file(tsv_path)
    assert table.weight("GCC") == pytest.approx(1.0)  # group max
    assert table.weight("GCT") == pytest.approx(2 / 5)


def test_write_table_provenance_sha256_matches_disk(tmp_path: Path) -> None:
    counts: dict[str, int] = dict.fromkeys(CODON_TABLE, 1)
    write_table(
        counts,
        organism="toy",
        path=tmp_path,
        source="unit-test CDS set",
        cds_count=7,
    )
    tsv_path = tmp_path / "toy.tsv"
    provenance = json.loads(
        (tmp_path / "toy.provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["sha256"] == sha256_hex(tsv_path.read_bytes())
    assert provenance["source"] == "unit-test CDS set"
    assert provenance["cds_count"] == 7
    assert isinstance(provenance["retrieved"], str)
    assert provenance["retrieved"]


def test_write_table_accepts_build_note_and_extra(tmp_path: Path) -> None:
    """A recount caller can override build/note and attach a re-derivation trail."""
    counts: dict[str, int] = dict.fromkeys(CODON_TABLE, 3)
    write_table(
        counts,
        organism="toy",
        path=tmp_path,
        source="pinned CDS set",
        cds_count=42,
        build="recounted from a pinned FASTA",
        note="Real genome-wide codon counts.",
        extra={"source_url": "https://example.org/toy.fa.gz", "source_sha256": "ab" * 32},
    )
    prov = json.loads((tmp_path / "toy.provenance.json").read_text(encoding="utf-8"))
    assert prov["build"] == "recounted from a pinned FASTA"
    assert prov["note"] == "Real genome-wide codon counts."
    assert prov["source_url"] == "https://example.org/toy.fa.gz"
    # Reserved keys still come from their own parameters, not from extra.
    assert prov["source"] == "pinned CDS set"
    assert prov["cds_count"] == 42


def test_write_table_extra_cannot_shadow_reserved_keys(tmp_path: Path) -> None:
    """`extra` must not be able to make the sidecar disagree with itself."""
    counts: dict[str, int] = dict.fromkeys(CODON_TABLE, 1)
    with pytest.raises(ValueError, match="reserved"):
        write_table(
            counts,
            organism="toy",
            path=tmp_path,
            source="s",
            extra={"sha256": "deadbeef", "note": "spoofed"},
        )
