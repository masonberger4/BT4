"""Tests for GenBank I/O (:mod:`bt4.io.genbank`).

The writer's job is to put BT4's honesty where a bench scientist actually looks:
residual violations become ``misc_feature`` spans on the map, positioned where
they occur. These tests pin that behaviour, the determinism the provenance story
needs, and round-tripping through both BT4's reader and (when installed) the
reference Biopython parser.
"""

from __future__ import annotations

import io

import pytest

from bt4 import api
from bt4.domain.context import CIRCULAR, ConstructContext
from bt4.io.genbank import (
    context_from_genbank,
    parse_genbank,
    write_genbank,
)

# A Trp run is a genuine feasibility floor: TGG is the only Trp codon, so the
# repeat it creates cannot be removed by ANY synonymous choice. That makes it the
# honest way to produce residual violations the writer must annotate.
_FLOOR_PROTEIN = "MWWWWWWWWKLDE"
_FLOOR_CONFIG = api.OptimizeConfig(max_homopolymer=None, max_repeat_length=6)


def _floor_result() -> api.Result:
    return api.optimize(_FLOOR_PROTEIN, _FLOOR_CONFIG)


def test_writer_emits_a_parsable_record() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    text = write_genbank(result, locus="demo")
    record = parse_genbank(text)
    assert record.locus == "demo"
    assert record.sequence == result.dna
    cds = record.feature("CDS")
    assert cds is not None
    assert (cds.start, cds.end) == (0, len(result.dna))
    assert text.endswith("//\n")


def test_residual_violations_are_annotated_on_the_map() -> None:
    """The headline: a defect the optimizer could not remove is visible as a feature."""
    result = _floor_result()
    assert result.metrics.hard_violations > 0, "premise: this protein has a floor"
    record = parse_genbank(write_genbank(result, locus="floor"))
    notes = [
        dict(f.qualifiers).get("note", "")
        for f in record.features
        if f.key == "misc_feature"
    ]
    assert any("BT4 max_repeat" in note for note in notes)
    assert any("HARD" in note for note in notes)


def test_overlapping_findings_merge_but_keep_their_count() -> None:
    """Merging is a rendering choice; the number of findings is still stated."""
    result = _floor_result()
    text = write_genbank(result, locus="floor")
    record = parse_genbank(text)
    repeat_features = [
        f
        for f in record.features
        if f.key == "misc_feature"
        and "max_repeat" in str(dict(f.qualifiers).get("note", ""))
    ]
    # Far fewer features than raw violations -- but the count is not lost.
    assert len(repeat_features) < result.metrics.hard_violations
    assert f"{result.metrics.hard_violations} hard" in text
    assert any(
        "overlapping findings" in str(dict(f.qualifiers).get("note", ""))
        for f in repeat_features
    )


def test_clean_run_annotates_no_violations() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    assert result.metrics.hard_violations == 0
    text = write_genbank(result)
    assert "BT4 max_repeat" not in text
    assert "0 hard, 0 soft" in text


def test_record_is_deterministic_and_carries_no_timestamp() -> None:
    # Invariant #7: identical input must give byte-identical output, so the LOCUS
    # date field (a clock reading) is deliberately absent.
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    first = write_genbank(result, locus="demo")
    second = write_genbank(result, locus="demo")
    assert first == second
    locus_line = first.splitlines()[0]
    assert "SYN" in locus_line
    # No dd-MMM-yyyy stamp anywhere in the header.
    assert not any(month in first for month in ("JAN-", "FEB-", "DEC-"))


def test_comment_carries_the_provenance_stamp() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    text = write_genbank(result)
    assert "config_hash" in text
    assert "Optimality: proven_optimal" in text
    assert "reference set: highly_expressed" in text


# --------------------------------------------------------------------------- #
# Construct context: the whole-construct map.
# --------------------------------------------------------------------------- #


def test_whole_construct_places_the_cds_after_the_leader() -> None:
    context = ConstructContext(upstream="GGCACCAGGTACC", downstream="TAAGCGGCCGC")
    result = api.optimize(
        "MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None, context=context)
    )
    record = parse_genbank(write_genbank(result, context=context, locus="construct"))
    assert record.sequence == context.assemble(result.dna)
    cds = record.feature("CDS")
    assert cds is not None
    assert cds.start == context.cds_offset
    assert cds.end == context.cds_offset + len(result.dna)
    # The supplied flanks are labelled as not-designed-by-BT4.
    notes = " ".join(
        str(dict(f.qualifiers).get("note", "")) for f in record.features
    )
    assert "supplied 5' construct context" in notes
    assert "supplied 3' construct context" in notes


def test_whole_construct_can_be_switched_off() -> None:
    context = ConstructContext(upstream="GGCACCAGGTACC")
    result = api.optimize(
        "MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None, context=context)
    )
    record = parse_genbank(
        write_genbank(result, context=context, whole_construct=False)
    )
    assert record.sequence == result.dna


def test_circular_topology_survives_the_round_trip() -> None:
    context = ConstructContext(upstream="GGCACCAGGTACC", topology=CIRCULAR)
    result = api.optimize(
        "MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None, context=context)
    )
    record = parse_genbank(write_genbank(result, context=context))
    assert record.topology == CIRCULAR


# --------------------------------------------------------------------------- #
# Reading: a vector map becomes the context a design is optimized inside.
# --------------------------------------------------------------------------- #


def test_context_from_genbank_splits_at_the_insertion_point() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    text = write_genbank(result, locus="vector")
    context = context_from_genbank(text, insertion_point=9)
    assert context.upstream == result.dna[:9]
    assert context.downstream == result.dna[9:]


def test_context_from_genbank_defaults_to_the_cds_feature() -> None:
    ctx = ConstructContext(upstream="GGCACCAGGTACC", downstream="TAAGCGGCCGC")
    result = api.optimize(
        "MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None, context=ctx)
    )
    text = write_genbank(result, context=ctx)
    # With no explicit insertion point the CDS feature's start is used, so the
    # recovered context is the original leader.
    recovered = context_from_genbank(text)
    assert recovered.upstream == ctx.upstream


def test_context_from_genbank_can_trim_the_flanks() -> None:
    result = api.optimize("MAALKHETQWY", api.OptimizeConfig(max_homopolymer=None))
    text = write_genbank(result, locus="vector")
    context = context_from_genbank(text, insertion_point=18, upstream_nt=6)
    assert len(context.upstream) == 6
    assert context.upstream == result.dna[12:18]


def test_parse_rejects_a_non_genbank_file() -> None:
    with pytest.raises(ValueError, match="LOCUS"):
        parse_genbank(">not a genbank\nACGT\n")


def test_parse_rejects_a_record_with_no_sequence() -> None:
    with pytest.raises(ValueError, match="ORIGIN"):
        parse_genbank("LOCUS       empty       0 bp    DNA     linear   SYN\n//\n")


# --------------------------------------------------------------------------- #
# Reference-parser compatibility (skipped unless Biopython is installed).
# --------------------------------------------------------------------------- #


def test_biopython_parses_the_record() -> None:
    SeqIO = pytest.importorskip("Bio.SeqIO")
    context = ConstructContext(upstream="GGCACCAGGTACC", downstream="TAAGCGGCCGC")
    result = api.optimize(
        _FLOOR_PROTEIN,
        api.OptimizeConfig(max_homopolymer=None, max_repeat_length=6, context=context),
    )
    text = write_genbank(result, context=context, locus="floor")
    record = SeqIO.read(io.StringIO(text), "genbank")
    assert str(record.seq).upper() == context.assemble(result.dna)
    kinds = {feature.type for feature in record.features}
    assert {"source", "CDS", "misc_feature"} <= kinds
    # The residual findings survive into the reference parser's feature table.
    notes = " ".join(
        " ".join(f.qualifiers.get("note", [])) for f in record.features
    )
    assert "BT4 max_repeat" in notes
