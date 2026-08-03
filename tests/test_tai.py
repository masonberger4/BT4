"""Tests for the tRNA-adaptation-index (tAI) table, term, and pipeline wiring.

These pin the faithfulness of the ``get.ws`` port (against the dos Reis E. coli
reference and structural properties), the ``delta == score`` honesty of the
:class:`~bt4.objectives.tai.TaiTerm` (invariant #4), and that a real, bundled,
provenance-hashed human tRNA table flows through ``bt4.api``.
"""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4 import api
from bt4.biomodels.codon.tai import (
    DOSREIS_2004_S,
    TaiTable,
    build_tai_weights,
    load_tai_provenance,
    load_tai_table,
)
from bt4.domain.genetic_code import AMINO_ACIDS, STOP, synonymous_codons, translate
from bt4.objectives.base import iter_codons
from bt4.objectives.tai import TaiTerm

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=40)


def _synonymous_backtranslate(protein: str) -> str:
    """First-codon-per-residue back-translation (for building test sequences)."""
    return "".join(synonymous_codons(aa)[0] for aa in [*protein, STOP])


# --------------------------------------------------------------------------- #
# Port faithfulness (get.ws) and table shape.
# --------------------------------------------------------------------------- #


def test_build_weights_shape_and_range() -> None:
    table = load_tai_table("homo_sapiens")
    w = table.relative_adaptiveness()
    # 60 scoreable codons: 61 sense codons minus Met (stops already excluded).
    assert len(w) == 60
    assert "ATG" not in w and "TAA" not in w and "TAG" not in w and "TGA" not in w
    assert max(w.values()) == pytest.approx(1.0)
    assert min(w.values()) > 0.0  # zero-W codons filled with the geometric mean


def test_ecoli_reference_optimal_codons_are_high() -> None:
    # dos Reis E. coli K-12 (sking=1): a tiny hand-built tGCN where Lys is read
    # only by a WC tRNA must give AAA the maximal w (== 1.0). This exercises the
    # bacterial path and the WC decoding.
    counts = {"TTT": 1}  # anticodon TTT reads Lys AAA by Watson-Crick (U:A)
    w = build_tai_weights(counts, sking=1)
    assert w["AAA"] == pytest.approx(1.0)


def test_bacterial_ata_scales_with_reader_copy_number() -> None:
    # sking=1 lysidine reading of Ile ATA must scale with its reader tRNA (anticodon
    # TAT) copy number, not be a bare constant (dos Reis W[ATA] = p * tRNA[ATA]).
    low = build_tai_weights({"TAT": 1, "AGC": 20}, sking=1)
    high = build_tai_weights({"TAT": 40, "AGC": 20}, sking=1)
    # More ATA-reading tRNA genes -> higher relative adaptiveness for ATA.
    assert high["ATA"] > low["ATA"]


def test_svalues_must_have_nine_entries() -> None:
    with pytest.raises(ValueError, match="9 entries"):
        build_tai_weights({"TTT": 1}, s_values=(0.0, 0.0))
    assert len(DOSREIS_2004_S) == 9


def test_wobble_penalty_lowers_weight() -> None:
    # A single Ala WC tRNA (anticodon AGC reads GCT). GCT (WC) should score above
    # a purely wobble-read synonymous codon in the same box.
    w = build_tai_weights({"AGC": 1})
    assert w["GCT"] > w["GCA"]  # GCT read WC; GCA only via I:A (heavy penalty)


# --------------------------------------------------------------------------- #
# TaiTerm: invariant #4 (delta == score) and orientation.
# --------------------------------------------------------------------------- #


@given(protein=_PROTEIN)
@settings(max_examples=50, deadline=None)
def test_tai_term_delta_equals_score(protein: str) -> None:
    term = TaiTerm(load_tai_table("homo_sapiens").relative_adaptiveness())
    dna = _synonymous_backtranslate(protein)
    summed = sum(term.delta("", codon, pos) for pos, codon in iter_codons(dna))
    assert term.score(dna) == pytest.approx(summed)


def test_tai_term_zero_for_met_and_stop() -> None:
    term = TaiTerm(load_tai_table("homo_sapiens").relative_adaptiveness())
    assert term.delta("", "ATG", 0) == 0.0  # Met excluded from tAI
    assert term.delta("", "TAA", 0) == 0.0  # stop excluded
    assert term.context_len() == 0


# --------------------------------------------------------------------------- #
# Table + provenance.
# --------------------------------------------------------------------------- #


def test_human_table_is_real_and_provenanced() -> None:
    table = load_tai_table("homo_sapiens")
    assert isinstance(table, TaiTable)
    assert sum(table.anticodon_counts.values()) == 431  # verified GtRNAdb hg38 count
    prov = load_tai_provenance("homo_sapiens")
    assert prov.total_genes == 431
    assert len(prov.sha256) == 64
    assert "GtRNAdb" in prov.source


def test_tai_of_sequence_in_unit_interval() -> None:
    table = load_tai_table("homo_sapiens")
    result = api.optimize("MAALKHETQWSNDECFGR")
    assert 0.0 < table.tai(result.dna) <= 1.0


def test_unknown_organism_has_no_tai_table() -> None:
    with pytest.raises(ValueError, match="no bundled tAI table"):
        load_tai_table("nonexistent_species")


def test_bundled_tai_organisms_load() -> None:
    # Human ships from Phase 2; mouse and yeast were added with real GtRNAdb data.
    for organism in ("homo_sapiens", "mus_musculus", "saccharomyces_cerevisiae"):
        table = load_tai_table(organism)
        w = table.relative_adaptiveness()
        assert len(w) == 60
        assert max(w.values()) == pytest.approx(1.0)
        assert min(w.values()) > 0.0


# --------------------------------------------------------------------------- #
# Pipeline / api wiring.
# --------------------------------------------------------------------------- #


def test_tai_weight_participates_and_reports() -> None:
    result = api.optimize("MAALKHETQWSNDECFGRPVIY", api.OptimizeConfig(tai_weight=2.0))
    assert "tai_logw" in result.metrics.objective.terms()
    assert "tai" in result.audit
    assert 0.0 < float(result.audit["tai"]) <= 1.0
    assert translate(result.dna) == result.protein + "*"
    assert result.certificate.is_proven_optimal


def test_tai_weight_hashes_trna_table_into_manifest() -> None:
    result = api.optimize("MAAL", api.OptimizeConfig(tai_weight=1.0))
    manifest = result.audit["manifest"]
    assert isinstance(manifest, dict)
    assert "trna_table_sha256" in manifest["inputs"]


def test_tai_weight_on_organism_without_data_raises() -> None:
    with pytest.raises(ValueError, match="no bundled tAI table"):
        api.optimize("MAAL", api.OptimizeConfig(organism="escherichia_coli", tai_weight=1.0))
