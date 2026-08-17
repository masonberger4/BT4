"""Tests for the tRNA-adaptation-index (tAI) table, term, and pipeline wiring.

These pin the faithfulness of the ``get.ws`` port (against the dos Reis E. coli
reference and structural properties), the ``delta == score`` honesty of the
:class:`~bt4.objectives.tai.TaiTerm` (invariant #4), and that a real, bundled,
provenance-hashed human tRNA table flows through ``bt4.api``.
"""

from __future__ import annotations

import math

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4 import api
from bt4.biomodels.codon.tables import available_organisms
from bt4.biomodels.codon.tai import (
    DOSREIS_2004_S,
    TaiTable,
    available_tai_organisms,
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


def test_bacterial_ata_uses_the_reference_constant_not_a_copy_number() -> None:
    """sking=1 sets W[ATA] to the bare lysidine contribution, per the reference.

    The reference get.ws is `if(sking == 1) W[35] = p[9]` -- a constant, with NO
    tRNA factor. That is deliberate: bacteria decode AUA with tRNA-Ile2, whose
    anticodon is CAU (lysidine-modified), *not* UAU. So the copy-number slot this
    codon would index (anticodon TAT) is empty in real bacteria -- E. coli K-12 has
    zero TAT genes -- and scaling by it would drive W[ATA] to 0, which the
    zero-filling step then replaces with the geometric mean: an arbitrary value.

    This test previously asserted the opposite (that ATA scales with TAT copies).
    That encoded a misreading of the reference and would have shipped a wrong Ile
    weight the moment a bacterial table existed, which it now does.
    """
    # W[ATA] does not move with the (biologically empty) TAT slot...
    none_ = build_tai_weights({"AGC": 20}, sking=1)
    some = build_tai_weights({"TAT": 40, "AGC": 20}, sking=1)
    assert none_["ATA"] == pytest.approx(some["ATA"])
    # ...and it equals the lysidine contribution p[8], normalized by max(W).
    p8 = 1.0 - DOSREIS_2004_S[8]
    assert some["ATA"] == pytest.approx(p8 / (1.0 * 20))
    # Under the eukaryotic model the same counts leave ATA with no reader at all,
    # so the two branches genuinely differ.
    assert build_tai_weights({"AGC": 20}, sking=0)["ATA"] != pytest.approx(none_["ATA"])


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
    # Every *bundled* organism now ships a tRNA table (E. coli was the last gap),
    # so the refusal is exercised with an organism that has no table at all rather
    # than one that merely used to lack one.
    with pytest.raises(ValueError, match="no bundled tAI table"):
        load_tai_table("nonexistent_organism")


def test_every_selectable_organism_has_a_tai_table() -> None:
    """tAI must be offered for every organism a user can actually pick.

    The inverse failure shipped for a long time: six organisms had authentic tRNA
    tables that were unreachable because they had no codon table. This asserts the
    other direction now that E. coli closes the last gap -- a codon table without
    tRNA data means tAI is silently unavailable exactly where a user asked for it.
    """
    assert set(available_organisms()) <= set(available_tai_organisms())


# --------------------------------------------------------------------------- #
# The bundled E. coli table: real GtRNAdb data, and the first bacterial one.
# --------------------------------------------------------------------------- #


def test_ecoli_table_matches_its_documented_totals() -> None:
    """The shipped counts must agree with the provenance that describes them."""
    table = load_tai_table("escherichia_coli")
    provenance = load_tai_provenance("escherichia_coli")
    assert sum(table.anticodon_counts.values()) == provenance.total_genes
    # 86 tRNA genes is the standard count for K-12 MG1655 (89 predictions less one
    # selenocysteine and two undetermined) -- an external cross-check, not a
    # restatement of our own arithmetic.
    assert provenance.total_genes == 86


def test_ecoli_is_scored_under_the_prokaryotic_model() -> None:
    """Super-kingdom comes from the table's provenance, never a hardcoded default."""
    assert load_tai_table("escherichia_coli").sking == 1
    assert load_tai_table("homo_sapiens").sking == 0


def test_ecoli_lysidine_path_gives_ata_a_real_weight() -> None:
    """The bacterial branch must not leave Ile ATA at the zero-filled geometric mean.

    E. coli has zero TAT-anticodon genes, so before the reference fix this codon's
    W was 0 and got replaced by the geometric mean of every other weight.
    """
    weights = load_tai_table("escherichia_coli").relative_adaptiveness()
    other = [w for codon, w in weights.items() if codon != "ATA"]
    geometric_mean = math.exp(sum(math.log(w) for w in other) / len(other))
    # The zero-filled value is exactly the geometric mean; the real lysidine weight
    # is not. (Its exact value, p[8]/max(W), is pinned in the unit test above, where
    # the inputs make max(W) known.)
    assert weights["ATA"] != pytest.approx(geometric_mean)
    assert 0.0 < weights["ATA"] < 1.0
    # Under the eukaryotic model the same counts would leave ATA unread entirely,
    # so the bacterial branch is genuinely doing something here.
    eukaryotic = build_tai_weights(load_tai_table("escherichia_coli").anticodon_counts)
    assert eukaryotic["ATA"] != pytest.approx(weights["ATA"])


def test_ecoli_provenance_pins_its_upstream_source() -> None:
    """Invariant #9: the table must be re-derivable from a pinned public source."""
    import json
    from importlib.resources import files

    raw = json.loads(
        files("bt4.biomodels.codon.data")
        .joinpath("escherichia_coli.trna.provenance.json")
        .read_text(encoding="utf-8")
    )
    assert raw["source_url"].startswith("https://gtrnadb.ucsc.edu/")
    assert len(raw["source_sha256"]) == 64
    assert raw["super_kingdom"] == "bacteria"
