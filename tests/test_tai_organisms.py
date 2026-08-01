"""Tests for the additional eukaryote tAI tables (mouse and yeast).

These pin that the two GtRNAdb-derived tRNA tables bundled alongside the human
one load through the same ``load_tai_table`` path, produce a well-formed
relative-adaptiveness vector (the 60 scoreable codons, every ``w`` in ``(0, 1]``
with a maximum of exactly ``1.0``), are discoverable via
``available_tai_organisms``, and carry a provenance sidecar whose declared
``total_genes`` equals the summed anticodon counts of the TSV it stamps
(invariant #9: provenance describes the actual bundled contents).
"""

from __future__ import annotations

import pytest

from bt4.biomodels.codon.tai import (
    available_tai_organisms,
    load_tai_provenance,
    load_tai_table,
)

_NEW_ORGANISMS = (
    "mus_musculus",
    "saccharomyces_cerevisiae",
    "rattus_norvegicus",
    "danio_rerio",
    "drosophila_melanogaster",
    "caenorhabditis_elegans",
    "arabidopsis_thaliana",
)


@pytest.mark.parametrize("organism", _NEW_ORGANISMS)
def test_new_organism_table_is_well_formed(organism: str) -> None:
    """Each new table yields 60 codons with every ``w`` in ``(0, 1]``, max 1.0."""
    table = load_tai_table(organism)
    w = table.relative_adaptiveness()
    assert len(w) == 60
    assert "ATG" not in w  # Met excluded from tAI
    assert all(v > 0.0 for v in w.values())  # zero-W codons filled, never annihilate
    assert max(w.values()) == pytest.approx(1.0)
    assert max(w.values()) <= 1.0
    assert all(v <= 1.0 for v in w.values())


@pytest.mark.parametrize("organism", _NEW_ORGANISMS)
def test_new_organism_is_discoverable(organism: str) -> None:
    """Both new organisms appear in the auto-discovered organism listing."""
    assert organism in available_tai_organisms()


@pytest.mark.parametrize("organism", _NEW_ORGANISMS)
def test_provenance_total_genes_matches_tsv_sum(organism: str) -> None:
    """The sidecar ``total_genes`` equals the summed anticodon counts it stamps."""
    table = load_tai_table(organism)
    prov = load_tai_provenance(organism)
    assert prov.total_genes == sum(table.anticodon_counts.values())
    assert len(prov.sha256) == 64
    assert "GtRNAdb" in prov.source


@pytest.mark.parametrize("organism", _NEW_ORGANISMS)
def test_provenance_is_citation_gated_not_public_domain(organism: str) -> None:
    """GtRNAdb has no explicit data license; the note must say so, honestly."""
    note = load_tai_provenance(organism).note
    assert "citation-gated" in note
    assert "tRNAscan-SE" in note


def test_known_new_organism_gene_totals() -> None:
    """Pin the re-counted GtRNAdb totals so a table swap changes this assertion."""
    totals = {
        "mus_musculus": 404,
        "saccharomyces_cerevisiae": 273,
        "rattus_norvegicus": 416,
        "danio_rerio": 8846,
        "drosophila_melanogaster": 289,
        "caenorhabditis_elegans": 580,
        "arabidopsis_thaliana": 584,
    }
    for organism, expected in totals.items():
        assert sum(load_tai_table(organism).anticodon_counts.values()) == expected, organism
