"""Tests for the **highly-expressed** codon reference tables.

These are the tables CAI is actually defined on: relative adaptiveness counted
over the 300 most abundant proteins per organism (PaxDb v6.1 whole-organism
integrated abundances, joined to the same release-pinned Ensembl CDS the
genome-wide tables are counted from), built by
``scripts/build_highly_expressed_tables.py``. They are the **default** reference
set for the eight organisms that have one.

Three kinds of claim are checked here, and only the third is really about BT4:

* **Completeness and refusal** -- every shipped table observes all 64 codons (so
  none needed smoothing, i.e. none contains an invented number), and the one
  organism without an abundance join, *A. thaliana*, has no table and says so
  rather than quietly falling back.
* **A re-derivation trail** -- pinned URLs and digests for all three sources
  (abundances, protein-to-gene mapping, CDS), the reference-set size, a digest of
  the ranked gene roster, and a join/filter tally that closes exactly.
* **External ground truth** (CLAUDE.md §8, "not just self-consistency"). Two
  independent checks that a fabricated or mis-joined table would fail:
  the classic *E. coli* and *S. cerevisiae* optimal codons must come out of the
  counts, and codon bias must be **stronger** in the highly-expressed set than
  genome-wide in every organism -- which is the entire premise of using a
  highly-expressed reference set at all.
"""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from bt4 import api
from bt4.biomodels.codon.tables import (
    GENOME_WIDE,
    HIGHLY_EXPRESSED,
    REFERENCE_SETS,
    CodonUsageTable,
    available_organisms,
    available_reference_sets,
    default_reference_set,
    load_provenance,
    load_table,
    sha256_hex,
)
from bt4.domain.genetic_code import CODON_TABLE, STOP

#: Organisms with a bundled highly-expressed table.
WITH_ABUNDANCE: tuple[str, ...] = (
    "caenorhabditis_elegans",
    "danio_rerio",
    "drosophila_melanogaster",
    "escherichia_coli",
    "homo_sapiens",
    "mus_musculus",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
)

#: Organisms deliberately shipped WITHOUT a highly-expressed table. Each entry
#: records the *measured* reason, because "we did not get to it" and "the
#: evidence does not support one" are different claims and only the second
#: justifies an absence:
#:
#: * *A. thaliana* -- PaxDb identifies its proteins by UniProt accession, which
#:   the pinned Ensembl Plants annotation does not carry, so the join would need
#:   an unpinned external mapping.
#: * *C. griseus* (CHO) -- PaxDb v6.1 holds **no whole-organism integrated**
#:   dataset for taxon 10029, only a single study
#:   (``10029-PXD014877_Mueller_Nature_2020``). Every table above is built from
#:   the integrated set; one built from a single study would carry the same
#:   ``highly_expressed`` label while meaning something materially different
#:   (one lab, one condition, no cross-study integration). Distinguishing them
#:   needs the reference-set label to carry its evidence class, which is a
#:   design change and not a data addition.
#: * *K. phaffii* -- PaxDb v6.1 has no dataset directory for taxon 644223 at all.
#:
#: *B. subtilis* is the near-miss worth recording: the integrated dataset **does**
#: exist (taxon 224308) and the join is available -- PaxDb writes ``BSU35360``
#: where Ensembl Bacteria writes ``BSU_35360``, and that declared ``^BSU`` ->
#: ``BSU_`` rewrite joins 4,042/4,052 = 99.8%. That is a locus-tag punctuation
#: difference derivable from the two pinned files alone, not a third-party
#: mapping, so it is admissible -- it just needs the builder to support a
#: declared per-spec rewrite, which is a separate change from adding the
#: genome-wide tables.
WITHOUT_ABUNDANCE: tuple[str, ...] = (
    "arabidopsis_thaliana",
    "bacillus_subtilis",
    "cricetulus_griseus_chok1gshd",
    "komagataella_phaffii",
)

#: The reference-set size every shipped table was built at.
TOP_N = 300


def _he(name: str) -> CodonUsageTable:
    return load_table(name, reference_set=HIGHLY_EXPRESSED)


def _raw_provenance(name: str) -> dict[str, object]:
    """The sidecar as raw JSON (``TableProvenance`` keeps only its own fields)."""
    resource = files("bt4.biomodels.codon.data").joinpath(
        f"{name}.highly_expressed.provenance.json"
    )
    return dict(json.loads(resource.read_text(encoding="utf-8")))


def _top_codon(table: CodonUsageTable, amino_acid: str) -> str:
    synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
    return max(synonyms, key=table.weight)


def _mean_dominance(table: CodonUsageTable) -> float:
    """Mean share taken by each degenerate amino acid's most-used codon.

    A simple, scale-free strength-of-codon-bias readout: 1/n for a perfectly
    even n-fold group, rising toward 1.0 as one synonym takes over.
    """
    by_aa: dict[str, list[float]] = {}
    for codon, value in table.frequency.items():
        by_aa.setdefault(CODON_TABLE[codon], []).append(value)
    shares = [
        max(values) / sum(values)
        for aa, values in by_aa.items()
        if aa != STOP and len(values) > 1
    ]
    return sum(shares) / len(shares)


# --------------------------------------------------------------------------- #
# What is bundled, what is default, and what is honestly absent.
# --------------------------------------------------------------------------- #


def test_every_organism_is_accounted_for() -> None:
    """The two lists above must together be exactly the bundled organisms."""
    assert set(WITH_ABUNDANCE) | set(WITHOUT_ABUNDANCE) == set(available_organisms())


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_highly_expressed_is_the_default(name: str) -> None:
    """Where it exists, the reference set CAI is defined on is what you get."""
    assert default_reference_set(name) == HIGHLY_EXPRESSED
    assert available_reference_sets(name) == (HIGHLY_EXPRESSED, GENOME_WIDE)
    assert load_table(name).reference_set == HIGHLY_EXPRESSED


@pytest.mark.parametrize("name", WITHOUT_ABUNDANCE)
def test_missing_reference_set_refuses_rather_than_substituting(name: str) -> None:
    """No silent fallback: the two tables answer different questions.

    Quietly handing back the genome-wide table would make a caller's CAI mean
    something other than what they asked for while still looking like success.
    """
    assert available_reference_sets(name) == (GENOME_WIDE,)
    assert default_reference_set(name) == GENOME_WIDE
    with pytest.raises(ValueError, match="no 'highly_expressed' codon table"):
        load_table(name, reference_set=HIGHLY_EXPRESSED)
    with pytest.raises(ValueError, match="no 'highly_expressed' codon table"):
        load_provenance(name, reference_set=HIGHLY_EXPRESSED)


def test_unknown_reference_set_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown reference set"):
        load_table("homo_sapiens", reference_set="tissue_specific")


def test_reference_set_tables_are_not_listed_as_organisms() -> None:
    """``<organism>.highly_expressed.tsv`` must not read as its own organism."""
    assert not any("." in name for name in available_organisms())


# --------------------------------------------------------------------------- #
# The counts themselves.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_all_64_codons_are_observed(name: str) -> None:
    """No shipped table needed smoothing -- so none contains an invented number.

    This is what fixes the reference-set size at 300: the builder refuses to
    write a table with an unobserved codon rather than papering over it, and 300
    is the smallest tested size at which no organism trips that refusal.
    """
    assert set(_he(name).frequency) == set(CODON_TABLE)


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_counts_are_positive_integers(name: str) -> None:
    for codon, value in _he(name).frequency.items():
        assert value == int(value), f"{name} {codon} is not a whole count"
        assert value >= 1


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_group_max_codon_has_unit_weight(name: str) -> None:
    table = _he(name)
    assert table.weight(_top_codon(table, "L")) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# Provenance: a complete, checkable re-derivation trail.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_provenance_sha256_matches_shipped_tsv(name: str) -> None:
    raw = (
        files("bt4.biomodels.codon.data")
        .joinpath(f"{name}.highly_expressed.tsv")
        .read_bytes()
    )
    prov = load_provenance(name, reference_set=HIGHLY_EXPRESSED)
    assert prov.sha256 == sha256_hex(raw)
    assert prov.reference_set == HIGHLY_EXPRESSED
    assert prov.cds_count == TOP_N


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_provenance_carries_all_three_pinned_sources(name: str) -> None:
    """Three inputs decide this table, so all three must be pinned and stamped.

    Abundances choose *which* genes; the peptide FASTA decides which gene each
    abundance belongs to; the CDS set supplies the sequence counted. Any one of
    them swapped silently would change the table, so a trail missing any one is
    not re-derivable (CLAUDE.md §8, invariant #9).
    """
    prov = _raw_provenance(name)
    for url_key, sha_key in (
        ("source_url", "source_sha256"),
        ("protein_to_gene_url", "protein_to_gene_sha256"),
        ("cds_url", "cds_sha256"),
    ):
        assert str(prov[url_key]).startswith("https://"), f"{name}: {url_key}"
        assert len(str(prov[sha_key])) == 64, f"{name}: {sha_key}"
    assert prov["reference_set"] == HIGHLY_EXPRESSED
    assert prov["top_n"] == TOP_N
    assert len(str(prov["gene_roster_sha256"])) == 64
    assert "pax-db.org" in str(prov["source_url"])
    # Pinned release, never a moving "latest" link.
    assert "/latest/" not in str(prov["source_url"])
    assert isinstance(prov["most_abundant_genes"], list)
    assert len(prov["most_abundant_genes"]) == 20


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_join_tally_accounts_for_every_abundance_row(name: str) -> None:
    """Every PaxDb row is matched, ambiguous, or unmatched -- exactly one.

    A join that does not close is describing a different run than the one that
    produced the shipped table, and "nothing is skipped silently" is not being
    kept. Ambiguous and unmatched are counted **separately** on purpose: an
    identifier the annotation resolves two ways is a mapping the builder refused
    to guess, not one the annotation lacks.
    """
    join = _raw_provenance(name)["join"]
    assert isinstance(join, dict)
    matched = sum(v for k, v in join.items() if k.startswith("rows_matched_via_"))
    unmatched = sum(v for k, v in join.items() if k.startswith("rows_unmatched_"))
    assert matched + unmatched == join["paxdb_rows"]
    # A subset of the matched rows, NOT a fourth part of the partition: the
    # identifier resolved, the gene just has no valid representative CDS.
    assert join["rows_matched_whose_gene_has_no_counted_cds"] <= matched
    assert join["genes_eligible_for_reference_set"] >= TOP_N
    assert (
        join["genes_joined"] - join["genes_excluded_organelle_encoded"]
        == join["genes_eligible_for_reference_set"]
    )
    # rows_* and genes_* are different units and must not be mixed: several
    # protein rows can collapse onto one gene, so the two families do not
    # reconcile by subtraction. Naming them apart is what stops a reader trying.
    assert all(key.startswith(("rows_", "genes_", "paxdb_", "organelle_")) for key in join)


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_stamped_totals_match_the_shipped_counts(name: str) -> None:
    prov = _raw_provenance(name)
    assert sum(_he(name).frequency.values()) == prov["total_codons_counted"]
    # The filter tally describes the whole CDS SOURCE, not this table, and its
    # key names must say so: the genome-wide sidecar's `cds_counted` means "the
    # number counted into this table", which here would be 300, not ~20,000.
    filters = prov["cds_source_filters"]
    assert isinstance(filters, dict)
    assert "cds_counted" not in filters
    assert filters["genes_with_a_representative_cds"] > TOP_N


def test_organelle_genes_are_actually_excluded_somewhere() -> None:
    """The organelle filter must be live, not decorative.

    Mitochondria translate with a different genetic code and their own tRNA
    pool, so their codon usage is not evidence about the nuclear machinery a BT4
    design will meet. Most organelle CDS never reach this filter -- under the
    standard code they read as having internal stops and are dropped as invalid
    first -- so the per-organism counts are small (1 each for mouse, rat and
    yeast). The assertion is therefore only that the filter still fires at all:
    if this total drops to zero it has stopped matching the region names
    upstream uses, and the guarantee is silently gone.
    """
    excluded = {}
    for name in WITH_ABUNDANCE:
        join = _raw_provenance(name)["join"]
        assert isinstance(join, dict)
        excluded[name] = join["genes_excluded_organelle_encoded"]
        assert excluded[name] >= 0
        # The CDS source's own organelle-record tally is stamped beside it, so a
        # zero above can never be misread as "this organism has no
        # organelle-encoded genes" -- human's mitochondrial genes are in the
        # annotation, they just fail the standard-code validity filter first.
        assert join["organelle_records_in_cds_source"] >= excluded[name]
    assert sum(excluded.values()) > 0, excluded
    assert _raw_provenance("homo_sapiens")["join"]["organelle_records_in_cds_source"] > 0


# --------------------------------------------------------------------------- #
# External ground truth.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("name", "amino_acid", "expected"),
    [
        # The classic E. coli optimal codons (Ikemura 1981; Sharp & Li 1987).
        # Most of these are ALSO where the genome-wide table disagrees, but not
        # all: E. coli Leu (CTG) and yeast Arg (AGA) / Leu (TTG) top both tables.
        # Asserting them anyway is the point -- a highly-expressed reference set
        # has to reproduce the known optimal codon whether or not it differs from
        # the genome-wide answer, and the contrast is asserted separately below.
        ("escherichia_coli", "F", "TTC"),
        ("escherichia_coli", "I", "ATC"),
        ("escherichia_coli", "V", "GTT"),
        ("escherichia_coli", "R", "CGT"),
        ("escherichia_coli", "G", "GGT"),
        ("escherichia_coli", "L", "CTG"),
        # ...and the classic S. cerevisiae optimal codons.
        ("saccharomyces_cerevisiae", "F", "TTC"),
        ("saccharomyces_cerevisiae", "Y", "TAC"),
        ("saccharomyces_cerevisiae", "N", "AAC"),
        ("saccharomyces_cerevisiae", "K", "AAG"),
        ("saccharomyces_cerevisiae", "H", "CAC"),
        ("saccharomyces_cerevisiae", "R", "AGA"),
        ("saccharomyces_cerevisiae", "L", "TTG"),
    ],
)
def test_classic_optimal_codons_come_out_of_the_counts(
    name: str, amino_acid: str, expected: str
) -> None:
    assert _top_codon(_he(name), amino_acid) == expected


@pytest.mark.parametrize(
    ("name", "amino_acid", "expected"),
    [
        # Five E. coli amino acids and two yeast ones, each of which points the
        # OTHER way genome-wide. (E. coli diverges at eight amino acids in total;
        # these five are the ones also asserted as classic optimal codons above.)
        # Asserting the contrast, not just the value, is what shows the two
        # tables were counted over different genes rather than duplicated.
        ("escherichia_coli", "F", "TTT"),
        ("escherichia_coli", "I", "ATT"),
        ("escherichia_coli", "V", "GTG"),
        ("escherichia_coli", "R", "CGC"),
        ("escherichia_coli", "G", "GGC"),
        ("saccharomyces_cerevisiae", "K", "AAA"),
        ("saccharomyces_cerevisiae", "N", "AAT"),
    ],
)
def test_the_genome_wide_table_disagrees_where_it_should(
    name: str, amino_acid: str, expected: str
) -> None:
    genome_wide = load_table(name, reference_set=GENOME_WIDE)
    assert _top_codon(genome_wide, amino_acid) == expected


@pytest.mark.parametrize("name", WITH_ABUNDANCE)
def test_codon_bias_is_stronger_in_the_highly_expressed_set(name: str) -> None:
    """The premise of the whole exercise, checked against the data.

    Translational selection acts hardest on the genes translated most, so a
    highly-expressed reference set must show *stronger* codon bias than a
    genome-wide count. This holds in all eight organisms; if a table were built
    from the wrong genes -- or from the same genes twice -- it would not.
    """
    strong = _mean_dominance(_he(name))
    weak = _mean_dominance(load_table(name, reference_set=GENOME_WIDE))
    assert strong > weak, f"{name}: highly_expressed={strong:.4f} genome_wide={weak:.4f}"


def test_translational_selection_is_strongest_where_the_literature_says() -> None:
    """Yeast and fly must separate from human and rat by a wide margin.

    dos Reis et al. (2004) find little evidence of translational selection in
    large vertebrate genomes, while yeast and *Drosophila* are textbook cases of
    strong selection. That contrast should fall straight out of how much the
    highly-expressed reference set moves each species' codon bias -- and it is
    the kind of ordering a mis-parsed table would get wrong.
    """

    def shift(name: str) -> float:
        return _mean_dominance(_he(name)) - _mean_dominance(
            load_table(name, reference_set=GENOME_WIDE)
        )

    strong = min(shift("saccharomyces_cerevisiae"), shift("drosophila_melanogaster"))
    weak = max(shift("homo_sapiens"), shift("rattus_norvegicus"))
    assert strong > 2 * weak, f"strong={strong:.4f} weak={weak:.4f}"


# --------------------------------------------------------------------------- #
# The reference set reaches the engine, the result, and the stamp.
# --------------------------------------------------------------------------- #


def test_the_choice_changes_the_delivered_sequence() -> None:
    """Not a label: picking a reference set picks different codons."""
    protein = "MKTAYIAKQRQISFVKSHFSRQ"
    strong = api.optimize(
        protein, api.OptimizeConfig(organism="escherichia_coli", reference_set=HIGHLY_EXPRESSED)
    )
    weak = api.optimize(
        protein, api.OptimizeConfig(organism="escherichia_coli", reference_set=GENOME_WIDE)
    )
    assert strong.dna != weak.dna
    # ...and the highly-expressed run uses the classic optimal codons.
    assert "CGT" in strong.dna and "CGC" in weak.dna


def test_the_result_says_which_reference_set_it_used() -> None:
    """A CAI without its reference set is a number with no question attached."""
    for reference_set in REFERENCE_SETS:
        result = api.optimize(
            "MAALKHETQW",
            api.OptimizeConfig(organism="homo_sapiens", reference_set=reference_set),
        )
        assert result.audit["codon_reference_set"] == reference_set


def test_the_reference_set_reaches_the_manifest() -> None:
    """Two reference sets must never stamp identically (invariant #9).

    They are different tables with different content hashes, so a run that used
    one cannot be reproduced from the other's stamp.
    """
    stamps = set()
    for reference_set in REFERENCE_SETS:
        result = api.optimize(
            "MAALKHETQW",
            api.OptimizeConfig(organism="escherichia_coli", reference_set=reference_set),
        )
        manifest = json.loads(api.result_to_json(result))["audit"]["manifest"]
        stamps.add(json.dumps(manifest, sort_keys=True))
    assert len(stamps) == len(REFERENCE_SETS)


def test_default_run_uses_the_highly_expressed_table() -> None:
    """Omitting the knob gives the organism's default, and says so."""
    explicit = api.optimize(
        "MAALKHETQW",
        api.OptimizeConfig(organism="homo_sapiens", reference_set=HIGHLY_EXPRESSED),
    )
    implicit = api.optimize("MAALKHETQW", api.OptimizeConfig(organism="homo_sapiens"))
    assert implicit.dna == explicit.dna
    assert implicit.audit["codon_reference_set"] == HIGHLY_EXPRESSED
