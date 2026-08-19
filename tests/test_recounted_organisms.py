"""Tests for the **genome-wide** codon tables recounted from pinned Ensembl CDS.

Every number in these nine tables is a codon occurrence count from a public,
release-pinned CDS FASTA, produced by ``scripts/build_organism_tables.py`` --
none is a hand-typed "representative published summary", which is what BT4 used
to ship for human, *E. coli* and yeast. So the tests here assert two things a
summary table cannot claim:

* **A complete re-derivation trail** -- source URL, the source file's own
  SHA-256, assembly, database release, and the filter tally -- so a third party
  can rebuild the exact shipped bytes (CLAUDE.md §8).
* **External ground truth** (§8, "not just self-consistency"): the tables must
  reproduce well-established, independently-published facts about each species'
  codon bias. A fabricated or mis-parsed table would fail these; internal
  consistency alone would not catch it.

They are also why the bundled GtRNAdb tRNA tables became reachable at all: tAI
needs an organism you can actually select, which needs a codon table.

**Every load here names ``GENOME_WIDE`` explicitly.** It is no longer the
default -- eight of the nine organisms default to their highly-expressed
reference set (see ``test_highly_expressed_tables.py``) -- and a test that
asserted "the top *E. coli* Phe codon is TTT" against whichever table happened to
be the default would silently start testing a different table.
"""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from bt4 import api
from bt4.biomodels.codon.tables import (
    GENOME_WIDE,
    CodonUsageTable,
    TableProvenance,
    available_organisms,
    default_reference_set,
    load_provenance,
    load_table,
    sha256_hex,
)
from bt4.biomodels.codon.tai import available_tai_organisms
from bt4.domain.genetic_code import CODON_TABLE

RECOUNTED: tuple[str, ...] = (
    "arabidopsis_thaliana",
    "bacillus_subtilis",
    "caenorhabditis_elegans",
    "cricetulus_griseus_chok1gshd",
    "danio_rerio",
    "drosophila_melanogaster",
    "escherichia_coli",
    "homo_sapiens",
    "komagataella_phaffii",
    "mus_musculus",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
)

_STOPS = ("TAA", "TAG", "TGA")


def _gw(name: str) -> CodonUsageTable:
    """The organism's genome-wide table -- the only one this file is about."""
    return load_table(name, reference_set=GENOME_WIDE)


def _gw_provenance(name: str) -> TableProvenance:
    """The genome-wide table's provenance sidecar."""
    return load_provenance(name, reference_set=GENOME_WIDE)


def _raw_provenance(name: str) -> dict[str, object]:
    """The provenance sidecar as raw JSON (``TableProvenance`` drops extra keys)."""
    text = files("bt4.biomodels.codon.data").joinpath(f"{name}.provenance.json")
    return dict(json.loads(text.read_text(encoding="utf-8")))


def _top_codon(table: CodonUsageTable, amino_acid: str) -> str:
    """The highest-frequency codon among ``amino_acid``'s synonyms."""
    synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
    return max(synonyms, key=table.weight)


def _gc3(table: CodonUsageTable) -> float:
    """Fraction of sense codons ending in G or C -- the classic GC3 readout."""
    sense = {c: v for c, v in table.frequency.items() if c not in _STOPS}
    return sum(v for c, v in sense.items() if c[2] in "GC") / sum(sense.values())


# --------------------------------------------------------------------------- #
# The tables load, are advertised, and are complete.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", RECOUNTED)
def test_table_loads_and_is_advertised(name: str) -> None:
    table = _gw(name)
    assert table.organism == name
    assert name in available_organisms()
    # A genome-scale CDS set observes every codon; the builder refuses to write
    # a table that does not, rather than smoothing an invented value into it.
    assert set(table.frequency) == set(CODON_TABLE)


@pytest.mark.parametrize("name", RECOUNTED)
def test_counts_are_positive_integers(name: str) -> None:
    table = _gw(name)
    for codon, value in table.frequency.items():
        assert value == int(value), f"{name} {codon} is not a whole count"
        assert value >= 1


@pytest.mark.parametrize("name", RECOUNTED)
def test_group_max_codon_has_unit_weight(name: str) -> None:
    # Leucine is a six-box amino acid; its most-used synonym must normalize to 1.
    assert _gw(name).weight(_top_codon(_gw(name), "L")) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# Provenance: a complete, checkable re-derivation trail.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_sha256_matches_shipped_tsv(name: str) -> None:
    raw = files("bt4.biomodels.codon.data").joinpath(f"{name}.tsv").read_bytes()
    assert _gw_provenance(name).sha256 == sha256_hex(raw)


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_is_a_real_recount_not_a_summary(name: str) -> None:
    prov = _gw_provenance(name)
    # The load-bearing distinction from the older bundled tables: a real CDS
    # count stands behind these numbers, and the note says so plainly.
    assert prov.cds_count is not None
    # A genome-wide set, not a sample. The bound accommodates the smallest
    # genome BT4 ships -- E. coli has only ~3,800 counted genes in total, which
    # IS its whole genome, so a mammal-sized floor would wrongly reject it.
    assert prov.cds_count > 3_000, "a genome-wide CDS set, not a sample"
    assert "Real genome-wide codon counts" in prov.note
    assert "REPRESENTATIVE" not in prov.note


@pytest.mark.parametrize("name", RECOUNTED)
def test_provenance_carries_the_rederivation_trail(name: str) -> None:
    prov = _raw_provenance(name)
    # Content-hashing our own output is not enough: without the SOURCE's hash and
    # URL nobody else can reproduce the table, which is what §8 asks for.
    for key in (
        "source_url",
        "source_sha256",
        "assembly",
        "genebuild",
        "database",
        "database_release",
        "total_codons_counted",
        "filters",
        "rebuild_command",
    ):
        assert key in prov, f"{name} provenance is missing {key!r}"
    assert str(prov["source_url"]).startswith("https://")
    assert len(str(prov["source_sha256"])) == 64
    assert f"release-{prov['database_release']}" in str(prov["source_url"])
    # Assembly and gene annotation are different things, and it is the ANNOTATION
    # that defines which CDS were counted -- so the sidecar must name both, and
    # the human-readable source line must not conflate them.
    assert str(prov["genebuild"]).strip()
    assert f"assembly {prov['assembly']}" in str(prov["source"])
    assert f"annotation {prov['genebuild']}" in str(prov["source"])


@pytest.mark.parametrize("name", RECOUNTED)
def test_stamped_totals_match_the_shipped_counts(name: str) -> None:
    prov = _raw_provenance(name)
    table = _gw(name)
    assert sum(table.frequency.values()) == prov["total_codons_counted"]
    filters = prov["filters"]
    assert isinstance(filters, dict)
    assert filters["cds_counted"] == _gw_provenance(name).cds_count
    assert filters["records_in_source"] >= filters["cds_counted"]


@pytest.mark.parametrize("name", RECOUNTED)
def test_filter_tally_accounts_for_every_source_record(name: str) -> None:
    """Kept + dropped must equal the source record count, exactly.

    This is the one substantive provenance claim checkable **offline**, with no
    network and no re-download: if the tally does not close, the sidecar is
    describing a different run than the one that produced the shipped TSV, and
    the "nothing is skipped silently" promise is not being kept.
    """
    filters = _raw_provenance(name)["filters"]
    assert isinstance(filters, dict)
    dropped = sum(v for k, v in filters.items() if k.startswith("dropped_"))
    assert dropped + filters["cds_counted"] == filters["records_in_source"]


# --------------------------------------------------------------------------- #
# External ground truth: independently-published facts about each species.
# --------------------------------------------------------------------------- #


def test_gc3_ordering_matches_known_genome_composition() -> None:
    """GC3 must order the species the way the literature does.

    Drosophila is strongly GC3-biased; mammals sit near 0.55-0.60; *C. elegans*
    and *Arabidopsis* are AT-rich coding genomes. Getting this ordering right is
    a real constraint on the data -- a mis-parsed or invented table would not.
    """
    gc3 = {name: _gc3(_gw(name)) for name in RECOUNTED}
    gc3["homo_sapiens"] = _gc3(_gw("homo_sapiens"))

    assert gc3["drosophila_melanogaster"] > gc3["danio_rerio"]
    assert gc3["danio_rerio"] > gc3["arabidopsis_thaliana"]
    assert gc3["arabidopsis_thaliana"] > gc3["caenorhabditis_elegans"]
    # AT-rich coding genomes sit well below half.
    for name in ("caenorhabditis_elegans", "arabidopsis_thaliana"):
        assert gc3[name] < 0.45
    # ...and the GC3-rich fly well above.
    assert gc3["drosophila_melanogaster"] > 0.60


def test_mammals_agree_with_each_other() -> None:
    """Mouse and rat were counted independently, yet must land near human.

    Three tables built from three separate CDS sets converging on the mammalian
    GC3 range is evidence the pipeline measured something real, not an artifact
    of one download.
    """
    human = _gc3(_gw("homo_sapiens"))
    for name in ("mus_musculus", "rattus_norvegicus"):
        assert abs(_gc3(_gw(name)) - human) < 0.05


@pytest.mark.parametrize(
    ("name", "amino_acid", "expected"),
    [
        # CTG is the preferred Leu codon across GC3-rich genomes...
        ("mus_musculus", "L", "CTG"),
        ("rattus_norvegicus", "L", "CTG"),
        ("danio_rerio", "L", "CTG"),
        ("drosophila_melanogaster", "L", "CTG"),
        # ...while AT-rich coding genomes prefer CTT.
        ("caenorhabditis_elegans", "L", "CTT"),
        ("arabidopsis_thaliana", "L", "CTT"),
    ],
)
def test_preferred_leucine_codon(name: str, amino_acid: str, expected: str) -> None:
    assert _top_codon(_gw(name), amino_acid) == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        # TGA dominates in GC-richer vertebrate/plant genomes, TAA in AT-rich ones.
        ("mus_musculus", "TGA"),
        ("rattus_norvegicus", "TGA"),
        ("danio_rerio", "TGA"),
        ("arabidopsis_thaliana", "TGA"),
        ("caenorhabditis_elegans", "TAA"),
        ("drosophila_melanogaster", "TAA"),
    ],
)
def test_preferred_stop_codon(name: str, expected: str) -> None:
    table = _gw(name)
    assert max(_STOPS, key=lambda c: table.frequency[c]) == expected


# --------------------------------------------------------------------------- #
# The point of the exercise: tAI is reachable for every organism that has it.
# --------------------------------------------------------------------------- #


def test_every_trna_table_has_a_selectable_organism() -> None:
    """No bundled tRNA table may be stranded without a codon table.

    tAI is only offered for an organism the user can actually select, and
    selection requires a codon-usage table. Before these six tables shipped, six
    of the eight bundled GtRNAdb tables were unreachable for exactly this reason
    -- this test is what keeps that from silently recurring.
    """
    stranded = sorted(set(available_tai_organisms()) - set(available_organisms()))
    assert stranded == [], f"tRNA data with no selectable organism: {stranded}"


@pytest.mark.parametrize("name", RECOUNTED)
def test_cai_in_unit_interval(name: str) -> None:
    table = _gw(name)
    dna = "ATG" + "GCA" + "CTT" + "CGG" + "AAA" + "TAA"
    assert 0.0 < table.cai(dna) <= 1.0


@pytest.mark.parametrize("name", RECOUNTED)
def test_table_content_hash_reaches_the_run_manifest(name: str) -> None:
    """Each new table's *content* hash must enter the provenance stamp (#9).

    Hashing config field names -- BT4's cautionary tale from BT3 (§10.10) -- would
    make two organisms stamp identically. Six new tables are six new chances to
    regress that, so assert the actual TSV digest appears in the manifest.
    """
    result = api.optimize(
        "MAALKHETQW", api.OptimizeConfig(organism=name, max_homopolymer=5)
    )
    manifest = json.loads(api.result_to_json(result))["audit"]["manifest"]
    # Deliberately the DEFAULT reference set, not the genome-wide one: the stamp
    # must track the table the run actually read, and which table that is now
    # varies by organism.
    assert manifest["inputs"]["codon_table_sha256"] == load_provenance(name).sha256
    assert result.audit["codon_reference_set"] == default_reference_set(name)


def test_each_organism_stamps_a_distinct_manifest() -> None:
    """Swapping the organism must change the stamp (invariant #9)."""
    stamps = set()
    for name in RECOUNTED:
        result = api.optimize(
            "MAALKHETQW", api.OptimizeConfig(organism=name, max_homopolymer=5)
        )
        manifest = json.loads(api.result_to_json(result))["audit"]["manifest"]
        stamps.add(json.dumps(manifest, sort_keys=True))
    assert len(stamps) == len(RECOUNTED)


def test_expected_strong_biases_in_the_microbial_tables() -> None:
    """Textbook biases must survive the switch from published values to recounts.

    *E. coli* strongly prefers CTG for Leu and *S. cerevisiae* prefers AGA for
    Arg. These held in the older hand-curated tables and must still hold now that
    the numbers are counted from the genome -- old and new agreeing here is what
    shows the recount reproduced the biology rather than merely replacing it.
    """
    assert _gw("escherichia_coli").weight("CTG") == pytest.approx(1.0)
    assert _gw("saccharomyces_cerevisiae").weight("AGA") == pytest.approx(1.0)


def test_every_bundled_organism_is_recounted() -> None:
    """No organism may quietly fall back to an undocumented table.

    BT4 used to ship three hand-typed "representative" tables, including the one
    for its DEFAULT organism -- so the most-used numbers were the least checkable.
    Every bundled organism is now a counted, re-derivable table, and this test is
    what keeps a future addition from reintroducing the old asymmetry.
    """
    assert set(available_organisms()) == set(RECOUNTED)


@pytest.mark.parametrize("name", RECOUNTED)
def test_no_alt_or_patch_region_survived_filtering(name: str) -> None:
    """No region that looks alternate/patch-like may be counted into a table.

    Ensembl publishes alternate haplotypes and patch scaffolds with their own
    gene IDs, so per-gene de-duplication does not collapse them and they inflate
    a table with duplicate copies of real genes. A name blocklist can only
    exclude conventions someone already knew about, and twice it did not: human's
    ``HG*_NOVEL_TEST`` patches leaked 12 genes (9 of them second copies of chr11
    olfactory receptors), and zebrafish's ``ALT_CTG*`` contigs leaked 4,127 genes
    -- **15.6% of that species' table** -- past a filter that looked complete.

    The builder therefore also counts anything it *kept* whose region name still
    looks alternate/patch-like, and stamps it. This test requires that count to
    be zero, so the next unknown naming variant fails here instead of quietly
    inflating a shipped table. It is checkable offline, with no re-download.
    """
    filters = _raw_provenance(name)["filters"]
    assert isinstance(filters, dict)
    assert filters["kept_suspicious_region"] == 0, (
        f"{name}: a region that looks like an alternate/patch locus was counted; "
        "the filter's naming list has fallen behind the source"
    )


# --------------------------------------------------------------------------- #
# The three industrial hosts, checked against what is known about them.
# --------------------------------------------------------------------------- #


def test_cho_lands_in_the_rodent_band() -> None:
    """CHO is a Chinese hamster line, so it must sit with mouse and rat.

    The strongest available external check for this table: it was counted from a
    separate download of a separate assembly, and if the pipeline had mis-parsed
    it there is no reason the answer would land within a point of two independently
    counted rodents. Deliberately a *band* check rather than a pinned value --
    the claim is "this is a rodent coding genome", which is what the data can
    support, not a target number.
    """
    cho = _gc3(_gw("cricetulus_griseus_chok1gshd"))
    for rodent in ("mus_musculus", "rattus_norvegicus"):
        assert abs(cho - _gc3(_gw(rodent))) < 0.02
    assert abs(cho - _gc3(_gw("homo_sapiens"))) < 0.03


def test_the_two_at_rich_industrial_hosts_are_at_rich() -> None:
    """*B. subtilis* and *K. phaffii* are low-GC genomes, unlike *E. coli*.

    A real discriminator: *E. coli* and *B. subtilis* are both bacteria counted
    through the identical Ensembl Bacteria path, so a pipeline artifact would move
    them together. Their genome GC differs by ~7 points (50.8% vs 43.5%) and the
    tables must reproduce that separation.
    """
    coli = _gc3(_gw("escherichia_coli"))
    subtilis = _gc3(_gw("bacillus_subtilis"))
    phaffii = _gc3(_gw("komagataella_phaffii"))

    assert subtilis < coli - 0.05
    assert subtilis < 0.50
    # The yeast is AT-richer still, in the band its relatives occupy.
    assert phaffii < 0.45
    assert phaffii > _gc3(_gw("saccharomyces_cerevisiae"))


@pytest.mark.parametrize(
    ("name", "amino_acid", "expected"),
    [
        # AT-rich bacteria and yeasts take the A/T-ending synonym for Lys and Glu;
        # this is the textbook signature of both hosts and is not a free parameter.
        ("bacillus_subtilis", "K", "AAA"),
        ("bacillus_subtilis", "E", "GAA"),
        ("komagataella_phaffii", "K", "AAA"),
        ("komagataella_phaffii", "E", "GAA"),
        # ...and the yeast prefers TTG for Leu, as S. cerevisiae does.
        ("komagataella_phaffii", "L", "TTG"),
        # CHO is mammalian and must take the G-ending synonyms instead.
        ("cricetulus_griseus_chok1gshd", "K", "AAG"),
        ("cricetulus_griseus_chok1gshd", "E", "GAG"),
        ("cricetulus_griseus_chok1gshd", "L", "CTG"),
    ],
)
def test_industrial_host_preferred_codons(name: str, amino_acid: str, expected: str) -> None:
    """The most-used synonym must be the one the literature reports."""
    table = _gw(name)
    synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
    assert max(synonyms, key=table.weight) == expected


def test_bacillus_start_codon_filter_is_a_gap_not_a_wrong_answer() -> None:
    """*B. subtilis* loses 22.5% of CDS to the ATG-start filter. Measured, not waved.

    The shared validity filter requires an ``ATG`` start, which drops 954 of 4,237
    *B. subtilis* records -- TTG (553), GTG (387), ATT (8), CTG (5), ATC (1). That
    is more than double the 9.6% it costs *E. coli*, because *B. subtilis* genuinely
    uses alternative starts, so the gap had to be measured for this organism rather
    than inherited from the *E. coli* finding.

    Counting the dropped genes back in (skipping the initiator, which is not a codon
    *choice* -- the ribosome uses fMet-tRNA whatever the triplet) moves **no** amino
    acid's most-used codon, and shifts relative adaptiveness by at most 0.023. So the
    filter costs precision, not correctness. This test pins the consequence that
    matters: the shipped table's preferred codons are the ones a complete count gives.
    """
    table = _gw("bacillus_subtilis")
    # The dropped genes are AT-rich alternative-start genes; if their absence had
    # skewed the table, the AT-ending preferences below are what would have moved.
    for amino_acid, expected in (("K", "AAA"), ("E", "GAA"), ("F", "TTT")):
        synonyms = [c for c, aa in CODON_TABLE.items() if aa == amino_acid]
        assert max(synonyms, key=table.weight) == expected
    # And the drop itself is pinned from the sidecar, so a future filter change
    # that silently alters what was counted fails here rather than passing quietly.
    # Read as raw JSON: `TableProvenance` keeps only its own fields, and an
    # attribute check against it would make this assertion dead code.
    sidecar = json.loads(
        files("bt4.biomodels.codon.data")
        .joinpath("bacillus_subtilis.provenance.json")
        .read_text(encoding="utf-8")
    )
    filters = sidecar["filters"]
    assert filters["dropped_no_atg_start"] == 954
    assert filters["records_in_source"] == 4237
    assert filters["cds_counted"] == 3283
