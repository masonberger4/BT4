"""Tests for the GENCODE -> BT4 splice-panel builder.

The builder's whole job is to get one thing right that is easy to get wrong silently:
which base is a splice site. A GTF is 1-based inclusive and always records
``start <= end`` regardless of strand, GENCODE emits minus-strand exons in *transcript*
order, and BT4 stores minus-strand windows already reverse-complemented -- three
conventions that have to compose correctly or every score in the gate is misaligned.

These tests pin that composition against a synthetic genome where the answer is known by
construction, plus the two traps that silently relabel true positives as negatives:

* a window contains **every** overlapping transcript's sites, not just its own;
* opposite-strand sites are skipped by default rather than quietly labelled negative.

They use a ``{chrom: sequence}`` mapping in place of a 3 GB FASTA, so the arithmetic is
checkable in CI without the download.
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
def builder() -> ModuleType:
    return _load("bt4_make_gencode_splice_panel", "make_gencode_splice_panel.py")


# A synthetic chromosome with two exons and one canonical intron, at known coordinates.
#   exon 1: 1-based [1, 30]      intron: [31, 70]      exon 2: [71, 100]
# The intron opens GT at 31-32 and closes AG at 69-70.
_INTRON = "GT" + "".join("ACTC"[i % 4] for i in range(36)) + "AG"
_CHROM = ("".join("ACTC"[i % 4] for i in range(30)) + _INTRON
          + "".join("ACTC"[i % 4] for i in range(30)))
_GENOME = {"chr1": _CHROM}
_EXONS = [(1, 30), (71, 100)]

# A real minus-strand gene carries its GT-AG on the MINUS strand, which reads CT...AC on
# the plus strand. Reverse-complementing the whole chromosome produces exactly that, at
# the mirrored coordinates -- so the same exon list is canonical here on `-`.
_CHROM_MINUS = "".join(
    {"A": "T", "C": "G", "G": "C", "T": "A"}[base] for base in reversed(_CHROM)
)


def test_the_fixture_really_carries_the_canonical_dinucleotides(builder: ModuleType) -> None:
    """Guard: if the synthetic genome is wrong, every assertion below is meaningless."""
    assert len(_CHROM) == 100
    assert _CHROM[30:32] == "GT"  # 1-based 31-32, the intron's first two bases
    assert _CHROM[68:70] == "AG"  # 1-based 69-70, its last two
    # The minus-strand fixture is the same intron on the other strand: CT...AC as read
    # on the plus strand, which is what a real minus-strand gene looks like there.
    assert len(_CHROM_MINUS) == 100
    assert _CHROM_MINUS[30:32] == "CT"
    assert _CHROM_MINUS[68:70] == "AC"


def test_plus_strand_sites_land_on_the_intronic_dinucleotides(builder: ModuleType) -> None:
    """Donor = the G of GT (first intronic base); acceptor = the G of AG (last)."""
    sites = builder.sites_from_exons(_EXONS, "+")
    assert sorted(sites) == [(31, "donor"), (70, "acceptor")]


def test_minus_strand_swaps_donor_and_acceptor(builder: ModuleType) -> None:
    """Transcription runs the other way, so the intron's first base is the high one."""
    sites = builder.sites_from_exons(_EXONS, "-")
    assert sorted(sites) == [(31, "acceptor"), (70, "donor")]


def test_exons_are_sorted_rather_than_trusted(builder: ModuleType) -> None:
    """GENCODE emits minus-strand exons in transcript order; a round-trip may not."""
    assert builder.sites_from_exons(reversed(_EXONS), "+") == builder.sites_from_exons(
        _EXONS, "+"
    )


def test_abutting_exons_produce_no_sites(builder: ModuleType) -> None:
    """No intron, no splice sites -- and no negative-width window either."""
    assert builder.sites_from_exons([(1, 30), (31, 60)], "+") == []


def test_one_site_pair_per_intron_not_per_exon(builder: ModuleType) -> None:
    """The spurious 'first acceptor' and 'last donor' are never generated at all."""
    three_exons = [(1, 30), (71, 100), (141, 170)]
    sites = builder.sites_from_exons(three_exons, "+")
    assert len(sites) == 4  # two introns
    assert sum(1 for _, kind in sites if kind == "donor") == 2


# --------------------------------------------------------------------------
# The composition: genomic coordinate -> index in the stored window


@pytest.mark.parametrize("strand", ["+", "-"])
def test_the_built_window_passes_bt4s_own_motif_check(
    builder: ModuleType, strand: str
) -> None:
    """The end-to-end property: BT4's reader accepts what this builder writes.

    This is the composition that has to be right -- GTF convention, strand handling,
    reverse-complementing, and index mapping -- and BT4's reader independently verifies
    the dinucleotide, so a passing panel is evidence rather than assertion.
    """
    from bt4.biomodels.splice.panel import SpliceWindow as BT4Window
    from bt4.biomodels.splice.panel import panel_from_windows

    reference = _CHROM if strand == "+" else _CHROM_MINUS
    transcripts = {
        "ENST1": builder.Transcript("ENST1", "GENE1", "chr1", strand, list(_EXONS)),
        "ENST2": builder.Transcript("ENST2", "GENE2", "chr3", strand, list(_EXONS)),
    }
    genome = {"chr1": reference, "chr3": reference}
    windows, counts = builder.build_windows(transcripts, genome, flank=0)
    assert len(windows) == 2, counts

    panel = panel_from_windows(
        [
            BT4Window(w.window_id, w.group, w.sequence, w.donors, w.acceptors, w.strand, w.note)
            for w in windows
        ],
        negative_construction="all other positions in the window",
    )
    assert panel.motif_consistency().fraction == 1.0
    assert panel.n_sites == 4


def test_minus_strand_indices_are_relative_to_the_stored_window(
    builder: ModuleType,
) -> None:
    """Stored reverse-complemented, so a genomic coordinate maps to ``w_end - g``."""
    assert builder.to_index(31, 1, 100, "+") == 30
    assert builder.to_index(31, 1, 100, "-") == 69
    assert builder.revcomp("ACGT") == "ACGT"
    assert builder.revcomp("GTAAGT") == "ACTTAC"


# --------------------------------------------------------------------------
# The two traps


def test_a_window_collects_every_overlapping_transcripts_sites(
    builder: ModuleType,
) -> None:
    """Trap 1: labelling only the centre transcript makes neighbours' sites negatives.

    A backend that correctly detects a neighbouring gene's real splice site would be
    punished for it -- a false positive that is actually a true one.
    """
    neighbour = [(101, 130), (171, 200)]
    chrom = _CHROM + _INTRON.join(("".join("ACTC"[i % 4] for i in range(30)),) * 2)
    genome = {"chr1": chrom}
    transcripts = {
        "ENST1": builder.Transcript("ENST1", "A", "chr1", "+", list(_EXONS)),
        "ENST2": builder.Transcript("ENST2", "B", "chr1", "+", list(neighbour)),
    }
    windows, _ = builder.build_windows(transcripts, genome, flank=200)
    # Both windows span both genes, so both must carry all four sites.
    for window in windows:
        assert len(window.donors) + len(window.acceptors) == 4, window.window_id


def test_antisense_overlapping_windows_are_skipped_by_default(
    builder: ModuleType,
) -> None:
    """Trap 2: the models are strand-specific, so an antisense site is not a site here.

    It is real sequence that looks exactly like one, though, so labelling it a negative
    is a claim. Skipping is the default because a silent false negative costs more than
    a smaller panel.
    """
    transcripts = {
        "ENST1": builder.Transcript("ENST1", "SENSE", "chr1", "+", list(_EXONS)),
        "ENST2": builder.Transcript("ENST2", "ANTI", "chr1", "-", list(_EXONS)),
    }
    windows, counts = builder.build_windows(transcripts, _GENOME, flank=0)
    assert windows == []
    assert counts["n_antisense"] == 2

    # Opting in keeps the window and records the disclosure in its note. (Only the sense
    # transcript survives here: this chromosome is canonical on `+`, so the antisense
    # one is not a real gene and the motif self-check correctly drops it -- which is
    # itself the self-check doing its job.)
    kept, counts = builder.build_windows(
        transcripts, _GENOME, flank=0, keep_antisense=True
    )
    assert [w.window_id for w in kept] == ["SENSE_ENST1"]
    assert "2 antisense site(s) present, scored as negatives" in kept[0].note


def test_a_window_with_an_assembly_gap_is_skipped(builder: ModuleType) -> None:
    """BT4's format forbids N; an unscoreable position is not a real negative."""
    gapped = {"chr1": _CHROM[:50] + "N" + _CHROM[51:]}
    transcripts = {"ENST1": builder.Transcript("ENST1", "G", "chr1", "+", list(_EXONS))}
    windows, counts = builder.build_windows(transcripts, gapped, flank=0)
    assert windows == []
    assert counts["n_gap"] == 1


# --------------------------------------------------------------------------
# Provenance and determinism


def test_only_mane_select_exons_are_read(builder: ModuleType, tmp_path: Path) -> None:
    """MANE Select is what keeps a panel stable across GENCODE releases.

    From v44 to v50 the protein-coding transcript count on the held-out chromosomes
    grows 4.1x while MANE Select grows 1.3%; an unfiltered panel's negative class fills
    with low-confidence transcript models.
    """
    gtf = tmp_path / "a.gtf"
    mane = 'tag "MANE_Select";'
    rows = [
        f'chr1\tH\texon\t1\t30\t.\t+\t.\tgene_name "KEEP"; transcript_id "T1"; {mane}',
        f'chr1\tH\texon\t71\t100\t.\t+\t.\tgene_name "KEEP"; transcript_id "T1"; {mane}',
        'chr1\tH\texon\t1\t30\t.\t+\t.\tgene_name "DROP"; transcript_id "T2"; tag "basic";',
        f'chr2\tH\texon\t1\t30\t.\t+\t.\tgene_name "OFF"; transcript_id "T3"; {mane}',
        "#comment",
    ]
    gtf.write_text("\n".join(rows) + "\n", encoding="utf-8")
    transcripts = builder.parse_gtf(gtf, ["chr1"])
    assert set(transcripts) == {"T1"}
    assert transcripts["T1"].gene_name == "KEEP"
    assert transcripts["T1"].exons == [(1, 30), (71, 100)]
    assert transcripts["T1"].span == (1, 100)


def test_the_builder_is_deterministic(builder: ModuleType) -> None:
    """Invariant #7: same inputs, byte-identical panel."""
    transcripts = {
        "ENST1": builder.Transcript("ENST1", "A", "chr1", "+", list(_EXONS)),
        "ENST2": builder.Transcript("ENST2", "B", "chr3", "+", list(_EXONS)),
    }
    genome = {"chr1": _CHROM, "chr3": _CHROM}
    first, _ = builder.build_windows(transcripts, genome, flank=0)
    second, _ = builder.build_windows(transcripts, genome, flank=0)
    assert first == second


def test_the_default_chromosomes_are_the_held_out_ones(builder: ModuleType) -> None:
    """Building from a training chromosome produces flattering nonsense."""
    assert builder.HELD_OUT_CHROMOSOMES == ("chr1", "chr3", "chr5", "chr7", "chr9")


def test_the_default_flank_gives_the_cnns_their_context(builder: ModuleType) -> None:
    """Both wrapped CNNs consume ~10 kb, so every interior site wants 5,000 nt a side."""
    assert builder.DEFAULT_FLANK == 5_000


# --------------------------------------------------------------------------
# Second review round


def test_one_non_canonical_intron_does_not_discard_a_small_gene(
    builder: ModuleType,
) -> None:
    """``len(bad) > len(same) * 0.1`` sounded like 10% slack and was not.

    A k-intron gene contributes only 2k sites, so the allowance was 0.2k -- below **1**
    for every gene with four introns or fewer. One legitimate GC-AG intron then discarded
    every site of the gene, the exact opposite of absorbing the minor spliceosome the
    tolerance exists for.
    """
    exon = "".join("ACTC"[i % 4] for i in range(30))
    canonical = "GT" + "".join("ACTC"[i % 4] for i in range(36)) + "AG"
    minor = "GC" + "".join("ACTC"[i % 4] for i in range(36)) + "AG"  # a real GC-AG intron
    chrom = exon + canonical + exon + minor + exon
    transcripts = {
        "T": builder.Transcript("T", "SMALL", "chr1", "+", [(1, 30), (71, 100), (141, 170)])
    }
    windows, counts = builder.build_windows(transcripts, {"chr1": chrom}, flank=0)
    assert counts["n_motif"] == 0
    assert len(windows) == 1
    assert len(windows[0].donors) + len(windows[0].acceptors) == 4


def test_a_truncated_flank_still_yields_correct_minus_strand_indices(
    builder: ModuleType,
) -> None:
    """Python slices truncate silently, and minus-strand indices are measured from w_end.

    So a fetch that ran off the contig used to shift every one of them. Real genes sit
    near contig ends, so this is reached on a real run. The window is shrunk to what the
    assembly actually has and the indices recomputed, rather than written wrong.
    """
    transcripts = {"T": builder.Transcript("T", "EDGE", "chr1", "-", list(_EXONS))}
    # flank=5000 asks for far more than this 100 nt contig can supply.
    windows, counts = builder.build_windows(
        transcripts, {"chr1": _CHROM_MINUS}, flank=5_000
    )
    assert counts["n_truncated"] == 0  # the span still fits, so it is usable
    assert len(windows) == 1
    window = windows[0]
    assert len(window.sequence) == len(_CHROM_MINUS)
    # The indices must still land on the canonical dinucleotides.
    for position in window.donors:
        assert window.sequence[position : position + 2] == "GT"
    for position in window.acceptors:
        assert window.sequence[position - 1 : position + 1] == "AG"


def test_a_window_whose_span_runs_off_the_contig_is_skipped(
    builder: ModuleType,
) -> None:
    """Truncation that cuts into the transcript itself cannot be salvaged."""
    transcripts = {"T": builder.Transcript("T", "EDGE", "chr1", "-", list(_EXONS))}
    windows, counts = builder.build_windows(
        transcripts, {"chr1": _CHROM_MINUS[:50]}, flank=0
    )
    assert windows == []
    assert counts["n_truncated"] == 1


def test_a_neighbours_site_with_no_room_for_its_motif_is_dropped(
    builder: ModuleType,
) -> None:
    """Otherwise the builder writes a panel BT4's own reader refuses.

    Sites are collected from every overlapping transcript, so one can land on the window
    boundary with no room for its own dinucleotide -- a legitimate site of a neighbouring
    gene that this particular window simply cannot carry.
    """
    from bt4.biomodels.splice.panel import SpliceWindow as BT4Window
    from bt4.biomodels.splice.panel import panel_from_windows

    # Two transcripts whose windows each clip the other's sites at an edge.
    neighbour = [(101, 130), (171, 200)]
    chrom = _CHROM + _INTRON.join(("".join("ACTC"[i % 4] for i in range(30)),) * 2)
    transcripts = {
        "A": builder.Transcript("A", "A", "chr1", "+", list(_EXONS)),
        "B": builder.Transcript("B", "B", "chr1", "+", list(neighbour)),
    }
    windows, _ = builder.build_windows(transcripts, {"chr1": chrom}, flank=0)
    # Whatever survives must be readable by BT4 -- that is the property under test.
    panel = panel_from_windows(
        [
            BT4Window(w.window_id, w.group, w.sequence, w.donors, w.acceptors, w.strand, w.note)
            for w in windows
        ],
        negative_construction="all other positions",
    )
    assert panel.motif_consistency().fraction >= 0.9
