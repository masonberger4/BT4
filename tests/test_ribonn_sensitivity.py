"""Tests for the zero-data RiboNN sensitivity checks (``scripts/ribonn_sensitivity.py``).

These run without torch, pandas or a RiboNN checkout: the checks are driven through
the ``null`` placeholder backend and through small fake predictors, so what is pinned
is the *reasoning* the script encodes -- grouping, spread arithmetic, tie handling,
confound correlation, and the honesty labelling -- not RiboNN's weights.

The load-bearing property: a backend that is **blind** to synonymous change must be
reported as blind, and a backend whose response is pure GC3 must be reported as a GC
detector. Those two verdicts are the whole reason the script exists.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import pathlib
import sys
from types import ModuleType

import pytest

from bt4.biomodels.expression import ExpressionResult

_SCRIPTS = pathlib.Path(__file__).resolve().parent.parent / "scripts"


def _load_script(name: str) -> ModuleType:
    """Import a maintainer script by path (they are not an installed package)."""
    path = _SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"bt4_script_{name}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


sens = _load_script("ribonn_sensitivity")


@dataclasses.dataclass(frozen=True)
class _FakePredictor:
    """A scriptable stand-in: score = f(dna), with a UTR context to swap."""

    scorer: object
    utr5: str = "GCCACC"
    utr3: str = "GCTAAT"
    label: str = "fake"

    @property
    def name(self) -> str:
        return self.label

    @property
    def calibrated(self) -> bool:
        return False

    def score_sequence(self, dna: str) -> ExpressionResult:
        return ExpressionResult(
            score=float(self.scorer(dna, self.utr5)),  # type: ignore[operator]
            model_name=self.name,
            calibrated=False,
            units="fake units",
        )


# --- grouping -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("KRas4B|Algorithm1|run3|acc=AF493917", "KRas4B"),
        ("Native", "Native"),
        ("  PDE3A |x", "PDE3A"),
    ],
)
def test_group_is_the_text_before_the_first_pipe(header: str, expected: str) -> None:
    # The grouping unit is the protein throughout this work: synonymous variants of
    # one protein are a dependent cluster, so they must land in one group.
    assert sens.group_of(header) == expected


def test_the_in_tree_panel_groups_into_three_proteins() -> None:
    from bt4.io.fasta import read_fasta

    records = read_fasta(sens.DEFAULT_PANEL)
    groups: dict[str, int] = {}
    for header, _dna in records:
        groups[sens.group_of(header)] = groups.get(sens.group_of(header), 0) + 1
    assert groups == {"KRas4B": 31, "Beclin1": 31, "PDE3A": 31}


# --- statistics ---------------------------------------------------------------


def test_sign_test_matches_hand_computed_values() -> None:
    assert sens.sign_test_p_value(0, 0) == 1.0
    assert sens.sign_test_p_value(5, 10) == pytest.approx(1.0)
    # 10/10 successes: two-sided exact p = 2 * (1/1024)
    assert sens.sign_test_p_value(10, 10) == pytest.approx(2.0 / 1024.0)
    assert sens.sign_test_p_value(0, 10) == pytest.approx(2.0 / 1024.0)
    # symmetric in successes vs failures
    assert sens.sign_test_p_value(8, 10) == pytest.approx(sens.sign_test_p_value(2, 10))


def test_sign_test_never_exceeds_one() -> None:
    # Doubling the smaller tail can overshoot for near-even splits; it must clamp.
    for trials in range(1, 12):
        for successes in range(trials + 1):
            assert 0.0 < sens.sign_test_p_value(successes, trials) <= 1.0


def test_gc3_reads_third_positions_only() -> None:
    assert sens.gc3_fraction("AAG AAC".replace(" ", "")) == 1.0  # G, C at wobble
    assert sens.gc3_fraction("AAAAAA") == 0.0
    assert sens.gc3_fraction("GGAGGA") == 0.0  # GC-rich but A at both third positions
    assert sens.gc3_fraction("") == 0.0


def test_read_flank_refuses_empty_and_non_acgt(tmp_path: pathlib.Path) -> None:
    # RiboNN refuses empty UTRs and the UTRs carry most of its signal, so an empty
    # or garbled flank must fail loudly here rather than deep inside the model.
    with pytest.raises(ValueError, match="empty"):
        sens.read_flank("   ")
    with pytest.raises(ValueError, match="non-ACGT"):
        sens.read_flank("ACGTX")
    fasta = tmp_path / "utr.fa"
    fasta.write_text(">u\nacg\ntaa\n", encoding="utf-8")
    assert sens.read_flank(str(fasta)) == "ACGTAA"  # path read, upper-cased, joined


# --- the two verdicts the script exists to produce ----------------------------


def _panel() -> list[tuple[str, str]]:
    """Two proteins x three synonymous variants, differing only in wobble bases."""
    return [
        ("P1|v1", "ATGAAAGGTTAA"),
        ("P1|v2", "ATGAAGGGCTAA"),
        ("P1|v3", "ATGAAGGGGTAA"),
        ("P2|v1", "ATGGATCCTTAA"),
        ("P2|v2", "ATGGACCCCTAA"),
        ("P2|v3", "ATGGACCCGTAA"),
    ]


def _tables() -> tuple[object, object]:
    from bt4.biomodels.codon.tables import load_table
    from bt4.biomodels.codon.tai import load_tai_table

    return load_table("homo_sapiens"), load_tai_table("homo_sapiens")


def test_a_blind_backend_is_reported_as_blind() -> None:
    # THE decisive verdict. A backend that returns the same number for every
    # synonymous variant must be reported as not responding -- that outcome closes
    # the calibration project honestly, so it must not be mistakable for a signal.
    table, tai = _tables()
    blind = _FakePredictor(scorer=lambda dna, utr5: 1.0)

    report = sens.check_cds_spread(blind, _panel(), table, tai)

    assert report["responds_to_synonymous_change"] is False
    assert report["median_within_group_sd"] == 0.0
    assert report["within_over_between"] is None  # zero denominator, not a ratio of 0


def test_a_gene_identity_backend_looks_responsive_between_but_not_within() -> None:
    # The failure mode the strict bar exists to catch: a backend that knows only
    # "which protein is this" has a large BETWEEN-group spread and zero WITHIN-group
    # spread. It must not be credited with synonymous-ranking skill.
    table, tai = _tables()
    by_length = _FakePredictor(scorer=lambda dna, utr5: 10.0 * len(dna))

    report = sens.check_cds_spread(by_length, _panel(), table, tai)

    assert report["responds_to_synonymous_change"] is False  # nothing within a protein
    assert report["between_group_sd"] == 0.0  # (both panel proteins share a length)
    for group in report["groups"]:
        assert group["score"]["sd"] == 0.0


def test_a_gc3_detector_is_exposed_by_the_confound_correlation() -> None:
    # A backend whose entire within-protein response IS GC3 does respond -- and the
    # report must say it is a GC detector, because BT4 already has GC for free.
    table, tai = _tables()
    gc3_only = _FakePredictor(scorer=lambda dna, utr5: sens.gc3_fraction(dna))

    report = sens.check_cds_spread(gc3_only, _panel(), table, tai)

    assert report["responds_to_synonymous_change"] is True
    assert report["median_abs_gc3_spearman"] == pytest.approx(1.0)


def test_within_over_between_ratio_is_computed_from_both_spreads() -> None:
    table, tai = _tables()
    # Score depends on BOTH the protein (a big offset) and the variant (a small term),
    # so the ratio is small but well defined -- the realistic shape.
    offsets = {"ATGA": 0.0, "ATGG": 100.0}

    def scorer(dna: str, utr5: str) -> float:
        return offsets[dna[:4]] + sens.gc3_fraction(dna)

    report = sens.check_cds_spread(_FakePredictor(scorer=scorer), _panel(), table, tai)
    ratio = report["within_over_between"]
    assert ratio is not None
    assert 0.0 < ratio < 0.1  # real within-protein movement, dwarfed by between


# --- utr-control --------------------------------------------------------------


def test_utr_control_detects_a_utr_blind_harness() -> None:
    # If swapping both UTRs leaves the score untouched, the harness is broken and
    # every null result elsewhere would be a wiring bug misread as biology.
    ignores_utr = _FakePredictor(scorer=lambda dna, utr5: float(len(dna)))
    report = sens.check_utr_control(
        ignores_utr, "ATGAAATAA", ("GCCACC", "GCTAAT"), ("AAACCC", "TTTGGG")
    )
    assert report["harness_ok"] is False


def test_utr_control_passes_when_the_utr_reaches_the_model() -> None:
    uses_utr = _FakePredictor(scorer=lambda dna, utr5: float(len(dna) + len(utr5) * 3))
    report = sens.check_utr_control(
        uses_utr, "ATGAAATAA", ("GCCACC", "GCTAAT"), ("AAACCCGGG", "TTTGGG")
    )
    assert report["harness_ok"] is True
    assert report["abs_difference"] == pytest.approx(9.0)  # (9 - 6) * 3


def test_utr_control_refuses_a_backend_with_no_utr_context() -> None:
    from bt4.biomodels.expression import NullExpressionModel

    with pytest.raises(ValueError, match="needs a UTR-aware backend"):
        sens.check_utr_control(
            NullExpressionModel(), "ATGAAATAA", ("A", "C"), ("G", "T")
        )


# --- direction ----------------------------------------------------------------


def test_direction_excludes_ties_instead_of_scoring_them_as_failures() -> None:
    # A blind backend ties every pair. Counting ties as failures would report
    # "0/N prefer the optimized design", which reads as a strong preference for the
    # DEOPTIMIZED one -- the opposite of the truth.
    table, _tai = _tables()
    blind = _FakePredictor(scorer=lambda dna, utr5: 1.0)
    proteins = ["MKV", "MGGA"]

    report = sens.check_direction(blind, proteins, table, "homo_sapiens", None, 0)

    assert report["pairs"] == 2
    assert report["ties"] == 2
    assert report["trials"] == 0
    assert report["successes"] == 0
    assert report["p_value"] == 1.0


def test_direction_counts_a_real_preference() -> None:
    table, _tai = _tables()
    # Score = CAI, so the max-CAI design wins every pair by construction.
    cai_backed = _FakePredictor(scorer=lambda dna, utr5: table.cai(dna))  # type: ignore[attr-defined]
    proteins = ["MKVGA", "MGGAKL", "MDEFKV"]

    report = sens.check_direction(cai_backed, proteins, table, "homo_sapiens", None, 0)

    assert report["ties"] == 0
    assert report["successes"] == report["trials"] == 3
    assert report["p_value"] == pytest.approx(0.25)  # 2 * (1/8), exact for n=3
    for row in report["proteins"]:
        assert row["cai_max_cai_design"] > row["cai_min_cai_design"]


# --- report shape / honesty ---------------------------------------------------


def test_report_always_carries_the_uncalibrated_labelling() -> None:
    report = sens.build_report(
        "cds-spread",
        _FakePredictor(scorer=lambda dna, utr5: sens.gc3_fraction(dna)),
        records=_panel(),
        utrs=("GCCACC", "GCTAAT"),
        alt_utrs=None,
        organism="homo_sapiens",
        reference_set=None,
        steps=3,
        seed=0,
        max_proteins=5,
    )
    assert report["backend"]["calibrated"] is False
    assert "no result here can promote a backend" in report["honesty"]
    # The reference set travels with any CAI number (CLAUDE.md §8).
    assert report["reference_set"] == "highly_expressed"
    # UTRs are identified by hash, never printed into the report.
    assert len(report["utr5_sha_prefix"]) == 12
    assert "GCCACC" not in report["utr5_sha_prefix"]


def test_unknown_check_and_missing_inputs_raise() -> None:
    kwargs = dict(
        records=_panel(),
        utrs=("GCCACC", "GCTAAT"),
        alt_utrs=None,
        organism="homo_sapiens",
        reference_set=None,
        steps=3,
        seed=0,
        max_proteins=5,
    )
    fake = _FakePredictor(scorer=lambda dna, utr5: 1.0)
    with pytest.raises(ValueError, match="unknown check"):
        sens.build_report("nope", fake, **kwargs)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="needs --utr5-alt"):
        sens.build_report("utr-control", fake, **kwargs)  # type: ignore[arg-type]


def test_cli_runs_the_null_backend_end_to_end(capsys: pytest.CaptureFixture[str]) -> None:
    # An end-to-end rehearsal that needs no weights: the null placeholder is the
    # reference for what "blind" looks like, and proves the plumbing before the
    # licensed model is ever pointed at it.
    code = sens.main(
        [
            "--check", "cds-spread", "--backend", "null",
            "--utr5", "GCCACC", "--utr3", "GCTAAT", "--json",
        ]
    )
    assert code == 0
    import json

    report = json.loads(capsys.readouterr().out)
    assert report["n_groups"] == 3  # the in-tree panel's three proteins
    assert report["n_records"] == 93
    assert report["responds_to_synonymous_change"] is False
    assert report["backend"]["calibrated"] is False


def test_cli_reports_a_missing_backend_without_a_traceback(
    capsys: pytest.CaptureFixture[str],
) -> None:
    code = sens.main(
        [
            "--check", "utr-control", "--backend", "null",
            "--utr5", "GCCACC", "--utr3", "GCTAAT",
            "--utr5-alt", "AAACCC", "--utr3-alt", "TTTGGG",
        ]
    )
    assert code == 2
    assert "needs a UTR-aware backend" in capsys.readouterr().err
