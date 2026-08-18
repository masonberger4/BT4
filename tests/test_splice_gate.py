"""Tests for the splice statistical-calibration gate.

Pangolin has passed its *integration-fidelity* gate -- BT4's adapter reproduces
upstream bit-for-bit. This gate asks the separate, still-unmet question: do those
numbers **mean** what they claim, and is the operating point BT4 thresholds at the
right one?

The tests pin the properties that make the gate hard to pass dishonestly:

* the two tasks are **kept apart** -- site prediction and variant effect answer
  different questions against different ground truth, and pooling them is refused;
* the verdict is **per stratum**, so a backend cannot certify on intronic strength
  while failing the exonic half where BT4 actually operates;
* a **single-class stratum is refused**, in both directions -- with no positives
  PR-AUC returns 0.0 (reading as incompetence), with no negatives 1.0 (reading as
  perfection), and both are artifacts of the panel rather than findings;
* the panel's **negative construction is mandatory**, because PR-AUC moves with
  prevalence and a threshold without a pinned denominator is passable by thinning
  negatives;
* the defaults **produce a report, not a pass**, so calling the gate with no
  thresholds can never certify anything;
* a **vacuous predictor** -- base rate everywhere, perfectly calibrated -- is caught
  by the skill scores even though its ECE looks excellent.
"""

from __future__ import annotations

import random

import pytest

from bt4.biomodels.splice.base import DEFAULT_SITE_PROBABILITY
from bt4.biomodels.splice.gate import (
    SITE_PREDICTION,
    VARIANT_EFFECT,
    SpliceSiteCase,
    SpliceVariantCase,
    verify_splice_gate,
)

_NEG = "all positions in test-chromosome gene bodies"


def _site_panel(
    *, n_per_chrom: int = 1000, prevalence: int = 100, quality: float = 0.08
) -> list[SpliceSiteCase]:
    """A synthetic site panel over two chromosomes with a separable-ish signal."""
    rng = random.Random(0)
    cases: list[SpliceSiteCase] = []
    for chrom in ("chr1", "chr3"):
        for i in range(n_per_chrom):
            is_site = i % prevalence == 0
            base = 0.85 if is_site else 0.05
            score = min(0.99, max(0.001, base + rng.gauss(0, quality)))
            # Assign the stratum from an index that is coprime with the prevalence
            # period, so positives land in BOTH kinds. (Deriving it from `i % 2`
            # would put every positive in one stratum, since the period is even.)
            kind = "donor" if (i // prevalence) % 2 == 0 else "acceptor"
            cases.append(
                SpliceSiteCase(
                    predicted=score, label=int(is_site), kind=kind, group=chrom
                )
            )
    return cases


def _variant_panel() -> list[SpliceVariantCase]:
    """A synthetic variant panel that is strong intronically and weak exonically."""
    rng = random.Random(1)
    cases: list[SpliceVariantCase] = []
    for gene in ("BRCA1", "FAS", "WT1"):
        for i in range(200):
            disruptive = i % 10 == 0
            # Intronic: well separated. Exonic: barely separated -- the real pattern.
            # Stratify on `i // 10` so disruptive variants land in BOTH regions;
            # `i % 2` would put every one of them in the same stratum.
            intronic = (i // 10) % 2 == 0
            spread = 0.10 if intronic else 0.30
            base = 0.80 if disruptive else 0.20
            score = min(0.99, max(0.001, base + rng.gauss(0, spread)))
            cases.append(
                SpliceVariantCase(
                    predicted=score,
                    label=int(disruptive),
                    region="intronic" if intronic else "exonic",
                    group=gene,
                )
            )
    return cases


# --------------------------------------------------------------------------
# Task separation


def test_site_panel_is_recognised_as_site_prediction() -> None:
    """The task is inferred from the case type, not passed in."""
    report = verify_splice_gate(_site_panel(), negative_construction=_NEG)
    assert report.task == SITE_PREDICTION
    assert {s.name for s in report.strata} == {"donor", "acceptor"}


def test_variant_panel_is_recognised_as_variant_effect() -> None:
    """Variant cases stratify by region, the split that matters for BT4."""
    report = verify_splice_gate(_variant_panel(), negative_construction=_NEG)
    assert report.task == VARIANT_EFFECT
    assert {s.name for s in report.strata} == {"exonic", "intronic"}


def test_mixing_the_two_case_types_is_refused() -> None:
    """They answer different questions against different ground truth."""
    mixed = [
        SpliceSiteCase(0.9, 1, "donor", "chr1"),
        SpliceVariantCase(0.9, 1, "exonic", "BRCA1"),
    ]
    with pytest.raises(ValueError, match="different tasks"):
        verify_splice_gate(mixed, negative_construction=_NEG)


# --------------------------------------------------------------------------
# The verdict is per stratum


def test_a_weak_stratum_fails_the_whole_gate() -> None:
    """The point of stratifying: intronic strength cannot carry an exonic failure.

    The synthetic panel is deliberately well separated intronically and barely
    separated exonically -- the pattern Smith & Kitzman measured (median prAUC 0.773
    vs 0.419). BT4 designs coding sequence, so the exonic stratum is its regime.
    """
    report = verify_splice_gate(
        _variant_panel(), negative_construction=_NEG, min_pr_auc=0.60
    )
    by_name = {s.name: s for s in report.strata}
    assert by_name["intronic"].pr_auc > by_name["exonic"].pr_auc
    assert report.passed is False
    assert any("exonic" in r for r in report.reasons)


def test_pooled_metrics_are_reported_but_never_the_verdict() -> None:
    """``overall`` exists for comparability and carries no pass authority."""
    report = verify_splice_gate(
        _variant_panel(), negative_construction=_NEG, min_pr_auc=0.60
    )
    # The pooled figure clears the bar the exonic stratum fails...
    assert report.overall.pr_auc > 0.60
    # ...and the gate still fails.
    assert report.passed is False


def test_passing_requires_every_stratum_to_clear() -> None:
    """A permissive bar every stratum clears does pass."""
    report = verify_splice_gate(
        _variant_panel(), negative_construction=_NEG, min_pr_auc=0.20
    )
    assert report.passed is True
    assert report.reasons == ()


# --------------------------------------------------------------------------
# Unscoreable strata are refused, not scored


def test_a_stratum_with_no_positives_is_refused() -> None:
    """PR-AUC would return 0.0, which reads as incompetence rather than no data."""
    cases = _site_panel()
    cases = [c for c in cases if not (c.kind == "donor" and c.label == 1)]
    cases.append(SpliceSiteCase(0.1, 0, "donor", "chr1"))
    report = verify_splice_gate(cases, negative_construction=_NEG)
    assert report.passed is False
    assert any("no positives" in r for r in report.reasons)
    assert "donor" not in {s.name for s in report.strata}


def test_a_stratum_with_no_negatives_is_refused() -> None:
    """The mirror case: PR-AUC would return 1.0, reading as perfection."""
    cases = [
        SpliceSiteCase(0.9, 1, "donor", "chr1"),
        SpliceSiteCase(0.8, 1, "donor", "chr3"),
        SpliceSiteCase(0.9, 1, "acceptor", "chr1"),
        SpliceSiteCase(0.1, 0, "acceptor", "chr3"),
    ]
    report = verify_splice_gate(cases, negative_construction=_NEG)
    assert any("no negatives" in r for r in report.reasons)
    assert "donor" not in {s.name for s in report.strata}


# --------------------------------------------------------------------------
# Provenance the numbers depend on


def test_negative_construction_is_mandatory() -> None:
    """Without a pinned denominator a PR-AUC threshold means nothing.

    PR-AUC's floor is the prevalence, and prevalence is a construction choice, so a
    gate that does not record how negatives were built is passable by sampling fewer
    of them.
    """
    with pytest.raises(ValueError, match="negative_construction is required"):
        verify_splice_gate(_site_panel(), negative_construction="   ")


def test_the_report_carries_its_panel_provenance() -> None:
    """Recorded verbatim, so a later reader can adjudicate the numbers."""
    report = verify_splice_gate(
        _site_panel(),
        negative_construction=_NEG,
        panel_note="GENCODE v44, chr1/chr3 held out",
    )
    assert report.negative_construction == _NEG
    assert "GENCODE v44" in report.panel_note


def test_prevalence_is_reported_beside_every_pr_auc() -> None:
    """A PR-AUC without its floor is not interpretable (section 8)."""
    report = verify_splice_gate(_site_panel(), negative_construction=_NEG)
    for stratum in (*report.strata, report.overall):
        assert 0.0 < stratum.prevalence < 1.0
        assert stratum.n_positive == pytest.approx(stratum.prevalence * stratum.n_cases)


def test_grouped_panel_needs_at_least_two_groups() -> None:
    """One chromosome cannot support a leakage-free evaluation."""
    single = [c for c in _site_panel() if c.group == "chr1"]
    with pytest.raises(ValueError, match="two leakage-control groups"):
        verify_splice_gate(single, negative_construction=_NEG)


# --------------------------------------------------------------------------
# The gate cannot certify by default, and catches vacuous models


def test_defaults_produce_a_report_not_a_pass() -> None:
    """Calling with no thresholds must never be how something gets certified.

    The defaults are permissive on purpose, so a bare call yields numbers; the
    thresholds are the maintainer's to set at gate time, on the evidence.
    """
    report = verify_splice_gate(_site_panel(), negative_construction=_NEG)
    assert report.min_pr_auc == 0.0
    assert report.max_ece == 1.0
    # It "passes" only in the sense that no bar was set -- which the report says.
    assert report.passed is True
    assert report.min_pr_auc_skill == 0.0


def test_a_vacuous_predictor_is_caught_by_the_skill_scores() -> None:
    """Base rate everywhere: perfectly calibrated, completely useless.

    Its ECE is excellent and its Brier tiny. Only the skill scores expose it, which
    is why the gate reports them alongside rather than trusting calibration alone.
    """
    cases = _site_panel()
    prevalence = sum(c.label for c in cases) / len(cases)
    vacuous = [
        SpliceSiteCase(prevalence, c.label, c.kind, c.group) for c in cases
    ]
    report = verify_splice_gate(
        vacuous, negative_construction=_NEG, min_pr_auc_skill=0.10
    )
    assert report.overall.ece < 0.02  # looks well calibrated
    assert report.overall.brier_skill == pytest.approx(0.0, abs=1e-9)
    assert report.overall.pr_auc_skill == pytest.approx(0.0, abs=1e-9)
    assert report.passed is False


def test_threshold_defaults_to_the_shared_operating_point() -> None:
    """The MCC is computed at BT4's own operating point, not an ad-hoc one."""
    report = verify_splice_gate(_site_panel(), negative_construction=_NEG)
    assert report.threshold == DEFAULT_SITE_PROBABILITY


@pytest.mark.parametrize("threshold", [0.0, 1.0, -0.5, 2.0])
def test_threshold_must_be_a_probability(threshold: float) -> None:
    """An operating point outside (0, 1) is not a cutoff."""
    with pytest.raises(ValueError, match="threshold must be in"):
        verify_splice_gate(
            _site_panel(), negative_construction=_NEG, threshold=threshold
        )


def test_empty_panel_is_refused() -> None:
    """An empty evaluation is a mistake worth surfacing, not a vacuous pass."""
    with pytest.raises(ValueError, match="at least one case"):
        verify_splice_gate([], negative_construction=_NEG)


def test_reliability_bins_are_reported_per_stratum() -> None:
    """The raw material of a reliability diagram, per stratum."""
    report = verify_splice_gate(_site_panel(), negative_construction=_NEG)
    for stratum in report.strata:
        assert stratum.reliability
        assert sum(c for _, _, c in stratum.reliability) == stratum.n_cases
