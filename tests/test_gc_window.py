"""Tests for :class:`~bt4.constraints.gc_window.GcWindowConstraint` (invariant #3).

The windowed GC rule is the one the synthesis vendors actually specify (bound GC
in *every* window, not just overall), and it is genuinely LOCAL -- so it must be
exact in the trellis and keep a proven-optimal certificate. The load-bearing
property is the same as every other local constraint: a sequence built respecting
``ok_suffix`` has zero hard violations under ``validate``, and the declared
``context_len`` must actually suffice, so a window straddling a codon boundary is
still caught.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from bt4._accel import gc_count
from bt4.constraints.gc_window import GcWindowConstraint
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.domain.result import Severity
from bt4.domain.scope import Scope

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=60)


def _build_feasible(protein: str, constraint: GcWindowConstraint) -> str:
    """Greedily first-fit a codon per residue that passes ``ok_suffix``."""
    dna = ""
    for aa in protein:
        for codon in synonymous_codons(aa):
            if constraint.ok_suffix(dna, codon):
                dna += codon
                break
        else:  # dead end -- stop with the valid prefix.
            break
    return dna


def _hard(constraint: GcWindowConstraint, dna: str) -> list[object]:
    return [v for v in constraint.validate(dna) if v.severity is Severity.HARD]


def _windows(dna: str, window: int) -> list[float]:
    """Reference GC fractions of every full window (the independent recompute)."""
    return [
        gc_count(dna[i : i + window]) / window for i in range(len(dna) - window + 1)
    ]


# --------------------------------------------------------------------------- #
# Invariant #3: ok_suffix-respecting builds have zero hard violations.
# --------------------------------------------------------------------------- #


@given(
    protein=_PROTEIN,
    window=st.integers(min_value=6, max_value=30),
    lo=st.floats(min_value=0.0, max_value=0.35),
    hi=st.floats(min_value=0.55, max_value=1.0),
)
def test_gc_window_ok_suffix_implies_validate_clean(
    protein: str, window: int, lo: float, hi: float
) -> None:
    constraint = GcWindowConstraint(window, lo, hi)
    dna = _build_feasible(protein, constraint)
    assume(len(dna) >= 3)
    assert _hard(constraint, dna) == []
    # The build is a genuine synonymous back-translation of a protein prefix.
    assert translate(dna) == protein[: len(dna) // 3]


@given(
    protein=_PROTEIN,
    window=st.integers(min_value=6, max_value=24),
    lo=st.floats(min_value=0.0, max_value=0.35),
    hi=st.floats(min_value=0.55, max_value=1.0),
)
def test_validate_agrees_with_an_independent_window_recompute(
    protein: str, window: int, lo: float, hi: float
) -> None:
    """`validate` flags exactly when an independent sliding recompute says so."""
    constraint = GcWindowConstraint(window, lo, hi)
    dna = _build_feasible(protein, constraint) or "ATG"
    fractions = _windows(dna, window)
    expected_bad = any(
        f < constraint.min_count / window or f > constraint.max_count / window
        for f in fractions
    )
    assert bool(_hard(constraint, dna)) == expected_bad


# --------------------------------------------------------------------------- #
# Positive detection.
# --------------------------------------------------------------------------- #


def test_validate_flags_a_gc_rich_stretch() -> None:
    constraint = GcWindowConstraint(10, 0.2, 0.6)
    violations = list(constraint.validate("AAAAAGGGGGGGGGGAAAAA"))
    assert violations
    assert all(v.constraint == "gc_window" for v in violations)
    assert all(v.severity is Severity.HARD for v in violations)
    assert all("above max" in v.detail for v in violations)


def test_validate_flags_a_gc_poor_stretch() -> None:
    constraint = GcWindowConstraint(10, 0.4, 1.0)
    violations = list(constraint.validate("AAAAAAAAAAAAAAA"))
    assert violations
    assert all("below min" in v.detail for v in violations)


def test_low_and_high_stretches_are_reported_separately() -> None:
    # A GC-poor run followed by a GC-rich run: two violations, never one span
    # described as both too low and too high.
    constraint = GcWindowConstraint(8, 0.3, 0.7)
    directions = [
        "below min" if "below min" in v.detail else "above max"
        for v in constraint.validate("A" * 16 + "GC" * 8)
    ]
    assert "below min" in directions
    assert "above max" in directions


def test_balanced_sequence_is_clean() -> None:
    assert list(GcWindowConstraint(10, 0.2, 0.6).validate("ACGTACGTACGTACGTACGT")) == []


def test_sequence_shorter_than_window_has_no_full_window() -> None:
    # No full-length window exists, so nothing is constrained -- and ok_suffix
    # agrees, which is what keeps the two sides consistent on short inputs.
    constraint = GcWindowConstraint(20, 0.4, 0.5)
    assert list(constraint.validate("GGGCCC")) == []


def test_violation_count_is_the_number_of_offending_windows() -> None:
    # One violation per offending window (not per merged region). The granularity
    # is load-bearing: the refinement pass drives this count down, and a merged
    # count would be a plateau it could not descend. Each violation spans exactly
    # one window.
    constraint = GcWindowConstraint(6, 0.0, 0.5)
    dna = "AAA" + "GC" * 12 + "AAA"
    violations = list(constraint.validate(dna))
    expected = sum(
        1
        for i in range(len(dna) - 6 + 1)
        if gc_count(dna[i : i + 6]) > constraint.max_count
    )
    assert len(violations) == expected > 1
    assert all(v.end - v.start == 6 for v in violations)
    # Ascending window-start order, so a viewer can merge adjacent spans itself.
    assert [v.start for v in violations] == sorted(v.start for v in violations)


def test_refinement_can_descend_the_per_window_count() -> None:
    # The end-to-end consequence of the granularity above: at a vendor-typical
    # 50-nt window (too wide for the trellis, so refinement-enforced) the annealer
    # actually reaches a compliant sequence instead of stalling on a plateau.
    from bt4 import api

    result = api.optimize(
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF",
        api.OptimizeConfig(
            gc_window_nt=50, gc_window_min=0.25, gc_window_max=0.65, max_homopolymer=None
        ),
    )
    assert result.audit["gc_window_enforced"] == "clean"
    assert result.audit["gc_window_residual"] == 0
    fractions = _windows(result.dna, 50)
    assert all(0.25 - 1e-9 <= f <= 0.65 + 1e-9 for f in fractions)
    # Refinement ran, so the certificate must say so rather than claim optimality.
    assert not result.certificate.is_proven_optimal


# --------------------------------------------------------------------------- #
# context_len suffices: boundary-crossing windows are vetoed by ok_suffix.
# --------------------------------------------------------------------------- #


def test_context_len_is_window_minus_one() -> None:
    assert GcWindowConstraint(50).context_len() == 49
    assert GcWindowConstraint(1).context_len() == 0


def test_ok_suffix_vetoes_window_completed_by_incoming_codon() -> None:
    constraint = GcWindowConstraint(6, 0.0, 0.5)  # max 3 GC per 6-nt window
    # prefix tail "GGG" + codon "GCA": the window GGGGCA has 5 GC > 3.
    assert constraint.ok_suffix("AAAGGG", "GCA") is False
    # A codon that keeps the window in bounds is fine.
    assert constraint.ok_suffix("AAAGGG", "AAA") is True


def test_ok_suffix_enforces_the_minimum_too() -> None:
    constraint = GcWindowConstraint(6, 0.5, 1.0)  # at least 3 GC per 6-nt window
    # Window "GCAAAA" has only 2 GC -- below the floor.
    assert constraint.ok_suffix("GCA", "AAA") is False
    assert constraint.ok_suffix("GCG", "CGC") is True


def test_ok_suffix_only_judges_windows_ending_in_the_new_codon() -> None:
    # A window lying entirely inside the prefix is not this codon's fault, so it
    # must not be re-vetoed. Here "GGGG" (4/4 GC, over the limit) sits at the very
    # start of the prefix, outside the context window; appending a benign codon
    # completes only in-bounds windows, so the extension is allowed.
    constraint = GcWindowConstraint(4, 0.0, 0.5)
    assert constraint.ok_suffix("GGGGAAAA", "AAA") is True
    # A window the incoming codon *does* complete is still judged: the codon's
    # first base closes "GGGA", which is 3/4 GC and over the limit.
    assert constraint.ok_suffix("GGGG", "AAA") is False


# --------------------------------------------------------------------------- #
# The rule stays exact in the trellis (this is the point of making it LOCAL).
# --------------------------------------------------------------------------- #


def test_narrow_window_stays_proven_optimal_end_to_end() -> None:
    """A window the trellis can carry is exact -- no honesty downgrade."""
    from bt4 import api
    from bt4.pipeline.optimize import _GC_WINDOW_TRELLIS_MAX_NT

    window = _GC_WINDOW_TRELLIS_MAX_NT
    config = api.OptimizeConfig(
        gc_window_nt=window,
        gc_window_min=0.20,
        gc_window_max=0.70,
        max_homopolymer=None,
    )
    result = api.optimize("MAALKHETQWYCDEFGHIKLMNPQRS", config)
    assert result.certificate.is_proven_optimal
    assert result.metrics.hard_violations == 0
    # And the delivered sequence really does satisfy the bound (invariant #2).
    fractions = _windows(result.dna, window)
    assert fractions
    assert all(0.20 - 1e-9 <= f <= 0.70 + 1e-9 for f in fractions)


def test_wide_window_is_refinement_enforced_and_says_so() -> None:
    """A window too wide for the trellis is enforced, reported, and never claims optimality.

    The vendor-typical 50 nt is past the trellis bound, so it is routed to the
    refinement pass. What must never happen is the BT3 failure mode: silently
    capping the context and still badging the result proven-optimal.
    """
    from bt4 import api
    from bt4.pipeline.optimize import _GC_WINDOW_TRELLIS_MAX_NT

    window = 50
    assert window > _GC_WINDOW_TRELLIS_MAX_NT, "premise: this window is refinement-routed"
    result = api.optimize(
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF",
        api.OptimizeConfig(
            gc_window_nt=window,
            gc_window_min=0.25,
            gc_window_max=0.65,
            max_homopolymer=None,
        ),
    )
    # Enforcement is reported per-rule, and optimality is not claimed.
    assert result.audit["gc_window_enforced"] in {"clean", "partial"}
    assert not result.certificate.is_proven_optimal
    # Reported residual == recomputed residual (invariant #2).
    fractions = _windows(result.dna, window)
    recomputed = sum(1 for f in fractions if not (0.25 - 1e-9 <= f <= 0.65 + 1e-9))
    assert result.audit["gc_window_residual"] == recomputed


# --------------------------------------------------------------------------- #
# Names, scope, and configuration guards.
# --------------------------------------------------------------------------- #


def test_name_scope_and_penalty() -> None:
    constraint = GcWindowConstraint(20, 0.3, 0.7)
    assert constraint.name == "gc_window"
    assert constraint.scope() is Scope.LOCAL
    assert constraint.penalty("GCG", "CGC") == 0.0


def test_count_thresholds_are_exact_at_representable_bounds() -> None:
    # 0.65 * 20 == 13.0 exactly in decimal but not in binary; the epsilon guard
    # must keep the threshold at 13, not 12.
    assert GcWindowConstraint(20, 0.0, 0.65).max_count == 13
    assert GcWindowConstraint(20, 0.65, 1.0).min_count == 13


def test_rejects_invalid_configuration() -> None:
    with pytest.raises(ValueError, match="window"):
        GcWindowConstraint(0)
    with pytest.raises(ValueError, match="gc_min"):
        GcWindowConstraint(10, -0.1, 1.0)
    with pytest.raises(ValueError, match="gc_max"):
        GcWindowConstraint(10, 0.0, 1.5)
    with pytest.raises(ValueError, match="must not exceed"):
        GcWindowConstraint(10, 0.8, 0.2)
