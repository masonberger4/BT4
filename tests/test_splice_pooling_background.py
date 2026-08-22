"""The pooling background is one explicit operating point, not two hidden copies.

``pool_log_odds`` summed ``max(0.0, logit(p))``. That expression has a threshold
inside it: ``logit(p) > 0`` exactly when ``p > 0.5``, so every position at or below
0.5 contributed **nothing** to pooled risk. The docstring described this as "only
above-background positions count" without saying that *background was 0.5*, and the
value was unparameterized.

That made two things silently coupled:

* the **visible** localization cutoff (``DEFAULT_SITE_THRESHOLD``, which the audit,
  the cross-check and ``api`` each spelled as their own literal ``0.5``), and
* the **invisible** pooling background, which decides ``pooled_risk`` and therefore
  ``delta_splicing``, the cross-backend agreement report, and the Studio banner.

Moving the first without the second would have made them disagree: a site flagged at
``p = 0.35`` would contribute exactly zero risk, and a designed variant introducing
five such sites would report ``delta_splicing == 0.0`` — indistinguishable from
introducing none. Since deriving a real operating point is precisely what the
statistical-calibration gate exists to do (published work puts SpliceAI's *delta*
cutoff at 0.2, not 0.5), the trap had to be closed before that number could move.

These tests pin: the default is byte-identical to the old behaviour, the background
is now movable, and every consumer reads one constant.

**The trap was then measured, and it was not a corner case.** Run against the
hash-verified Pangolin weights on the designed-CDS panel — BT4's own regime, a native
CDS and thirty synonymous redesigns for each of three proteins (93 sequences) — *only
6 of those 93* carried any position above 0.5, and all six were designs of a single
protein. Peak scores were 0.323 / 0.435 and varied more than twofold between the
native and its designs, and every sub-background sequence pooled to a risk of exactly
``0.0``. So ``delta_splicing`` was identically zero for those candidates, the rank
agreements computed from those deltas were Spearman correlations of constants (they
printed ``+0.000``), and the CLI reported it as the backend being
unable to rank the candidates — a statement about BT4's pooling, read as a statement
about the model.

Lowering the background is *not* the fix; that is the same uncalibrated knob pointed
somewhere more flattering. The second half of this file pins the fix that was made
instead: :class:`PooledRisk` carries the counts that tell a floored zero from a
measured one, and :func:`pool_top_k_logit` is a **background-free** statistic that
still separates sequences the hinge has flattened.
"""

from __future__ import annotations

import math

import pytest

from bt4.biomodels.splice import DEFAULT_SITE_PROBABILITY, SpliceResult
from bt4.biomodels.splice.audit import DEFAULT_SITE_THRESHOLD
from bt4.biomodels.splice.base import (
    logit,
    pool_log_odds,
    pool_top_k_logit,
    pooled_risk,
    pooled_risk_detail,
)
from bt4.pipeline.splice_crosscheck import DEFAULT_CROSSCHECK_THRESHOLD


def _legacy_pool(probs: list[float], top_k: int = 3) -> float:
    """The exact expression the implementation used before it was parameterized."""
    logits = sorted((logit(p) for p in probs), reverse=True)
    return math.fsum(max(0.0, value) for value in logits[:top_k])


# --------------------------------------------------------------------------
# The default must not have changed


@pytest.mark.parametrize(
    "probs",
    [
        [],
        [0.5],
        [0.1, 0.9, 0.5],
        [0.99, 0.98, 0.97, 0.96],
        [0.2, 0.3],
        [0.0001, 0.9999],
        [0.5, 0.5, 0.5],
    ],
)
def test_default_background_reproduces_the_legacy_expression(probs: list[float]) -> None:
    """Parameterizing the background changed no shipped number."""
    assert pool_log_odds(probs) == _legacy_pool(probs)


def test_default_background_is_the_zero_crossing_of_logit() -> None:
    """0.5 is where the old ``max(0.0, logit(p))`` switched on -- now stated."""
    assert DEFAULT_SITE_PROBABILITY == 0.5
    assert logit(DEFAULT_SITE_PROBABILITY) == pytest.approx(0.0)


# --------------------------------------------------------------------------
# The trap this closes


def test_sub_background_sites_pool_to_exactly_zero() -> None:
    """The defect, pinned: five genuine sites below background register as none.

    This is why the localization cutoff could not be lowered on its own.
    """
    sites = [0.25, 0.30, 0.35, 0.40, 0.45]
    assert pool_log_odds(sites, 3) == 0.0
    # ...so a variant introducing all five is indistinguishable from one introducing
    # nothing, under the default background.
    clean = [0.001] * 5
    assert pool_log_odds(clean, 3) - pool_log_odds(sites, 3) == 0.0


def test_a_lower_background_makes_those_same_sites_register() -> None:
    """With the background moved, the identical scores carry real risk."""
    sites = [0.25, 0.30, 0.35, 0.40, 0.45]
    assert pool_log_odds(sites, 3, background=0.2) > 2.0
    clean = [0.001] * 5
    delta = pool_log_odds(clean, 3, background=0.2) - pool_log_odds(sites, 3, background=0.2)
    assert delta < -2.0  # the variant is now correctly scored as worse


def test_background_is_monotone_in_the_expected_direction() -> None:
    """Lowering the background can only add risk, never remove it."""
    probs = [0.15, 0.35, 0.55, 0.75]
    risks = [pool_log_odds(probs, 4, background=b) for b in (0.8, 0.5, 0.3, 0.1)]
    assert risks == sorted(risks)


def test_pooled_risk_threads_the_background_through() -> None:
    """The `SpliceResult` entry point exposes the same knob, not a fixed 0.5."""
    result = SpliceResult(
        donor=(0.3, 0.4), acceptor=(0.35, 0.0), model_name="stub", calibrated=False
    )
    assert pooled_risk(result) == 0.0
    assert pooled_risk(result, background=0.2) > 0.0


@pytest.mark.parametrize("background", [0.0, 1.0, -0.1, 1.5])
def test_background_must_be_a_probability(background: float) -> None:
    """A background outside (0, 1) has no finite log-odds -- refuse it."""
    with pytest.raises(ValueError, match="background must be in"):
        pool_log_odds([0.5], 3, background=background)


# --------------------------------------------------------------------------
# One constant, not four literals


def test_every_consumer_reads_the_same_operating_point() -> None:
    """The audit, the cross-check and ``api`` no longer carry their own copies.

    Four independent ``0.5`` literals could be moved one at a time; a single
    constant cannot drift against itself.
    """
    assert DEFAULT_SITE_THRESHOLD == DEFAULT_SITE_PROBABILITY
    assert DEFAULT_CROSSCHECK_THRESHOLD == DEFAULT_SITE_PROBABILITY


def test_api_defaults_track_the_constant() -> None:
    """``api.splice_audit`` / ``api.splice_crosscheck`` default to the constant."""
    import inspect

    from bt4 import api

    for fn in (api.splice_audit, api.splice_crosscheck):
        default = inspect.signature(fn).parameters["threshold"].default
        assert default == DEFAULT_SITE_PROBABILITY, fn.__name__


def test_no_bare_threshold_literals_remain_in_source() -> None:
    """Structural: a new ``0.5`` default would reintroduce the drift."""
    from pathlib import Path

    src = Path(__file__).resolve().parent.parent / "src" / "bt4"
    offenders = []
    for path in src.rglob("*.py"):
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if "threshold: float = 0.5" in stripped or "THRESHOLD: float = 0.5" in stripped:
                offenders.append(f"{path.relative_to(src)}:{n}")
    assert not offenders, f"bare 0.5 threshold defaults reintroduced: {offenders}"


# --------------------------------------------------------------------------
# The fix: a floored zero is distinguishable from a measured one


# Pangolin's actual peak per-position scores on the designed-CDS panel, measured with
# the hash-verified weights. Every one is below the 0.5 background; they differ from
# each other by more than twofold. Kept as literals so this regression does not need
# the licensed weights to run.
_MEASURED_PEAKS = (0.128, 0.276, 0.323, 0.435, 0.445)


def test_pool_top_k_logit_takes_no_background_at_all() -> None:
    """The response statistic has no operating point to choose, honestly or otherwise."""
    import inspect

    params = inspect.signature(pool_top_k_logit).parameters
    assert "background" not in params
    assert set(params) == {"probs", "top_k"}


def test_pool_top_k_logit_is_monotone_below_the_background_where_risk_is_flat() -> None:
    """The property the fix turns on: response separates what risk has flattened."""
    weak = [0.10, 0.12, 0.14]
    strong = [0.40, 0.42, 0.44]  # every value still below 0.5
    assert pool_log_odds(weak, 3) == pool_log_odds(strong, 3) == 0.0  # risk: identical
    assert pool_top_k_logit(strong, 3) > pool_top_k_logit(weak, 3)  # response: not


def test_pool_top_k_logit_is_negative_below_the_background() -> None:
    """It is not a risk, and the sign says so rather than a docstring alone."""
    assert pool_top_k_logit(list(_MEASURED_PEAKS), 3) < 0.0


def test_response_equals_risk_exactly_when_the_hinge_does_not_bind() -> None:
    """``logit(0.5) == 0``, so above background the two poolings coincide.

    This is why adding the response changed no shipped number: where BT4's pooled risk
    was ever informative, the response *is* it.
    """
    above = [0.9, 0.8, 0.7, 0.2]
    assert pool_top_k_logit(above, 3) == pytest.approx(pool_log_odds(above, 3))


def test_pooled_risk_detail_reproduces_pooled_risk_exactly() -> None:
    """The attribution rides along; it does not move the number."""
    for donor, acceptor in (
        ((0.9, 0.1), (0.3, 0.8)),
        (_MEASURED_PEAKS, (0.0,) * 5),
        ((), ()),
    ):
        result = SpliceResult(
            donor=tuple(donor), acceptor=tuple(acceptor), model_name="s", calibrated=False
        )
        assert pooled_risk_detail(result).risk == pooled_risk(result)


def test_below_background_separates_the_two_meanings_of_zero() -> None:
    """The defect, now attributable: floored-to-zero versus genuinely-zero."""
    floored = SpliceResult(
        donor=_MEASURED_PEAKS, acceptor=(0.0,) * 5, model_name="s", calibrated=False
    )
    detail = pooled_risk_detail(floored)
    assert detail.risk == 0.0
    assert detail.below_background is True  # ...because nothing cleared background
    assert detail.n_above_background == 0
    assert detail.max_score == 0.445  # ...and this is the magnitude discarded

    # A sequence that genuinely has no signal reports the same risk and a different why.
    quiet = SpliceResult(
        donor=(0.001, 0.002), acceptor=(0.0,), model_name="s", calibrated=False
    )
    assert pooled_risk_detail(quiet).risk == 0.0
    assert pooled_risk_detail(quiet).max_score == 0.002


def test_below_background_is_false_when_there_is_nothing_to_score() -> None:
    """An empty result is not evidence the pooling floored anything."""
    empty = SpliceResult(donor=(), acceptor=(), model_name="s", calibrated=False)
    detail = pooled_risk_detail(empty)
    assert detail.risk == 0.0
    assert detail.n_scored == 0
    assert detail.below_background is False


def test_the_measured_regression_delta_is_zero_while_response_is_not() -> None:
    """The exact shape of what was found, pinned: a real difference reported as none.

    Two sequences whose peak scores differ more than twofold, both entirely below
    background. ``delta_splicing`` cannot tell them apart; the response can.
    """

    def _result(peaks: tuple[float, ...]) -> SpliceResult:
        return SpliceResult(
            donor=peaks, acceptor=(0.0,) * len(peaks), model_name="s", calibrated=False
        )

    native = _result((0.128, 0.110, 0.098))
    design = _result((0.276, 0.240, 0.201))
    native_d = pooled_risk_detail(native)
    design_d = pooled_risk_detail(design)

    assert native_d.risk - design_d.risk == 0.0  # delta_splicing: exactly nothing
    assert native_d.below_background and design_d.below_background
    assert native_d.response - design_d.response < -1.0  # the design scores higher


# --------------------------------------------------------------------------
# The localization cutoff and the pooling background now move together


def test_the_audit_pools_against_the_threshold_it_localizes_at() -> None:
    """Lowering the visible cutoff must lower the hidden one, or they disagree.

    The trap this file opens with, closed at the consumer: a site flagged at 0.35 used
    to contribute exactly zero to the pooled risk that is meant to summarize the flags,
    because the audit passed the caller's ``threshold`` to localization and let pooling
    keep the default 0.5.
    """
    from bt4.biomodels.splice import ConsensusPwmSplicePredictor
    from bt4.biomodels.splice.audit import audit_splice

    predictor = ConsensusPwmSplicePredictor()
    candidate = "ATGGTAAGTACCGGCGTGAGTGCCAAATTTGGC"
    reference = "ATGAAACCCTTTAAACCCTTTAAACCCTTTGGC"

    low = audit_splice([predictor], [candidate], reference, threshold=0.2)
    flags = low.candidates[0].by_backend[0].flags
    if any(flag.score <= DEFAULT_SITE_PROBABILITY for flag in flags):
        # A site flagged below the default background must now carry pooled risk.
        assert low.candidates[0].by_backend[0].pooled_risk > 0.0


def test_the_audit_default_threshold_is_unchanged() -> None:
    """Coupling the two must not have moved the shipped default path."""
    from bt4.biomodels.splice import ConsensusPwmSplicePredictor
    from bt4.biomodels.splice.audit import DEFAULT_SITE_THRESHOLD, audit_splice

    predictor = ConsensusPwmSplicePredictor()
    seq = "ATGGTAAGTACCGGCGTGAGTGCCAAATTTGGC"
    ref = "ATGAAACCCTTTAAACCCTTTAAACCCTTTGGC"
    explicit = audit_splice([predictor], [seq], ref, threshold=DEFAULT_SITE_THRESHOLD)
    assert DEFAULT_SITE_THRESHOLD == DEFAULT_SITE_PROBABILITY
    # The default background and the default threshold are the same number, so the
    # audit's pooled risk equals the bare `pooled_risk` it used to call.
    assert explicit.candidates[0].by_backend[0].pooled_risk == pooled_risk(
        predictor.score_sequence(seq), explicit.top_k
    )
