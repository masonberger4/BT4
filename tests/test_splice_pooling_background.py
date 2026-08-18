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
"""

from __future__ import annotations

import math

import pytest

from bt4.biomodels.splice import DEFAULT_SITE_PROBABILITY, SpliceResult
from bt4.biomodels.splice.audit import DEFAULT_SITE_THRESHOLD
from bt4.biomodels.splice.base import logit, pool_log_odds, pooled_risk
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
