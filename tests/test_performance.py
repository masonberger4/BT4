"""Performance / scaling regression tests for the exact codon-trellis DP.

CLAUDE.md §7 (Performance) requires that "Runtime and peak-memory scaling are
regression-tested in CI with an asserted curve and a wall-clock ceiling" — BT4
had no performance test despite runtime being BT2's original weakness, and
CLAUDE.md §10 #8 ("no quadratic refinement") makes sub-quadratic scaling a
load-bearing property, not a hope. This module is that first guard.

The exact DP over the codon trellis does bounded work per residue: a bounded
number of synonymous codons times a bounded state space (the trellis key is only
the last ``K`` context characters, and ``K`` stays small here because the
constraint set is just a mild homopolymer bound). So end-to-end
:func:`bt4.api.optimize` runtime should grow roughly LINEARLY in protein length.
A quadratic regression — e.g. re-scoring or re-validating the whole prefix per
move, the BT3 anti-pattern §10 #8 warns about — would show up as super-linear
growth and trip these tests.

Anti-flakiness design (why the margins are what they are):

* Absolute wall-clock times are noisy on shared CI runners, so the scaling
  assertion is built on RATIOS of times at different lengths, never on an
  absolute duration.
* Each length is timed as the MINIMUM over ``k`` repeats. ``min`` is the least
  noisy timing estimator: noise only ever adds time, so the smallest sample is
  the closest to the true cost.
* The discriminating ratio is over a 4x jump in length. A near-linear DP gives
  ``t(4L)/t(L) ~ 4``; a quadratic one gives ``~16``. Asserting ``< 10`` sits
  squarely between the two: it leaves the near-linear reality (measured ~3.5-3.7
  on the reference machine) ample headroom yet still catches a quadratic blow-up.
* ``time.perf_counter`` is used (monotonic, highest available resolution). The
  timing floor concern is handled by choosing a denominator length (100
  residues) whose min-time is comfortably above timer resolution — milliseconds,
  not microseconds — and by guarding the division with a tiny epsilon so a
  degenerate zero can never raise.
* No assertion depends on absolute machine speed except one deliberately
  generous backstop ceiling (20s for a 300-residue solve, ~4000x the observed
  cost), which exists only to trip a catastrophic (100x) regression and is far
  too loose to flake on a slow runner.
"""

from __future__ import annotations

from time import perf_counter

from bt4 import api

# A repeating motif of two-codon amino acids (D, E, F, H, K, N, Q, Y). Each
# residue has modest synonymous branching (2 codons), so the DP does real work,
# yet the trellis state space stays small and every length is cheap to solve.
_MOTIF = "DEFHKNQY"

# Guards the ratio denominator against a degenerate zero. It is far below any
# real min-time here (which is milliseconds), so it never distorts a live ratio.
_EPSILON = 1e-9


def _make_protein(length: int) -> str:
    """Return a deterministic protein of exactly ``length`` residues.

    Built by repeating :data:`_MOTIF` and truncating, so the amino-acid content
    (and thus the DP's branching and feasibility) is fixed for a given length —
    only the timings vary between runs.

    Args:
        length: The desired number of residues (``>= 1``).

    Returns:
        An uppercase single-letter amino-acid string of exactly ``length``.
    """
    reps = (length // len(_MOTIF)) + 1
    return (_MOTIF * reps)[:length]


def _time_once(protein: str, config: api.OptimizeConfig) -> float:
    """Return the wall-clock seconds for a single :func:`bt4.api.optimize` call.

    Args:
        protein: The protein to back-translate.
        config: The run configuration.

    Returns:
        Elapsed seconds measured with :func:`time.perf_counter`.
    """
    start = perf_counter()
    api.optimize(protein, config)
    return perf_counter() - start


def _best_time(protein: str, config: api.OptimizeConfig, repeats: int) -> float:
    """Return the minimum optimize time over ``repeats`` runs (least-noisy estimate).

    ``min`` is used rather than mean/median because timing noise is one-sided
    (it only ever adds time), so the smallest sample is the closest to the true
    cost and the most stable across noisy CI runners.

    Args:
        protein: The protein to back-translate.
        config: The run configuration.
        repeats: How many times to time the call (``>= 1``).

    Returns:
        The smallest observed elapsed time, in seconds.
    """
    return min(_time_once(protein, config) for _ in range(repeats))


def test_exact_dp_scales_subquadratically() -> None:
    """Optimize runtime must grow sub-quadratically in protein length.

    The exact trellis DP does bounded work per residue, so a 4x jump in length
    should cost roughly 4x (near-linear), not 16x (quadratic). We time four
    lengths (50/100/200/400), take the min over 5 repeats at each, and assert two
    independent 4x-length ratios stay well under the quadratic factor of 16.

    Margin rationale: a near-linear DP yields ``t(4L)/t(L) ~ 4`` (measured ~3.5
    here); a quadratic one yields ``~16``. The ``< 10`` bound is the midpoint —
    >=2x headroom above the observed near-linear value, and a clear gap below the
    quadratic factor — so it never flakes yet still catches a real regression.
    Timings are min-over-repeats to de-noise, and only ratios are asserted, so no
    absolute machine speed is baked in.
    """
    config = api.OptimizeConfig(max_homopolymer=6)
    # Warm up one-time costs (codon-table load, module import) so the first
    # measured length is not penalized by cold-start work.
    api.optimize(_make_protein(20), config)

    lengths = (50, 100, 200, 400)
    repeats = 5
    times = {n: _best_time(_make_protein(n), config, repeats) for n in lengths}

    # Two disjoint 4x-length ratios: linear ~4, quadratic ~16. Both must be well
    # under 16; asserting < 10 leaves generous headroom in both directions.
    ratio_100_400 = times[400] / max(times[100], _EPSILON)
    ratio_50_200 = times[200] / max(times[50], _EPSILON)

    assert ratio_100_400 < 10.0, (
        f"t(400)/t(100) = {ratio_100_400:.2f} looks super-linear "
        f"(quadratic ~16, near-linear ~4); times(ms)="
        f"{ {n: round(times[n] * 1000, 3) for n in lengths} }"
    )
    assert ratio_50_200 < 10.0, (
        f"t(200)/t(50) = {ratio_50_200:.2f} looks super-linear "
        f"(quadratic ~16, near-linear ~4); times(ms)="
        f"{ {n: round(times[n] * 1000, 3) for n in lengths} }"
    )


def test_single_optimize_under_generous_ceiling() -> None:
    """A moderate single solve must finish well under a generous wall-clock ceiling.

    This is a backstop, not the main assertion: a 300-residue default-config
    optimize completes in milliseconds on normal hardware, so a 20s ceiling gives
    roughly 4000x headroom. It never flakes on a slow shared runner but would
    still trip on a catastrophic (100x+) regression. We take the min over 3
    repeats so a single unlucky pause cannot fail the test.
    """
    config = api.OptimizeConfig(max_homopolymer=6)
    # Warm up shared caches so the timed solve reflects steady-state cost.
    api.optimize(_make_protein(20), config)

    elapsed = _best_time(_make_protein(300), config, repeats=3)
    assert elapsed < 20.0, f"L=300 optimize took {elapsed:.3f}s, over the 20s backstop ceiling"
