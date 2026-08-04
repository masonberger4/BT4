"""Cross-backend splice agreement -- a first-class uncertainty signal.

CLAUDE.md section 6 makes backend *agreement* a design goal, not redundancy:
"Running *both* [SpliceAI and Pangolin] and reporting their agreement /
disagreement is a first-class uncertainty signal (section 8), not redundancy."
This module computes that signal for any set of
:class:`~bt4.biomodels.splice.base.SplicePredictor` backends -- the labeled PWM
baseline, Pangolin, and (later) SpliceAI -- over a panel of candidate sequences
scored against one reference.

The comparison is done at the **backend-agnostic** level -- pooled
Delta-splicing per candidate -- because backends differ in what per-position
scores mean (Pangolin reports one combined splice-site probability; the PWM
baseline reports separated donor / acceptor). Delta-splicing (larger is better)
and its ranking across candidates are directly comparable across backends, so
this module reports:

* each backend's Delta-splicing vector over the candidates;
* pairwise **Spearman rank correlation** of those vectors -- do the backends
  agree on *which* redesigns lower splice risk?
* the fraction of candidates on which all backends agree on the *sign* of
  Delta-splicing (did the redesign help or hurt?).

It **reports, it does not judge**: no backend is declared right, and a low
correlation is surfaced as uncertainty, not resolved. Pure standard library (no
numpy); it depends only on the :class:`SplicePredictor` protocol and
:mod:`bt4.domain` (transitively via the protocol).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels.splice.base import DEFAULT_TOP_K, SplicePredictor, pooled_risk

__all__ = [
    "AgreementReport",
    "backend_agreement",
    "pearson",
    "spearman",
]


def _ranks(values: Sequence[float]) -> list[float]:
    """Return fractional (tie-averaged) ranks of ``values``.

    Ties share the average of the ranks they span, so Spearman on the ranks is
    well defined in the presence of equal scores.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0  # 1-based, averaged over the tie block
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    """Return the Pearson correlation of ``x`` and ``y``.

    Args:
        x: First series.
        y: Second series (same length as ``x``).

    Returns:
        The Pearson correlation in ``[-1, 1]``, or ``0.0`` when either series has
        zero variance (correlation is undefined; ``0.0`` is the honest "no
        detectable linear relationship" default).

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    if len(x) != len(y):
        raise ValueError(f"series differ in length: {len(x)} vs {len(y)}")
    n = len(x)
    if n < 2:
        raise ValueError("correlation needs at least two points")
    mean_x = math.fsum(x) / n
    mean_y = math.fsum(y) / n
    dx = [xi - mean_x for xi in x]
    dy = [yi - mean_y for yi in y]
    cov = math.fsum(a * b for a, b in zip(dx, dy, strict=True))
    var_x = math.fsum(a * a for a in dx)
    var_y = math.fsum(b * b for b in dy)
    if var_x == 0.0 or var_y == 0.0:
        return 0.0
    return cov / math.sqrt(var_x * var_y)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Return the Spearman rank correlation of ``x`` and ``y``.

    Computed as the Pearson correlation of the tie-averaged ranks, so it measures
    monotonic (rank) agreement rather than linear agreement.

    Raises:
        ValueError: If the series differ in length or are shorter than 2.
    """
    return pearson(_ranks(x), _ranks(y))


@dataclass(frozen=True, slots=True)
class AgreementReport:
    """Cross-backend splice-risk agreement over a candidate panel.

    Attributes:
        backends: The backend names, in input order.
        delta_by_backend: For each backend, its Delta-splicing (larger is better)
            for each candidate vs the reference, aligned to the candidate order.
        rank_correlations: Pairwise Spearman rank correlation of the
            Delta-splicing vectors, keyed by ``(name_a, name_b)`` for each
            unordered backend pair. Empty when fewer than two candidates (rank
            correlation is undefined) or fewer than two backends.
        sign_agreement: Fraction of candidates on which *all* backends agree on
            the sign of Delta-splicing (treating a value within ``sign_epsilon``
            of zero as "no change"). ``1.0`` when there is only one backend.
        sign_epsilon: The magnitude below which a Delta-splicing is treated as
            zero for the sign-agreement count.
        n_candidates: Number of candidates scored.
    """

    backends: tuple[str, ...]
    delta_by_backend: dict[str, tuple[float, ...]]
    rank_correlations: dict[tuple[str, str], float]
    sign_agreement: float
    sign_epsilon: float
    n_candidates: int


def _sign(value: float, epsilon: float) -> int:
    """Return -1 / 0 / +1 for ``value``, with a dead-band of ``epsilon`` at zero."""
    if value > epsilon:
        return 1
    if value < -epsilon:
        return -1
    return 0


def backend_agreement(
    predictors: Sequence[SplicePredictor],
    candidates: Sequence[str],
    reference: str,
    *,
    top_k: int | None = None,
    sign_epsilon: float = 1e-9,
) -> AgreementReport:
    """Compare splice backends over ``candidates`` scored against ``reference``.

    For each backend and each candidate, computes Delta-splicing (larger is
    better) as the contract's ``pooled_risk(reference) - pooled_risk(candidate)``
    -- equal to :meth:`SplicePredictor.delta_splicing`, but scoring the reference
    once per backend rather than per candidate -- then reports pairwise Spearman
    rank agreement and sign agreement across backends.

    Args:
        predictors: The splice backends to compare. Names must be distinct.
        candidates: Candidate coding sequences (e.g. frontier members) to rank.
        reference: The reference sequence each candidate is compared against.
        top_k: Optional pooling depth override. When ``None`` each backend is
            pooled at its own depth (its ``top_k`` attribute, else
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`); when set, every
            backend is pooled at this ``top_k`` uniformly. Either way the reference
            is scored once per backend and reused across candidates.
        sign_epsilon: Dead-band around zero for sign agreement.

    Returns:
        An :class:`AgreementReport`.

    Raises:
        ValueError: If ``predictors`` or ``candidates`` is empty, or two
            predictors share a name.
    """
    if not predictors:
        raise ValueError("need at least one splice backend to compare")
    if not candidates:
        raise ValueError("need at least one candidate sequence to compare")
    names = tuple(p.name for p in predictors)
    if len(set(names)) != len(names):
        raise ValueError(f"backend names must be distinct, got {names}")

    delta_by_backend: dict[str, tuple[float, ...]] = {}
    for predictor in predictors:
        # Pool at `top_k` when overridden, else at the backend's own depth (its
        # `top_k` attribute, or DEFAULT_TOP_K if it exposes none). The contract
        # defines delta_splicing as pooled_risk(reference) - pooled_risk(designed),
        # so scoring the reference ONCE per backend and reusing it equals calling
        # delta_splicing per candidate -- while avoiding re-running the reference
        # through a heavy CNN C times (CLAUDE.md section 7, "everything incremental").
        depth = top_k if top_k is not None else getattr(predictor, "top_k", DEFAULT_TOP_K)
        ref_risk = pooled_risk(predictor.score_sequence(reference), depth)
        deltas = [
            ref_risk - pooled_risk(predictor.score_sequence(cand), depth)
            for cand in candidates
        ]
        delta_by_backend[predictor.name] = tuple(deltas)

    rank_correlations: dict[tuple[str, str], float] = {}
    if len(names) >= 2 and len(candidates) >= 2:
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                rank_correlations[(names[i], names[j])] = spearman(
                    delta_by_backend[names[i]], delta_by_backend[names[j]]
                )

    if len(names) == 1:
        sign_agreement = 1.0
    else:
        agree = 0
        for idx in range(len(candidates)):
            signs = {
                _sign(delta_by_backend[name][idx], sign_epsilon) for name in names
            }
            if len(signs) == 1:
                agree += 1
        sign_agreement = agree / len(candidates)

    return AgreementReport(
        backends=names,
        delta_by_backend=delta_by_backend,
        rank_correlations=rank_correlations,
        sign_agreement=sign_agreement,
        sign_epsilon=sign_epsilon,
        n_candidates=len(candidates),
    )
