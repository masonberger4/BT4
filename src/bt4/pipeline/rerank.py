"""Post-solve frontier reranking by a learned expression head (a scaffold hook).

The learned expression model is non-local: it cannot decompose over a bounded
codon-trellis context, so it never runs inside the per-move optimizer loop
(CLAUDE.md §15). Instead it reranks the *already-computed* Pareto frontier as a
post-solve pass -- exactly as ``--refine`` scores non-local folding on top of the
exact-DP seed. This module is that hook.

The load-bearing honesty rule (CLAUDE.md §10.5, and the guardrail from the
expression scaffold): **an uncalibrated score never steers delivery.** The hook
always *annotates* each frontier point with the predictor's score for reporting,
but it only re-picks the delivered point (``frontier.chosen``) when the predictor
reports ``calibrated is True``. With the default neutral placeholder
(``calibrated is False``, every score ``0.0``) this is a pure reporting no-op: it
stamps ``expression_calibrated=False`` on each result and changes nothing about
which sequence ships.

This layer composes :mod:`bt4.domain`, :mod:`bt4.biomodels.expression`, and the
:class:`~bt4.pipeline.optimize.FrontierResult` it reranks; it imports nothing
above ``pipeline``.
"""

from __future__ import annotations

import dataclasses

from bt4.biomodels.expression import ExpressionPredictor
from bt4.biomodels.expression import default as expression_default
from bt4.domain import Frontier
from bt4.pipeline.optimize import FrontierResult

__all__ = ["rerank_by_expression"]


def rerank_by_expression(
    result: FrontierResult, predictor: ExpressionPredictor | None = None
) -> FrontierResult:
    """Annotate a frontier with expression scores; re-pick delivery only if calibrated.

    Every result gets ``expression_score`` / ``expression_model`` /
    ``expression_calibrated`` / ``expression_units`` added to its audit. When
    ``predictor.calibrated`` is ``True`` the delivered point (``frontier.chosen``)
    is moved to the highest-scoring result; when it is ``False`` -- as with the
    default placeholder -- the delivered point is left exactly as the solver
    chose it (an uncalibrated score must not steer delivery, CLAUDE.md §10.5).

    Args:
        result: The frontier to rerank (from :func:`bt4.api.frontier`).
        predictor: The expression backend; defaults to
            :func:`bt4.biomodels.expression.default` (the neutral placeholder).

    Returns:
        A new :class:`~bt4.pipeline.optimize.FrontierResult` with annotated
        results and a possibly-updated ``chosen`` index (unchanged unless the
        predictor is calibrated).
    """
    backend = predictor or expression_default()
    annotated = []
    scores: list[float] = []
    for r in result.results:
        er = backend.score_sequence(r.dna)
        scores.append(er.score)
        audit = {
            **dict(r.audit),
            "expression_score": er.score,
            "expression_model": er.model_name,
            "expression_calibrated": er.calibrated,
            "expression_units": er.units,
        }
        annotated.append(dataclasses.replace(r, audit=audit))

    chosen = result.frontier.chosen
    manifest = result.manifest
    if backend.calibrated and scores:
        # A calibrated head may steer delivery to the highest predicted expression.
        chosen = max(range(len(scores)), key=lambda i: scores[i])
        # Because it chose WHICH sequence ships, its identity has to enter the
        # provenance stamp (invariant #9: a stamp must not map to two different
        # delivered sequences). Without this, reranking the same frontier with two
        # different calibrated heads yields different DNA under byte-identical
        # manifests -- and identical to the un-reranked run as well. `_refine`
        # already stamps the folding backend this way, and
        # `assemble_and_rank_candidates` stamps `predictor` / `predictor_calibrated`;
        # this is the one steering site that did not.
        manifest = dataclasses.replace(
            result.manifest,
            extra={
                **dict(result.manifest.extra),
                "expression_model": backend.name,
                "expression_calibrated": "True",
                # Which attestation authorized the calibration that steered this pick;
                # two heads promoted by different claims must stamp differently.
                "expression_attestation": getattr(backend, "attestation_sha256", ""),
            },
        )

    frontier = Frontier(points=result.frontier.points, chosen=chosen)
    return FrontierResult(
        frontier=frontier, results=tuple(annotated), manifest=manifest
    )
