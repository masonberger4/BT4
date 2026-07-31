"""Optional headless HTTP API over :mod:`bt4.api`.

This is BT4's ``service`` layer: a thin, async FastAPI wrapper that exposes the
same stable engine the CLI and desktop app call. It imports **only** the public
``bt4.api`` surface (plus the top-level ``bt4`` package for its version string),
never the optimizer, pipeline, biomodels, or domain layers — the layering
contract is enforced by import-linter.

The service adds no science of its own: every route delegates to ``bt4.api`` and
serializes the result. Engine failures surface as ``ValueError`` (the infeasible
case is a subclass) and are translated into HTTP 400 responses.

Run it with ``uvicorn bt4.service.api:app``.
"""

from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from bt4 import __version__, api

__all__ = ["ConfigModel", "OptimizeRequest", "ValidateRequest", "app", "create_app"]


class ConfigModel(BaseModel):
    """Request-body mirror of :class:`bt4.api.OptimizeConfig`.

    Fields and defaults match the dataclass one-for-one so a JSON body can omit
    anything it does not care about. :meth:`to_config` converts it back into the
    immutable engine config.

    Attributes:
        organism: Codon-usage table key or alias.
        gc_target: Desired GC fraction in ``[0, 1]``.
        cai_weight: Weight on the CAI objective in a single solve.
        gc_weight: Weight on the GC-proximity objective in a single solve.
        max_homopolymer: Longest allowed single-base run, or ``None`` to disable.
        forbidden_motifs: Substrings that may not appear in the sequence.
        avoid_reverse_complement: Also forbid each motif's reverse complement.
        beam: ``None`` for exact DP, or an int beam-width cap.
        seed: Master seed recorded in the manifest.
    """

    organism: str = "homo_sapiens"
    gc_target: float = 0.55
    cai_weight: float = 1.0
    gc_weight: float = 0.0
    max_homopolymer: int | None = 6
    forbidden_motifs: list[str] = []
    avoid_reverse_complement: bool = True
    beam: int | None = None
    seed: int = 0

    def to_config(self) -> api.OptimizeConfig:
        """Build an :class:`bt4.api.OptimizeConfig` from this request model.

        Returns:
            An immutable engine config with ``forbidden_motifs`` frozen to a
            tuple.
        """
        return api.OptimizeConfig(
            organism=self.organism,
            gc_target=self.gc_target,
            cai_weight=self.cai_weight,
            gc_weight=self.gc_weight,
            max_homopolymer=self.max_homopolymer,
            forbidden_motifs=tuple(self.forbidden_motifs),
            avoid_reverse_complement=self.avoid_reverse_complement,
            beam=self.beam,
            seed=self.seed,
        )


class OptimizeRequest(BaseModel):
    """Body for ``POST /optimize`` and ``POST /frontier``.

    Attributes:
        protein: A stop-free single-letter amino-acid string.
        config: Run configuration; defaults to the engine defaults.
        steps: Number of scalarization weights swept for the frontier route.
    """

    protein: str
    config: ConfigModel = ConfigModel()
    steps: int = 11


class ValidateRequest(BaseModel):
    """Body for ``POST /validate``.

    Attributes:
        dna: An ACGT coding sequence to audit.
        config: Run configuration whose constraints define the audit.
    """

    dna: str
    config: ConfigModel = ConfigModel()


def create_app() -> FastAPI:
    """Build the BT4 HTTP service application.

    Returns:
        A :class:`fastapi.FastAPI` app exposing health, organism listing,
        optimize, frontier, and validate routes, each delegating to
        :mod:`bt4.api`.
    """
    app = FastAPI(title="BT4 service", version=__version__)

    @app.get("/health")
    def health() -> dict[str, str]:
        """Report liveness and the engine version."""
        return {"status": "ok", "version": __version__}

    @app.get("/organisms")
    def organisms() -> dict[str, list[str]]:
        """List the codon-usage organisms the engine can target."""
        return {"organisms": list(api.available_organisms())}

    @app.post("/optimize")
    def optimize(request: OptimizeRequest) -> dict[str, object]:
        """Back-translate a protein into an optimized coding sequence."""
        try:
            result = api.optimize(request.protein, request.config.to_config())
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return api.result_to_dict(result)

    @app.post("/frontier")
    def frontier(request: OptimizeRequest) -> dict[str, object]:
        """Compute the CAI/GC Pareto frontier for a protein."""
        try:
            frontier_result = api.frontier(
                request.protein, request.config.to_config(), request.steps
            )
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return {
            "points": [api.result_to_dict(r) for r in frontier_result.results],
            "delivered_index": frontier_result.frontier.chosen,
        }

    @app.post("/validate")
    def validate(request: ValidateRequest) -> dict[str, object]:
        """Audit a caller-supplied coding sequence (no optimization)."""
        try:
            report = api.validate(request.dna, request.config.to_config())
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return {
            "dna": report.dna,
            "is_feasible": report.is_feasible,
            "violations": [
                {
                    "constraint": v.constraint,
                    "severity": v.severity.value,
                    "start": v.start,
                    "end": v.end,
                    "detail": v.detail,
                }
                for v in report.violations
            ],
            "metrics": {
                "gc": report.metrics.gc,
                "length_nt": report.metrics.length_nt,
                "hard_violations": report.metrics.hard_violations,
                "soft_violations": report.metrics.soft_violations,
            },
        }

    return app


app = create_app()
