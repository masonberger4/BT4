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
from pydantic import BaseModel, ConfigDict

from bt4 import __version__, api

__all__ = ["ConfigModel", "OptimizeRequest", "ValidateRequest", "app", "create_app"]


class ConfigModel(BaseModel):
    """Request-body mirror of :class:`bt4.api.OptimizeConfig`.

    Every field of :class:`bt4.api.OptimizeConfig` is mirrored here with the same
    default, so a JSON body can set any knob the engine supports and omit the
    rest. ``extra="forbid"`` rejects unknown keys with a 422 rather than silently
    dropping them, so a caller never believes a mistyped option took effect.
    :meth:`to_config` converts it back into the immutable engine config.
    """

    model_config = ConfigDict(extra="forbid")

    organism: str = "homo_sapiens"
    gc_target: float = 0.55
    cai_weight: float = 1.0
    tai_weight: float = 0.0
    gc_weight: float = 0.0
    cpb_weight: float = 0.0
    cpb_reference_cds: list[str] = []
    max_homopolymer: int | None = 6
    max_gc_run: int | None = None
    max_repeat_length: int | None = None
    forbidden_motifs: list[str] = []
    forbidden_presets: list[str] = []
    avoid_reverse_complement: bool = True
    restriction_enzymes: list[str] = []
    ramp_weight: float = 0.0
    ramp_codons: int = 35
    cpg_weight: float = 0.0
    cpg_mode: str = "deplete"
    minmax_weight: float = 0.0
    minmax_direction: str = "max"
    tandem_unit: int | None = None
    tandem_copies: int = 3
    inverted_stem: int | None = None
    inverted_loop: int = 0
    avoid_internal_start: bool = False
    avoid_uorf: bool = False
    uorf_region_nt: int = 100
    refine: bool = False
    refine_iterations: int = 2000
    folding_weight: float = 1.0
    gc_min: int | None = None
    gc_max: int | None = None
    beam: int | None = None
    seed: int = 0

    def to_config(self) -> api.OptimizeConfig:
        """Build an :class:`bt4.api.OptimizeConfig` from this request model.

        Returns:
            An immutable engine config with list fields frozen to tuples.
        """
        return api.OptimizeConfig(
            organism=self.organism,
            gc_target=self.gc_target,
            cai_weight=self.cai_weight,
            tai_weight=self.tai_weight,
            gc_weight=self.gc_weight,
            cpb_weight=self.cpb_weight,
            cpb_reference_cds=tuple(self.cpb_reference_cds),
            max_homopolymer=self.max_homopolymer,
            max_gc_run=self.max_gc_run,
            max_repeat_length=self.max_repeat_length,
            forbidden_motifs=tuple(self.forbidden_motifs),
            forbidden_presets=tuple(self.forbidden_presets),
            avoid_reverse_complement=self.avoid_reverse_complement,
            restriction_enzymes=tuple(self.restriction_enzymes),
            ramp_weight=self.ramp_weight,
            ramp_codons=self.ramp_codons,
            cpg_weight=self.cpg_weight,
            cpg_mode=self.cpg_mode,
            minmax_weight=self.minmax_weight,
            minmax_direction=self.minmax_direction,
            tandem_unit=self.tandem_unit,
            tandem_copies=self.tandem_copies,
            inverted_stem=self.inverted_stem,
            inverted_loop=self.inverted_loop,
            avoid_internal_start=self.avoid_internal_start,
            avoid_uorf=self.avoid_uorf,
            uorf_region_nt=self.uorf_region_nt,
            refine=self.refine,
            refine_iterations=self.refine_iterations,
            folding_weight=self.folding_weight,
            gc_min=self.gc_min,
            gc_max=self.gc_max,
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
        """Compute the multi-objective Pareto frontier for a protein."""
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
