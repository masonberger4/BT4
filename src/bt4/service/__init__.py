"""BT4 service layer: an optional headless HTTP API over :mod:`bt4.api`.

This layer is gated behind the ``bt4[service]`` extra (FastAPI, uvicorn). It
imports only the public ``bt4.api`` surface, so it can be scripted over HTTP
without depending on the optimizer, pipeline, or biomodels. Import
:data:`app` (or call :func:`create_app`) to serve it, e.g.
``uvicorn bt4.service.api:app``.
"""

from __future__ import annotations

from bt4.service.api import app, create_app

__all__ = ["app", "create_app"]
