"""Bundled fidelity attestations, and the opt-in seam that honors them.

:mod:`bt4.biomodels.splice.attestation` defines *what* an attestation is and the
one function that can flip a wrapped CNN to ``calibrated=True``
(:func:`~bt4.biomodels.splice.attestation.verified_predictor`). This module is the
other half: *where the committed attestations live* and *when BT4 applies them*.

**Promotion is opt-in, and deliberately so.** A committed attestation records that
BT4's adapter reproduces the published model bit-for-bit -- a fact about the
wrapper. It is not a statement that the user wants to run a licensed,
noncommercial model, nor that the model's scores are calibrated *probabilities*
for coding-sequence design (they are not; see below). So the attestation alone
changes nothing: a caller must ask, via ``BT4_SPLICE_USE_ATTESTED=1`` or an
explicit ``enabled=True``. With the switch off -- the default --
:func:`bt4.biomodels.splice.default` still returns the honest PWM baseline and
every audit stays banner-led "UNCALIBRATED (advisory)".

**What a promotion does and does not assert.** It asserts *integration fidelity*:
the adapter is faithful to upstream. It does **not** assert that a
``P(splice)`` of 0.5 means a 50% chance of splicing in BT4's regime. Across eight
predictors including this one, median prAUC is **0.419 on exonic variants** versus
0.773 intronic (Smith & Kitzman, *Genome Biol* 2023) -- and BT4 designs coding
sequence, so its entire regime is the weaker half. A statistical-calibration gate
is separate and still unmet (Part B of
``docs/DESIGN_splice_cnn_calibration.md``).

**Regime scoping is enforced elsewhere, not here.** An attestation is captured on
the bare-CDS path, where the adapter pads its ~10 kb window with literal ``N``.
Scoring with real flanks is a different input regime the gate never exercised, so
:func:`~bt4.biomodels.splice.base.score_in_context` and
``pipeline.splice_audit._FlankedPredictor`` clear ``calibrated`` themselves. A
predictor promoted here can still report uncalibrated downstream, and that is
correct.
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any, TypeVar, cast

from bt4.biomodels.splice.attestation import (
    AttestationError,
    FidelityAttestation,
    load_attestation,
)

__all__ = [
    "USE_ATTESTED_ENV_VAR",
    "attested_promotion_enabled",
    "bundled_attestation",
    "bundled_attestation_path",
    "promote_if_attested",
]

USE_ATTESTED_ENV_VAR = "BT4_SPLICE_USE_ATTESTED"
"""Set to a truthy value to honor bundled attestations (default: off)."""

_TRUTHY = frozenset({"1", "true", "yes", "on"})

_DATA_DIR = Path(__file__).with_name("data")

_P = TypeVar("_P")
"""A promoted predictor is the *same* type as the one passed in.

``verified_predictor`` returns ``dataclasses.replace(predictor, ...)``, so the
class never changes -- only ``fidelity_verified``. Preserving that in the
signature keeps callers strictly typed rather than degrading to ``Any``.
"""


def bundled_attestation_path(backend: str) -> Path:
    """Return the path a committed attestation for ``backend`` would occupy."""
    return _DATA_DIR / f"{backend}.attestation.json"


@lru_cache(maxsize=4)
def bundled_attestation(backend: str) -> FidelityAttestation | None:
    """Return the committed attestation for ``backend``, or ``None`` if none ships.

    Args:
        backend: ``"pangolin"`` or ``"spliceai"``.

    Returns:
        The parsed attestation, or ``None`` when no file is bundled (which is the
        honest state for any backend whose gate has not been run).

    Raises:
        AttestationError: If a file *is* bundled but is malformed -- a corrupt
            attestation is refused loudly rather than treated as absent, because
            silently falling back would hide a packaging error.
    """
    path = bundled_attestation_path(backend)
    if not path.is_file():
        return None
    return load_attestation(path)


def attested_promotion_enabled() -> bool:
    """Return whether the environment opts in to honoring bundled attestations."""
    return os.environ.get(USE_ATTESTED_ENV_VAR, "").strip().lower() in _TRUTHY


def promote_if_attested(predictor: _P, *, enabled: bool | None = None) -> _P:
    """Promote ``predictor`` to ``calibrated=True`` when opted in and attested.

    The single place BT4 itself applies a committed attestation. Returns the
    predictor **unchanged** when the caller has not opted in, when no attestation
    is bundled for this backend, or when the predictor is not an attestable
    wrapped CNN -- so it is safe to call on any :class:`SplicePredictor`.

    It never *downgrades* silently: a bundled attestation that does not match the
    adapter's pinned weights, or a configuration the attestation does not cover,
    raises :class:`AttestationError` from
    :func:`~bt4.biomodels.splice.attestation.verified_predictor` rather than
    quietly returning an uncalibrated predictor. A mismatch is a packaging or
    configuration error worth surfacing.

    Args:
        predictor: Any splice predictor. Non-CNN backends are returned as-is.
        enabled: Override the opt-in. ``None`` (the default) consults
            :data:`USE_ATTESTED_ENV_VAR`.

    Returns:
        The promoted predictor, or ``predictor`` unchanged.

    Raises:
        AttestationError: If a bundled attestation exists and is opted into but
            does not match this predictor.
    """
    if not (attested_promotion_enabled() if enabled is None else enabled):
        return predictor

    backend = _attestable_backend(predictor)
    if backend is None:
        return predictor

    attestation = bundled_attestation(backend)
    if attestation is None:
        return predictor

    from bt4.biomodels.splice.attestation import verified_predictor

    return cast(_P, verified_predictor(predictor, attestation))


def _attestable_backend(predictor: Any) -> str | None:
    """Return the backend id for an attestable predictor, else ``None``.

    Mirrors ``attestation._backend_name`` but *reports* rather than raises: this
    function is called on every backend in a mixed list (baseline included), where
    "not attestable" is the normal case and not an error.
    """
    cls = type(predictor).__name__
    if cls == "PangolinSplicePredictor":
        return "pangolin"
    if cls == "SpliceAiSplicePredictor":
        return "spliceai"
    return None


# Re-exported so callers handling promotion failures need only this module.
_ = AttestationError
