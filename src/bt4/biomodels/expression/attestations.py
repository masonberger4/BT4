"""Where an expression attestation is found, and the opt-in seam that honors it.

:mod:`bt4.biomodels.expression.attestation` defines *what* an attestation is and the one
function that can flip a learned expression head to ``calibrated=True``
(:func:`~bt4.biomodels.expression.attestation.verified_predictor`). This module is the
other half, mirroring :mod:`bt4.biomodels.splice.attestations`: *where BT4 looks for an
attestation* and *when it applies one*.

**Why this module had to exist.** Until it did, ``verified_predictor`` had **no caller
anywhere in** ``src/`` -- it was reachable only from its own tests. A maintainer could
run the acceptance gate, earn a claim and commit a record, and nothing would change for
a single user: the head stayed ``calibrated=False``, the candidate set stayed in
discovery order, and the only way to use the result was to write the promotion call by
hand. The gate was built, the record format was built, the calibrated-gating on the
consuming side (:mod:`bt4.pipeline.candidates`, :mod:`bt4.pipeline.rerank`) was built --
and the seam between them was missing. This is that seam.

**Promotion is opt-in, and stays opt-in.** An attestation records that a head passed a
gate *on one panel, in one scope*. It is not a statement that the user wants to run a
licensed non-commercial model, nor that this panel's scope is theirs. So an attestation
alone changes nothing: a caller must ask, via ``BT4_EXPRESSION_USE_ATTESTED=1`` or an
explicit ``enabled=True``. With the switch off -- the default --
:func:`bt4.biomodels.expression.default` still returns the neutral
:class:`~bt4.biomodels.expression.baseline.NullExpressionModel`, every shipped head
reports ``calibrated is False``, and the candidate set stays labelled *discovery order,
not a ranking*.

**Nothing is bundled today, and that is the honest state.** No
``data/ribonn.attestation.json`` ships, because no gate has been run -- and fabricating
one to light up the calibrated path is exactly what CLAUDE.md §10.6 forbids. So the
switch resolves nothing by default and every caller is told so rather than left with a
control that silently does nothing.

**A local attestation is first-class, not a workaround.** The splice side can bundle its
attestations because they are earned against *published* model weights that anyone with
the same files can re-verify. An expression attestation is earned against a
**maintainer's own measured panel**, which is frequently unpublished data -- and
committing one publishes the panel's hash and (see the UTR-hash note in
:mod:`~bt4.biomodels.expression.attestation`) its UTR context. So BT4 also reads a
maintainer-supplied path from ``$BT4_EXPRESSION_ATTESTATION``, letting someone who ran
the gate use the result on their own machine without committing anything. Resolution
order is explicit argument, then that path, then a bundled file.
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any, TypeVar, cast

from bt4.biomodels.expression.attestation import (
    ExpressionAttestation,
    ExpressionAttestationError,
    load_expression_attestation,
)

__all__ = [
    "ATTESTATION_PATH_ENV_VAR",
    "USE_ATTESTED_ENV_VAR",
    "attested_expression_backends",
    "attested_promotion_enabled",
    "bundled_expression_attestation",
    "bundled_expression_attestation_path",
    "promote_if_attested",
    "resolve_expression_attestation",
]

USE_ATTESTED_ENV_VAR = "BT4_EXPRESSION_USE_ATTESTED"
"""Set to a truthy value to honor a resolvable attestation (default: off)."""

ATTESTATION_PATH_ENV_VAR = "BT4_EXPRESSION_ATTESTATION"
"""Optional path to a maintainer's own attestation JSON, preferred over a bundled one."""

_TRUTHY = frozenset({"1", "true", "yes", "on"})

_DATA_DIR = Path(__file__).with_name("data")

_P = TypeVar("_P")
"""A promoted predictor is the *same* type as the one passed in.

``verified_predictor`` returns ``dataclasses.replace(predictor, ...)``, so the class
never changes -- only ``fidelity_verified``. Preserving that in the signature keeps
callers strictly typed rather than degrading to ``Any``.
"""


def bundled_expression_attestation_path(backend: str = "ribonn") -> Path:
    """Return the path a committed attestation for ``backend`` would occupy."""
    return _DATA_DIR / f"{backend}.attestation.json"


@lru_cache(maxsize=4)
def bundled_expression_attestation(backend: str = "ribonn") -> ExpressionAttestation | None:
    """Return the committed attestation for ``backend``, or ``None`` if none ships.

    ``None`` is the shipped state and the honest one: no expression head has passed its
    acceptance gate.

    Raises:
        ExpressionAttestationError: If a file *is* bundled but malformed -- refused
            loudly rather than treated as absent, because falling back silently would
            hide a packaging error.
    """
    path = bundled_expression_attestation_path(backend)
    if not path.is_file():
        return None
    return load_expression_attestation(path)


def resolve_expression_attestation(
    backend: str = "ribonn", *, attestation: ExpressionAttestation | None = None
) -> ExpressionAttestation | None:
    """Return the attestation BT4 would use for ``backend``, or ``None``.

    Resolution order: the explicit ``attestation`` argument, then
    ``$BT4_EXPRESSION_ATTESTATION`` (a maintainer's own record, which need never be
    committed), then a bundled file. Never raises for a *missing* attestation -- that is
    the normal state -- but does raise for one that is present and unusable, since a
    typo'd path silently behaving like "no attestation" is how a user ends up believing
    a promotion happened when it did not.

    Args:
        backend: Which head to resolve for (``"ribonn"``).
        attestation: An already-loaded record, which wins over both env and bundle.

    Returns:
        The attestation, or ``None`` when none resolves.

    Raises:
        ExpressionAttestationError: If the resolved record is malformed, or if
            ``$BT4_EXPRESSION_ATTESTATION`` names a file that does not exist.
    """
    if attestation is not None:
        return attestation
    raw = os.environ.get(ATTESTATION_PATH_ENV_VAR, "").strip()
    if raw:
        path = Path(raw).expanduser()
        if not path.is_file():
            raise ExpressionAttestationError(
                f"${ATTESTATION_PATH_ENV_VAR} points at {path}, which is not a readable "
                "file. Refusing to fall back silently: a typo here would look exactly "
                "like 'no attestation', and the head would stay uncalibrated for a "
                "reason nobody could see."
            )
        return load_expression_attestation(path)
    return bundled_expression_attestation(backend)


def attested_promotion_enabled() -> bool:
    """Return whether the environment opts in to honoring an expression attestation."""
    return os.environ.get(USE_ATTESTED_ENV_VAR, "").strip().lower() in _TRUTHY


def attested_expression_backends() -> tuple[str, ...]:
    """Return backends that can run here *and* have a resolvable attestation.

    What a frontend needs to decide whether offering the promotion opt-in is meaningful:
    an empty tuple means the switch would do nothing, and a control that silently does
    nothing is worse than one that says why it is unavailable (CLAUDE.md §6.6). Never
    raises -- a malformed or mis-pointed attestation means "not offerable here", which is
    reported by absence and explained by :func:`resolve_expression_attestation` when the
    user actually asks for it.
    """
    from bt4.biomodels.expression.ribonn import RiboNNExpressionModel

    try:
        runnable = RiboNNExpressionModel().available()
    except (OSError, ValueError, ImportError):
        runnable = False
    if not runnable:
        return ()
    try:
        resolved = resolve_expression_attestation("ribonn")
    except (ExpressionAttestationError, OSError):
        return ()
    return ("ribonn",) if resolved is not None else ()


def promote_if_attested(
    predictor: _P,
    *,
    enabled: bool | None = None,
    attestation: ExpressionAttestation | None = None,
) -> _P:
    """Promote ``predictor`` to ``calibrated=True`` when opted in and attested.

    The single place BT4 itself applies an expression attestation. Returns the predictor
    **unchanged** when the caller has not opted in, when no attestation resolves, or when
    the predictor is not an attestable head -- so it is safe to call on any
    :class:`~bt4.biomodels.expression.base.ExpressionPredictor`, including the neutral
    placeholder.

    It never *downgrades* silently: an attestation that resolves but does not match the
    predictor's scope -- species, cell-type selection, ``top_k``, UTR context, or the
    adapter's pinned weights -- raises
    :class:`~bt4.biomodels.expression.attestation.ExpressionAttestationError` from
    :func:`~bt4.biomodels.expression.attestation.verified_predictor` rather than quietly
    returning an uncalibrated head. A caller who asked for a calibrated ranking and got
    an uncalibrated one without being told is precisely the failure this layer exists to
    prevent.

    Args:
        predictor: Any expression predictor. Non-attestable heads are returned as-is.
        enabled: Override the opt-in. ``None`` (the default) consults
            :data:`USE_ATTESTED_ENV_VAR`, so a GUI can offer the choice per run without
            mutating the process environment by passing ``True``/``False``.
        attestation: Use this record instead of resolving one.

    Returns:
        The promoted predictor, or ``predictor`` unchanged.

    Raises:
        ExpressionAttestationError: If an attestation resolves and is opted into but does
            not match this predictor.
    """
    if not (attested_promotion_enabled() if enabled is None else enabled):
        return predictor
    if _attestable_backend(predictor) is None:
        return predictor
    resolved = resolve_expression_attestation("ribonn", attestation=attestation)
    if resolved is None:
        return predictor

    from bt4.biomodels.expression.attestation import verified_predictor

    return cast(_P, verified_predictor(predictor, resolved))


def _attestable_backend(predictor: Any) -> str | None:
    """Return the backend id for an attestable head, else ``None``.

    Mirrors ``attestation.verified_predictor``'s type check but *reports* rather than
    raises: this runs on every head a caller might hand over, where "not attestable"
    (the placeholder) is the normal case and not an error.
    """
    return "ribonn" if type(predictor).__name__ == "RiboNNExpressionModel" else None


# Re-exported so callers handling promotion failures need only this module.
_ = ExpressionAttestationError
