"""License-clean fidelity attestations for the wrapped splice CNN backends.

A wrapped CNN backend (:class:`~bt4.biomodels.splice.PangolinSplicePredictor`,
:class:`~bt4.biomodels.splice.SpliceAiSplicePredictor`) may only report
``calibrated is True`` after its **integration-fidelity gate** passes -- proof the
adapter reproduces the published model's own per-position scores on a captured
reference panel (CLAUDE.md sections 6, 10.6). The problem this module solves is
*where a passing gate result may live*: the raw per-position panel scores are
license-encumbered (GPL-derived for Pangolin, CC BY-NC for SpliceAI) and must
never enter MIT-licensed BT4, yet the **facts a pass produces** -- a boolean, an
aggregate deviation, the tolerance, the pinned weight SHA-256s, the tool version
-- are not licensed outputs and are safe to record.

A :class:`FidelityAttestation` captures exactly those license-clean scalars and
nothing else. It is the record-of-truth for a passing gate, and it layers four
ways (the "layer them" decision):

* **(b) committed attestation** -- a maintainer who ran the gate against the
  licensed weights commits the attestation JSON (scalars + weight SHAs + version).
  It is re-verifiable by anyone who holds the same weights, and its
  :meth:`~FidelityAttestation.content_hash` can enter the provenance manifest so a
  calibrated splice audit stays reproducible-from-stamp.
* **(a) private execution** -- the gate itself runs where the weights live (a
  private CI/secret env); only the attestation leaves that venue.
* **(c) user opt-in** -- a user with their own install + panel runs
  ``verify_*_fidelity``, builds an attestation, and calls
  :func:`verified_predictor` to flip ``calibrated`` for their session.
* **(d) fallback** -- with no attestation the backends stay ``calibrated=False``
  and :func:`bt4.biomodels.splice.default` returns the honest PWM baseline.

**The load-bearing honesty rule:** an attestation stores *only* the
:class:`~bt4.biomodels.splice.FidelityReport` scalars plus the public weight
SHA-256s. It must **never** carry a ``FidelityCase``/``SpliceAiFidelityCase`` raw
``expected_*`` per-position array -- those *are* the licensed model outputs. This
is enforced structurally: :data:`_ALLOWED_FIELDS` pins the dataclass shape and a
test asserts no raw-score field is ever added.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any, Protocol

__all__ = [
    "MAX_ATTESTATION_TOLERANCE",
    "AttestationError",
    "FidelityAttestation",
    "attest_backend",
    "load_attestation",
    "verified_predictor",
]

# Backend identifiers the attestation layer recognizes.
_BACKENDS = ("pangolin", "spliceai")

# The loosest tolerance an attestation may claim and still be honored. The
# integration-fidelity gate proves *bit-for-bit* reproduction of the published
# model, so a legitimate pass deviates by ~1e-3 or less; a looser claim is not a
# real fidelity assertion and is refused (a floor against a bogus attestation).
MAX_ATTESTATION_TOLERANCE = 1e-3

# The exact, license-clean field set of a FidelityAttestation. Pinned here so the
# honesty test can assert the dataclass never grows a raw per-position-score field
# (which would leak the licensed model outputs into MIT BT4).
_ALLOWED_FIELDS = frozenset(
    {
        "backend",
        "passed",
        "max_abs_deviation",
        "n_cases",
        "tolerance",
        "weight_sha256",
        "bt4_version",
        "schema_version",
    }
)


class AttestationError(ValueError):
    """Raised when an attestation is missing, malformed, or does not match.

    A refusal, never a silent downgrade: a backend is only ever promoted to
    ``calibrated=True`` against an attestation that structurally matches the
    adapter's pinned weights and clears the tolerance floor.
    """


class _ScalarReport(Protocol):
    """The scalar fields shared by every ``*FidelityReport`` (duck-typed).

    Both :class:`~bt4.biomodels.splice.FidelityReport` and
    :class:`~bt4.biomodels.splice.SpliceAiFidelityReport` expose exactly these,
    and only these, so an attestation built from either is license-clean by
    construction (the raw per-position arrays live on the *case*, never the report).
    """

    passed: bool
    max_abs_deviation: float
    n_cases: int
    tolerance: float


@dataclass(frozen=True, slots=True)
class FidelityAttestation:
    """A license-clean record that a splice backend's fidelity gate passed.

    Holds only the derived scalars of a passing gate plus the public pinned weight
    SHA-256s -- never a raw per-position model score. Two attestations for the same
    weights and tool version are byte-identical (no timestamp, per the manifest's
    no-wall-clock rule, CLAUDE.md section 9 #7), so :meth:`content_hash` is a
    stable provenance stamp.

    Attributes:
        backend: Which wrapped backend this attests (``"pangolin"`` or
            ``"spliceai"``).
        passed: Whether the gate passed (always ``True`` for a stored attestation;
            :func:`attest_backend` refuses to record a failing gate).
        max_abs_deviation: Largest absolute per-position deviation the gate saw.
        n_cases: How many reference cases were checked.
        tolerance: The absolute tolerance the gate used.
        weight_sha256: The pinned ``{weight_filename: sha256}`` map the gate ran
            against, as a sorted tuple of pairs -- these are public content hashes
            of the weight files, not licensed outputs.
        bt4_version: The BT4 version that produced the attestation (provenance).
        schema_version: Attestation schema version (for forward compatibility).
    """

    backend: str
    passed: bool
    max_abs_deviation: float
    n_cases: int
    tolerance: float
    weight_sha256: tuple[tuple[str, str], ...]
    bt4_version: str
    schema_version: int = 1

    def to_dict(self) -> dict[str, Any]:
        """Return the canonical, JSON-serializable form (license-clean scalars)."""
        return {
            "backend": self.backend,
            "passed": self.passed,
            "max_abs_deviation": self.max_abs_deviation,
            "n_cases": self.n_cases,
            "tolerance": self.tolerance,
            "weight_sha256": {name: sha for name, sha in self.weight_sha256},
            "bt4_version": self.bt4_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> FidelityAttestation:
        """Rebuild an attestation from :meth:`to_dict` output.

        Raises:
            AttestationError: If a required field is missing or the backend is
                unknown.
        """
        extra = set(data) - _ALLOWED_FIELDS
        if extra:
            # An unexpected field could smuggle in a raw-score array; refuse it.
            raise AttestationError(f"unexpected attestation field(s): {sorted(extra)}")
        try:
            weights = data["weight_sha256"]
            att = cls(
                backend=data["backend"],
                passed=data["passed"],
                max_abs_deviation=float(data["max_abs_deviation"]),
                n_cases=int(data["n_cases"]),
                tolerance=float(data["tolerance"]),
                weight_sha256=tuple(sorted((str(k), str(v)) for k, v in weights.items())),
                bt4_version=str(data["bt4_version"]),
                schema_version=int(data.get("schema_version", 1)),
            )
        except (KeyError, TypeError, AttributeError) as exc:
            raise AttestationError(f"malformed attestation: {exc}") from exc
        if att.backend not in _BACKENDS:
            raise AttestationError(
                f"unknown backend {att.backend!r} (known: {sorted(_BACKENDS)})"
            )
        return att

    def content_hash(self) -> str:
        """Return a deterministic SHA-256 over the canonical form.

        Timestamp-free and key-sorted, so the same passing gate always yields the
        same hash -- suitable as a provenance-manifest stamp for a calibrated
        splice audit (CLAUDE.md section 9 #9).
        """
        canonical = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def attest_backend(
    backend: str,
    report: _ScalarReport,
    weight_sha256: dict[str, str],
    *,
    bt4_version: str,
) -> FidelityAttestation:
    """Build an attestation from a **passing** fidelity report.

    Args:
        backend: ``"pangolin"`` or ``"spliceai"``.
        report: A :class:`~bt4.biomodels.splice.FidelityReport` /
            :class:`~bt4.biomodels.splice.SpliceAiFidelityReport` -- only its
            license-clean scalar fields are read.
        weight_sha256: The ``{weight_filename: sha256}`` the gate ran against
            (the adapter's :data:`PINNED_WEIGHT_SHA256`).
        bt4_version: The producing BT4 version.

    Returns:
        A :class:`FidelityAttestation`.

    Raises:
        AttestationError: If ``backend`` is unknown, the report did **not** pass
            (a failing gate is never recorded as an attestation), or the report's
            tolerance is looser than :data:`MAX_ATTESTATION_TOLERANCE`.
    """
    if backend not in _BACKENDS:
        raise AttestationError(
            f"unknown backend {backend!r} (known: {sorted(_BACKENDS)})"
        )
    if not report.passed:
        raise AttestationError("refusing to attest a failing fidelity gate")
    if report.tolerance > MAX_ATTESTATION_TOLERANCE:
        raise AttestationError(
            f"gate tolerance {report.tolerance} looser than the attestation floor "
            f"{MAX_ATTESTATION_TOLERANCE}; not a bit-for-bit fidelity claim"
        )
    return FidelityAttestation(
        backend=backend,
        passed=True,
        max_abs_deviation=float(report.max_abs_deviation),
        n_cases=int(report.n_cases),
        tolerance=float(report.tolerance),
        weight_sha256=tuple(sorted((str(k), str(v)) for k, v in weight_sha256.items())),
        bt4_version=bt4_version,
        schema_version=1,
    )


def load_attestation(path: str | Path) -> FidelityAttestation:
    """Load a committed / user-supplied attestation JSON.

    Raises:
        AttestationError: If the file is missing, unreadable, or malformed.
    """
    try:
        raw = Path(path).read_text(encoding="utf-8")
        data = json.loads(raw)
    except (OSError, json.JSONDecodeError) as exc:
        raise AttestationError(f"cannot load attestation {path!r}: {exc}") from exc
    if not isinstance(data, dict):
        raise AttestationError(f"attestation {path!r} is not a JSON object")
    return FidelityAttestation.from_dict(data)


def _pinned_weights(backend: str) -> dict[str, str]:
    """Return the adapter's pinned weight SHA-256 map (lazy import, keeps core light)."""
    if backend == "pangolin":
        from bt4.biomodels.splice.pangolin import PINNED_WEIGHT_SHA256

        return dict(PINNED_WEIGHT_SHA256)
    if backend == "spliceai":
        from bt4.biomodels.splice.spliceai import PINNED_WEIGHT_SHA256

        return dict(PINNED_WEIGHT_SHA256)
    raise AttestationError(f"unknown backend {backend!r} (known: {sorted(_BACKENDS)})")


def _backend_name(predictor: Any) -> str:
    """Map a predictor instance to its attestation backend id."""
    cls = type(predictor).__name__
    if cls == "PangolinSplicePredictor":
        return "pangolin"
    if cls == "SpliceAiSplicePredictor":
        return "spliceai"
    raise AttestationError(
        f"{cls} is not an attestable wrapped CNN backend "
        "(only Pangolin / SpliceAI can be calibrated via a fidelity attestation)"
    )


def _check_full_weight_coverage(predictor: Any, backend: str) -> None:
    """Refuse to promote a configuration that does not load every pinned weight.

    An honored attestation is *required* to claim the adapter's full pinned map --
    a subset fails the equality check above and is refused. But Pangolin's tissue
    set selects both which weight files load and which output channel is read
    (``TISSUE_OUTPUTS`` maps each tissue to a ``weight_index`` / ``output_channel``),
    so a gate run at ``tissues=("heart",)`` touches 3 of the 12 pinned files and
    one of the four channels, while still recording the full 12-file map. Promoting
    the default four-tissue predictor against that attestation would claim
    calibration for nine weight files and three channels the gate never executed.

    ``verify_pangolin_fidelity``'s docstring already states the requirement ("its
    tissue set must match how the panel was captured"); this makes it enforced
    rather than advisory, and it lives at the promotion seam so it holds for every
    caller. SpliceAI has no equivalent axis -- its adapter always loads all five
    ensemble members -- so only Pangolin is checked.

    Raises:
        AttestationError: If the predictor's configuration would load only part of
            the pinned weight set.
    """
    if backend != "pangolin":
        return
    from bt4.biomodels.splice.pangolin import DEFAULT_TISSUES

    tissues = tuple(getattr(predictor, "tissues", DEFAULT_TISSUES))
    if set(tissues) != set(DEFAULT_TISSUES):
        raise AttestationError(
            f"attestation covers all {len(DEFAULT_TISSUES)} Pangolin tissue heads "
            f"({', '.join(DEFAULT_TISSUES)}), but this predictor is configured for "
            f"{', '.join(tissues)}; a partial configuration loads only part of the "
            "attested weight set, so the attestation does not cover it"
        )


def verified_predictor(predictor: Any, attestation: FidelityAttestation) -> Any:
    """Return ``predictor`` promoted to ``calibrated=True`` iff ``attestation`` matches.

    This is the single seam that flips a wrapped CNN backend to calibrated. It
    refuses -- never silently downgrades -- unless **all** of the following hold:

    * the attestation ``passed``;
    * its ``backend`` matches the predictor's type;
    * its ``tolerance`` is no looser than :data:`MAX_ATTESTATION_TOLERANCE`;
    * its ``weight_sha256`` exactly matches the adapter's own
      :data:`PINNED_WEIGHT_SHA256` (so the attestation is bound to the same weights
      the adapter will hash-verify before loading -- a mismatch means the
      attestation is for different weights and is not trustworthy here).

    On success it returns ``dataclasses.replace(predictor, fidelity_verified=True)``
    so the backend reports ``calibrated is True``; the scores themselves are still
    computed live from the user's own hash-pinned weights, never from the
    attestation (which carries no scores).

    Args:
        predictor: A :class:`~bt4.biomodels.splice.PangolinSplicePredictor` or
            :class:`~bt4.biomodels.splice.SpliceAiSplicePredictor` instance.
        attestation: A passing :class:`FidelityAttestation` for that backend.

    Returns:
        The promoted predictor (``calibrated is True``).

    Raises:
        AttestationError: On any mismatch or a non-attestable predictor type.
    """
    import dataclasses

    backend = _backend_name(predictor)
    if attestation.backend != backend:
        raise AttestationError(
            f"attestation is for {attestation.backend!r}, predictor is {backend!r}"
        )
    if not attestation.passed:
        raise AttestationError("attestation did not pass; refusing to calibrate")
    if attestation.tolerance > MAX_ATTESTATION_TOLERANCE:
        raise AttestationError(
            f"attestation tolerance {attestation.tolerance} looser than floor "
            f"{MAX_ATTESTATION_TOLERANCE}"
        )
    pinned = tuple(sorted((str(k), str(v)) for k, v in _pinned_weights(backend).items()))
    if attestation.weight_sha256 != pinned:
        raise AttestationError(
            "attestation weight SHA-256 set does not match the adapter's pinned "
            "weights; the attestation is for different weights"
        )
    _check_full_weight_coverage(predictor, backend)
    return dataclasses.replace(predictor, fidelity_verified=True)


# The dataclass shape is honesty-load-bearing: assert at import that it carries
# exactly the license-clean fields and nothing that could hold a raw model score.
assert {f.name for f in fields(FidelityAttestation)} == _ALLOWED_FIELDS, (
    "FidelityAttestation fields drifted from the license-clean allowed set; a new "
    "field must never carry raw per-position model scores (licensed outputs)."
)
