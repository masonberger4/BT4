"""SpliceAI-backed splice predictor -- the second wrapped published CNN.

:class:`SpliceAiSplicePredictor` is BT4's second *wrapped published* splice
backend (CLAUDE.md section 6: "wrap published SpliceAI + Pangolin -- no
self-training"), the cross-check to :class:`~bt4.biomodels.splice.pangolin.PangolinSplicePredictor`.
It runs the already-validated **SpliceAI** deep-learning splice model
(Jaganathan et al., *Cell* 2019, doi:10.1016/j.cell.2018.12.015) as an
*inference-only* backend behind the existing
:class:`~bt4.biomodels.splice.base.SplicePredictor` contract, feeding SpliceAI's
per-nucleotide acceptor / donor probabilities into BT4's already-shipped
Delta-splicing framing and top-k / log-odds pooling
(:mod:`bt4.biomodels.splice.base`). Running both SpliceAI and Pangolin and
reporting their agreement (:func:`bt4.biomodels.splice.backend_agreement`) is the
first-class uncertainty signal of CLAUDE.md sections 6 and 8.

The honesty rules mirror the Pangolin adapter, with three SpliceAI specifics:

* **SpliceAI is *more* restrictively licensed than Pangolin -- NOT MIT, NOT even
  GPL.** Its **source code is PolyForm Strict License 1.0.0** (noncommercial;
  forbids redistribution and derivative works) and its **model weights are
  CC BY-NC 4.0** (noncommercial; require attribution). (The stray ``GPLv3`` in
  the upstream ``setup.py`` is contradicted by the authoritative ``LICENSE`` /
  ``spliceai/models/LICENSE`` files and is not authoritative.) Both are
  incompatible with bundling into MIT-licensed BT4, so -- exactly as BT4 wraps
  GPL ViennaRNA and GPL Pangolin -- this adapter **lazily imports the user's own
  installed ``spliceai`` package and its weights and bundles neither code nor
  weights**; BT4 ships only the SHA-256 content hashes (facts), the wrapping
  glue, and attribution. **Honest scope note:** because the weights are
  noncommercial (CC BY-NC), this optional backend is for academic / noncommercial
  use only -- a caveat that does not touch BT4's MIT core (the backend is opt-in,
  installed by the user).

* **``calibrated`` starts ``False`` and is gated on an integration-fidelity
  check** (:func:`verify_spliceai_fidelity`) -- the adapter must be shown to
  reproduce upstream SpliceAI's own scores on a captured reference panel, *not* a
  from-scratch training gate. That gate has since PASSED on a maintainer machine
  holding the weights, and its
  license-clean attestation ships (the reference *panel* still does not -- it is the
  licensed model's own output). So the shipped adapter is ``calibrated is False``
  **by default**, and reports ``True`` only under the explicit
  ``BT4_SPLICE_USE_ATTESTED=1`` opt-in; :func:`bt4.biomodels.splice.default` keeps
  returning the labeled PWM baseline either way. Even promoted, the flag asserts
  *integration fidelity* -- the adapter reproduces the published model -- never that
  the scores are calibrated probabilities for designed coding sequence.

* **Weights are hash-pinned and verified before they are loaded.** Every
  ``spliceai{1..5}.h5`` file is checked against :data:`PINNED_WEIGHT_SHA256`
  *before* Keras ``load_model`` touches it, so the adapter stays
  reproducible-from-manifest (unlike BT3's live ASSP scrape, CLAUDE.md section
  10.15) and never loads bytes it has not content-verified.

**Clean donor / acceptor split.** Unlike Pangolin (one combined ``P(splice)``),
SpliceAI outputs a 3-way per-position softmax -- ``[null, acceptor, donor]`` --
so :attr:`SpliceResult.acceptor` and :attr:`SpliceResult.donor` are *both*
populated from channels 1 and 2. :func:`~bt4.biomodels.splice.base.pooled_risk`
pools ``donor`` union ``acceptor``, counting a donor site and an acceptor site as
distinct real sites (correct, no double-count).

The heavy dependency (TensorFlow / Keras) is imported **only inside methods**,
never at module load, so importing this module -- and ``import bt4`` -- stays
lightweight (CLAUDE.md section 3). This module depends only on :mod:`bt4.domain`,
the standard library, and -- lazily, inside methods -- the optional
``tensorflow`` / ``keras`` and ``spliceai`` packages.
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SpliceResult,
    pooled_risk,
)
from bt4.domain.sequence import validate_dna

__all__ = [
    "ACCEPTOR_CHANNEL",
    "DONOR_CHANNEL",
    "PINNED_WEIGHT_SHA256",
    "SPLICEAI_CONTEXT",
    "SPLICEAI_ENSEMBLE_SIZE",
    "SPLICEAI_FLANK",
    "SpliceAiFidelityCase",
    "SpliceAiFidelityReport",
    "SpliceAiSplicePredictor",
    "verify_spliceai_fidelity",
]


SPLICEAI_CONTEXT: int = 10_000
"""Total receptive-field context SpliceAI-10k crops from its input.

SpliceAI outputs per-position scores for the middle ``N - 10000`` bases of an
``N``-base input. To score every position of a coding sequence, the adapter pads
it with :data:`SPLICEAI_FLANK` ``N`` bases on each side (SpliceAI's documented
convention for custom sequences).
"""

SPLICEAI_FLANK: int = SPLICEAI_CONTEXT // 2
"""Number of ``N`` bases padded on each side so output aligns 1:1 with the CDS."""

SPLICEAI_ENSEMBLE_SIZE: int = 5
"""SpliceAI's published ensemble: ``spliceai1.h5`` .. ``spliceai5.h5``, averaged."""

ACCEPTOR_CHANNEL: int = 1
"""Output channel holding the per-position acceptor probability (``y[:, 1]``)."""

DONOR_CHANNEL: int = 2
"""Output channel holding the per-position donor probability (``y[:, 2]``)."""

# SHA-256 of SpliceAI's published weight files (github.com/Illumina/SpliceAI,
# spliceai/models/spliceai{1..5}.h5). Weights are CC BY-NC 4.0 artifacts kept OUT
# of this repository; only these content hashes ship, and each file is verified
# against its pin *before* it is loaded (CLAUDE.md sections 6, 10.15).
PINNED_WEIGHT_SHA256: dict[str, str] = {
    "spliceai1.h5": "e1fd5adcef7489d604b10e79c40078ef790d51ef048c4ce3869c9119ac5de42b",
    "spliceai2.h5": "6ab042b82ab966b6d3582cb31b96f0859ea08a864f168d69e83aa14450a3b66e",
    "spliceai3.h5": "e2e790bde53dfdf410c6dc434a86122a7d12f3f38dc2ef45d85986e9ecf22fad",
    "spliceai4.h5": "ca88ac9e58e69ba6fdeed319b72f063f164c9abf7392eaccef903e94c1d99dd6",
    "spliceai5.h5": "791cd22c62a80a08d2ca674615a93ce8159d7b55bd157cfef2983b1bd6b41391",
}

_WEIGHT_FILE_TEMPLATE: str = "spliceai{index}.h5"
"""Upstream weight-file name pattern (``index`` in 1..5)."""

_WEIGHTS_ENV_VAR: str = "BT4_SPLICEAI_MODEL_DIR"
"""Optional env var overriding where the pinned SpliceAI weight files are found."""

# One-hot channel order A, C, G, T -- matching SpliceAI's map so the pinned
# weights receive the input they were trained on. 'N' (padding) maps to all-zero.
# SpliceAI uses a POSITION-MAJOR [L][4] layout (Keras channels-last), unlike
# Pangolin's channel-major [4][L].
_BASE_ROW: dict[str, tuple[float, float, float, float]] = {
    "A": (1.0, 0.0, 0.0, 0.0),
    "C": (0.0, 1.0, 0.0, 0.0),
    "G": (0.0, 0.0, 1.0, 0.0),
    "T": (0.0, 0.0, 0.0, 1.0),
    "N": (0.0, 0.0, 0.0, 0.0),
}


def _one_hot_rows(seq: str) -> list[list[float]]:
    """One-hot encode ``seq`` to position-major ``[len][4]`` (A, C, G, T order).

    Pure standard-library encoding (no numpy, no TensorFlow), so it is
    unit-testable without the optional deps. Matches SpliceAI's ``one_hot_encode``
    base ordering and channels-last layout; ``N`` maps to an all-zero row.

    Args:
        seq: An upper-case string over ``{A,C,G,T,N}``.

    Returns:
        ``len(seq)`` rows, each a 4-element ``[A,C,G,T]`` indicator.
    """
    return [list(_BASE_ROW[base]) for base in seq]


def _split_acceptor_donor(
    per_position: Sequence[Sequence[float]],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Split SpliceAI's per-position 3-way output into (acceptor, donor) tracks.

    Args:
        per_position: One ``[null, acceptor, donor]`` triple per nucleotide.

    Returns:
        ``(acceptor, donor)`` where ``acceptor[i] = row[ACCEPTOR_CHANNEL]`` and
        ``donor[i] = row[DONOR_CHANNEL]``. Pure; testable without TensorFlow.
    """
    acceptor = tuple(float(row[ACCEPTOR_CHANNEL]) for row in per_position)
    donor = tuple(float(row[DONOR_CHANNEL]) for row in per_position)
    return acceptor, donor


def _sha256_file(path: Path) -> str:
    """Return the lowercase hex SHA-256 of the file at ``path`` (streamed)."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_weight_file(path: Path) -> None:
    """Verify ``path``'s bytes match the pinned SHA-256 for its file name.

    Args:
        path: A SpliceAI ``spliceai{n}.h5`` weight file.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        KeyError: If the file name is not in :data:`PINNED_WEIGHT_SHA256`.
        ValueError: If the file's SHA-256 does not match its pin -- the weights
            are unknown or altered and are refused, not loaded.
    """
    if not path.is_file():
        raise FileNotFoundError(f"SpliceAI weight file not found: {path}")
    expected = PINNED_WEIGHT_SHA256.get(path.name)
    if expected is None:
        raise KeyError(
            f"SpliceAI weight file '{path.name}' is not a pinned weight "
            f"(known: {sorted(PINNED_WEIGHT_SHA256)})"
        )
    actual = _sha256_file(path)
    if actual != expected:
        raise ValueError(
            f"SpliceAI weight file '{path.name}' failed its SHA-256 pin "
            f"(expected {expected}, got {actual}); refusing to load unverified weights"
        )


def _ambient_keras_is_v3() -> bool:
    """Whether the Keras that ``import tensorflow.keras`` would yield is 3.x.

    Reported from the installed packages rather than the TensorFlow version, because a
    user can pin either side: ``TF_USE_LEGACY_KERAS`` and an explicit ``keras<3`` both
    make a modern TensorFlow serve Keras 2. Any failure to determine it answers
    ``False``, which keeps the historical order for every environment this cannot read.
    """
    for module_name in ("tensorflow.keras", "keras"):
        try:
            module = __import__(module_name, fromlist=["__version__"])
        except ImportError:
            continue
        version = getattr(module, "__version__", "")
        if isinstance(version, str) and version:
            return version.split(".")[0] == "3"
    return False


def _import_keras() -> Any:
    """Import and return a Keras ``load_model`` callable (lazy, guarded).

    SpliceAI's weights are 2019 **Keras 2** ``.h5`` graphs, and Keras 3 cannot load
    them. From TensorFlow 2.16 ``tensorflow.keras`` *is* Keras 3, so preference order
    is decided by the installed Keras **generation**, not by module availability:

    1. ``tf_keras`` -- the legacy shim -- whenever the ambient Keras is 3.x. It exists
       for exactly this case, and is what ``TF_USE_LEGACY_KERAS=1`` selects.
    2. ``tensorflow.keras`` / ``keras`` otherwise (TF <= 2.15, where they are Keras 2).
    3. ``tf_keras`` as a last resort.

    Ordering by availability alone -- trying ``tensorflow.keras`` first and falling
    back on ``ImportError`` -- makes the shim **unreachable**, because under TF >= 2.16
    ``tensorflow.keras`` imports perfectly well; it is simply the wrong Keras. The
    failure then surfaces at ``load_model`` as an opaque deserialization error about a
    file that is not corrupt.

    Returns:
        The resolved ``load_model`` function.

    Raises:
        ModuleNotFoundError: If none of the Keras entry points import. Install the
            ``bt4[splice-spliceai]`` extra (TensorFlow) plus the ``spliceai``
            package (see the module docstring), or use the baseline predictor.
    """
    order = ("tensorflow.keras.models", "keras.models", "tf_keras.models")
    if _ambient_keras_is_v3():
        order = ("tf_keras.models", "tensorflow.keras.models", "keras.models")
    for module_name in order:
        try:
            module = __import__(module_name, fromlist=["load_model"])
        except ImportError:
            continue
        return module.load_model
    raise ModuleNotFoundError(
        "TensorFlow / Keras is not installed; install the 'bt4[splice-spliceai]' "
        "extra and the CC BY-NC 'spliceai' package, or use the baseline splice predictor"
    )


def _resolve_model_dir(explicit: str | None = None) -> Path | None:
    """Resolve the directory holding SpliceAI's weight files, or ``None``.

    Resolution order: an explicit path, then the :data:`_WEIGHTS_ENV_VAR`
    environment variable, then the ``models`` directory bundled with the installed
    ``spliceai`` package. Uses :func:`importlib.util.find_spec` rather than
    importing ``spliceai`` (whose ``__init__`` installs a SIGINT handler and pulls
    heavy deps). Never raises: returns ``None`` when nothing resolves.

    Args:
        explicit: An explicit weights directory, or ``None``.

    Returns:
        The resolved directory, or ``None`` if none could be found.
    """
    for candidate in (explicit, os.environ.get(_WEIGHTS_ENV_VAR)):
        if candidate:
            path = Path(candidate)
            return path if path.is_dir() else None
    try:
        spec = importlib.util.find_spec("spliceai")
    except (ImportError, ValueError):
        return None
    if spec is None or spec.origin is None:
        return None
    models = Path(spec.origin).resolve().parent / "models"
    return models if models.is_dir() else None


@lru_cache(maxsize=2)
def _load_ensemble(model_dir: str) -> list[Any]:
    """Load and cache the 5-model SpliceAI ensemble from ``model_dir``.

    Every weight file's SHA-256 pin is verified *before* Keras loads it. Cached
    per ``model_dir`` so the expensive load happens once per process.

    Args:
        model_dir: Directory holding ``spliceai1.h5`` .. ``spliceai5.h5``.

    Returns:
        The five loaded Keras models, in order.

    Raises:
        ModuleNotFoundError: If TensorFlow / Keras is not importable.
        FileNotFoundError / KeyError / ValueError: From :func:`_verify_weight_file`.
    """
    load_model = _import_keras()
    base = Path(model_dir)
    models: list[Any] = []
    for index in range(1, SPLICEAI_ENSEMBLE_SIZE + 1):
        path = base / _WEIGHT_FILE_TEMPLATE.format(index=index)
        _verify_weight_file(path)
        models.append(load_model(str(path)))
    return models


@dataclass(frozen=True, slots=True)
class SpliceAiSplicePredictor:
    """Wrapped published SpliceAI model behind the ``SplicePredictor`` contract.

    Runs the installed CC BY-NC ``spliceai`` package's 5-model Keras ensemble as
    an inference-only backend, reporting per-position acceptor and donor
    probabilities. See the module docstring for the licensing (no-bundle,
    noncommercial), the hash-pinning, and why ``calibrated`` starts ``False``.

    Attributes:
        top_k: Number of strongest sites summed by :meth:`delta_splicing`'s
            top-k / log-odds pooling. Defaults to
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`.
        model_dir: Explicit weights directory, or ``None`` to resolve from the
            :data:`_WEIGHTS_ENV_VAR` env var or the installed ``spliceai`` package.
        fidelity_verified: Whether this instance has passed
            :func:`verify_spliceai_fidelity` against a captured upstream panel.
            Mirrored by :attr:`calibrated`; ``False`` by default. A committed
            attestation ships and can flip it, but only under the explicit
            ``BT4_SPLICE_USE_ATTESTED=1`` opt-in.
    """

    top_k: int = DEFAULT_TOP_K
    model_dir: str | None = None
    fidelity_verified: bool = field(default=False)

    def __post_init__(self) -> None:
        """Validate the pooling depth.

        Raises:
            ValueError: If ``top_k`` is not a positive integer.
        """
        if self.top_k <= 0:
            raise ValueError(f"top_k must be a positive integer, got {self.top_k}")

    @property
    def name(self) -> str:
        """Backend identifier."""
        return "spliceai"

    @property
    def calibrated(self) -> bool:
        """``True`` only once the integration-fidelity gate has passed.

        Mirrors :attr:`fidelity_verified`. Stays ``False`` until
        :func:`verify_spliceai_fidelity` confirms the adapter reproduces upstream
        SpliceAI's own scores -- which it has, on a maintainer machine holding the
        weights. Shipped instances are uncalibrated by default and report ``True``
        only under the ``BT4_SPLICE_USE_ATTESTED=1`` opt-in, where it means the
        adapter is faithful, not that the scores are calibrated for designed coding
        sequence (CLAUDE.md sections 6 and 10.6).
        """
        return self.fidelity_verified

    def weights_dir(self) -> Path | None:
        """Return the weights directory this predictor would load from, or ``None``.

        Resolution order: this predictor's ``model_dir``, then the
        ``BT4_SPLICEAI_MODEL_DIR`` environment variable, then the ``models``
        directory inside the installed ``spliceai`` package. Never raises and
        imports no heavy dependency (the package is located via
        :func:`importlib.util.find_spec`, so ``spliceai``'s import-time SIGINT
        handler is never installed).

        Public because "which files is BT4 actually about to hash?" is a question
        users and maintenance scripts legitimately need answered -- notably
        ``scripts/check_splice_weights.py`` (step A3 of
        ``docs/DESIGN_splice_cnn_calibration.md``), which must not reach across a
        layer for a private resolver (CLAUDE.md section 10.9).

        Returns:
            The resolved directory, or ``None`` when nothing resolves.
        """
        return _resolve_model_dir(self.model_dir)

    def available(self) -> bool:
        """Return whether this backend can run (never raises).

        ``True`` iff Keras / TensorFlow imports and a weights directory containing
        all five ``spliceai{1..5}.h5`` files resolves. Does not hash-verify the
        files (that happens at load); a lightweight existence check keeps
        ``available`` cheap and non-raising, mirroring the Pangolin adapter.
        """
        try:
            _import_keras()
        except ImportError:
            return False
        model_dir = _resolve_model_dir(self.model_dir)
        if model_dir is None:
            return False
        for index in range(1, SPLICEAI_ENSEMBLE_SIZE + 1):
            if not (model_dir / _WEIGHT_FILE_TEMPLATE.format(index=index)).is_file():
                return False
        return True

    def _acceptor_donor(self, dna: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return SpliceAI's per-position ``(acceptor, donor)`` tracks for ``dna``.

        Raises:
            ModuleNotFoundError: If TensorFlow / Keras is not importable.
            RuntimeError: If no weights directory resolves, or the model's output
                length does not align to the input (an architecture-mismatch guard).
            FileNotFoundError / KeyError / ValueError: From weight verification.
        """
        seq = validate_dna(dna)
        model_dir = _resolve_model_dir(self.model_dir)
        if model_dir is None:
            raise RuntimeError(
                "no SpliceAI weights directory resolved; set "
                f"${_WEIGHTS_ENV_VAR}, pass model_dir=, or install the 'spliceai' package"
            )
        import numpy as np  # lazy: TF ships numpy

        ensemble = _load_ensemble(str(model_dir))
        padded = "N" * SPLICEAI_FLANK + seq + "N" * SPLICEAI_FLANK
        x = np.asarray([_one_hot_rows(padded)], dtype="float32")
        y = np.mean([model.predict(x, verbose=0) for model in ensemble], axis=0)[0]

        acceptor, donor = _split_acceptor_donor(y.tolist())
        if len(acceptor) != len(seq):
            raise RuntimeError(
                f"SpliceAI output length {len(acceptor)} != input length {len(seq)}; "
                f"the model's context crop may differ from the expected {SPLICEAI_CONTEXT}"
            )
        return acceptor, donor

    def site_scores(self, dna: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return SpliceAI's per-position ``(acceptor, donor)`` probabilities.

        Both tracks are the ensemble-averaged SpliceAI probabilities at each
        nucleotide of ``dna``, aligned 1:1.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            An ``(acceptor, donor)`` pair, each one probability per input
            nucleotide.

        Raises:
            ModuleNotFoundError / RuntimeError / FileNotFoundError / KeyError /
            ValueError: As documented on :meth:`_acceptor_donor`.
        """
        return self._acceptor_donor(dna)

    def score_sequence(self, dna: str) -> SpliceResult:
        """Return per-position acceptor / donor site scores as a :class:`SpliceResult`.

        SpliceAI's 3-way softmax maps cleanly onto the contract: channel 1 ->
        ``acceptor``, channel 2 -> ``donor``, both populated (unlike Pangolin's
        single combined track). :func:`pooled_risk` pools ``donor`` union
        ``acceptor``, counting a donor and an acceptor site as distinct sites.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            A :class:`SpliceResult`; scores are calibrated probabilities only when
            :attr:`calibrated` is ``True`` (i.e. after the fidelity gate passes).

        Raises:
            ModuleNotFoundError / RuntimeError / FileNotFoundError / KeyError /
            ValueError: As documented on :meth:`_acceptor_donor`.
        """
        acceptor, donor = self._acceptor_donor(dna)
        return SpliceResult(
            donor=donor,
            acceptor=acceptor,
            model_name=self.name,
            calibrated=self.calibrated,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Return the negated added splice risk of ``designed`` vs ``reference``.

        See :meth:`bt4.biomodels.splice.base.SplicePredictor.delta_splicing` for
        the fixed *larger-is-better* orientation. Concretely returns
        ``pooled_risk(reference) - pooled_risk(designed)`` using top-k / log-odds
        pooling over the combined acceptor+donor sites, so it is ``0.0`` for
        identical sequences, positive when the redesign lowers predicted splice
        risk, and negative when it raises it.

        Raises:
            ModuleNotFoundError / RuntimeError / FileNotFoundError / KeyError /
            ValueError: As documented on :meth:`_acceptor_donor`.
        """
        designed_risk = pooled_risk(self.score_sequence(designed_dna), self.top_k)
        reference_risk = pooled_risk(self.score_sequence(reference_dna), self.top_k)
        return reference_risk - designed_risk


@dataclass(frozen=True, slots=True)
class SpliceAiFidelityCase:
    """One captured-from-upstream reference case for the integration-fidelity gate.

    A maintainer with the CC BY-NC SpliceAI install and weights runs upstream
    SpliceAI's custom-sequence scoring on ``sequence`` and records the resulting
    per-position acceptor and donor probabilities here.
    :func:`verify_spliceai_fidelity` then asserts the adapter reproduces them.

    Attributes:
        sequence: The coding sequence scored.
        expected_acceptor: Upstream SpliceAI's per-position acceptor probability.
        expected_donor: Upstream SpliceAI's per-position donor probability.
    """

    sequence: str
    expected_acceptor: tuple[float, ...]
    expected_donor: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class SpliceAiFidelityReport:
    """Result of a SpliceAI integration-fidelity check.

    Attributes:
        passed: Whether every case's both tracks matched within tolerance.
        max_abs_deviation: The largest absolute per-position deviation observed
            across both the acceptor and donor tracks.
        n_cases: How many cases were checked.
        tolerance: The absolute tolerance used.
    """

    passed: bool
    max_abs_deviation: float
    n_cases: int
    tolerance: float


def verify_spliceai_fidelity(
    predictor: SpliceAiSplicePredictor,
    panel: Sequence[SpliceAiFidelityCase],
    *,
    tolerance: float = 1e-3,
) -> SpliceAiFidelityReport:
    """Check that ``predictor`` reproduces upstream SpliceAI's own scores.

    The integration-fidelity gate from CLAUDE.md section 6: proof the BT4 adapter
    faithfully reproduces the published model's outputs (both acceptor and donor
    tracks) on a captured reference panel -- *not* a from-scratch training gate.
    On success a maintainer may promote via
    ``dataclasses.replace(predictor, fidelity_verified=True)``. No reference panel
    ships (capturing one needs the CC BY-NC weights and yields licensed outputs),
    so shipped predictors never pass this gate.

    Args:
        predictor: The adapter to verify.
        panel: Captured reference cases. Must be non-empty.
        tolerance: Maximum allowed absolute per-position deviation.

    Returns:
        A :class:`SpliceAiFidelityReport`. ``passed`` is ``True`` iff every case's
        acceptor and donor scores match within ``tolerance``.

    Raises:
        ValueError: If ``panel`` is empty, or a case's expected-score length does
            not match the sequence length.
        ModuleNotFoundError / RuntimeError: If the adapter cannot run (deps or
            weights missing) -- fidelity cannot be asserted without running it.
    """
    if not panel:
        raise ValueError("fidelity panel must contain at least one case")
    max_dev = 0.0
    for case in panel:
        length = len(validate_dna(case.sequence))
        if len(case.expected_acceptor) != length or len(case.expected_donor) != length:
            raise ValueError(
                "fidelity case expected-score length does not match sequence length "
                f"{length} (acceptor={len(case.expected_acceptor)}, "
                f"donor={len(case.expected_donor)})"
            )
        acceptor, donor = predictor.site_scores(case.sequence)
        for got, want in zip(acceptor, case.expected_acceptor, strict=True):
            max_dev = max(max_dev, abs(got - want))
        for got, want in zip(donor, case.expected_donor, strict=True):
            max_dev = max(max_dev, abs(got - want))
    return SpliceAiFidelityReport(
        passed=max_dev <= tolerance,
        max_abs_deviation=max_dev,
        n_cases=len(panel),
        tolerance=tolerance,
    )
