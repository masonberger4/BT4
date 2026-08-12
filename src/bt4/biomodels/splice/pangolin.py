"""Pangolin-backed splice predictor -- a wrapped published CNN, lazily imported.

:class:`PangolinSplicePredictor` is the first of BT4's *wrapped published* splice
backends (CLAUDE.md section 6: "wrap published SpliceAI + Pangolin -- no
self-training"). It runs the already-validated **Pangolin** deep-learning splice
model (Zeng & Li, *Genome Biology* 2022, doi:10.1186/s13059-022-02664-4) as an
*inference-only* backend behind the existing
:class:`~bt4.biomodels.splice.base.SplicePredictor` contract, feeding Pangolin's
per-nucleotide splice-site probabilities into BT4's already-shipped Delta-splicing
framing and top-k / log-odds pooling (:mod:`bt4.biomodels.splice.base`).

Three facts shape this module, each an honesty rule made structural:

* **Pangolin is GPL-3.0 -- NOT MIT.** The constitution's roadmap called Pangolin
  "MIT"; the upstream repository (https://github.com/tkzeng/Pangolin) is licensed
  **GNU GPL v3**, and its bundled model weights are GPL artifacts too. BT4 is MIT,
  so it **neither bundles Pangolin's weights nor copies its source**. Instead --
  exactly as :class:`~bt4.biomodels.folding.vienna.ViennaFoldingModel` wraps an
  *installed* ViennaRNA rather than reimplementing it -- this adapter **lazily
  imports the user's own installed ``pangolin`` package** (its ``nn.Module`` and
  trained weights) and drives it. The GPL obligation stays with the user's
  install; BT4 ships only the glue (encoding, ensemble averaging, the honest
  Delta-splicing wiring) and the *SHA-256 content hashes* of the published
  weights (facts, not code). See :mod:`bt4.biomodels.splice` and CLAUDE.md
  section 6 ("verify its weight license before bundling anything").

* **``calibrated`` starts ``False`` and is gated on an integration-fidelity
  check.** Wrapping a validated model is not the same as having verified the
  wrapping. :attr:`PangolinSplicePredictor.calibrated` mirrors
  :attr:`fidelity_verified`, which is ``False`` unless the adapter has been shown
  to reproduce upstream Pangolin's own scores on a captured reference panel
  (:func:`verify_pangolin_fidelity`) -- *not* a from-scratch training gate, just
  proof the adapter faithfully reproduces the published model (CLAUDE.md section
  6). No reference panel ships (capturing it needs the GPL weights and produces
  GPL-derived numbers), so the shipped adapter is always ``calibrated is False``
  and :func:`bt4.biomodels.splice.default` keeps returning the labeled PWM
  baseline. Until the gate passes, Pangolin's numbers are a strong prior, never a
  validated result.

* **Weights are hash-pinned and verified before they are unpickled.** Every
  weight file is checked against :data:`PINNED_WEIGHT_SHA256` (the SHA-256 of the
  production ``.3.v2`` weights the installed ``pangolin`` CLI loads -- folds 1-3,
  the retrained ensemble, *not* the older ``.3`` demo weights) *before*
  ``torch.load`` touches it, so the adapter (a) stays reproducible-from-manifest
  -- unlike BT3's live ASSP scrape (CLAUDE.md section 10.15) -- and (b) never
  unpickles bytes it has not content-verified. Unknown or altered weights are
  refused, not run.

**Honest scope (unchanged from the contract).** Pangolin predicts splice-*site
presence*; a lower pooled Delta-splicing means lower *predicted cryptic-splice
risk* -- a strong prior, but not validated expression gain (the same
CAI-as-weak-proxy caution, CLAUDE.md sections 1 and 6).

**Pangolin reports a single combined splice-site probability per position** (a
2-class ``P(splice) / P(not)`` softmax per tissue), not a separated donor /
acceptor pair like SpliceAI. The adapter therefore reports that combined
per-position probability in :attr:`SpliceResult.donor` and leaves
:attr:`SpliceResult.acceptor` all-zero, so :func:`pooled_risk` (which pools
``donor`` union ``acceptor``) counts each site exactly once. Use
:meth:`PangolinSplicePredictor.site_scores` for the combined track without
relying on that field convention.

The heavy dependencies (``torch`` and the GPL ``pangolin`` package) are imported
**only inside methods**, never at module load, so importing this module -- and
``import bt4`` -- stays lightweight (CLAUDE.md section 3). This module depends
only on :mod:`bt4.domain`, the standard library, and -- lazily, inside methods --
the optional ``torch`` / ``pangolin`` packages.
"""

from __future__ import annotations

import hashlib
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
    "DEFAULT_TISSUES",
    "PANGOLIN_CONTEXT",
    "PANGOLIN_FLANK",
    "PINNED_WEIGHT_SHA256",
    "TISSUE_OUTPUTS",
    "FidelityCase",
    "FidelityReport",
    "PangolinSplicePredictor",
    "TissueOutput",
    "verify_pangolin_fidelity",
]


PANGOLIN_CONTEXT: int = 10_000
"""Total receptive-field context Pangolin crops from its input (``2*sum(AR*(W-1))``).

Pangolin's ``forward`` pass pads its skip tensor by ``-CL/2`` on each side, where
``CL = 2 * sum(AR * (W - 1)) = 10000`` for the published architecture. An input
of length ``N`` therefore yields per-position scores for the middle ``N - 10000``
bases. To score every position of a coding sequence, the adapter pads it with
:data:`PANGOLIN_FLANK` ``N`` bases on each side (Pangolin's documented convention
for sequences shorter than the context).
"""

PANGOLIN_FLANK: int = PANGOLIN_CONTEXT // 2
"""Number of ``N`` bases padded on each side so output aligns 1:1 with the CDS."""


@dataclass(frozen=True, slots=True)
class TissueOutput:
    """How to read one tissue's P(splice) head from the Pangolin ensemble.

    Attributes:
        weight_index: The ``i`` in Pangolin's ``final.{j}.{i}.3.v2`` weight files
            (``j`` in 1..3 are the production folds, see :data:`_CV_FOLDS`); this
            ``i`` selects the trained output group.
        output_channel: The channel of the model's 12-way output tensor holding
            this tissue's ``P(splice)`` probability (Pangolin's channel map).
    """

    weight_index: int
    output_channel: int


# Pangolin's four tissue P(splice) heads. The weight_index / output_channel pairs
# are Pangolin's own published usage constants: the production CLI
# (pangolin/pangolin.py) loads weight indices ``[0,2,4,6]`` and reads output
# channels ``[1,4,7,10]`` for Heart / Liver / Brain / Testis P(splice). These are
# interface facts about how to read the published model, not GPL source code.
TISSUE_OUTPUTS: dict[str, TissueOutput] = {
    "heart": TissueOutput(weight_index=0, output_channel=1),
    "liver": TissueOutput(weight_index=2, output_channel=4),
    "brain": TissueOutput(weight_index=4, output_channel=7),
    "testis": TissueOutput(weight_index=6, output_channel=10),
}

DEFAULT_TISSUES: tuple[str, ...] = ("heart", "liver", "brain", "testis")
"""Default tissue set: average P(splice) across all four for a tissue-agnostic risk.

Averaging the four tissue heads gives a tissue-agnostic cryptic-splice-site
probability, appropriate for coding-sequence redesign (which is not
tissue-targeted). Callers may pass a single tissue for speed.
"""

_CV_FOLDS: tuple[int, ...] = (1, 2, 3)
"""The three cross-validation folds the production Pangolin CLI averages.

Pangolin's command-line tool (`pangolin/pangolin.py`) loads the **retrained
``.3.v2`` weights over folds 1-3** (``for j in range(1,4)``), and its
``compute_score`` averages those three folds per tissue. (Its older demo script
``custom_usage.py`` used the pre-``v2`` ``.3`` weights over five folds; the CLI is
the model users actually run, so the adapter tracks the CLI.)
"""

# SHA-256 of Pangolin's production P(splice) weight files -- the retrained
# ``final.{fold}.{weight_index}.3.v2`` files (folds 1-3) that the installed
# `pangolin` CLI loads (github.com/tkzeng/Pangolin). Weights are GPL artifacts
# kept OUT of this repository; only these content hashes ship, and each file is
# verified against its pin *before* it is unpickled (CLAUDE.md sections 6, 10.15).
PINNED_WEIGHT_SHA256: dict[str, str] = {
    "final.1.0.3.v2": "f0478fab173b75f7f7e9fe96688bad6c50fa4a46d70557f423b110caaf565501",
    "final.2.0.3.v2": "c4c6bb4880fa6fb28b14182ae3ea0600edb07056158f55325b5e6e6e48fc9f26",
    "final.3.0.3.v2": "ec685a6e7105a4486c1f89a005458a13deb3fe7171f13d434f4877e386d10676",
    "final.1.2.3.v2": "559c05de3e1ce65c2515ca3e92ef85edb0ec2e47686ca58060e25891ce06eb3a",
    "final.2.2.3.v2": "48758ba8b95eee9aa9feea52672ef06ca1b34111299c27f8a710f734d8b9aae5",
    "final.3.2.3.v2": "7cb576c2b24db4fdd6970c4ca4fb7c20ae1b1d8ae80645ebbe689848b5743129",
    "final.1.4.3.v2": "c50b12e0c0af776d5674ca5e346493f8265783494d4df383364de9c1136657f6",
    "final.2.4.3.v2": "e03303bed4fd6f135ec0f6c1b192cce954ea42d0646f44d17b4a6fbb2b1f610e",
    "final.3.4.3.v2": "9476d2e25520d7ff15bece0cd5d3b657e3b1dd3cc5fcab1d9c3b62bea7a0c5b6",
    "final.1.6.3.v2": "2aae563fa18a8a9b6699c6c96e0d32b8ec7543f8f805fb3bc9de77302cc9f66e",
    "final.2.6.3.v2": "7d3c0b1b2a60067b940dec315567874fbc8bcd322f1b7c76bf969f51f0f53f7f",
    "final.3.6.3.v2": "756e7721a382cace24e9bfea5b543af5623f2487d9a3efe7385e9c76367005fd",
}

_WEIGHT_FILE_TEMPLATE: str = "final.{fold}.{index}.3.v2"
"""Upstream weight-file name pattern (the production ``.v2`` files, folds 1-3)."""

# One-hot channel order A, C, G, T -- matching Pangolin's IN_MAP so the pinned
# weights receive the input they were trained on. 'N' (padding) maps to all-zero.
_BASE_ROW: dict[str, tuple[float, float, float, float]] = {
    "A": (1.0, 0.0, 0.0, 0.0),
    "C": (0.0, 1.0, 0.0, 0.0),
    "G": (0.0, 0.0, 1.0, 0.0),
    "T": (0.0, 0.0, 0.0, 1.0),
    "N": (0.0, 0.0, 0.0, 0.0),
}

_WEIGHTS_ENV_VAR: str = "BT4_PANGOLIN_MODEL_DIR"
"""Optional env var overriding where the pinned Pangolin weight files are found."""


def _one_hot_channels(seq: str) -> list[list[float]]:
    """One-hot encode ``seq`` to channel-major ``[4][len]`` (A, C, G, T order).

    Pure standard-library encoding (no numpy, no torch), so it is unit-testable
    without the optional deps. Matches Pangolin's ``IN_MAP`` base ordering; ``N``
    maps to an all-zero column.

    Args:
        seq: An upper-case string over ``{A,C,G,T,N}``.

    Returns:
        Four lists (one per base channel), each ``len(seq)`` long.
    """
    channels: list[list[float]] = [[0.0] * len(seq) for _ in range(4)]
    for pos, base in enumerate(seq):
        row = _BASE_ROW[base]
        for chan in range(4):
            channels[chan][pos] = row[chan]
    return channels


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
        path: A Pangolin ``final.{fold}.{i}.3`` weight file.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        KeyError: If the file name is not in :data:`PINNED_WEIGHT_SHA256`.
        ValueError: If the file's SHA-256 does not match its pin -- the weights
            are unknown or altered and are refused, not run.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Pangolin weight file not found: {path}")
    expected = PINNED_WEIGHT_SHA256.get(path.name)
    if expected is None:
        raise KeyError(
            f"Pangolin weight file '{path.name}' is not a pinned P(splice) weight "
            f"(known: {sorted(PINNED_WEIGHT_SHA256)})"
        )
    actual = _sha256_file(path)
    if actual != expected:
        raise ValueError(
            f"Pangolin weight file '{path.name}' failed its SHA-256 pin "
            f"(expected {expected}, got {actual}); refusing to load unverified weights"
        )


def _import_torch() -> Any:
    """Import and return the ``torch`` module (lazy, guarded).

    Raises:
        ModuleNotFoundError: If ``torch`` is not installed. Install the
            ``bt4[splice-pangolin]`` extra (PyTorch) plus the GPL ``pangolin``
            package (see the module docstring), or use the baseline predictor.
    """
    try:
        import torch

        return torch
    except ImportError as exc:
        raise ModuleNotFoundError(
            "PyTorch is not installed; install the 'bt4[splice-pangolin]' extra "
            "and the GPL 'pangolin' package, or use the baseline splice predictor"
        ) from exc


def _import_pangolin_model() -> tuple[Any, Any, Any, Any]:
    """Import the installed ``pangolin`` model class and its architecture constants.

    BT4 does not reimplement or bundle Pangolin's GPL ``nn.Module``; it drives the
    user's installed ``pangolin`` package. Returns ``(Pangolin, L, W, AR)``.

    Raises:
        ModuleNotFoundError: If the GPL ``pangolin`` package is not importable.
            Install it from https://github.com/tkzeng/Pangolin (GPL-3.0), which
            provides both the module and its weights.
    """
    try:
        from pangolin.model import AR, L, Pangolin, W

        return Pangolin, L, W, AR
    except ImportError as exc:
        raise ModuleNotFoundError(
            "The 'pangolin' package (GPL-3.0) is not installed; install it from "
            "https://github.com/tkzeng/Pangolin to use the Pangolin splice backend, "
            "or use the baseline splice predictor"
        ) from exc


def _resolve_model_dir(explicit: str | None = None) -> Path | None:
    """Resolve the directory holding Pangolin's weight files, or ``None``.

    Resolution order: an explicit path, then the :data:`_WEIGHTS_ENV_VAR`
    environment variable, then the ``models`` directory bundled with the installed
    ``pangolin`` package. Never raises: returns ``None`` when nothing resolves so
    :meth:`PangolinSplicePredictor.available` can report unavailability cleanly.

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
        import pangolin
    except ImportError:
        return None
    pkg_file = getattr(pangolin, "__file__", None)
    if not pkg_file:
        return None
    models = Path(pkg_file).resolve().parent / "models"
    return models if models.is_dir() else None


@lru_cache(maxsize=4)
def _load_ensemble(model_dir: str, tissues: tuple[str, ...]) -> dict[str, list[Any]]:
    """Load and cache the Pangolin CV ensemble for ``tissues`` from ``model_dir``.

    Each tissue loads its three cross-validation weight files (folds 1-3, per
    :data:`_CV_FOLDS`), verifying every file's SHA-256 pin *before* unpickling it.
    Cached per ``(model_dir, tissues)`` so the expensive load happens once per
    process.

    Args:
        model_dir: Directory holding the ``final.{fold}.{i}.3.v2`` weight files.
        tissues: The tissue keys to load (subset of :data:`TISSUE_OUTPUTS`).

    Returns:
        A map ``tissue -> [3 loaded, eval-mode Pangolin models]``.

    Raises:
        ModuleNotFoundError: If ``torch`` or ``pangolin`` is not importable.
        FileNotFoundError / KeyError / ValueError: From :func:`_verify_weight_file`.
    """
    torch = _import_torch()
    pangolin_cls, dim_l, widths, rates = _import_pangolin_model()
    base = Path(model_dir)
    ensemble: dict[str, list[Any]] = {}
    for tissue in tissues:
        spec = TISSUE_OUTPUTS[tissue]
        models: list[Any] = []
        for fold in _CV_FOLDS:
            path = base / _WEIGHT_FILE_TEMPLATE.format(fold=fold, index=spec.weight_index)
            _verify_weight_file(path)
            state = torch.load(path, map_location=torch.device("cpu"))
            model = pangolin_cls(dim_l, widths, rates)
            model.load_state_dict(state)
            model.eval()
            models.append(model)
        ensemble[tissue] = models
    return ensemble


@dataclass(frozen=True, slots=True)
class PangolinSplicePredictor:
    """Wrapped published Pangolin splice model behind the ``SplicePredictor`` contract.

    Runs the installed GPL ``pangolin`` package's CNN ensemble as an
    inference-only backend, reporting a combined per-position splice-site
    probability. See the module docstring for the GPL boundary, the
    hash-pinning, and why ``calibrated`` starts ``False``.

    Attributes:
        tissues: Which Pangolin tissue P(splice) heads to average. Defaults to
            all four (:data:`DEFAULT_TISSUES`) for a tissue-agnostic risk.
        top_k: Number of strongest sites summed by :meth:`delta_splicing`'s
            top-k / log-odds pooling. Defaults to
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`.
        model_dir: Explicit weights directory, or ``None`` to resolve from the
            :data:`_WEIGHTS_ENV_VAR` env var or the installed ``pangolin`` package.
        fidelity_verified: Whether this instance has passed
            :func:`verify_pangolin_fidelity` against a captured upstream panel.
            Mirrored by :attr:`calibrated`; ``False`` by default and in every
            shipped configuration (no reference panel ships).
    """

    tissues: tuple[str, ...] = DEFAULT_TISSUES
    top_k: int = DEFAULT_TOP_K
    model_dir: str | None = None
    fidelity_verified: bool = field(default=False)

    def __post_init__(self) -> None:
        """Validate the tissue set and pooling depth.

        Raises:
            ValueError: If ``tissues`` is empty, names an unknown tissue, or
                ``top_k`` is not a positive integer.
        """
        if not self.tissues:
            raise ValueError("tissues must name at least one Pangolin tissue head")
        unknown = [t for t in self.tissues if t not in TISSUE_OUTPUTS]
        if unknown:
            raise ValueError(
                f"unknown Pangolin tissue(s) {unknown}; known: {sorted(TISSUE_OUTPUTS)}"
            )
        if self.top_k <= 0:
            raise ValueError(f"top_k must be a positive integer, got {self.top_k}")

    @property
    def name(self) -> str:
        """Backend identifier (includes the averaged tissue set for provenance)."""
        return f"pangolin[{'+'.join(self.tissues)}]"

    @property
    def calibrated(self) -> bool:
        """``True`` only once the integration-fidelity gate has passed.

        Mirrors :attr:`fidelity_verified`. Wrapping a validated model is not the
        same as verifying the wrapping, so this stays ``False`` until
        :func:`verify_pangolin_fidelity` confirms the adapter reproduces upstream
        Pangolin's own scores. No reference panel ships, so shipped instances are
        always uncalibrated (CLAUDE.md sections 6 and 10.6).
        """
        return self.fidelity_verified

    def available(self) -> bool:
        """Return whether this backend can run (never raises).

        ``True`` iff ``torch`` and the ``pangolin`` package import and a weights
        directory containing this predictor's required files resolves. Does not
        hash-verify the files (that happens at load); a lightweight existence
        check keeps ``available`` cheap and non-raising, mirroring
        :meth:`ViennaFoldingModel.available`.
        """
        try:
            _import_torch()
            _import_pangolin_model()
        except ImportError:
            return False
        model_dir = _resolve_model_dir(self.model_dir)
        if model_dir is None:
            return False
        for tissue in self.tissues:
            index = TISSUE_OUTPUTS[tissue].weight_index
            for fold in _CV_FOLDS:
                name = _WEIGHT_FILE_TEMPLATE.format(fold=fold, index=index)
                if not (model_dir / name).is_file():
                    return False
        return True

    def site_scores(self, dna: str) -> tuple[float, ...]:
        """Return Pangolin's combined per-position splice-site probability.

        This is the tissue-averaged, CV-ensembled ``P(splice)`` at each
        nucleotide of ``dna`` -- Pangolin's single combined site probability, not
        a separated donor / acceptor pair. Prefer this over
        :attr:`SpliceResult.donor` when you want the combined track explicitly.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            One probability in ``[0, 1]`` per input nucleotide, aligned 1:1.

        Raises:
            ModuleNotFoundError: If ``torch`` / ``pangolin`` is not importable.
            RuntimeError: If no weights directory resolves, or the model's output
                length does not align to the input (an architecture-mismatch guard).
            FileNotFoundError / KeyError / ValueError: From weight verification.
        """
        seq = validate_dna(dna)
        model_dir = _resolve_model_dir(self.model_dir)
        if model_dir is None:
            raise RuntimeError(
                "no Pangolin weights directory resolved; set "
                f"${_WEIGHTS_ENV_VAR}, pass model_dir=, or install the 'pangolin' package"
            )
        torch = _import_torch()
        ensemble = _load_ensemble(str(model_dir), self.tissues)

        padded = "N" * PANGOLIN_FLANK + seq + "N" * PANGOLIN_FLANK
        tensor = torch.tensor([_one_hot_channels(padded)], dtype=torch.float32)

        n = len(seq)
        per_tissue: list[Any] = []
        with torch.no_grad():
            for tissue in self.tissues:
                channel = TISSUE_OUTPUTS[tissue].output_channel
                fold_scores = [model(tensor)[0, channel, :] for model in ensemble[tissue]]
                per_tissue.append(torch.stack(fold_scores).mean(dim=0))
            combined = torch.stack(per_tissue).mean(dim=0)

        scores = [float(value) for value in combined]
        if len(scores) != n:
            raise RuntimeError(
                f"Pangolin output length {len(scores)} != input length {n}; "
                f"the model's context crop may differ from the expected {PANGOLIN_CONTEXT}"
            )
        return tuple(scores)

    def score_sequence(self, dna: str) -> SpliceResult:
        """Return per-position site scores as a :class:`SpliceResult`.

        Pangolin reports one combined splice-site probability per position (it
        does not separate donor from acceptor), so the combined track is placed
        in ``donor`` and ``acceptor`` is all-zero -- keeping
        :func:`pooled_risk` (which pools ``donor`` union ``acceptor``) counting
        each site exactly once. See :meth:`site_scores` for the combined track by
        itself.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            A :class:`SpliceResult`; scores are calibrated probabilities only when
            :attr:`calibrated` is ``True`` (i.e. after the fidelity gate passes).

        Raises:
            ModuleNotFoundError / RuntimeError / FileNotFoundError / KeyError /
            ValueError: As documented on :meth:`site_scores`.
        """
        sites = self.site_scores(dna)
        return SpliceResult(
            donor=sites,
            acceptor=(0.0,) * len(sites),
            model_name=self.name,
            calibrated=self.calibrated,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Return the negated added splice risk of ``designed`` vs ``reference``.

        See :meth:`bt4.biomodels.splice.base.SplicePredictor.delta_splicing` for
        the fixed *larger-is-better* orientation. Concretely returns
        ``pooled_risk(reference) - pooled_risk(designed)`` using top-k / log-odds
        pooling, so it is ``0.0`` for identical sequences, positive when the
        redesign lowers predicted splice risk, and negative when it raises it.

        Raises:
            ModuleNotFoundError / RuntimeError / FileNotFoundError / KeyError /
            ValueError: As documented on :meth:`site_scores`.
        """
        designed_risk = pooled_risk(self.score_sequence(designed_dna), self.top_k)
        reference_risk = pooled_risk(self.score_sequence(reference_dna), self.top_k)
        return reference_risk - designed_risk


@dataclass(frozen=True, slots=True)
class FidelityCase:
    """One captured-from-upstream reference case for the integration-fidelity gate.

    A maintainer with the GPL Pangolin install and weights runs upstream
    Pangolin's production per-position scoring (the ``.v2`` / folds-1-3 ensemble
    the ``pangolin`` CLI uses, reading channels ``[1,4,7,10]``) on ``sequence``
    for the predictor's configured tissue set and records the resulting combined
    per-position probabilities here. :func:`verify_pangolin_fidelity` then asserts
    the adapter reproduces them.

    Attributes:
        sequence: The coding sequence scored.
        expected_site_scores: Upstream Pangolin's combined per-position
            ``P(splice)`` for ``sequence``, aligned 1:1.
    """

    sequence: str
    expected_site_scores: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class FidelityReport:
    """Result of an integration-fidelity check.

    Attributes:
        passed: Whether every case matched within tolerance.
        max_abs_deviation: The largest absolute per-position deviation observed.
        n_cases: How many cases were checked.
        tolerance: The absolute tolerance used.
    """

    passed: bool
    max_abs_deviation: float
    n_cases: int
    tolerance: float


def verify_pangolin_fidelity(
    predictor: PangolinSplicePredictor,
    panel: Sequence[FidelityCase],
    *,
    tolerance: float = 1e-3,
) -> FidelityReport:
    """Check that ``predictor`` reproduces upstream Pangolin's own scores.

    This is the integration-fidelity gate from CLAUDE.md section 6: it is *not* a
    from-scratch training gate, only proof the BT4 adapter faithfully reproduces
    the published model's outputs on a captured reference panel. On success a
    maintainer may construct a promoted predictor via
    ``dataclasses.replace(predictor, fidelity_verified=True)`` so its
    :attr:`~PangolinSplicePredictor.calibrated` reports ``True``. No reference
    panel ships (capturing one needs the GPL weights and yields GPL-derived
    numbers), so shipped predictors never pass this gate.

    Args:
        predictor: The adapter to verify (its tissue set must match how the panel
            was captured).
        panel: Captured reference cases. Must be non-empty.
        tolerance: Maximum allowed absolute per-position deviation.

    Returns:
        A :class:`FidelityReport`. ``passed`` is ``True`` iff every case's
        per-position scores match within ``tolerance``.

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
        expected = case.expected_site_scores
        if len(expected) != len(validate_dna(case.sequence)):
            raise ValueError(
                "fidelity case expected_site_scores length "
                f"{len(expected)} != sequence length {len(case.sequence)}"
            )
        actual = predictor.site_scores(case.sequence)
        for got, want in zip(actual, expected, strict=True):
            max_dev = max(max_dev, abs(got - want))
    return FidelityReport(
        passed=max_dev <= tolerance,
        max_abs_deviation=max_dev,
        n_cases=len(panel),
        tolerance=tolerance,
    )
