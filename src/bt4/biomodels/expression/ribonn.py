"""Wrapped RiboNN translation-efficiency model behind the ``ExpressionPredictor``.

RiboNN (Zheng, Persyn, Wang et al., *Nat Biotechnol* 2025; Sanofi / Cenik Lab) is
a published multitask convolutional network that predicts mRNA **translation
efficiency (TE)** from transcript sequence across >140 mammalian cell types. This
adapter exposes it behind BT4's :class:`~bt4.biomodels.expression.base.ExpressionPredictor`
contract as the Phase-4 learned expression head (CLAUDE.md §6/§9).

**It wraps, it does not reimplement.** Unlike a pip package, RiboNN ships as a
repository plus a large multitask architecture (``TECNNMultiTasks``). Faithfully
re-deriving that network by hand would be exactly the fidelity risk the honesty
gates exist to catch, so this adapter **drives the user's own RiboNN checkout**:
it lazily imports the repo's ``src`` (its ``predict``/``data``/``model`` modules)
and loads the user's own weights, exactly as :mod:`bt4.biomodels.splice.pangolin`
drives the user's installed Pangolin. Nothing is vendored. Point it at the clone
via the ``BT4_RIBONN_DIR`` environment variable (or the ``repo_dir`` field), with
weights laid out as ``<weights_dir>/<species>/<run_id>/state_dict.pth`` plus a
``<species>/runs.csv`` (the layout of the published Zenodo ``weights.zip``);
``weights_dir`` defaults to ``<repo_dir>/models`` and is overridable via
``BT4_RIBONN_WEIGHTS``.

**License.** RiboNN's code and weights are each under a **Sanofi non-commercial**
license (academic / non-commercial use only). BT4 is open-source and
non-commercial, so this is compatible, and — like SpliceAI's CC BY-NC weights —
the material is **never bundled**: the user supplies their own clone and weights.
The bundled file :data:`ribonn_sha256.json <PINNED_WEIGHT_SHA256>` contains only
the *public content hashes* of the released weight files (no licensed bytes); the
adapter verifies every weight it loads against those pins *before* ``torch.load``
runs, keeping runs reproducible-from-manifest.

**Honest scope / units.** RiboNN's target is a **centered-log-ratio (CLR)
compositional residual** of TE (Methods; the model applies no log/exp itself), so
the score is reported in native ``"RiboNN CLR-residual TE"`` units — *never*
exponentiated. Orientation is larger-is-better (higher predicted TE). The score
is a *single validated model output*, not a hand-weighted composite, so it does
not trip the §10.5 magic-scalar trap. Nonetheless ``calibrated`` starts **False**
and mirrors :attr:`fidelity_verified`: reproducing upstream RiboNN faithfully is
necessary but not sufficient, because RiboNN was validated on *natural-transcript*
TE while BT4 optimizes *CDS variants with fixed UTRs* (RiboNN's own ablation puts
only ~31% of per-nucleotide signal in the CDS). Promotion to ``calibrated=True``
is earned only by passing the acceptance gate
(:func:`bt4.biomodels.expression.verify_expression_gate`) on a CDS-variant panel
from the deployment regime — a maintainer step, never assigned here.

Because it varies only the CDS with UTRs held fixed, the decision-relevant quantity
is :meth:`~RiboNNExpressionModel.delta_logte` — ``logTE(designed) - logTE(reference)``
— which isolates the CDS-attributable signal (a negative Δ flags a CDS change
predicted to *reduce* expression), analogous to
:meth:`bt4.biomodels.splice.PangolinSplicePredictor.delta_splicing`.

Importing this module stays light: ``torch``, ``pandas`` and the user's RiboNN
``src`` are imported only inside methods, never at module load.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

from bt4.biomodels.expression.base import ExpressionResult
from bt4.domain.sequence import validate_dna

__all__ = [
    "PINNED_WEIGHT_SHA256",
    "RiboNNExpressionModel",
    "load_pinned_sha256",
]

# Environment variables that point the adapter at the user's own RiboNN clone and
# (optionally) a separate weights directory. Nothing is bundled; the user supplies
# both the code and the Sanofi-non-commercial weights.
_ENV_REPO_DIR = "BT4_RIBONN_DIR"
_ENV_WEIGHTS_DIR = "BT4_RIBONN_WEIGHTS"

_SPECIES = ("human", "mouse")

# RiboNN's data loader hard-asserts that every CDS is length-3N and ends in a stop
# codon (it crashes the *whole* batch otherwise), so the adapter validates up front.
_STOP_CODONS = ("TAA", "TGA", "TAG")

_MANIFEST_PATH = Path(__file__).with_name("data") / "ribonn_sha256.json"


@lru_cache(maxsize=1)
def load_pinned_sha256() -> dict[str, str]:
    """Return the bundled ``{relative_path: sha256}`` pin map for RiboNN weights.

    These are the public SHA-256 content hashes of the released Zenodo
    ``weights.zip`` state-dict files (90 human + 90 mouse cross-validation runs) --
    *not* licensed model bytes. Keys are ``"<species>/<run_id>/state_dict.pth"``.
    """
    with _MANIFEST_PATH.open(encoding="utf-8") as handle:
        data: dict[str, str] = json.load(handle)
    return data


# The public pinned content hashes (loaded from the bundled manifest).
PINNED_WEIGHT_SHA256: dict[str, str] = load_pinned_sha256()


def _sha256_file(path: Path) -> str:
    """Return the SHA-256 hex digest of ``path`` (streamed, constant memory)."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class RiboNNExpressionModel:
    """Wrapped RiboNN TE predictor (drives the user's own checkout + weights).

    Attributes:
        species: ``"human"`` (default) or ``"mouse"`` -- selects the weight set.
        top_k: Number of top cross-validation runs (by ``val_r2``) to ensemble,
            matching RiboNN's ``predict_using_nested_cross_validation_models``.
        utr5: Fixed 5' UTR context prepended to every scored CDS. Empty by default,
            but **scoring requires it non-empty**: RiboNN's data loader reads an
            all-empty UTR column as NaN and cannot preprocess it, and the UTRs carry
            most of RiboNN's signal in any case. Because BT4 varies only the CDS,
            holding the real UTRs constant is the intended use (a modeling choice,
            documented, never silently zero-padded into the model's own channels).
        utr3: Fixed 3' UTR context (as ``utr5``; likewise required non-empty to score).
        repo_dir: Path to the RiboNN clone (its ``src`` package). ``None`` resolves
            from ``$BT4_RIBONN_DIR``.
        weights_dir: Directory holding ``<species>/<run_id>/state_dict.pth`` and
            ``<species>/runs.csv``. ``None`` resolves from ``$BT4_RIBONN_WEIGHTS``
            then ``<repo_dir>/models``.
        fidelity_verified: Whether this instance passed the CDS-variant acceptance
            gate; mirrored by :attr:`calibrated`. ``False`` by default and in every
            shipped configuration.
    """

    species: str = "human"
    top_k: int = 5
    utr5: str = ""
    utr3: str = ""
    repo_dir: str | None = None
    weights_dir: str | None = None
    fidelity_verified: bool = field(default=False)

    def __post_init__(self) -> None:
        if self.species not in _SPECIES:
            raise ValueError(f"species must be one of {_SPECIES}, got {self.species!r}")
        if self.top_k < 1:
            raise ValueError(f"top_k must be >= 1, got {self.top_k}")

    @property
    def name(self) -> str:
        """Stable backend identifier, e.g. ``"ribonn[human]"``."""
        return f"ribonn[{self.species}]"

    @property
    def calibrated(self) -> bool:
        """``True`` only after the CDS-variant acceptance gate has passed.

        Mirrors :attr:`fidelity_verified`. Wrapping a validated model is not the
        same as validating it *for BT4's regime*, so the shipped adapter is always
        ``calibrated is False`` (CLAUDE.md §6/§8/§10.6).
        """
        return self.fidelity_verified

    def _resolve_repo_dir(self) -> Path | None:
        raw = self.repo_dir or os.environ.get(_ENV_REPO_DIR)
        if not raw:
            return None
        path = Path(raw).expanduser()
        return path if (path / "src").is_dir() else None

    def _resolve_weights_dir(self, repo_dir: Path | None) -> Path | None:
        raw = self.weights_dir or os.environ.get(_ENV_WEIGHTS_DIR)
        if raw:
            candidate = Path(raw).expanduser()
        elif repo_dir is not None:
            candidate = repo_dir / "models"
        else:
            return None
        return candidate if (candidate / self.species).is_dir() else None

    def available(self) -> bool:
        """Return whether the adapter can actually run (deps + checkout + weights).

        ``True`` only when the RiboNN clone is resolvable, its ``<species>`` weight
        directory exists, and ``torch``/``pandas`` import. Never raises; a ``False``
        return is the honest "cannot score here" signal that keeps
        :func:`bt4.biomodels.expression.default` on the neutral placeholder.
        """
        import importlib.util

        repo = self._resolve_repo_dir()
        if repo is None or self._resolve_weights_dir(repo) is None:
            return False
        # Probe the optional heavy deps without importing them (keeps this light and
        # avoids typed-import fragility): find_spec returns None when absent.
        return all(importlib.util.find_spec(mod) is not None for mod in ("pandas", "torch"))

    def _verify_weights(self, weights_dir: Path) -> None:
        """Hash-verify every ``<species>`` state dict against the bundled pins.

        Raises before any weight is loaded if a file is missing, unpinned, or its
        digest does not match :data:`PINNED_WEIGHT_SHA256` -- so unverified or
        tampered bytes are never handed to ``torch.load`` (CLAUDE.md §4.3).
        """
        prefix = f"{self.species}/"
        expected = {k: v for k, v in PINNED_WEIGHT_SHA256.items() if k.startswith(prefix)}
        for rel, want in expected.items():
            path = weights_dir / rel
            if not path.is_file():
                raise FileNotFoundError(f"RiboNN weight file missing: {path}")
            got = _sha256_file(path)
            if got != want:
                raise ValueError(
                    f"RiboNN weight {rel} sha256 {got} != pinned {want}; refusing to load"
                )

    def _predict_te(self, dnas: list[str]) -> list[float]:
        """Run RiboNN on ``dnas`` and return one mean-TE (over cell types) per input.

        Drives the user's RiboNN checkout: verifies weights, writes a temporary
        input table (``tx_id, utr5_sequence, cds_sequence, utr3_sequence``) using the
        fixed UTR context, calls the repo's
        ``predict_using_nested_cross_validation_models``, and averages the
        ``predicted_TE_*`` cell-type columns. Heavy deps and the repo's ``src`` are
        imported here, never at module load.
        """
        import csv
        import importlib
        import sys
        import tempfile

        # RiboNN hard-asserts CDS is 3N ending in a stop codon; a bad row crashes the
        # whole batch, so refuse up front with a clear, per-input message (a pure input
        # check, done before touching the checkout).
        for i, dna in enumerate(dnas):
            if len(dna) % 3 != 0 or dna[-3:].upper() not in _STOP_CODONS:
                raise ValueError(
                    f"RiboNN requires each CDS to be length-3N ending in a stop codon; "
                    f"input {i} (len {len(dna)}, ends {dna[-3:]!r}) does not"
                )

        # RiboNN's data loader runs pandas ``.str`` preprocessing on the UTR columns;
        # an all-empty UTR column is read back as NaN (float) and crashes deep inside
        # that loader. Refuse up front with a clear message -- and the UTRs carry most
        # of RiboNN's signal anyway, so an empty-UTR score would not be meaningful.
        if not self.utr5 or not self.utr3:
            raise ValueError(
                "RiboNN requires non-empty 5' and 3' UTR context (an all-empty UTR "
                "column is read as NaN and breaks its .str preprocessing). Set "
                "utr5=/utr3= to the transcript's real UTRs; they are held fixed while "
                "the CDS varies."
            )

        repo = self._resolve_repo_dir()
        if repo is None:
            raise RuntimeError(
                f"RiboNN clone not found; set ${_ENV_REPO_DIR} to the repository path"
            )
        weights = self._resolve_weights_dir(repo)
        if weights is None:
            raise RuntimeError(
                f"RiboNN weights for {self.species!r} not found; set ${_ENV_WEIGHTS_DIR}"
            )
        self._verify_weights(weights)

        repo_str = str(repo)
        if repo_str not in sys.path:
            sys.path.insert(0, repo_str)
        predict_mod = importlib.import_module("src.predict")
        predict = predict_mod.predict_using_nested_cross_validation_models

        pd = importlib.import_module("pandas")  # Any-typed; keeps mypy env-independent

        runs_csv = weights / self.species / "runs.csv"
        run_df = pd.read_csv(runs_csv)

        tx_ids = [f"bt4_{i}" for i in range(len(dnas))]
        with tempfile.TemporaryDirectory() as tmp:
            in_path = Path(tmp) / "bt4_ribonn_input.tsv"
            with in_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence"])
                for tx_id, dna in zip(tx_ids, dnas, strict=True):
                    writer.writerow([tx_id, self.utr5, dna, self.utr3])
            out_df = _run_predict_with_models_layout(
                predict, weights, str(in_path), self.species, run_df, self.top_k
            )

        return _reduce_te_by_tx_id(out_df, tx_ids)

    def score_sequence(self, dna: str) -> ExpressionResult:
        """Return RiboNN's predicted mean TE for ``dna`` (CLR-residual units).

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (its start/stop codons
                included, per RiboNN's input contract).

        Returns:
            An :class:`ExpressionResult`; ``score`` is the mean predicted TE across
            cell types in native CLR-residual units (larger is better), and
            ``calibrated`` mirrors :attr:`calibrated` (``False`` until the gate
            passes -- never read it as a validated expression prediction otherwise).
        """
        validate_dna(dna)
        (te,) = self._predict_te([dna])
        return ExpressionResult(
            score=te,
            model_name=self.name,
            calibrated=self.calibrated,
            units=f"RiboNN CLR-residual TE (mean over {self.species} cell types)",
        )

    def delta_logte(self, designed_dna: str, reference_dna: str) -> float:
        """Return ``TE(designed) - TE(reference)`` -- the CDS-attributable signal.

        Holding the UTRs fixed, this isolates what the CDS redesign changes: a
        positive value predicts higher expression, a negative value flags a CDS
        change RiboNN predicts will *reduce* expression (the maintainer's
        limiting-sequence framing). Analogous to
        :meth:`bt4.biomodels.splice.PangolinSplicePredictor.delta_splicing`.
        """
        validate_dna(designed_dna)
        validate_dna(reference_dna)
        designed, reference = self._predict_te([designed_dna, reference_dna])
        return designed - reference


def _reduce_te_by_tx_id(out_df: Any, tx_ids: list[str]) -> list[float]:
    """Reduce RiboNN's raw prediction table to one mean-TE per input, in input order.

    RiboNN returns the ensemble's predictions as **multiple rows per input** -- one per
    cross-validation model kept by ``top_k`` -- each carrying the ``predicted_TE_*``
    per-cell-type columns. The decision-relevant scalar is the mean over cell types
    *and* over the ensemble, so we average the per-cell-type columns within each row and
    then average those rows per ``tx_id``. Grouping collapses the ensemble (and is a
    no-op when the table already has one row per input); a plain ``set_index`` left
    duplicate ``tx_id`` labels, so ``float(ordered[tx_id])`` saw a Series and raised.

    A ``tx_id`` absent from the output was dropped by RiboNN's length cap (5'UTR > 1381
    or CDS+3'UTR > 11937 nt), so it is surfaced honestly rather than as a ``KeyError``.
    """
    te_cols = [c for c in out_df.columns if str(c).startswith("predicted_TE_")]
    if not te_cols:
        raise RuntimeError("RiboNN output had no predicted_TE_* columns")
    means = out_df[te_cols].mean(axis=1)
    ordered = out_df.assign(_te=means).groupby("tx_id")["_te"].mean()
    results: list[float] = []
    for i, tx_id in enumerate(tx_ids):
        if tx_id not in ordered.index:
            raise ValueError(
                f"RiboNN produced no prediction for input {i}; its CDS likely "
                "exceeds RiboNN's length cap (CDS+3'UTR <= 11937 nt)"
            )
        results.append(float(ordered[tx_id]))
    return results


def _run_predict_with_models_layout(
    predict: Any, weights_dir: Path, input_path: str, species: str, run_df: Any, top_k: int
) -> Any:
    """Call RiboNN's predict fn with the ``models/<species>/...`` layout it hard-codes.

    RiboNN loads weights from the **literal** relative path
    ``models/<species>/<run_id>/state_dict.pth``, so a ``models`` directory must
    resolve from the process working directory. ``weights_dir`` is the directory
    that *contains* ``<species>/...``. We satisfy the contract without moving the
    user's (large, licensed) weights: if ``weights_dir`` is already named ``models``
    we run from its parent; otherwise we run from a temp dir holding a ``models``
    symlink to ``weights_dir``. The working directory is always restored.

    Raises:
        RuntimeError: If neither layout can be arranged (e.g. the platform refuses
            the symlink and the directory is not named ``models``).
    """
    import os as _os
    import tempfile

    prev = _os.getcwd()
    if weights_dir.name == "models":
        run_root = weights_dir.parent
        try:
            _os.chdir(run_root)
            return predict(
                input_path=input_path, species=species, run_df=run_df,
                top_k_models_to_use=top_k,
            )
        finally:
            _os.chdir(prev)

    with tempfile.TemporaryDirectory() as tmp:
        link = Path(tmp) / "models"
        try:
            _os.symlink(weights_dir, link, target_is_directory=True)
        except OSError as exc:
            raise RuntimeError(
                "RiboNN loads weights from a hard-coded 'models/' path; could not "
                f"link it to {weights_dir} ({exc}). Point $BT4_RIBONN_WEIGHTS at a "
                "directory named 'models' (as the Zenodo weights.zip extracts)."
            ) from exc
        try:
            _os.chdir(tmp)
            return predict(
                input_path=input_path, species=species, run_df=run_df,
                top_k_models_to_use=top_k,
            )
        finally:
            _os.chdir(prev)
