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
only ~31% of per-*nucleotide* signal in the CDS -- though note its
length-integrated total attribution is 22/73/5 for 5'UTR/CDS/3'UTR, so the CDS is
the majority of the total signal; both numbers are real and quoting only one
misleads). The load-bearing gap is different and sharper: RiboNN has never been
shown to discriminate *synonymous CDS variants of the same protein under a fixed
UTR*, which is exactly BT4's regime. Promotion to ``calibrated=True``
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
import math
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
    "RiboNNFoldPrediction",
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


@dataclass(frozen=True, slots=True)
class RiboNNFoldPrediction:
    """One RiboNN ensemble-fold prediction for one input (immutable).

    RiboNN's ``predict`` returns **one row per input per outer fold** (10 folds,
    each row already the mean of that fold's ``top_k`` models), tagged with a
    ``fold`` column. :meth:`RiboNNExpressionModel.score_many` averages those folds,
    which is correct for a *novel designed* sequence -- no fold ever saw it.

    It is **wrong** for a natural transcript, because nine of the ten folds had that
    transcript's own label in training. Reproducing RiboNN's published held-out
    accuracy therefore requires keeping only the fold that held the transcript out,
    which needs the per-fold values rather than their mean -- hence this record and
    :meth:`RiboNNExpressionModel.predict_folds`.

    Attributes:
        index: Position of the sequence in the ``dnas`` list that was scored.
        fold: RiboNN's outer-fold id for this prediction.
        te: Mean predicted TE over the selected cell types, in native CLR-residual
            units (larger is better).
    """

    index: int
    fold: int
    te: float


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
        batch_size: Inference batch size handed to RiboNN's ``predict`` (its own
            default is 1024). **Purely a memory/speed knob -- it cannot change a
            score.** RiboNN pads every transcript to a *fixed* width
            (``max_utr5_len + max_cds_utr3_len`` = 13318), not to the longest member
            of a batch, and its predict dataloader is built with ``shuffle=False``,
            so batch composition affects neither the encoding nor the row order.
            BT4 defaults to 64 because a 1024-row batch of
            ``(channels, 13318)`` float32 tensors is hundreds of MB *per batch*
            before worker prefetch, which OOMs an ordinary CPU box. Raise it on a
            GPU.
        cell_types: Which of RiboNN's per-cell-type outputs to average. Empty (the
            default) means **all** of them -- 78 columns for human, 68 for mouse.
            Averaging all 78 is the right summary for a generic design, but it is a
            **scope error** when comparing against a measurement from one cell line:
            a HEK293T ribosome-load panel should be scored against
            ``cell_types=("HEK293T",)``, not against the mean of 78 tissues. Names
            are matched against RiboNN's ``predicted_TE_<name>`` output columns, and
            an unmatched name **raises** (listing what is available) rather than
            being quietly dropped.
        num_workers: DataLoader worker processes handed to RiboNN's ``predict``
            (its own default is 4). BT4 defaults to **0**, which is a correctness
            requirement rather than a tuning choice: this adapter scores from a
            mutated ``sys.path`` and a temporary working directory (see
            :func:`_run_predict_with_models_layout`), and on any platform whose
            multiprocessing start method is *spawn* (Windows, macOS) each worker
            re-imports the module and does not inherit either, so workers hang or
            fail. RiboNN also rebuilds the dataloader once per ensemble member
            (``top_k`` x folds = up to 50 times), paying the spawn cost every
            time. Worker count never affects a score.
        fidelity_verified: Whether this instance passed the CDS-variant acceptance
            gate; mirrored by :attr:`calibrated`. ``False`` by default and in every
            shipped configuration.
        attestation_sha256: Content hash of the attestation that promoted this head,
            or ``""`` when it is uncalibrated. Set **only** by
            :func:`~bt4.biomodels.expression.attestation.verified_predictor`, and folded
            into the run manifest so a result steered by a calibrated head records
            *which* claim authorized it -- two runs promoted by different attestations
            must not share a provenance stamp (invariant #9).
    """

    species: str = "human"
    top_k: int = 5
    utr5: str = ""
    utr3: str = ""
    repo_dir: str | None = None
    weights_dir: str | None = None
    batch_size: int = 64
    num_workers: int = 0
    cell_types: tuple[str, ...] = ()
    fidelity_verified: bool = field(default=False)
    attestation_sha256: str = ""

    def __post_init__(self) -> None:
        if self.species not in _SPECIES:
            raise ValueError(f"species must be one of {_SPECIES}, got {self.species!r}")
        if self.top_k < 1:
            raise ValueError(f"top_k must be >= 1, got {self.top_k}")
        if self.batch_size < 1:
            raise ValueError(f"batch_size must be >= 1, got {self.batch_size}")
        if self.num_workers < 0:
            raise ValueError(f"num_workers must be >= 0, got {self.num_workers}")
        if any(not name.strip() for name in self.cell_types):
            raise ValueError("cell_types must not contain blank names")
        if len(set(self.cell_types)) != len(self.cell_types):
            raise ValueError(f"cell_types contains duplicates: {self.cell_types}")

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

    def _run_predict(self, dnas: list[str]) -> tuple[Any, list[str]]:
        """Run RiboNN on ``dnas`` and return its ``(raw prediction table, tx_ids)``.

        Drives the user's RiboNN checkout: validates the inputs, verifies weights,
        writes a temporary input table
        (``tx_id, utr5_sequence, cds_sequence, utr3_sequence``) using the fixed UTR
        context, and calls the repo's
        ``predict_using_nested_cross_validation_models``. Heavy deps and the repo's
        ``src`` are imported here, never at module load.

        The *raw* table is returned (one row per input per outer fold) so both the
        fold-averaged summary (:meth:`_predict_te`) and the fold-resolved view
        (:meth:`predict_folds`) are served by one invocation and one code path --
        there is no second scoring route to drift.
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
                predict,
                weights,
                str(in_path),
                self.species,
                run_df,
                self.top_k,
                batch_size=self.batch_size,
                num_workers=self.num_workers,
            )

        return out_df, tx_ids

    def _predict_te(self, dnas: list[str]) -> list[float]:
        """Return one fold-averaged mean-TE per input, in input order."""
        out_df, tx_ids = self._run_predict(dnas)
        return _reduce_te_by_tx_id(out_df, tx_ids, self.cell_types)

    def predict_folds(self, dnas: list[str]) -> list[RiboNNFoldPrediction]:
        """Return RiboNN's **per-outer-fold** predictions, not their mean.

        :meth:`score_many` averages RiboNN's ten outer folds, which is right for the
        novel designed sequences BT4 produces -- no fold ever saw them. It is wrong
        for a *natural* transcript, where nine of the ten folds trained on that
        transcript's own label, so the fold-averaged number is optimistic and cannot
        be compared against RiboNN's published held-out accuracy.

        This exposes the fold identity so a caller can keep only the fold that held a
        given transcript out. That is what makes the free adapter-validation check
        possible: score RiboNN's own published labels, keep the matching fold, and
        confirm the held-out r^2 lands near the published value while the other nine
        folds sit visibly higher. If they are indistinguishable, the fold semantics
        are wrong and every downstream number is uninterpretable.

        Args:
            dnas: Coding sequences over ``{A,C,G,T}`` (start/stop codons included).
                An empty list returns ``[]`` without touching the checkout.

        Returns:
            One :class:`RiboNNFoldPrediction` per ``(input, fold)`` pair, sorted by
            ``(index, fold)``. Values are in native CLR-residual units and carry **no
            calibration claim** -- this is a diagnostic surface, not a prediction API.
        """
        if not dnas:
            return []
        for dna in dnas:
            validate_dna(dna)
        out_df, tx_ids = self._run_predict(dnas)
        return _reduce_te_by_tx_id_and_fold(out_df, tx_ids, self.cell_types)

    def score_sequence(self, dna: str) -> ExpressionResult:
        """Return RiboNN's predicted mean TE for ``dna`` (CLR-residual units).

        A thin single-sequence wrapper over :meth:`score_many` (one RiboNN
        invocation), kept for the :class:`ExpressionPredictor` contract. When
        scoring a whole candidate set, call :meth:`score_many` directly so the
        large fixed per-call overhead is paid once, not once per sequence.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (its start/stop codons
                included, per RiboNN's input contract).

        Returns:
            An :class:`ExpressionResult`; ``score`` is the mean predicted TE across
            cell types in native CLR-residual units (larger is better), and
            ``calibrated`` mirrors :attr:`calibrated` (``False`` until the gate
            passes -- never read it as a validated expression prediction otherwise).
        """
        (result,) = self.score_many([dna])
        return result

    def score_many(self, dnas: list[str]) -> list[ExpressionResult]:
        """Score a whole candidate set in a **single** RiboNN invocation.

        The batched counterpart to :meth:`score_sequence`. RiboNN's cost is
        dominated by fixed per-call overhead -- weight hashing, model load, and
        its DataLoader worker spawn (heavy on Windows) -- almost all of it *per
        invocation*, not per sequence (CLAUDE.md §6;
        ``docs/DESIGN_expression_splice_flow.md`` step 1). Routing the whole set
        through one :meth:`_predict_te` call (one temporary TSV, one ``predict``
        invocation -- RiboNN's ``top_k``-model ensemble runs inside that single
        call) amortizes that overhead, so scoring a frontier costs roughly the
        wall-clock of scoring a single sequence. This reuses the existing batched
        ``_predict_te`` path; it does not add a second scoring route.

        Per-input validation is preserved: every sequence must be valid DNA and
        (inside ``_predict_te``) length-3N ending in a stop codon, and ``utr5`` /
        ``utr3`` must be non-empty -- the same guards as :meth:`score_sequence`.

        Args:
            dnas: Coding sequences over ``{A,C,G,T}`` (start/stop codons included).
                An empty list returns ``[]`` without touching the checkout.

        Returns:
            One :class:`ExpressionResult` per input, **in input order**. Each
            ``score`` is the mean predicted TE across cell types in native
            CLR-residual units (larger is better); ``calibrated`` mirrors
            :attr:`calibrated` (``False`` until the gate passes).
        """
        if not dnas:
            return []
        for dna in dnas:
            validate_dna(dna)
        units = self._units()
        return [
            ExpressionResult(
                score=te,
                model_name=self.name,
                calibrated=self.calibrated,
                units=units,
            )
            for te in self._predict_te(dnas)
        ]

    def _units(self) -> str:
        """Return the units label, naming the cell-type selection when explicit.

        The label travels with every score into reports and the manifest, so it must
        say *which* cell types were averaged -- "mean over human cell types" and
        "mean over HEK293T" are different numbers and must not share a label.
        """
        if self.cell_types:
            return (
                "RiboNN CLR-residual TE (mean over "
                + ", ".join(self.cell_types)
                + ")"
            )
        return f"RiboNN CLR-residual TE (mean over all {self.species} cell types)"

    def delta_logte(self, designed_dna: str, reference_dna: str) -> float:
        """Return ``TE(designed) - TE(reference)`` -- the CDS-attributable signal.

        Holding the UTRs fixed, this isolates what the CDS redesign changes: a
        positive value predicts higher expression, a negative value flags a CDS
        change RiboNN predicts will *reduce* expression (the maintainer's
        limiting-sequence framing). Analogous to
        :meth:`bt4.biomodels.splice.PangolinSplicePredictor.delta_splicing`. A
        thin single-design wrapper over :meth:`delta_logte_many`; batch a whole
        candidate set through that method to score the reference only once.
        """
        (delta,) = self.delta_logte_many([designed_dna], reference_dna)
        return delta

    def delta_logte_many(self, designed: list[str], reference: str) -> list[float]:
        """Return ``TE(designed_i) - TE(reference)`` for every design, batched.

        Because the UTRs are held fixed, the reference is one shared baseline: it
        is scored **once**, appended to the batch, not once per design. Every
        design *and* the reference go through a single RiboNN invocation (one
        :meth:`_predict_te` call), amortizing the fixed per-call overhead across
        the whole set -- so ranking a frontier by ΔTE costs roughly one RiboNN
        call, not one per candidate (CLAUDE.md §6;
        ``docs/DESIGN_expression_splice_flow.md`` step 1). The same per-input
        validation as :meth:`delta_logte` is preserved.

        Args:
            designed: Candidate coding sequences (each with start/stop codons). An
                empty list returns ``[]`` without scoring the reference.
            reference: The baseline CDS each design is compared against, scored once.

        Returns:
            One CDS-attributable Δ per design, **in input order** (a positive value
            predicts higher expression; a negative value flags a predicted
            reduction).
        """
        if not designed:
            return []
        # Validate the designs before the reference, so the single-design
        # ``delta_logte`` delegation reports errors in the same order it used to.
        for dna in designed:
            validate_dna(dna)
        validate_dna(reference)
        # Score every design plus the single shared reference in ONE invocation;
        # the reference is the last element, so it is scored exactly once.
        *designed_te, reference_te = self._predict_te([*designed, reference])
        return [te - reference_te for te in designed_te]


_TE_PREFIX = "predicted_TE_"


def _select_te_columns(out_df: Any, cell_types: tuple[str, ...]) -> list[str]:
    """Return the ``predicted_TE_*`` columns to average, honouring ``cell_types``.

    With ``cell_types`` empty, every per-cell-type column is used (RiboNN's own
    summary: 78 human / 68 mouse). With names given, only those are used -- and an
    unmatched name **raises** rather than being silently ignored, because quietly
    averaging the wrong set of tissues is precisely the kind of scope error that
    would show up later as an unexplained calibration failure.

    Raises:
        RuntimeError: If the table carries no ``predicted_TE_*`` columns at all.
        ValueError: If a requested cell type has no matching column.
    """
    available = [c for c in out_df.columns if str(c).startswith(_TE_PREFIX)]
    if not available:
        raise RuntimeError("RiboNN output had no predicted_TE_* columns")
    if not cell_types:
        return available
    by_name = {str(c)[len(_TE_PREFIX) :]: c for c in available}
    selected: list[str] = []
    missing: list[str] = []
    for name in cell_types:
        column = by_name.get(name)
        if column is None:
            missing.append(name)
        else:
            selected.append(column)
    if missing:
        raise ValueError(
            f"RiboNN has no output for cell type(s) {sorted(missing)}; "
            f"available: {sorted(by_name)}"
        )
    return selected


def _reduce_te_by_tx_id_and_fold(
    out_df: Any, tx_ids: list[str], cell_types: tuple[str, ...] = ()
) -> list[RiboNNFoldPrediction]:
    """Reduce RiboNN's table to one :class:`RiboNNFoldPrediction` per (input, fold).

    RiboNN emits one row per input per outer fold (each row already that fold's
    ``top_k``-model mean) tagged with a ``fold`` column. This averages the selected
    per-cell-type columns within each row and keeps the fold identity, so a caller
    can either average the folds (a novel sequence: no fold saw it) or keep only the
    holdout fold (a natural transcript: nine folds trained on its label).

    Rows are returned sorted by ``(index, fold)``. A table without a ``fold`` column
    is treated as a single unlabelled fold (``-1``), which keeps synthetic
    single-fold tables usable.

    Raises:
        ValueError: If an input has no prediction row -- RiboNN's length caps
            (5'UTR > 1381 or CDS+3'UTR > 11937 nt) drop such rows.
    """
    te_cols = _select_te_columns(out_df, cell_types)
    frame = out_df.assign(_te=out_df[te_cols].mean(axis=1))
    position = {tx_id: i for i, tx_id in enumerate(tx_ids)}
    row_tx_ids = [str(value) for value in frame["tx_id"]]
    row_tes = [float(value) for value in frame["_te"]]
    # A table with no ``fold`` column is one unlabelled fold; that keeps synthetic
    # single-fold tables (and any future single-fold upstream mode) usable.
    row_folds = (
        [int(value) for value in frame["fold"]]
        if "fold" in frame.columns
        else [-1] * len(row_tx_ids)
    )

    grouped: dict[tuple[int, int], list[float]] = {}
    for tx_id, fold, te in zip(row_tx_ids, row_folds, row_tes, strict=True):
        index = position.get(tx_id)
        if index is None:
            continue  # a row for something we did not ask about; ignore it
        grouped.setdefault((index, fold), []).append(te)

    seen = {index for index, _ in grouped}
    for i, _tx_id in enumerate(tx_ids):
        if i not in seen:
            raise ValueError(
                f"RiboNN produced no prediction for input {i}; its CDS likely "
                "exceeds RiboNN's length cap (CDS+3'UTR <= 11937 nt)"
            )
    return [
        RiboNNFoldPrediction(index=index, fold=fold, te=math.fsum(values) / len(values))
        for (index, fold), values in sorted(grouped.items())
    ]


def _reduce_te_by_tx_id(
    out_df: Any, tx_ids: list[str], cell_types: tuple[str, ...] = ()
) -> list[float]:
    """Reduce RiboNN's raw prediction table to one mean-TE per input, in input order.

    The fold-averaging summary: :func:`_reduce_te_by_tx_id_and_fold` then a mean over
    folds per input. Correct for the *novel designed* sequences BT4 produces, since no
    ensemble fold has seen them; see :class:`RiboNNFoldPrediction` for why it is wrong
    for a natural transcript.
    """
    per_fold = _reduce_te_by_tx_id_and_fold(out_df, tx_ids, cell_types)
    by_index: dict[int, list[float]] = {}
    for record in per_fold:
        by_index.setdefault(record.index, []).append(record.te)
    return [
        math.fsum(by_index[i]) / len(by_index[i]) for i in range(len(tx_ids))
    ]


def _run_predict_with_models_layout(
    predict: Any,
    weights_dir: Path,
    input_path: str,
    species: str,
    run_df: Any,
    top_k: int,
    *,
    batch_size: int,
    num_workers: int,
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
    # ``predict_using_nested_cross_validation_models`` DOES expose ``batch_size``
    # (default 1024) and ``num_workers`` (default 4) -- an earlier version of this
    # comment asserted it did not, and that was wrong. Both are forwarded, because
    # both upstream defaults are actively hostile here:
    #
    # * ``num_workers=4`` spawns dataloader workers that re-import the module. This
    #   adapter has just mutated ``sys.path`` and is about to ``chdir`` into a
    #   temporary directory, and a *spawned* worker (Windows, macOS) inherits
    #   neither -- so the workers hang or fail. RiboNN also rebuilds the dataloader
    #   once per ensemble member (up to 50 times), paying that cost every time.
    # * ``batch_size=1024`` allocates 1024 fixed-width ``(channels, 13318)`` float32
    #   tensors at once -- hundreds of MB per batch before prefetch -- which OOMs an
    #   ordinary CPU box.
    #
    # Neither knob can change a score: RiboNN pads to a fixed width rather than to
    # the longest member of a batch, and its predict dataloader is constructed with
    # ``shuffle=False`` (``reorder = stage == "train"``), so batch composition
    # affects neither encoding nor row order. They are memory/throughput only.
    import os as _os
    import tempfile

    prev = _os.getcwd()
    if weights_dir.name == "models":
        run_root = weights_dir.parent
        try:
            _os.chdir(run_root)
            return predict(
                input_path=input_path, species=species, run_df=run_df,
                top_k_models_to_use=top_k, batch_size=batch_size,
                num_workers=num_workers,
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
                top_k_models_to_use=top_k, batch_size=batch_size,
                num_workers=num_workers,
            )
        finally:
            _os.chdir(prev)
