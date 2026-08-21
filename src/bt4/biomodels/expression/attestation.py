"""License-clean attestations that a learned expression head passed its gate.

The wrapped RiboNN head may report ``calibrated is True`` only after an acceptance gate
passes **on data from the regime it serves** -- synonymous CDS variants of one protein
under a fixed UTR (CLAUDE.md §6/§8/§10.6). Today promotion is a bare
``dataclasses.replace(model, fidelity_verified=True)``, which is exactly the
assign-it-by-hand move the constitution forbids. This module is the seam that replaces
it, mirroring :mod:`bt4.biomodels.splice.attestation`.

**What may be recorded.** RiboNN's weights are Sanofi *non-commercial*, so its raw
per-sequence outputs are licence-encumbered and must never enter MIT-licensed BT4. The
*facts a passing gate produces* -- a boolean, rank and coverage scalars, the pinned
weight SHA-256s, the panel's content hash, the tool version -- are not model outputs and
are safe to commit. :data:`_ALLOWED_FIELDS` pins that shape and a test asserts it never
grows a raw-score field.

**What makes this different from the splice attestation.** A splice attestation asserts
*fidelity* -- that the adapter reproduces a published model bit-for-bit -- and is
therefore gated on a tolerance floor. An expression attestation asserts something
strictly harder: that the head is **useful for a job**, on a panel, against baselines. So
it is gated on four floors that a fidelity claim does not need, and it carries its
**scope**:

* the run must have been **within-group** -- a pooled run credits between-protein skill,
  which is not the regime BT4 deploys in, so a pooled result cannot certify anything
  here no matter how good it looks;
* the head must have **beaten every baseline** -- BT4 already computes CAI inside the
  optimizer loop, so a head that cannot beat it has earned nothing;
* the interval must be **informative** (width below the label IQR), because split
  conformal is valid for any score function and a constant predictor passes coverage;
* ``min_spearman`` must clear :data:`MIN_ATTESTATION_SPEARMAN`, so a self-serving
  threshold cannot self-certify.

**Scope is part of the claim, not a footnote.** An attestation earned on HEK293T ribosome
load, human weights, one UTR context, certifies *that*. :func:`verified_predictor`
refuses a predictor whose species or cell-type selection differs, because the mean of 78
tissues is a different quantity from one cell line and a claim about one is not a claim
about the other.

**The scope is taken from the run, not from the caller's word for it.** An earlier
version accepted ``species`` and ``cell_types`` as free text and copied them straight
into the record, so a gate run averaging all 78 cell types could be filed as a HEK293T
result and every later check would accept it. :func:`attest_expression` now *derives*
the scope from the comparison's own :class:`GateScope
<bt4.pipeline.expression_gate.GateScope>` -- what the gate actually scored -- and treats
anything the caller declares as a **cross-check that refuses on mismatch**, never an
override. Where the panel itself declares a fact (its ``species`` / ``cell_type`` /
``readout`` columns) the declaration is additionally checked against the panel bytes;
:attr:`ExpressionAttestation.verified_against_panel` records exactly which fields got
that second check, so a reader can tell a verified scope from a merely-declared one.

**Two more things the scope binds, because they change the number.**
:attr:`~ExpressionAttestation.top_k` (how many cross-validation runs were ensembled) and
:attr:`~ExpressionAttestation.utr_context_sha256` (the UTR contexts the panel was
measured in) are part of what produced the score, so a head configured differently is
not the head that was gated and :func:`verified_predictor` refuses it. ``batch_size`` and
``num_workers`` are deliberately **not** bound: RiboNN pads to a fixed width and does not
shuffle when predicting, so neither can change a score, and binding them would be false
precision.

**On the UTR hash.** ``utr_context_sha256`` is a plain SHA-256 of each ``(utr5, utr3)``
pair, so a *short* UTR is recoverable from it by brute force. That is stated rather than
papered over: committing an attestation publishes the panel's UTR context. It is not
licensed material (RiboNN's weights and outputs are; a panel's UTRs are the maintainer's
own data, usually from the paper the panel came from) -- but a maintainer whose panel is
not public should keep the attestation local (see
:mod:`bt4.biomodels.expression.attestations`) rather than commit it.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any, Protocol, cast

__all__ = [
    "MAX_ATTESTATION_COVERAGE_TOLERANCE",
    "MAX_ATTESTATION_WIDTH_OVER_IQR",
    "MIN_ATTESTATION_SPEARMAN",
    "ExpressionAttestation",
    "ExpressionAttestationError",
    "attest_expression",
    "load_expression_attestation",
    "utr_context_sha256",
    "verified_predictor",
]

_BACKENDS = ("ribonn",)

# Floors. These are not thresholds for a *run* -- a maintainer sets those, and records
# them -- but the loosest claim the attestation layer will honour. Their job is to stop a
# gate run configured to pass from certifying anything.
MIN_ATTESTATION_SPEARMAN = 0.30
"""Loosest ``min_spearman`` an attestation may have been gated at.

Not a community standard -- there is none -- but a pre-committed floor, so a run
configured with ``min_spearman=0.05`` cannot produce a certificate."""

MAX_ATTESTATION_COVERAGE_TOLERANCE = 0.10
"""Loosest coverage tolerance an attestation may claim. A gate allowed to miss its
target by more than this is not making an uncertainty claim worth recording."""

MAX_ATTESTATION_WIDTH_OVER_IQR = 1.0
"""An interval as wide as the spread of the labels says nothing, however valid its
coverage. This is the vacuity floor."""

_ALLOWED_FIELDS = frozenset(
    {
        "backend",
        "species",
        "cell_types",
        "readout",
        "within_group",
        "recalibrate",
        "passed",
        "spearman",
        "spearman_ci_low",
        "n_test",
        "n_groups_ranked",
        "target_coverage",
        "empirical_coverage",
        "width_over_iqr",
        "min_spearman",
        "coverage_tolerance",
        "best_baseline",
        "best_baseline_spearman",
        "panel_sha256",
        "weight_sha256",
        "top_k",
        "utr_context_sha256",
        "verified_against_panel",
        "scoring_source",
        "bt4_version",
        "schema_version",
    }
)

# Domain separation, so a UTR hash can never collide with any other SHA-256 in this
# record and a reader can tell what was hashed from the tag alone.
_UTR_HASH_TAG = "bt4-expression-utr-context-v1"


def utr_context_sha256(utr5: str, utr3: str) -> str:
    """Return the stable content hash of one ``(utr5, utr3)`` context.

    The gate records one of these per distinct context in the panel, and
    :func:`verified_predictor` recomputes it from the predictor's own fixed UTRs, so an
    attestation earned under one transcript context cannot promote a head configured for
    another. Case-insensitive and whitespace-free, matching how both the panel reader and
    the adapter normalise a flank.

    See the module docstring for why this hash is not a secret: a short UTR is
    brute-forceable from it.
    """
    payload = "\t".join(
        (_UTR_HASH_TAG, "".join(utr5.split()).upper(), "".join(utr3.split()).upper())
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class ExpressionAttestationError(ValueError):
    """Raised when an attestation is missing, malformed, out of scope, or unearned.

    Always a refusal, never a silent downgrade: a head is promoted only against an
    attestation that structurally matches the adapter's pinned weights *and* its scope.
    """


class _ComparisonLike(Protocol):
    """The scalars an attestation needs from a gate comparison (duck-typed).

    Declared as a Protocol so this module never imports :mod:`bt4.pipeline` -- the
    layering forbids it (CLAUDE.md §3), and duck-typing keeps the attestation buildable
    from any equivalently-shaped result.
    """

    @property
    def panel_hash(self) -> str: ...

    @property
    def head(self) -> Any: ...

    @property
    def best_baseline(self) -> str: ...

    @property
    def best_baseline_spearman(self) -> float: ...

    @property
    def beats_every_baseline(self) -> bool: ...

    @property
    def interval_is_informative(self) -> bool: ...

    @property
    def promotable(self) -> bool: ...

    @property
    def scope(self) -> Any: ...


class _ScopeLike(Protocol):
    """How a gate run was actually configured (duck-typed, see :class:`_ComparisonLike`).

    This is the record that makes a declared scope checkable instead of decorative.
    """

    @property
    def species(self) -> str: ...

    @property
    def cell_types(self) -> tuple[str, ...]: ...

    @property
    def top_k(self) -> int: ...

    @property
    def readout(self) -> str: ...

    @property
    def utr_context_sha256(self) -> tuple[str, ...]: ...

    @property
    def panel_species(self) -> tuple[str, ...]: ...

    @property
    def panel_cell_types(self) -> tuple[str, ...]: ...

    @property
    def panel_readouts(self) -> tuple[str, ...]: ...

    @property
    def scoring_source(self) -> str: ...


@dataclass(frozen=True, slots=True)
class ExpressionAttestation:
    """A licence-clean record that an expression head passed its gate, and for what.

    Holds only derived scalars plus public content hashes -- never a raw model score.
    Timestamp-free, so :meth:`content_hash` is a stable provenance stamp (invariant #7).

    Attributes:
        backend: Which head this attests (``"ribonn"``).
        species: The weight set the gate ran against.
        cell_types: The cell-type selection scored, sorted. Empty means "all of them",
            which is a *different* quantity from any single cell line.
        readout: What the panel measured (e.g. ``"mean_ribosome_load"``), so the claim
            names the assay it was earned on.
        within_group: Always ``True`` for a stored attestation -- a pooled run cannot
            certify BT4's regime and :func:`attest_expression` refuses to record one.
        recalibrate: Whether the affine link was fitted on the calibration fold.
        passed: Whether the gate passed (always ``True`` when stored).
        spearman: The primary rank metric on the head's raw predictions.
        spearman_ci_low: Lower bound of the cluster-bootstrap CI -- the number that had
            to exceed every baseline.
        n_test: Test-fold size.
        n_groups_ranked: Test groups that contributed a rank correlation. This, not
            ``n_test``, is the effective sample size for a cross-protein claim.
        target_coverage: The conformal level requested.
        empirical_coverage: The realized coverage.
        width_over_iqr: Median interval width over the label IQR (the vacuity check).
        min_spearman: The threshold the run was gated at, recorded so a reader can see
            what was pre-committed.
        coverage_tolerance: The coverage tolerance the run allowed.
        best_baseline: The strongest baseline on that panel.
        best_baseline_spearman: Its primary metric -- what the head had to beat.
        panel_sha256: The panel's content hash, binding the claim to exact bytes.
        weight_sha256: The pinned ``{weight_path: sha256}`` for this species, as a
            sorted tuple of pairs. Public content hashes, not licensed outputs.
        top_k: How many cross-validation runs the gated head ensembled. Part of the
            model that produced the number, so a differently-configured head is a
            different head.
        utr_context_sha256: One :func:`utr_context_sha256` per distinct ``(utr5, utr3)``
            context in the panel, sorted. A predictor is promoted only when its own fixed
            UTR context is one of these -- an attestation earned under one transcript
            context does not certify another. (Not a secret; see the module docstring.)
        verified_against_panel: Which scope fields were cross-checked against the panel's
            own columns rather than merely declared -- a subset of ``("cell_types",
            "readout", "species")``, sorted. A field absent here was taken on the
            maintainer's word because the panel did not declare it.
        scoring_source: ``"gate"`` when the gate invoked the backend itself, or
            ``"caller_supplied"`` when it was handed pre-computed values. Recorded
            because that is exactly the step at which the link between the named backend
            and the numbers stops being mechanical.
        bt4_version: The BT4 version that produced the attestation.
        schema_version: Attestation schema version.
    """

    backend: str
    species: str
    cell_types: tuple[str, ...]
    readout: str
    within_group: bool
    recalibrate: bool
    passed: bool
    spearman: float
    spearman_ci_low: float
    n_test: int
    n_groups_ranked: int
    target_coverage: float
    empirical_coverage: float
    width_over_iqr: float
    min_spearman: float
    coverage_tolerance: float
    best_baseline: str
    best_baseline_spearman: float
    panel_sha256: str
    weight_sha256: tuple[tuple[str, str], ...]
    top_k: int
    utr_context_sha256: tuple[str, ...]
    verified_against_panel: tuple[str, ...]
    scoring_source: str
    bt4_version: str
    schema_version: int = 2

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready mapping (weight and cell-type tuples become lists)."""
        return {
            "backend": self.backend,
            "species": self.species,
            "cell_types": list(self.cell_types),
            "readout": self.readout,
            "within_group": self.within_group,
            "recalibrate": self.recalibrate,
            "passed": self.passed,
            "spearman": self.spearman,
            "spearman_ci_low": self.spearman_ci_low,
            "n_test": self.n_test,
            "n_groups_ranked": self.n_groups_ranked,
            "target_coverage": self.target_coverage,
            "empirical_coverage": self.empirical_coverage,
            "width_over_iqr": self.width_over_iqr,
            "min_spearman": self.min_spearman,
            "coverage_tolerance": self.coverage_tolerance,
            "best_baseline": self.best_baseline,
            "best_baseline_spearman": self.best_baseline_spearman,
            "panel_sha256": self.panel_sha256,
            "weight_sha256": [list(pair) for pair in self.weight_sha256],
            "top_k": self.top_k,
            "utr_context_sha256": list(self.utr_context_sha256),
            "verified_against_panel": list(self.verified_against_panel),
            "scoring_source": self.scoring_source,
            "bt4_version": self.bt4_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ExpressionAttestation:
        """Rebuild an attestation, refusing any unexpected field.

        The refusal is the point: an extra key is how a raw per-sequence score array
        would be smuggled into a licence-clean record.

        Raises:
            ExpressionAttestationError: On an unexpected or missing field.
        """
        extra = set(data) - _ALLOWED_FIELDS
        if extra:
            raise ExpressionAttestationError(
                f"unexpected attestation field(s): {sorted(extra)}"
            )
        missing = _ALLOWED_FIELDS - set(data) - {"schema_version"}
        if missing:
            raise ExpressionAttestationError(
                f"attestation is missing field(s): {sorted(missing)}"
            )
        return cls(
            backend=str(data["backend"]),
            species=str(data["species"]),
            cell_types=tuple(str(name) for name in data["cell_types"]),
            readout=str(data["readout"]),
            within_group=bool(data["within_group"]),
            recalibrate=bool(data["recalibrate"]),
            passed=bool(data["passed"]),
            spearman=float(data["spearman"]),
            spearman_ci_low=float(data["spearman_ci_low"]),
            n_test=int(data["n_test"]),
            n_groups_ranked=int(data["n_groups_ranked"]),
            target_coverage=float(data["target_coverage"]),
            empirical_coverage=float(data["empirical_coverage"]),
            width_over_iqr=float(data["width_over_iqr"]),
            min_spearman=float(data["min_spearman"]),
            coverage_tolerance=float(data["coverage_tolerance"]),
            best_baseline=str(data["best_baseline"]),
            best_baseline_spearman=float(data["best_baseline_spearman"]),
            panel_sha256=str(data["panel_sha256"]),
            weight_sha256=tuple(
                (str(pair[0]), str(pair[1])) for pair in data["weight_sha256"]
            ),
            top_k=int(data["top_k"]),
            utr_context_sha256=tuple(str(h) for h in data["utr_context_sha256"]),
            verified_against_panel=tuple(
                str(name) for name in data["verified_against_panel"]
            ),
            scoring_source=str(data["scoring_source"]),
            bt4_version=str(data["bt4_version"]),
            schema_version=int(data.get("schema_version", 2)),
        )

    def content_hash(self) -> str:
        """Return a stable SHA-256 over the attestation's content (no wall-clock)."""
        payload = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _pinned_weights(species: str) -> dict[str, str]:
    """Return the adapter's pinned weight hashes for one species.

    Filtered by species prefix exactly as
    :meth:`~bt4.biomodels.expression.ribonn.RiboNNExpressionModel._verify_weights` does,
    so an attestation is bound to the same 90 files the adapter will hash-verify before
    loading.
    """
    from bt4.biomodels.expression.ribonn import PINNED_WEIGHT_SHA256

    prefix = f"{species}/"
    return {k: v for k, v in PINNED_WEIGHT_SHA256.items() if k.startswith(prefix)}


def _require_scope(comparison: _ComparisonLike) -> _ScopeLike:
    """Return the comparison's run scope, refusing a comparison that carries none.

    A comparison without a scope cannot say how it was scored, so nothing it reports can
    be bound to a configuration. Refusing is the point: the whole purpose of this layer
    is that the recorded scope is the run's, not the caller's.
    """
    scope: Any = getattr(comparison, "scope", None)
    if scope is None:
        raise ExpressionAttestationError(
            "this gate comparison carries no scope, so how it was scored cannot be "
            "checked. Re-run it through bt4.pipeline.expression_gate.run_panel_gate."
        )
    return cast(_ScopeLike, scope)


def attest_expression(
    comparison: _ComparisonLike,
    *,
    backend: str = "ribonn",
    species: str | None = None,
    cell_types: tuple[str, ...] | None = None,
    readout: str | None = None,
    bt4_version: str,
) -> ExpressionAttestation:
    """Record a passing gate comparison as a licence-clean attestation.

    Refuses -- never downgrades -- unless the run actually earned a claim:

    * the gate ``passed`` and the comparison is ``promotable``;
    * the run was **within-group**, because a pooled run credits between-protein skill
      and cannot certify the regime BT4 deploys in;
    * the head **beat every baseline**;
    * the interval is informative (``width_over_iqr`` below
      :data:`MAX_ATTESTATION_WIDTH_OVER_IQR`);
    * the run's ``min_spearman`` clears :data:`MIN_ATTESTATION_SPEARMAN` and its
      coverage tolerance clears :data:`MAX_ATTESTATION_COVERAGE_TOLERANCE`, so a
      threshold set to pass cannot self-certify.

    **The scope comes from the run.** ``species``, ``cell_types`` and ``readout`` are
    taken from ``comparison.scope`` -- what the gate actually scored. Passing them
    explicitly asks for a **cross-check**: a value that disagrees with the run is a
    refusal, not an override. This is the hole that made the previous version's scope
    decorative: a run averaging all 78 cell types could be filed as a HEK293T result and
    nothing downstream could tell.

    Where the panel declares the same fact in its own ``species`` / ``cell_type`` /
    ``readout`` columns, that is checked too, and the field is listed in
    :attr:`~ExpressionAttestation.verified_against_panel`. A field the panel does not
    declare is recorded as merely declared -- an honest gap, not a silent pass.

    Args:
        comparison: A gate comparison (see :mod:`bt4.pipeline.expression_gate`).
        backend: Which head is attested.
        species: Optional cross-check against the species the gate ran with.
        cell_types: Optional cross-check against the cell-type selection it scored.
        readout: What the panel measured. Optional when the panel's own ``readout``
            column declares exactly one; required otherwise, because a number must name
            the question it answers.
        bt4_version: The BT4 version producing the record.

    Returns:
        The :class:`ExpressionAttestation`.

    Raises:
        ExpressionAttestationError: On an unknown backend, any unearned claim, or any
            declared scope value that disagrees with the run or the panel.
    """
    if backend not in _BACKENDS:
        raise ExpressionAttestationError(
            f"unknown expression backend {backend!r}; known: {list(_BACKENDS)}"
        )
    scope = _require_scope(comparison)
    head = comparison.head
    if not head.passed:
        raise ExpressionAttestationError(
            "refusing to attest a failing gate"
        )
    if not head.within_group:
        raise ExpressionAttestationError(
            "refusing to attest a POOLED run: pooled scoring credits between-protein "
            "skill, which is not the regime BT4 deploys in. Re-run with within_group."
        )
    if not comparison.beats_every_baseline:
        raise ExpressionAttestationError(
            f"refusing to attest a head that does not beat every baseline (best was "
            f"{comparison.best_baseline!r} at {comparison.best_baseline_spearman:.3f}, "
            f"head CI lower bound {head.spearman_ci_low:.3f})"
        )
    if head.width_over_iqr >= MAX_ATTESTATION_WIDTH_OVER_IQR:
        raise ExpressionAttestationError(
            f"refusing to attest a vacuous interval: width/IQR "
            f"{head.width_over_iqr:.3f} >= {MAX_ATTESTATION_WIDTH_OVER_IQR}"
        )
    if head.min_spearman < MIN_ATTESTATION_SPEARMAN:
        raise ExpressionAttestationError(
            f"refusing to attest a run gated at min_spearman={head.min_spearman}, "
            f"below the floor {MIN_ATTESTATION_SPEARMAN}"
        )
    if head.coverage_tolerance > MAX_ATTESTATION_COVERAGE_TOLERANCE:
        raise ExpressionAttestationError(
            f"refusing to attest a run allowed to miss coverage by "
            f"{head.coverage_tolerance}, looser than the floor "
            f"{MAX_ATTESTATION_COVERAGE_TOLERANCE}"
        )

    run_species, run_cell_types, run_readout, verified = _resolved_scope(
        scope, species=species, cell_types=cell_types, readout=readout
    )

    pinned = _pinned_weights(run_species)
    if not pinned:
        raise ExpressionAttestationError(
            f"no pinned weights for species {run_species!r}; cannot bind an attestation"
        )
    if not scope.utr_context_sha256:
        raise ExpressionAttestationError(
            "the gate run recorded no UTR context, so the claim cannot be bound to the "
            "transcript context it was earned in"
        )
    return ExpressionAttestation(
        backend=backend,
        species=run_species,
        cell_types=run_cell_types,
        readout=run_readout,
        within_group=head.within_group,
        recalibrate=head.recalibrate,
        passed=True,
        spearman=head.spearman,
        spearman_ci_low=head.spearman_ci_low,
        n_test=head.n_test,
        n_groups_ranked=head.n_groups_ranked,
        target_coverage=head.target_coverage,
        empirical_coverage=head.empirical_coverage,
        width_over_iqr=head.width_over_iqr,
        min_spearman=head.min_spearman,
        coverage_tolerance=head.coverage_tolerance,
        best_baseline=comparison.best_baseline,
        best_baseline_spearman=comparison.best_baseline_spearman,
        panel_sha256=comparison.panel_hash,
        weight_sha256=tuple(sorted((str(k), str(v)) for k, v in pinned.items())),
        top_k=scope.top_k,
        utr_context_sha256=tuple(sorted(scope.utr_context_sha256)),
        verified_against_panel=verified,
        scoring_source=scope.scoring_source,
        bt4_version=bt4_version,
    )


def _resolved_scope(
    scope: _ScopeLike,
    *,
    species: str | None,
    cell_types: tuple[str, ...] | None,
    readout: str | None,
) -> tuple[str, tuple[str, ...], str, tuple[str, ...]]:
    """Resolve the scope from the run, cross-checking anything the caller declared.

    Returns ``(species, cell_types, readout, verified_against_panel)``. Every disagreement
    -- caller vs run, or run vs the panel's own columns -- is a refusal. ``verified``
    lists the fields the panel independently confirmed, so a record can distinguish a
    checked scope from one taken on trust.
    """
    verified: list[str] = []

    run_species = scope.species
    if species is not None and species != run_species:
        raise ExpressionAttestationError(
            f"declared species {species!r} but the gate ran against {run_species!r}; "
            "the scope is the run's, not the caller's"
        )
    declared_by_panel = tuple(sorted({name for name in scope.panel_species if name}))
    if declared_by_panel:
        if declared_by_panel != (run_species,):
            raise ExpressionAttestationError(
                f"the panel declares species {list(declared_by_panel)} but the gate ran "
                f"against {run_species!r}"
            )
        verified.append("species")

    run_cell_types = tuple(sorted(scope.cell_types))
    if cell_types is not None and tuple(sorted(cell_types)) != run_cell_types:
        shown = list(run_cell_types) or "every cell type (no selection)"
        raise ExpressionAttestationError(
            f"declared cell types {sorted(cell_types)} but the gate scored {shown}. "
            "Averaging a different set of cell types is a different quantity, so the "
            "declaration is refused rather than recorded."
        )
    panel_cells = tuple(sorted({name for name in scope.panel_cell_types if name}))
    if panel_cells:
        if panel_cells != run_cell_types:
            shown = list(run_cell_types) or "every cell type (no selection)"
            raise ExpressionAttestationError(
                f"the panel was measured in {list(panel_cells)} but the gate scored "
                f"{shown}"
            )
        verified.append("cell_types")

    run_readout = readout if readout is not None else scope.readout
    if not run_readout:
        raise ExpressionAttestationError(
            "no readout: the panel declares none (or declares several) and none was "
            "given. A measured number must name the question it answers."
        )
    panel_readouts = tuple(sorted({name for name in scope.panel_readouts if name}))
    if panel_readouts:
        if run_readout not in panel_readouts:
            raise ExpressionAttestationError(
                f"declared readout {run_readout!r} is not one the panel measures "
                f"({list(panel_readouts)})"
            )
        verified.append("readout")

    return run_species, run_cell_types, run_readout, tuple(sorted(verified))


def load_expression_attestation(path: str | Path) -> ExpressionAttestation:
    """Load an attestation from JSON on disk.

    Raises:
        ExpressionAttestationError: If the file is malformed or carries an unexpected
            field.
        OSError: If the file cannot be read.
    """
    with Path(path).open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ExpressionAttestationError("attestation file must contain a JSON object")
    return ExpressionAttestation.from_dict(data)


def verified_predictor(predictor: Any, attestation: ExpressionAttestation) -> Any:
    """Return ``predictor`` promoted to ``calibrated=True`` iff ``attestation`` matches.

    **The single seam that flips an expression head.** It refuses unless all hold:

    * the attestation ``passed`` and was earned **within-group**;
    * its ``backend`` matches the predictor's type;
    * its ``species`` matches the predictor's -- human and mouse are different models;
    * its ``cell_types`` match the predictor's **exactly** -- an attestation earned on
      one cell line does not certify a head averaging all 78, because those are
      different quantities;
    * its ``top_k`` matches -- a differently-sized ensemble is a different model, so it
      is not the model that was gated;
    * the predictor's own fixed ``(utr5, utr3)`` context is one the attestation covers --
      the gate measured ranking *inside* a transcript context, and nothing was measured
      about another one;
    * its floors hold (in case the file was hand-edited after
      :func:`attest_expression` built it);
    * its ``weight_sha256`` exactly matches the adapter's pins for that species, so the
      claim is bound to the same weights the adapter will hash-verify before loading.

    ``batch_size`` and ``num_workers`` are deliberately **not** checked: RiboNN pads to a
    fixed width and builds its predict dataloader with ``shuffle=False``, so neither can
    change a score. Binding them would refuse a head that is provably the gated one.

    Scores are still computed live from the user's own weights; the attestation carries
    none.

    Args:
        predictor: A :class:`~bt4.biomodels.expression.RiboNNExpressionModel`.
        attestation: A passing attestation for that head and scope.

    Returns:
        The promoted predictor (``calibrated is True``), carrying the attestation's
        :meth:`~ExpressionAttestation.content_hash` so a run stamps which claim
        authorized it. Nothing about what the head *computes* changes.

    Raises:
        ExpressionAttestationError: On any mismatch, unearned claim, or a
            non-attestable predictor type.
    """
    import dataclasses

    from bt4.biomodels.expression.ribonn import RiboNNExpressionModel

    if not isinstance(predictor, RiboNNExpressionModel):
        raise ExpressionAttestationError(
            f"{type(predictor).__name__} is not an attestable expression head"
        )
    if attestation.backend != "ribonn":
        raise ExpressionAttestationError(
            f"attestation is for {attestation.backend!r}, predictor is 'ribonn'"
        )
    if not attestation.passed:
        raise ExpressionAttestationError("attestation did not pass; refusing to calibrate")
    if not attestation.within_group:
        raise ExpressionAttestationError(
            "attestation was earned on a POOLED run, which cannot certify BT4's "
            "within-protein regime; refusing to calibrate"
        )
    if attestation.species != predictor.species:
        raise ExpressionAttestationError(
            f"attestation covers species {attestation.species!r}, predictor is "
            f"{predictor.species!r}"
        )
    if attestation.cell_types != tuple(sorted(predictor.cell_types)):
        raise ExpressionAttestationError(
            f"attestation covers cell types {list(attestation.cell_types)}, predictor "
            f"scores {sorted(predictor.cell_types)}. Averaging a different set of cell "
            "types is a different quantity, so the claim does not transfer."
        )
    if attestation.min_spearman < MIN_ATTESTATION_SPEARMAN:
        raise ExpressionAttestationError(
            f"attestation was gated at min_spearman={attestation.min_spearman}, below "
            f"the floor {MIN_ATTESTATION_SPEARMAN}"
        )
    if attestation.coverage_tolerance > MAX_ATTESTATION_COVERAGE_TOLERANCE:
        raise ExpressionAttestationError(
            f"attestation coverage tolerance {attestation.coverage_tolerance} is looser "
            f"than the floor {MAX_ATTESTATION_COVERAGE_TOLERANCE}"
        )
    if attestation.width_over_iqr >= MAX_ATTESTATION_WIDTH_OVER_IQR:
        raise ExpressionAttestationError(
            f"attestation interval is vacuous: width/IQR {attestation.width_over_iqr} "
            f">= {MAX_ATTESTATION_WIDTH_OVER_IQR}"
        )
    if attestation.top_k != predictor.top_k:
        raise ExpressionAttestationError(
            f"attestation was earned at top_k={attestation.top_k}, predictor ensembles "
            f"top_k={predictor.top_k}. A different ensemble size is a different model."
        )
    context = utr_context_sha256(predictor.utr5, predictor.utr3)
    if context not in attestation.utr_context_sha256:
        raise ExpressionAttestationError(
            "the predictor's fixed 5'/3' UTR context is not one this attestation covers "
            f"(it covers {len(attestation.utr_context_sha256)} context(s), none matching "
            f"{context[:16]}...). The gate measured ranking inside the panel's own "
            "transcript context; nothing was measured about another one."
        )
    pinned = tuple(
        sorted((str(k), str(v)) for k, v in _pinned_weights(predictor.species).items())
    )
    if attestation.weight_sha256 != pinned:
        raise ExpressionAttestationError(
            "attestation weight hashes do not match this adapter's pinned weights; "
            "it was earned against different weights and is not trustworthy here"
        )
    return dataclasses.replace(
        predictor,
        fidelity_verified=True,
        # Record WHICH claim authorized this, so a manifest can tell two promotions
        # apart (invariant #9). The attestation is timestamp-free, so this is stable.
        attestation_sha256=attestation.content_hash(),
    )


# Structural guard: the dataclass shape IS the licence-clean contract, so a drift
# between it and _ALLOWED_FIELDS must fail at import rather than at review time.
assert {f.name for f in fields(ExpressionAttestation)} == _ALLOWED_FIELDS, (
    "ExpressionAttestation fields drifted from _ALLOWED_FIELDS; a new field must be "
    "reviewed for licence-cleanliness (never a raw per-sequence model score)"
)
