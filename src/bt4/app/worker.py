"""Background engine workers for BT4 Studio.

Every long-running engine call the desktop app makes runs here, on a
``QThread``, so the window never blocks. Each worker wraps one :mod:`bt4.api`
entry point, communicates purely through Qt signals, and never touches widgets.

They share :class:`_EngineWorker`, which owns the signal trio and the
never-raise contract: any engine error is caught and reported (as the exception
object) through :attr:`~_EngineWorker.failed`, so bad input or a missing
optional backend can never crash the UI thread and the window can translate the
failure into plain language keyed on its type.
"""

from __future__ import annotations

from PySide6 import QtCore

from bt4 import api

__all__ = [
    "CandidatesResult",
    "CandidatesWorker",
    "CrossCheckWorker",
    "LibraryWorker",
    "OptimizeWorker",
]


class CandidatesResult:
    """The bundle a :class:`CandidatesWorker` delivers on success.

    Groups the step-3 candidate set with its optional step-4 splice audit so the
    window receives both in one ``finished`` signal. ``audit`` is ``None`` only if
    the candidate set was empty (the audit is skipped rather than raising).
    """

    __slots__ = ("audit", "candidate_set")

    def __init__(
        self,
        candidate_set: api.CandidateSet,
        audit: api.SpliceAuditReport | None,
    ) -> None:
        self.candidate_set = candidate_set
        self.audit = audit


class _EngineWorker(QtCore.QObject):
    """Base for the app's background workers: signals plus a never-raise ``run``.

    Subclasses implement :meth:`compute` (the synchronous engine call) and
    inherit the signal trio and the error contract. ``compute`` is deliberately
    public and thread-agnostic so tests -- and any synchronous caller -- can drive
    the engine without a thread or an event loop.

    Signals:
        finished(object): Emitted with the computed result on success.
        failed(object): Emitted with the raised exception on any error, so the
            window can translate it by type (e.g. ``api.InfeasibleError``).
        progress(int, str): Emitted with a percentage and a short status label.
    """

    finished = QtCore.Signal(object)
    failed = QtCore.Signal(object)
    progress = QtCore.Signal(int, str)

    def compute(self) -> object:
        """Run the engine call synchronously and return its result."""
        raise NotImplementedError

    @QtCore.Slot()
    def run(self) -> None:
        """Execute :meth:`compute` and emit progress/finished/failed signals.

        This is the slot connected to ``QThread.started``. It never raises: an
        engine ``ValueError`` (including ``InfeasibleError``) or any other
        exception is converted into a :attr:`failed` signal carrying the
        exception object -- a worker must never crash the UI thread.
        """
        try:
            self.progress.emit(0, "starting")
            self.finished.emit(self.compute())
        except Exception as exc:
            self.failed.emit(exc)


class OptimizeWorker(_EngineWorker):
    """Runs one frontier optimization on a background thread."""

    def __init__(self, protein: str, config: api.OptimizeConfig, steps: int) -> None:
        """Store the run parameters (no Qt parent, so it can move to a thread).

        Args:
            protein: The single-letter amino-acid string to back-translate.
            config: The optimization configuration.
            steps: Number of frontier scalarization steps to sweep.
        """
        super().__init__()
        self._protein = protein
        self._config = config
        self._steps = steps
        self._cancelled = False

    def cancel(self) -> None:
        """Request cooperative cancellation (safe to call from the GUI thread).

        Sets a plain flag the compute loop polls between frontier points. Setting
        a bool across threads is atomic under CPython, and this must NOT go
        through a queued signal: while a solve is running the worker thread is
        busy and never processes its own event loop, so a queued slot would not
        run until the whole optimization had already finished.
        """
        self._cancelled = True

    def compute(self) -> api.FrontierResult:
        """Run the frontier optimization synchronously and return the result.

        Reports per-point progress and honors :meth:`cancel`.

        Returns:
            The computed :class:`~bt4.api.FrontierResult`.
        """
        return api.frontier(
            self._protein,
            self._config,
            self._steps,
            on_progress=self._emit_progress,
            should_cancel=lambda: self._cancelled,
        )

    def _emit_progress(self, done: int, total: int) -> None:
        """Translate a (done, total) grid position into a progress signal."""
        pct = round(100 * done / total) if total else 0
        self.progress.emit(pct, f"solving frontier point {done} of {total}")


class CandidatesWorker(_EngineWorker):
    """Assembles + expression-ranks a candidate set and splice-audits it off-thread.

    Runs design-flow steps 3-4 (:func:`bt4.api.candidates` then
    :func:`bt4.api.splice_audit`) on a background thread so the window never
    blocks -- the expression head and the splice audit in particular can invoke
    wrapped CNN backends, which are far too slow for the GUI thread.
    """

    def __init__(
        self,
        protein: str,
        config: api.OptimizeConfig,
        *,
        steps: int,
        n: int,
        repeat_variants: int,
        include_cnns: bool,
        predictor: api.ExpressionPredictor | None = None,
    ) -> None:
        """Store the run parameters (no Qt parent, so it can move to a thread).

        Args:
            protein: The single-letter amino-acid string to back-translate.
            config: The optimization configuration.
            steps: Frontier scalarization grid resolution.
            n: Maximum candidates to keep after scoring.
            repeat_variants: Repeat-refined variants to attempt for a GLOBAL-rule
                violating seed.
            include_cnns: When ``True``, run every backend
                :func:`bt4.api.available_splice_backends` reports (adding the
                wrapped SpliceAI / Pangolin CNNs if installed); otherwise the
                audit runs the honest PWM baseline only.
            predictor: Expression head to annotate the set with; ``None`` uses
                :func:`bt4.api.expression_model` (the neutral placeholder). An
                uncalibrated head -- which is every shipped head today -- only
                annotates: :func:`bt4.api.candidates` keeps the set in discovery
                order and leaves the solver's pick delivered (CLAUDE.md §10.6).
        """
        super().__init__()
        self._protein = protein
        self._config = config
        self._steps = steps
        self._n = n
        self._repeat_variants = repeat_variants
        self._include_cnns = include_cnns
        self._predictor = predictor

    def compute(self) -> CandidatesResult:
        """Assemble + rank candidates, then splice-audit them, synchronously.

        Returns:
            A :class:`CandidatesResult` pairing the candidate set with its splice
            audit (the audit is ``None`` only when the set is empty).
        """
        self.progress.emit(10, "assembling candidate set")
        candidate_set = api.candidates(
            self._protein,
            self._config,
            steps=self._steps,
            n=self._n,
            repeat_variants=self._repeat_variants,
            predictor=self._predictor,
        )
        if not candidate_set.candidates:
            return CandidatesResult(candidate_set, None)
        self.progress.emit(60, "auditing splice sites")
        predictors = api.available_splice_backends() if self._include_cnns else None
        audit = api.splice_audit(candidate_set, predictors=predictors)
        self.progress.emit(100, "done")
        return CandidatesResult(candidate_set, audit)


class CrossCheckWorker(_EngineWorker):
    """Runs an opt-in, out-of-loop splice cross-check on one delivered sequence.

    Wraps :func:`bt4.api.splice_crosscheck`. For the ASSP backend this is the one
    place BT4 Studio touches the network, so it is explicitly user-triggered and
    runs here rather than on the GUI thread (the service is rate-limited and
    retried with backoff, so a call can take seconds). The cross-check is
    *never blocking*: an outage comes back as ``available is False`` with a
    reason rather than an exception, so the run it audits is never failed
    (CLAUDE.md §10.15).
    """

    def __init__(self, dna: str, *, backend: str = "assp") -> None:
        """Store the sequence and backend name (no Qt parent, so it can move)."""
        super().__init__()
        self._dna = dna
        self._backend = backend

    def compute(self) -> api.SpliceCrossCheck:
        """Run the cross-check synchronously and return its report."""
        self.progress.emit(30, f"cross-checking with {self._backend}")
        report = api.splice_crosscheck(self._dna, backend=self._backend)
        self.progress.emit(100, "done")
        return report


class LibraryWorker(_EngineWorker):
    """Samples a sequence library off-thread (:func:`bt4.api.library`).

    Library mode is a stochastic **sampler, not an optimizer**: every member
    carries the ``SAMPLED`` certificate and makes no optimality or expression
    claim (CLAUDE.md §9, Phase 5). The window labels it accordingly.
    """

    def __init__(
        self,
        protein: str,
        config: api.OptimizeConfig,
        *,
        n: int,
        temperature: float,
        seed: int,
    ) -> None:
        """Store the sampling parameters (no Qt parent, so it can move)."""
        super().__init__()
        self._protein = protein
        self._config = config
        self._n = n
        self._temperature = temperature
        self._seed = seed

    def compute(self) -> api.LibraryResult:
        """Draw the library synchronously and return it."""
        self.progress.emit(20, f"sampling {self._n} sequence(s)")
        result = api.library(
            self._protein,
            self._config,
            self._n,
            seed=self._seed,
            temperature=self._temperature,
        )
        self.progress.emit(100, "done")
        return result
