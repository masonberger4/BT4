"""Background optimization worker for BT4 Studio.

The worker wraps :func:`bt4.api.frontier` so the Pareto sweep runs off the GUI
thread and the window never blocks. It communicates purely through Qt signals
and never touches widgets, which makes it safe to move onto a ``QThread``. Any
engine error is caught and reported (as the exception object) through
:attr:`failed`, so bad input can never crash the UI thread and the window can
translate the failure into a plain-language message keyed on its type.
"""

from __future__ import annotations

from PySide6 import QtCore

from bt4 import api

__all__ = ["CandidatesResult", "CandidatesWorker", "OptimizeWorker"]


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


class OptimizeWorker(QtCore.QObject):
    """Runs one frontier optimization on a background thread.

    Signals:
        finished(object): Emitted with the :class:`~bt4.api.FrontierResult` on
            success.
        failed(object): Emitted with the raised exception on any error, so the
            window can translate it by type (e.g. ``api.InfeasibleError``).
        progress(int, str): Emitted with a percentage and a short status label.
    """

    finished = QtCore.Signal(object)
    failed = QtCore.Signal(object)
    progress = QtCore.Signal(int, str)

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

        Factored out of :meth:`run` so tests (and any synchronous caller) can
        drive the engine without a thread or an event loop. Reports per-point
        progress and honors :meth:`cancel`.

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

    @QtCore.Slot()
    def run(self) -> None:
        """Execute the optimization and emit progress/finished/failed signals.

        This is the slot connected to ``QThread.started``. It never raises: an
        engine ``ValueError`` (including ``InfeasibleError``) or any other
        exception is converted into a :attr:`failed` signal carrying the
        exception object.
        """
        try:
            self.progress.emit(0, "starting")
            result = self.compute()
            self.finished.emit(result)
        except ValueError as exc:
            # Bad protein / dna / organism, InfeasibleError, or a cancel-with-no-
            # points -- all ValueErrors the window turns into friendly text.
            self.failed.emit(exc)
        except Exception as exc:  # a worker must never crash the UI thread
            self.failed.emit(exc)


class CandidatesWorker(QtCore.QObject):
    """Assembles + expression-ranks a candidate set and splice-audits it off-thread.

    Runs design-flow steps 3-4 (:func:`bt4.api.candidates` then
    :func:`bt4.api.splice_audit`) on a background thread so the window never
    blocks -- the splice audit in particular can invoke wrapped CNN backends,
    which are far too slow for the GUI thread. Communicates purely through Qt
    signals; any engine error is caught and reported through :attr:`failed`.

    Signals:
        finished(object): Emitted with a :class:`CandidatesResult` on success.
        failed(object): Emitted with the raised exception on any error.
        progress(int, str): Emitted with a percentage and a short status label.
    """

    finished = QtCore.Signal(object)
    failed = QtCore.Signal(object)
    progress = QtCore.Signal(int, str)

    def __init__(
        self,
        protein: str,
        config: api.OptimizeConfig,
        *,
        steps: int,
        n: int,
        repeat_variants: int,
        include_cnns: bool,
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
        """
        super().__init__()
        self._protein = protein
        self._config = config
        self._steps = steps
        self._n = n
        self._repeat_variants = repeat_variants
        self._include_cnns = include_cnns

    def compute(self) -> CandidatesResult:
        """Assemble + rank candidates, then splice-audit them, synchronously.

        Factored out of :meth:`run` so tests (and any synchronous caller) can
        drive the flow without a thread or an event loop.

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
        )
        if not candidate_set.candidates:
            return CandidatesResult(candidate_set, None)
        self.progress.emit(60, "auditing splice sites")
        predictors = api.available_splice_backends() if self._include_cnns else None
        audit = api.splice_audit(candidate_set, predictors=predictors)
        self.progress.emit(100, "done")
        return CandidatesResult(candidate_set, audit)

    @QtCore.Slot()
    def run(self) -> None:
        """Execute steps 3-4 and emit progress/finished/failed signals.

        Connected to ``QThread.started``. It never raises: an engine
        ``ValueError`` (including ``InfeasibleError``) or any other exception is
        converted into a :attr:`failed` signal carrying the exception object.
        """
        try:
            self.progress.emit(0, "starting")
            result = self.compute()
            self.finished.emit(result)
        except ValueError as exc:
            self.failed.emit(exc)
        except Exception as exc:  # a worker must never crash the UI thread
            self.failed.emit(exc)
