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

__all__ = ["OptimizeWorker"]


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
