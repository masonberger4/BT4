"""Background optimization worker for BT4 Studio.

The worker wraps :func:`bt4.api.frontier` so the Pareto sweep runs off the GUI
thread and the window never blocks. It communicates purely through Qt signals
and never touches widgets, which makes it safe to move onto a ``QThread``. Any
engine error (``ValueError`` -- which ``bt4.optimize.InfeasibleError`` subclasses
-- or anything unexpected) is caught and reported through :attr:`failed`, so bad
input can never crash the UI thread.
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
        failed(str): Emitted with a human-readable message on any error.
        progress(int, str): Emitted with a percentage and a short status label.
    """

    finished = QtCore.Signal(object)
    failed = QtCore.Signal(str)
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

    def compute(self) -> api.FrontierResult:
        """Run the frontier optimization synchronously and return the result.

        Factored out of :meth:`run` so tests (and any synchronous caller) can
        drive the engine without a thread or an event loop.

        Returns:
            The computed :class:`~bt4.api.FrontierResult`.
        """
        return api.frontier(self._protein, self._config, self._steps)

    @QtCore.Slot()
    def run(self) -> None:
        """Execute the optimization and emit progress/finished/failed signals.

        This is the slot connected to ``QThread.started``. It never raises: an
        engine ``ValueError`` (including ``InfeasibleError``) or any other
        exception is converted into a :attr:`failed` signal.
        """
        try:
            self.progress.emit(0, "starting")
            result = self.compute()
            self.progress.emit(100, "done")
            self.finished.emit(result)
        except ValueError as exc:
            # Bad protein / dna / organism, or InfeasibleError (a ValueError).
            self.failed.emit(str(exc))
        except Exception as exc:
            # A worker must never crash the UI thread; report anything else too.
            self.failed.emit(str(exc))
