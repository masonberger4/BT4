"""BT4 Studio -- the native PySide6 desktop app for BT4 (CLAUDE.md section 6.6).

This package is the offline-first, in-tree desktop surface for BT4. It calls the
stable :mod:`bt4.api` and nothing below it, and runs each optimization on a
background thread so the window never blocks.

Importing ``bt4.app`` is deliberately cheap: this module never imports PySide6 or
pyqtgraph at package-import time. The Qt UI is pulled in lazily, only when
:func:`main` is actually called, so tools that merely ``import bt4.app`` (e.g.
entry-point discovery) do not pay the cost of loading Qt.
"""

from __future__ import annotations

__all__ = ["main"]


def main() -> int:
    """Launch BT4 Studio and return its process exit code.

    The Qt UI is imported here, inside the function, so importing this package
    stays lightweight (no PySide6 at module top).

    Returns:
        The exit code from the Qt event loop.
    """
    from bt4.app import studio

    return studio.main()
