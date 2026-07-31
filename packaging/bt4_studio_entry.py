"""PyInstaller entry point for the packaged BT4 Studio desktop app.

Kept as a tiny standalone script (not an ``-m`` invocation) so PyInstaller has a
concrete file to analyze. It defers entirely to :func:`bt4.app.main`.
"""

from __future__ import annotations

import sys

from bt4.app import main

if __name__ == "__main__":
    sys.exit(main())
