"""BT4 — constrained, multi-objective protein to mRNA back-translation optimizer.

The top-level package is intentionally lightweight: importing ``bt4`` pulls in
only the pure, stdlib-only ``domain`` layer. Heavy optional dependencies
(torch, ViennaRNA, OR-Tools, FastAPI) live behind extras and are imported
lazily by the layers that need them.
"""

from __future__ import annotations

__version__ = "0.0.0"
"""Single source of truth for the BT4 version (surfaced via ``bt4 --version``)."""
