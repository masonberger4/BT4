"""BT4 cli layer -- the only layer permitted to print.

The command-line entry point is a thin, print-oriented shell over
:mod:`bt4.api`. All computation happens behind the api; the CLI only parses
arguments and renders results.
"""

from __future__ import annotations
