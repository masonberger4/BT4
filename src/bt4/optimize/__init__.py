"""BT4 optimize layer - solver backends over the codon trellis.

This layer owns the solvers that turn a protein plus an objective and a set of
constraints into a coding sequence with an honest optimality certificate. The
first backend is :func:`solve_exact`, an exact (or explicitly beam-truncated)
dynamic program.
"""

from __future__ import annotations

from bt4.optimize.anneal_refine import anneal_refine
from bt4.optimize.exact_dp import InfeasibleError, SolveResult, solve_exact
from bt4.optimize.sample import sample_sequences

__all__ = [
    "InfeasibleError",
    "SolveResult",
    "anneal_refine",
    "sample_sequences",
    "solve_exact",
]
