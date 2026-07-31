"""Locality scope shared by objective terms and constraints.

BT4 partitions every objective term and constraint by *locality* so the planner
knows which can be solved exactly in the codon trellis and which are genuinely
non-local (and must live in the refinement layer). Declaring scope is mandatory:
BT3's sin was smuggling a non-local splice term into the DP with a ``window=0``
hack — BT4 makes locality an explicit, typed contract instead.
"""

from __future__ import annotations

from enum import Enum

__all__ = ["Scope"]


class Scope(Enum):
    """How much sequence context a term or constraint truly depends on.

    Ordered from most local to least local:

    * ``LOCAL`` — depends only on a bounded trailing window of fixed length
      (``context_len``). Exactly solvable per-site in the codon trellis.
    * ``PAIRWISE`` — depends on the immediately preceding codon (e.g. codon-pair
      bias). Solvable by extending the trellis state with the previous codon.
    * ``POSITIONAL`` — LOCAL, but with position-dependent weights (e.g. the 5'
      ramp over the first ~30-50 codons).
    * ``GLOBAL`` — depends on the whole sequence (e.g. total CpG budget, 5'
      folding ΔG, whole-sequence splice risk). Not exactly per-site; handled by
      budgets/relaxation or the refinement layer.
    """

    LOCAL = "local"
    PAIRWISE = "pairwise"
    POSITIONAL = "positional"
    GLOBAL = "global"
