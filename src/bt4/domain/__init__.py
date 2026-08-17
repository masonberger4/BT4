"""BT4 domain layer — pure, stdlib-only. Depends on nothing else in bt4.

This package holds the vocabulary the whole system is written in: the genetic
code, sequence validation, the multi-objective primitives, the optimality
certificate, and the immutable result types. Nothing here imports optimization
logic, biological models, or heavy third-party dependencies.
"""

from __future__ import annotations

from .certificate import OptimalityCertificate, OptimalityStatus
from .context import CIRCULAR, LINEAR, ConstructContext
from .contracts import Constraint, ObjectiveTerm
from .genetic_code import (
    AA_TO_CODONS,
    AMINO_ACIDS,
    CODON_TABLE,
    STOP,
    is_stop,
    synonymous_codons,
    translate,
)
from .objective import Frontier, ObjectiveVector, dominates, pareto_front
from .relax import SoftConstraint, is_relaxable, relax_constraint
from .result import Metrics, Result, Severity, Violation
from .scope import Scope
from .sequence import DNA_BASES, gc_fraction, validate_dna, validate_protein

__all__ = [
    "AA_TO_CODONS",
    "AMINO_ACIDS",
    "CIRCULAR",
    "CODON_TABLE",
    "DNA_BASES",
    "LINEAR",
    "STOP",
    "Constraint",
    "ConstructContext",
    "Frontier",
    "Metrics",
    "ObjectiveTerm",
    "ObjectiveVector",
    "OptimalityCertificate",
    "OptimalityStatus",
    "Result",
    "Scope",
    "Severity",
    "SoftConstraint",
    "Violation",
    "dominates",
    "gc_fraction",
    "is_relaxable",
    "is_stop",
    "pareto_front",
    "relax_constraint",
    "synonymous_codons",
    "translate",
    "validate_dna",
    "validate_protein",
]
