"""Golden (characterization) regression tests for the ``bt4.api`` optimizer.

This suite pins the *current* delivered output of a single-solve optimization
over a fixed protein panel and a fixed :class:`~bt4.api.OptimizeConfig`. It is a
tripwire: any future change that alters the emitted DNA, the reported CAI/GC, or
the optimality certificate will fail here loudly instead of drifting silently.

Golden values are characterization snapshots of the current engine; regenerate
by re-running :func:`bt4.api.optimize` over the panel and updating the constants
below if (and only if) an intended behavior change lands.

The panel spans a short peptide, a medium sequence, one rich in rare residues
(tryptophan and methionine, which have a single codon each), an all-hydrophobic
stretch, and a longer ~50-aa protein.
"""

from __future__ import annotations

import pytest

from bt4 import api
from bt4.domain.genetic_code import STOP, translate

# The fixed panel: name -> protein. Kept stable so the golden constants below
# stay meaningful across runs.
PANEL: dict[str, str] = {
    "short": "MAAL",
    "medium": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIE",
    "rare_residues": "MWCWMHWQWYWMW",
    "hydrophobic": "AVLIFMAVLIFMAVLIF",
    "long": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEKQANGTWPADEFHLMNCYVGR",
}

# The one configuration under which every golden value below was captured.
CONFIG: api.OptimizeConfig = api.OptimizeConfig(
    organism="homo_sapiens",
    gc_target=0.55,
    max_homopolymer=6,
    forbidden_motifs=("GAATTC",),
    seed=0,
)

# Pinned delivered coding sequences (includes the appended stop codon).
EXPECTED_DNA: dict[str, str] = {
    "short": "ATGGCCGCCCTGTAA",
    "medium": (
        "ATGAAGACCGCCTACATCGCCAAGCAGAGACAGATCAGCTTCGTGAAGAGCCACTTCAGC"
        "AGACAGCTGGAGGAGAGACTGGGCCTGATCGAGTAA"
    ),
    "rare_residues": "ATGTGGTGCTGGATGCACTGGCAGTGGTACTGGATGTGGTAA",
    "hydrophobic": (
        "GCCGTGCTGATCTTCATGGCCGTGCTGATCTTCATGGCCGTGCTGATCTTCTAA"
    ),
    "long": (
        "ATGAAGACCGCCTACATCGCCAAGCAGAGACAGATCAGCTTCGTGAAGAGCCACTTCAGC"
        "AGACAGCTGGAGGAGAGACTGGGCCTGATCGAGAAGCAGGCCAACGGCACCTGGCCCGCC"
        "GACGAGTTCCACCTGATGAACTGCTACGTGGGCAGATAA"
    ),
}

# Pinned CAI (Sharp & Li) of the delivered sequence, rounded to 6 places.
EXPECTED_CAI: dict[str, float] = {
    "short": 1.0,
    "medium": 1.0,
    "rare_residues": 1.0,
    "hydrophobic": 1.0,
    "long": 1.0,
}

# Pinned GC fraction of the delivered sequence, rounded to 6 places.
EXPECTED_GC: dict[str, float] = {
    "short": 0.6,
    "medium": 0.541667,
    "rare_residues": 0.52381,
    "hydrophobic": 0.537037,
    "long": 0.578616,
}

_NAMES: list[str] = list(PANEL)


@pytest.mark.parametrize("name", _NAMES)
def test_golden_dna_matches_snapshot(name: str) -> None:
    """The delivered DNA is byte-identical to the pinned golden snapshot."""
    result = api.optimize(PANEL[name], CONFIG)
    assert result.dna == EXPECTED_DNA[name]


@pytest.mark.parametrize("name", _NAMES)
def test_golden_metrics_match_snapshot(name: str) -> None:
    """Reported CAI and GC match the pinned golden values within tolerance."""
    result = api.optimize(PANEL[name], CONFIG)
    assert float(result.audit["cai"]) == pytest.approx(EXPECTED_CAI[name], abs=1e-6)
    assert result.metrics.gc == pytest.approx(EXPECTED_GC[name], abs=1e-6)


@pytest.mark.parametrize("name", _NAMES)
def test_golden_certificate_is_proven_optimal(name: str) -> None:
    """Every panel entry is solved to certified global optimality by exact DP."""
    result = api.optimize(PANEL[name], CONFIG)
    assert result.certificate.status.value == "proven_optimal"
    assert result.certificate.solver == "exact_dp"


@pytest.mark.parametrize("name", _NAMES)
def test_golden_roundtrips_and_is_feasible(name: str) -> None:
    """Round-trip (#1) holds and the delivered sequence has no hard violations."""
    protein = PANEL[name]
    result = api.optimize(protein, CONFIG)
    assert translate(result.dna) == protein + STOP
    assert result.metrics.hard_violations == 0


@pytest.mark.parametrize("name", _NAMES)
def test_golden_is_deterministic(name: str) -> None:
    """Determinism (#7): two solves of the same input yield identical DNA."""
    first = api.optimize(PANEL[name], CONFIG)
    second = api.optimize(PANEL[name], CONFIG)
    assert first.dna == second.dna
