"""Tests for library / degenerate-design mode (CLAUDE.md §9, Phase 5).

Library mode is a **stochastic sampler**, not an optimizer, so these tests pin
the honesty properties that must survive that: every sampled sequence still
round-trips (#1), still respects the LOCAL hard constraints (#3), still carries
metrics recomputed from its own DNA (#2), and carries a certificate that says it
was *sampled* -- never optimal, never an expression prediction. Determinism (#7)
is exercised too: same inputs, same library.
"""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from bt4 import api
from bt4._accel import max_homopolymer_run
from bt4.biomodels.codon.tables import load_table
from bt4.domain import OptimalityStatus, gc_fraction
from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons, translate
from bt4.optimize import InfeasibleError

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=30)
_RICH = "MAALKHETQWSNDECFGRPVIY"


@given(protein=_PROTEIN)
@settings(max_examples=40, deadline=None)
def test_library_roundtrips(protein: str) -> None:
    # Invariant #1: every sampled member translates back to protein + stop.
    lib = api.library(protein, api.OptimizeConfig(max_homopolymer=None), n=4, seed=7)
    assert len(lib.results) == 4
    for res in lib.results:
        assert translate(res.dna) == res.protein + "*"


def test_library_is_deterministic() -> None:
    # Invariant #7: identical (protein, config, n, seed, temperature) => identical
    # library, byte-for-byte, including the shared manifest stamp.
    config = api.OptimizeConfig(max_homopolymer=5)
    first = api.library(_RICH, config, n=8, seed=123, temperature=1.3)
    second = api.library(_RICH, config, n=8, seed=123, temperature=1.3)
    assert [r.dna for r in first.results] == [r.dna for r in second.results]
    assert first.manifest.stamp == second.manifest.stamp
    assert first.results[0].audit["manifest"] == second.results[0].audit["manifest"]


def test_library_respects_local_constraints() -> None:
    # Invariant #3: sampling only ever picks codons passing ok_suffix, so a
    # homopolymer + GC-run + forbidden-motif set yields zero hard violations.
    config = api.OptimizeConfig(
        max_homopolymer=4,
        max_gc_run=5,
        forbidden_motifs=("GAATTC",),
        avoid_reverse_complement=False,
    )
    lib = api.library(_RICH * 2, config, n=12, seed=99)
    for res in lib.results:
        assert res.metrics.hard_violations == 0
        assert max_homopolymer_run(res.dna) <= 4
        assert "GAATTC" not in res.dna
        assert not any(
            v.constraint in {"homopolymer", "gc_run", "forbidden_motif"}
            for v in res.hard_violations
        )


def test_library_is_diverse() -> None:
    # A protein rich in high-degeneracy residues sampled at T=1 must not collapse
    # to a single sequence.
    lib = api.library(_RICH, api.OptimizeConfig(max_homopolymer=None), n=30, seed=5)
    assert lib.distinct > 1
    # Distinct count and Hamming diversity are honest, bounded readouts.
    assert 0.0 <= lib.mean_pairwise_hamming <= 1.0


def test_library_distribution_tracks_table_frequencies() -> None:
    # Distribution sanity (deterministic via the seed): for a single-residue
    # protein sampled at T=1, the empirical per-codon frequencies roughly track
    # the table frequencies. Loose tolerance -- a sanity check, not a strict test.
    table = load_table("homo_sapiens")
    leu_codons = synonymous_codons("L")
    total_freq = sum(table.frequency[c] for c in leu_codons)
    expected = {c: table.frequency[c] / total_freq for c in leu_codons}

    protein = "L" * 30
    lib = api.library(protein, api.OptimizeConfig(max_homopolymer=None), n=120, seed=2024)
    counts = {c: 0 for c in leu_codons}
    total = 0
    for res in lib.results:
        # Every codon but the trailing stop encodes L.
        for i in range(0, len(res.dna) - 3, 3):
            counts[res.dna[i : i + 3]] += 1
            total += 1
    assert total == 120 * 30
    for codon in leu_codons:
        empirical = counts[codon] / total
        assert abs(empirical - expected[codon]) < 0.06


def test_library_certificate_is_sampled() -> None:
    lib = api.library(_RICH, api.OptimizeConfig(max_homopolymer=5), n=6, seed=1)
    for res in lib.results:
        cert = res.certificate
        assert cert.status is OptimalityStatus.SAMPLED
        assert cert.status.value == "sampled"
        assert not cert.is_proven_optimal
        assert cert.gap is None
        assert cert.solver == "library_sampler"


def test_library_metrics_are_recomputed() -> None:
    # Invariant #2: reported CAI/GC equal a fresh recomputation from the DNA.
    table = load_table("homo_sapiens")
    lib = api.library(_RICH, api.OptimizeConfig(max_homopolymer=5), n=5, seed=42)
    for res in lib.results:
        assert res.audit["cai"] == pytest.approx(table.cai(res.dna))
        assert res.metrics.gc == pytest.approx(gc_fraction(res.dna))
        assert res.metrics.length_nt == len(res.dna)


def test_library_infeasible_raises() -> None:
    # Trp (W) has the single codon TGG; forbidding it leaves that residue with no
    # feasible codon, so the sampler raises rather than looping forever.
    config = api.OptimizeConfig(
        forbidden_motifs=("TGG",),
        avoid_reverse_complement=False,
        max_homopolymer=None,
    )
    with pytest.raises(InfeasibleError):
        api.library("AWA", config, n=3, seed=0)


def test_library_global_constraints_not_enforced_but_validated() -> None:
    # GLOBAL rules (max_repeat_length) are not enforced during sampling -- no
    # refinement runs, so the certificate stays SAMPLED -- but the full local +
    # global set is validated on each member, so any residual repeat is honestly
    # reported in that member's violations rather than silently dropped.
    config = api.OptimizeConfig(max_repeat_length=8, max_homopolymer=None)
    lib = api.library("LLLLRRRRSSSS", config, n=6, seed=3)
    for res in lib.results:
        assert res.certificate.status is OptimalityStatus.SAMPLED
        assert translate(res.dna) == res.protein + "*"
        # Any max-repeat violations present are reported (never hidden); the count
        # in metrics matches the reported hard violations.
        assert res.metrics.hard_violations == len(res.hard_violations)


def test_library_config_seed_is_honored_and_stamped() -> None:
    # When no explicit seed is passed, config.seed drives the draw and enters the
    # manifest; two different config seeds give different libraries and stamps.
    a = api.library(_RICH, api.OptimizeConfig(max_homopolymer=5, seed=1), n=6)
    b = api.library(_RICH, api.OptimizeConfig(max_homopolymer=5, seed=2), n=6)
    assert [r.dna for r in a.results] != [r.dna for r in b.results]
    assert a.manifest.stamp != b.manifest.stamp
