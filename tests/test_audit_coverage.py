"""Regression tests for the adversarial-audit fixes and coverage gaps.

Each test pins a CONFIRMED audit finding's fix or fills a coverage gap the audit
flagged (CLAUDE.md invariants #6/#8/#9 and the data-honesty rules).
"""

from __future__ import annotations

import subprocess
import sys
from importlib.resources import files

import pytest

from bt4 import api
from bt4.biomodels.codon.tables import CodonUsageTable
from bt4.biomodels.codon.tai import available_tai_organisms, load_tai_provenance
from bt4.domain.genetic_code import CODON_TABLE, STOP, translate
from bt4.provenance import content_hash, resolve_git_commit

_DATA = "bt4.biomodels.codon.data"


def test_core_import_is_light() -> None:
    # Invariant: `import bt4` (and the api) must not pull heavy optional deps.
    code = (
        "import bt4, bt4.api, sys;"
        "bad=[m for m in ('torch','ortools','RNA','fastapi','PySide6') if m in sys.modules];"
        "print(bad); assert not bad, bad"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_tai_provenance_sha_matches_tsv() -> None:
    # C-minor: the bundled tRNA table content must match its provenance SHA-256.
    for organism in available_tai_organisms():
        data = files(_DATA).joinpath(f"{organism}.trna.tsv").read_bytes()
        assert content_hash(data) == load_tai_provenance(organism).sha256, organism


def test_manifest_records_git_commit() -> None:
    # C4: the manifest must carry the source git commit (or None outside a repo).
    result = api.optimize("MAAL")
    manifest = result.audit["manifest"]
    assert isinstance(manifest, dict)
    assert manifest["git_commit"] == resolve_git_commit()


def test_folding_backend_enters_manifest_stamp() -> None:
    # C3: the folding backend identity must change the stamp (invariant #9), so a
    # stamp can never map to sequences from two different folding backends.
    from bt4.pipeline.optimize import OptimizeConfig, _manifest

    cfg = OptimizeConfig()
    m_a = _manifest(cfg, {"folding_model": "a", "folding_calibrated": False})
    m_b = _manifest(cfg, {"folding_model": "b", "folding_calibrated": False})
    assert m_a.config_hash != m_b.config_hash


def test_cai_missing_codon_raises_valueerror() -> None:
    # C-minor: cai() must raise the documented ValueError (not a bare KeyError)
    # for a valid codon absent from a sparse-but-valid table.
    from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons

    freq = {synonymous_codons(aa)[0]: 1.0 for aa in AMINO_ACIDS}
    freq[synonymous_codons(STOP)[0]] = 1.0
    table = CodonUsageTable("sparse", freq)
    missing = next(
        c for c in CODON_TABLE if c not in freq and CODON_TABLE[c] not in (STOP, "M", "W")
    )
    with pytest.raises(ValueError, match="no weight"):
        table.cai(missing)


def test_all_stops_forbidden_is_infeasible() -> None:
    # Invariant #8 (negative side): forbidding all three stops leaves no feasible
    # appended stop codon, so the solve is genuinely infeasible.
    config = api.OptimizeConfig(
        forbidden_motifs=("TAA", "TAG", "TGA"),
        avoid_reverse_complement=False,
        max_homopolymer=None,
    )
    with pytest.raises(api.InfeasibleError):
        api.optimize("MA", config)


@pytest.mark.parametrize("protein", ["M", "W", "MW", "MWM"])
def test_degenerate_single_codon_proteins(protein: str) -> None:
    # All-single-codon proteins (no synonymous choice) must optimize and trace a
    # (degenerate) frontier without error.
    result = api.optimize(protein)
    assert translate(result.dna) == protein + "*"
    frontier = api.frontier(protein, steps=5)
    assert frontier.delivered() is not None


def test_validate_non_multiple_of_three() -> None:
    # run_validate's non-multiple-of-3 branch: audits without computing objectives.
    report = api.validate("ATGA")
    assert report.metrics.length_nt == 4
    assert report.metrics.objective.terms() == frozenset()
