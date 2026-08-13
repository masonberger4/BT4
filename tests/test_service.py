"""Tests for the optional headless HTTP service layer."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient

from bt4.biomodels.codon.tables import GENOME_WIDE, HIGHLY_EXPRESSED
from bt4.service.api import create_app

client = TestClient(create_app())

# An organism whose two reference sets disagree at eight amino acids.
ORGANISM = "escherichia_coli"


def test_health() -> None:
    """The health route reports ok and a version string."""
    response = client.get("/health")
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "ok"
    assert isinstance(body["version"], str)
    assert body["version"]


def test_organisms() -> None:
    """The organisms route lists targetable codon tables including human."""
    response = client.get("/organisms")
    assert response.status_code == 200
    assert "homo_sapiens" in response.json()["organisms"]


def test_optimize() -> None:
    """Optimize returns a proven-optimal, well-formed coding sequence."""
    response = client.post(
        "/optimize",
        json={"protein": "MAALKHETQW", "config": {"max_homopolymer": 5}},
    )
    assert response.status_code == 200
    body = response.json()
    dna = body["dna"]
    assert body["certificate"]["status"] == "proven_optimal"
    assert dna.startswith("ATG")
    assert dna[-3:] in {"TAA", "TAG", "TGA"}


def test_frontier() -> None:
    """Frontier returns a non-empty point set with a valid delivered index."""
    response = client.post("/frontier", json={"protein": "MAALKHETQW", "steps": 5})
    assert response.status_code == 200
    body = response.json()
    assert body["points"]
    assert body["delivered_index"] >= 0


def test_validate() -> None:
    """Validate flags a homopolymer violation and reports infeasibility."""
    response = client.post(
        "/validate",
        json={"dna": "GGGCCCAAAAAA", "config": {"max_homopolymer": 3}},
    )
    assert response.status_code == 200
    body = response.json()
    assert body["is_feasible"] is False
    assert any("homopolymer" in v["constraint"].lower() for v in body["violations"])


def test_optimize_invalid_protein() -> None:
    """An invalid amino-acid string yields HTTP 400."""
    response = client.post("/optimize", json={"protein": "MAZX"})
    assert response.status_code == 400


def test_optimize_dinuc_budget() -> None:
    """A CpG count budget routes through the exact bucketed DP over the service."""
    response = client.post(
        "/optimize",
        json={
            "protein": "MARPGARSTKLE",
            "config": {"dinuc_budget": "CG", "dinuc_max": 3},
        },
    )
    assert response.status_code == 200
    body = response.json()
    assert body["certificate"]["status"] == "proven_optimal"
    dna = body["dna"]
    cpg = sum(1 for i in range(len(dna) - 1) if dna[i : i + 2] == "CG")
    assert cpg <= 3


def test_config_accepts_new_fields_and_rejects_unknown() -> None:
    """ConfigModel mirrors the full engine config and forbids unknown keys.

    Audit C2: the request model previously exposed only 9 of the config's fields
    and silently dropped the rest. It now accepts every knob (e.g. tai_weight,
    restriction_enzymes) and rejects a mistyped/unknown key with a 422 rather than
    silently ignoring it.
    """
    ok = client.post(
        "/optimize",
        json={
            "protein": "MAALKHETQWSNDECFGR",
            "config": {"tai_weight": 1.0, "restriction_enzymes": ["EcoRI"], "cpg_weight": 2.0},
        },
    )
    assert ok.status_code == 200
    bad = client.post(
        "/optimize",
        json={"protein": "MAAL", "config": {"not_a_real_field": 1}},
    )
    assert bad.status_code == 422


def test_service_reference_set_reaches_the_engine() -> None:
    body = {
        "protein": "MAALKHETQW",
        "config": {"organism": ORGANISM, "reference_set": GENOME_WIDE},
    }
    response = client.post("/optimize", json=body)
    assert response.status_code == 200, response.text
    assert response.json()["audit"]["codon_reference_set"] == GENOME_WIDE


def test_service_rejects_an_unknown_reference_set() -> None:
    """An unusable value must fail loudly, not fall back to a default."""
    response = client.post(
        "/optimize",
        json={"protein": "MAAL", "config": {"reference_set": "tissue_specific"}},
    )
    assert response.status_code >= 400


def test_service_organisms_lists_reference_sets() -> None:
    payload = client.get("/organisms").json()
    assert "homo_sapiens" in payload["organisms"]
    assert payload["reference_sets"]["homo_sapiens"] == [HIGHLY_EXPRESSED, GENOME_WIDE]
    assert payload["reference_sets"]["arabidopsis_thaliana"] == [GENOME_WIDE]
