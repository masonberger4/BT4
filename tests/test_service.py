"""Tests for the optional headless HTTP service layer."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient

from bt4.service.api import create_app

client = TestClient(create_app())


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
