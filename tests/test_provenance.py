"""Tests for deterministic, content-addressed provenance."""

from __future__ import annotations

from bt4.provenance import build_manifest, config_hash, content_hash


def test_content_hash_matches_for_str_and_bytes() -> None:
    assert content_hash("ACGT") == content_hash(b"ACGT")


def test_config_hash_is_order_independent_but_value_sensitive() -> None:
    a = config_hash({"organism": "homo_sapiens", "beam": 8})
    b = config_hash({"beam": 8, "organism": "homo_sapiens"})
    assert a == b  # key order does not matter
    c = config_hash({"organism": "homo_sapiens", "beam": 16})
    assert a != c  # a value change does


def test_swapped_input_contents_change_the_stamp() -> None:
    base = dict(bt4_version="0.0.0", config={"organism": "custom"}, seed=7)
    m1 = build_manifest(**base, inputs={"table": content_hash("table-v1")})
    m2 = build_manifest(**base, inputs={"table": content_hash("table-v2")})
    # Same config, but different table *contents* => different stamp.
    assert m1.config_hash == m2.config_hash
    assert m1.stamp != m2.stamp


def test_manifest_is_reproducible_from_identical_inputs() -> None:
    kw = dict(
        bt4_version="0.0.0",
        config={"organism": "homo_sapiens", "beam": 8},
        inputs={"table": content_hash("t")},
        seed=42,
        git_commit="abc123",
    )
    assert build_manifest(**kw).stamp == build_manifest(**kw).stamp
