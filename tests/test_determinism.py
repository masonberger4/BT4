"""End-to-end determinism guard (the CI `-k determinism` job targets this).

As BT4 grows, this module asserts that identical inputs produce byte-identical
outputs across the surfaces that exist. Today that is the provenance stamp and
the naive back-translation; every new deterministic surface should add a case
here rather than trusting determinism by inspection.
"""

from __future__ import annotations

from bt4.domain.genetic_code import AMINO_ACIDS, synonymous_codons
from bt4.provenance import build_manifest, content_hash


def test_determinism_manifest_stamp_is_stable() -> None:
    def make() -> str:
        return build_manifest(
            bt4_version="0.0.0",
            config={"organism": "homo_sapiens", "beam": 8, "seed": 1},
            inputs={"table": content_hash("homo_sapiens")},
            seed=1,
        ).stamp

    assert make() == make()


def test_determinism_backtranslation_is_stable() -> None:
    protein = "".join(sorted(AMINO_ACIDS)) * 3

    def make() -> str:
        return "".join(synonymous_codons(aa)[0] for aa in protein)

    assert make() == make()
