"""The round-trip invariant, property-tested.

BT4's bedrock guarantee: any synonymous back-translation of a protein translates
back to that protein (plus the stop). There is no optimizer yet, so we exercise
the genetic code directly with a naive deterministic back-translator — the same
invariant the real Engine will have to preserve.
"""

from __future__ import annotations

from hypothesis import given
from hypothesis import strategies as st

from bt4.domain.genetic_code import AMINO_ACIDS, STOP, synonymous_codons, translate

_PROTEIN = st.text(alphabet=sorted(AMINO_ACIDS), min_size=1, max_size=200)


def _naive_backtranslate(protein: str, *, append_stop: bool = True) -> str:
    """Pick the first (sorted) synonymous codon per residue — deterministic."""
    dna = "".join(synonymous_codons(aa)[0] for aa in protein)
    if append_stop:
        dna += synonymous_codons(STOP)[0]
    return dna


@given(protein=_PROTEIN)
def test_roundtrip_with_stop(protein: str) -> None:
    dna = _naive_backtranslate(protein, append_stop=True)
    assert translate(dna) == protein + STOP


@given(protein=_PROTEIN)
def test_roundtrip_without_stop(protein: str) -> None:
    dna = _naive_backtranslate(protein, append_stop=False)
    assert translate(dna) == protein


@given(protein=_PROTEIN)
def test_backtranslation_is_deterministic(protein: str) -> None:
    assert _naive_backtranslate(protein) == _naive_backtranslate(protein)
