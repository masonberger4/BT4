"""The Rust accelerator and its pure-Python fallback must agree exactly."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4 import _accel
from bt4._accel import (
    _py_gc_count,
    _py_max_homopolymer_run,
    _py_reverse_complement,
)

_DNA = st.text(alphabet="ACGT", min_size=1, max_size=200)


def test_reports_which_path_is_active() -> None:
    assert isinstance(_accel.ACCELERATED, bool)


def test_reverse_complement_known() -> None:
    assert _py_reverse_complement("ATGC") == "GCAT"


def test_max_homopolymer_run_known() -> None:
    assert _py_max_homopolymer_run("AAATTGCCCC") == 4
    assert _py_max_homopolymer_run("ACGT") == 1


@given(seq=_DNA)
def test_active_impl_matches_pure_python(seq: str) -> None:
    # Whether or not the Rust extension is loaded, the public functions must
    # equal the reference pure-Python implementations.
    assert _accel.reverse_complement(seq) == _py_reverse_complement(seq)
    assert _accel.gc_count(seq) == _py_gc_count(seq)
    assert _accel.max_homopolymer_run(seq) == _py_max_homopolymer_run(seq)


def test_rejects_non_acgt() -> None:
    with pytest.raises(ValueError):
        _accel.reverse_complement("ACGX")
