"""The Rust accelerator and its pure-Python fallback must agree exactly."""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from bt4 import _accel
from bt4._accel import (
    _py_gc_count,
    _py_longest_repeat,
    _py_max_gc_run,
    _py_max_homopolymer_run,
    _py_reverse_complement,
    _py_trellis_solve,
)
from bt4.constraints.gc_run import _max_gc_run
from bt4.constraints.max_repeat import MaxRepeatConstraint
from bt4.constraints.rules import HomopolymerConstraint
from bt4.domain.result import Severity
from bt4.optimize.exact_dp import _delta_tables, _precompute_structure

_DNA = st.text(alphabet="ACGT", min_size=1, max_size=200)

# Crafted sequences pinning each repeat flavour the primitives must agree on: a
# palindrome, an inverted (reverse-complement) repeat, a far-apart direct repeat,
# an overlapping direct repeat, a GC-run edge, and a repeat-free control.
_CRAFTED = (
    "AATT",  # palindrome (equals its own reverse complement)
    "AAGGACCCTT",  # inverted repeat: AAGG at 0, revcomp CCTT at 6
    "AAACCCAAAC",  # direct repeat: AAAC at 0 and 6
    "AAAAAA",  # overlapping direct repeats of A
    "GCGCGCAT",  # a long mixed GC run then a break
    "ACGGATTCAG",  # no offending 4-mer (clean at max_length=3)
    "ATGC",  # reverse-complement palindrome ACGT-style (ATGC -> GCAT)
    "ACGT",  # palindromic 4-mer (equals its own reverse complement)
)


def test_reports_which_path_is_active() -> None:
    assert isinstance(_accel.ACCELERATED, bool)


def test_reverse_complement_known() -> None:
    assert _py_reverse_complement("ATGC") == "GCAT"


def test_max_homopolymer_run_known() -> None:
    assert _py_max_homopolymer_run("AAATTGCCCC") == 4
    assert _py_max_homopolymer_run("ACGT") == 1


def test_max_gc_run_known() -> None:
    assert _py_max_gc_run("") == 0
    assert _py_max_gc_run("ATAT") == 0  # no G or C at all
    assert _py_max_gc_run("AC") == 1  # a single strong base is a run of one
    assert _py_max_gc_run("GCGCGC") == 6  # mixed G/C is ONE run, not a homopolymer
    assert _py_max_gc_run("ATGCGCAAA") == 4  # positions 2..5 (G,C,G,C)
    assert _py_max_gc_run("GGGATCCCC") == 4  # the longer of GGG (3) and CCCC (4)


def test_longest_repeat_known() -> None:
    assert _py_longest_repeat("") == 0
    assert _py_longest_repeat("A") == 0  # a single base is not a repeat
    assert _py_longest_repeat("AC") == 0  # no repeat and no revcomp coincidence
    assert _py_longest_repeat("AAACCCAAAC") == 4  # direct repeat AAAC (0 and 6)
    assert _py_longest_repeat("AAAAAA") == 5  # overlapping direct repeat of A's
    assert _py_longest_repeat("AATT") == 4  # palindrome (== own reverse complement)
    assert _py_longest_repeat("AAGGACCCTT") == 4  # inverted: AAGG vs revcomp CCTT
    assert _py_longest_repeat("ACGT") == 4  # ACGT is its own reverse complement


@given(seq=_DNA)
def test_active_impl_matches_pure_python(seq: str) -> None:
    # Whether or not the Rust extension is loaded, the public functions must
    # equal the reference pure-Python implementations.
    assert _accel.reverse_complement(seq) == _py_reverse_complement(seq)
    assert _accel.gc_count(seq) == _py_gc_count(seq)
    assert _accel.max_homopolymer_run(seq) == _py_max_homopolymer_run(seq)
    assert _accel.max_gc_run(seq) == _py_max_gc_run(seq)
    assert _accel.longest_repeat(seq) == _py_longest_repeat(seq)


@pytest.mark.parametrize("seq", _CRAFTED)
def test_active_impl_matches_pure_python_crafted(seq: str) -> None:
    # The crafted palindrome / inverted / direct-repeat cases exercise every
    # branch of the RC-aware repeat notion on both the active and reference paths.
    assert _accel.max_gc_run(seq) == _py_max_gc_run(seq)
    assert _accel.longest_repeat(seq) == _py_longest_repeat(seq)


@given(seq=_DNA)
def test_max_gc_run_matches_constraint_reference(seq: str) -> None:
    # The accelerator must equal the GcRunConstraint's own pure-Python scan, the
    # value its ok_suffix veto depends on (invariant #2, §7 wiring).
    assert _accel.max_gc_run(seq) == _max_gc_run(seq)


@given(seq=_DNA, m=st.integers(min_value=1, max_value=8))
def test_longest_repeat_matches_max_repeat_validate(seq: str, m: int) -> None:
    # The load-bearing equivalence the fast-path short-circuit relies on:
    # longest_repeat(seq) > m  iff  MaxRepeatConstraint(m).validate(seq) flags a
    # hard violation. This is what lets validate return early when
    # longest_repeat(seq) <= m without changing any observable behavior.
    has_violation = any(
        v.severity is Severity.HARD for v in MaxRepeatConstraint(m).validate(seq)
    )
    assert (_accel.longest_repeat(seq) > m) == has_violation


@pytest.mark.parametrize("seq", _CRAFTED)
@pytest.mark.parametrize("m", [1, 2, 3, 4, 5])
def test_longest_repeat_matches_max_repeat_validate_crafted(seq: str, m: int) -> None:
    has_violation = any(
        v.severity is Severity.HARD for v in MaxRepeatConstraint(m).validate(seq)
    )
    assert (_accel.longest_repeat(seq) > m) == has_violation


_CODON_W = {
    "ATG": 0.0, "GCT": 0.1, "GCC": 0.9, "GCA": 0.3, "GCG": 0.2,
    "AAA": 0.4, "AAG": 0.8, "TGT": 0.5, "TGC": 0.7,
    "TAA": 0.6, "TAG": 0.2, "TGA": 0.1,
}


def _codon_delta(prefix: str, codon: str, pos: int) -> float:
    return _CODON_W.get(codon, 0.0)


@given(
    protein=st.text(alphabet="MAKC", min_size=0, max_size=8),
    max_run=st.integers(min_value=2, max_value=4),
    beam=st.one_of(st.none(), st.integers(min_value=1, max_value=5)),
)
def test_trellis_solve_native_matches_pure_python(
    protein: str, max_run: int, beam: int | None
) -> None:
    # The Rust trellis_solve and its pure-Python twin must return byte-identical
    # (dna, scalar, pruned) over the same precomputed transition tables.
    residues = [*protein, "*"]
    cons = (HomopolymerConstraint(max_run),)
    context_len = max(c.context_len() for c in cons)
    structure = _precompute_structure(residues, cons, context_len)
    assert structure is not None
    ld = _delta_tables(structure, _codon_delta)
    args = (
        list(structure.codons),
        [list(x) for x in structure.layer_from],
        [list(x) for x in structure.layer_to],
        [list(x) for x in structure.layer_codon],
        ld,
    )
    native = _accel.trellis_solve(*args, beam)
    pure = _py_trellis_solve(*args, beam)
    assert native == pure


def test_trellis_solve_infeasible_returns_none_both_paths() -> None:
    # A layer with no allowed transitions -> None on both the native and pure path.
    tables = ([("ATG")], [[]], [[]], [[]], [[]])  # one layer, zero transitions
    assert _accel.trellis_solve(*tables, None) is None
    assert _py_trellis_solve(*tables, None) is None


def test_rejects_non_acgt() -> None:
    with pytest.raises(ValueError):
        _accel.reverse_complement("ACGX")
    with pytest.raises(ValueError):
        _accel.max_gc_run("ACGX")
    with pytest.raises(ValueError):
        _accel.longest_repeat("ACGX")
