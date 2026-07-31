//! BT4 native accelerator.
//!
//! Small, correct, deterministic hot-loop primitives operating on uppercase
//! `ACGT` DNA strings. A pure-Python fallback exists separately, so this module
//! keeps its surface intentionally minimal. No `unsafe`, no external crates
//! beyond `pyo3`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Complement a single uppercase base, or `None` if it is not `A`, `C`, `G`, `T`.
#[inline]
fn complement_base(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        _ => None,
    }
}

/// Raise `ValueError` if `seq` contains any character other than `A`, `C`, `G`, `T`.
#[inline]
fn ensure_acgt(seq: &str) -> PyResult<()> {
    for &b in seq.as_bytes() {
        if !matches!(b, b'A' | b'C' | b'G' | b'T') {
            return Err(PyValueError::new_err(
                "sequence must contain only uppercase A, C, G, T",
            ));
        }
    }
    Ok(())
}

/// Return the reverse complement of an uppercase `ACGT` sequence.
///
/// Raises `ValueError` on any non-`ACGT` character.
#[pyfunction]
fn reverse_complement(seq: &str) -> PyResult<String> {
    let bytes = seq.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        match complement_base(b) {
            Some(c) => out.push(c),
            None => {
                return Err(PyValueError::new_err(
                    "sequence must contain only uppercase A, C, G, T",
                ))
            }
        }
    }
    // SAFETY-free: complements of ASCII ACGT are ASCII, so this is valid UTF-8.
    Ok(String::from_utf8(out).expect("complemented ACGT is valid ASCII"))
}

/// Return the number of `G` or `C` bases in an uppercase `ACGT` sequence.
///
/// Raises `ValueError` on any non-`ACGT` character.
#[pyfunction]
fn gc_count(seq: &str) -> PyResult<usize> {
    ensure_acgt(seq)?;
    Ok(seq
        .as_bytes()
        .iter()
        .filter(|&&b| b == b'G' || b == b'C')
        .count())
}

/// Return the length of the longest run of a single repeated base.
///
/// Returns `0` for the empty sequence. Raises `ValueError` on any non-`ACGT`
/// character.
#[pyfunction]
fn max_homopolymer_run(seq: &str) -> PyResult<usize> {
    ensure_acgt(seq)?;
    let bytes = seq.as_bytes();
    let mut best: usize = 0;
    let mut run: usize = 0;
    let mut prev: u8 = 0;
    for &b in bytes {
        if b == prev {
            run += 1;
        } else {
            run = 1;
            prev = b;
        }
        if run > best {
            best = run;
        }
    }
    Ok(best)
}

/// The `bt4_native` extension module.
#[pymodule]
fn bt4_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(gc_count, m)?)?;
    m.add_function(wrap_pyfunction!(max_homopolymer_run, m)?)?;
    Ok(())
}
