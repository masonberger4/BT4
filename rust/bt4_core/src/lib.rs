//! BT4 native accelerator.
//!
//! Small, correct, deterministic hot-loop primitives operating on uppercase
//! `ACGT` DNA strings. A pure-Python fallback exists separately, so this module
//! keeps its surface intentionally minimal. No `unsafe`, no external crates
//! beyond `pyo3`.

use std::collections::HashMap;

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

/// Return the length of the longest run of consecutive `{G, C}` bases.
///
/// A "GC run" is a maximal stretch of positions whose base is `G` or `C`; the
/// bases may be mixed, so `GCGC` counts as a run of four (this is the "GC
/// length" semantics, distinct from a single-base homopolymer). Returns `0` for
/// the empty sequence. Raises `ValueError` on any non-`ACGT` character.
///
/// Mirrors `bt4.constraints.gc_run._max_gc_run` exactly on `ACGT` input.
#[pyfunction]
fn max_gc_run(seq: &str) -> PyResult<usize> {
    ensure_acgt(seq)?;
    let mut best: usize = 0;
    let mut run: usize = 0;
    for &b in seq.as_bytes() {
        if b == b'G' || b == b'C' {
            run += 1;
            if run > best {
                best = run;
            }
        } else {
            run = 0;
        }
    }
    Ok(best)
}

/// Return the length of the longest reverse-complement-aware repeat in `seq`.
///
/// This is the largest `L` such that some length-`L` substring occurs at two
/// distinct start positions (a *direct* repeat) OR some length-`L` substring's
/// reverse complement occurs anywhere in `seq` (an *inverted* repeat; a
/// *palindrome* when a substring equals its own reverse complement). Returns `0`
/// when no substring of length `>= 1` repeats in that sense.
///
/// It is computed as the maximum of two longest-common-substring dynamic
/// programs: `seq` against itself off the main diagonal (direct repeats, overlaps
/// allowed) and `seq` against its own reverse complement (inverted repeats /
/// palindromes). Because "offending" is monotone in length, this matches
/// `MaxRepeatConstraint`'s notion exactly: `longest_repeat(seq) > max_length`
/// iff `MaxRepeatConstraint(max_length).validate(seq)` yields a hard violation.
///
/// O(n^2) time, O(n) extra space. Raises `ValueError` on any non-`ACGT`
/// character.
#[pyfunction]
fn longest_repeat(seq: &str) -> PyResult<usize> {
    ensure_acgt(seq)?;
    let s = seq.as_bytes();
    let n = s.len();
    if n == 0 {
        return Ok(0);
    }
    // Reverse complement of `s` (valid because `ensure_acgt` passed).
    let mut rc = vec![0u8; n];
    for i in 0..n {
        rc[i] = complement_base(s[n - 1 - i]).expect("ensure_acgt guarantees ACGT");
    }
    let mut best: usize = 0;
    let mut prev = vec![0usize; n];
    let mut cur = vec![0usize; n];
    // Direct repeats: longest common substring of `s` with itself, excluding the
    // `i == j` main diagonal so the two copies sit at distinct positions.
    for i in 0..n {
        for j in 0..n {
            if i != j && s[i] == s[j] {
                let d = if j > 0 { prev[j - 1] } else { 0 };
                let v = d + 1;
                cur[j] = v;
                if v > best {
                    best = v;
                }
            } else {
                cur[j] = 0;
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    // Inverted / palindromic repeats: longest common substring of `s` and `rc`.
    // X occurs in both iff X is in `s` and reverse_complement(X) is in `s`.
    for x in prev.iter_mut() {
        *x = 0;
    }
    for i in 0..n {
        for j in 0..n {
            if s[i] == rc[j] {
                let d = if j > 0 { prev[j - 1] } else { 0 };
                let v = d + 1;
                cur[j] = v;
                if v > best {
                    best = v;
                }
            } else {
                cur[j] = 0;
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    Ok(best)
}

/// Run the exact codon-trellis layer DP over precomputed, position-independent
/// transition tables and return the maximum-scalar coding sequence.
///
/// This is the Rust port of the inner loop of `bt4.optimize.exact_dp.solve_exact`.
/// Because the DP is callback-driven in Python (`scalar_delta` and each
/// constraint's `ok_suffix`) and Rust must never call back into Python, the
/// caller *precomputes* every transition in Python — but only in the regime where
/// the objective is position-independent, so each `(context, codon)` transition
/// has a single well-defined `(allowed, delta, new_context)`. The caller then
/// passes, per trellis layer, the flat lists of the **allowed** transitions:
///
/// * `codons` — the global codon-string table, indexed by codon id.
/// * `layer_from[l]` — source context id of each allowed transition in layer `l`.
/// * `layer_to[l]` — destination context id (the merged trailing-context key).
/// * `layer_codon[l]` — codon id placed by the transition (index into `codons`).
/// * `layer_delta[l]` — the pre-summed scalar objective delta of the transition,
///   computed in Python so the float summation order is fixed and the result is
///   bit-for-bit identical to the pure-Python DP.
///
/// Context id `0` is the empty start context. The DP keeps, per context id, the
/// best `(score, dna)` — highest scalar, ties broken toward the lexicographically
/// smaller DNA, exactly as `solve_exact._wins` does. `beam`, when set, keeps at
/// most that many contexts per layer (highest scalar first, ties by smaller DNA)
/// and the returned `pruned` flag reports whether any layer was truncated.
///
/// Returns `None` when some layer has no reachable context (infeasible), so the
/// caller can raise `InfeasibleError` with the offending constraint names; else
/// `(best_dna, best_scalar, pruned)`.
#[pyfunction]
#[pyo3(signature = (codons, layer_from, layer_to, layer_codon, layer_delta, beam=None))]
fn trellis_solve(
    codons: Vec<String>,
    layer_from: Vec<Vec<u32>>,
    layer_to: Vec<Vec<u32>>,
    layer_codon: Vec<Vec<u32>>,
    layer_delta: Vec<Vec<f64>>,
    beam: Option<usize>,
) -> PyResult<Option<(String, f64, bool)>> {
    let n_layers = layer_from.len();
    if layer_to.len() != n_layers
        || layer_codon.len() != n_layers
        || layer_delta.len() != n_layers
    {
        return Err(PyValueError::new_err(
            "layer_from/layer_to/layer_codon/layer_delta must have equal length",
        ));
    }

    // The trellis starts with a single empty-context state (id 0, DNA "").
    let mut cur: HashMap<u32, (f64, String)> = HashMap::new();
    cur.insert(0u32, (0.0f64, String::new()));
    let mut pruned = false;

    for li in 0..n_layers {
        let froms = &layer_from[li];
        let tos = &layer_to[li];
        let cods = &layer_codon[li];
        let dels = &layer_delta[li];
        if tos.len() != froms.len() || cods.len() != froms.len() || dels.len() != froms.len() {
            return Err(PyValueError::new_err(
                "per-layer transition lists must have equal length",
            ));
        }
        let mut next: HashMap<u32, (f64, String)> = HashMap::new();
        for t in 0..froms.len() {
            if let Some((score, dna)) = cur.get(&froms[t]) {
                let ns = score + dels[t];
                let codon = &codons[cods[t] as usize];
                let mut nd = String::with_capacity(dna.len() + codon.len());
                nd.push_str(dna);
                nd.push_str(codon);
                let to = tos[t];
                // Merge rule identical to solve_exact._wins: keep the strictly
                // higher scalar, break exact ties toward the smaller DNA.
                let replace = match next.get(&to) {
                    None => true,
                    Some((cs, cd)) => ns > *cs || (ns == *cs && nd < *cd),
                };
                if replace {
                    next.insert(to, (ns, nd));
                }
            }
        }
        if next.is_empty() {
            return Ok(None);
        }
        if let Some(b) = beam {
            if next.len() > b {
                // Keep the `b` highest-scalar states, ties toward smaller DNA -
                // the same key as solve_exact's `sorted(..., key=(-score, dna))`.
                let mut items: Vec<(u32, (f64, String))> = next.into_iter().collect();
                items.sort_by(|a, c| {
                    match c.1 .0.partial_cmp(&a.1 .0).unwrap_or(std::cmp::Ordering::Equal) {
                        std::cmp::Ordering::Equal => a.1 .1.cmp(&c.1 .1),
                        other => other,
                    }
                });
                items.truncate(b);
                next = items.into_iter().collect();
                pruned = true;
            }
        }
        cur = next;
    }

    // Best over the final layer: highest scalar, ties toward smaller DNA.
    let mut best: Option<(f64, String)> = None;
    for (_k, (s, d)) in cur.into_iter() {
        let replace = match &best {
            None => true,
            Some((bs, bd)) => s > *bs || (s == *bs && d < *bd),
        };
        if replace {
            best = Some((s, d));
        }
    }
    Ok(best.map(|(s, d)| (d, s, pruned)))
}

/// The `bt4_native` extension module.
#[pymodule]
fn bt4_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(gc_count, m)?)?;
    m.add_function(wrap_pyfunction!(max_homopolymer_run, m)?)?;
    m.add_function(wrap_pyfunction!(max_gc_run, m)?)?;
    m.add_function(wrap_pyfunction!(longest_repeat, m)?)?;
    m.add_function(wrap_pyfunction!(trellis_solve, m)?)?;
    Ok(())
}
