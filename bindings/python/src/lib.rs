//! Python bindings for snpick, built on the `snpick` library crate.
//!
//! Build a wheel with maturin:
//!   cd bindings/python && maturin develop --release
//!
//! Then in Python:
//!   >>> import snpick
//!   >>> snpick.stats("alignment.fasta")
//!   {'sequences': 250, 'alignment_length': 4411532, 'variable_sites': 158934, ...}

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use snpick_core::fasta::{get_ref_seq, index_fasta};
use snpick_core::input::map_input;
use snpick_core::scan::{analyze, pass1_scan};
use snpick_core::types::build_lookup;

/// Compute variable-site statistics for a FASTA alignment (gzip and '-' for stdin are accepted).
///
/// Returns a dict with `sequences`, `alignment_length`, `variable_sites`, `constant_sites`,
/// `ambiguous_sites`, and the `fconst` tuple (A, C, G, T).
#[pyfunction]
#[pyo3(signature = (path, include_gaps = false))]
fn stats(py: Python<'_>, path: &str, include_gaps: bool) -> PyResult<Py<PyDict>> {
    let input = map_input(path).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let data = &input.mmap[..];
    let lookup = build_lookup(include_gaps);
    let (records, seq_len, layout) =
        index_fasta(data).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let bm = pass1_scan(data, &records, seq_len, layout, &lookup);
    let rs = get_ref_seq(data, &records[0], seq_len, layout);
    let (v, sc) = analyze(&bm, &rs, &lookup, include_gaps);

    let d = PyDict::new(py);
    d.set_item("sequences", records.len())?;
    d.set_item("alignment_length", seq_len)?;
    d.set_item("variable_sites", v.len())?;
    d.set_item("constant_sites", sc.constant.total())?;
    d.set_item("ambiguous_sites", sc.ambiguous)?;
    d.set_item("fconst", (sc.constant.a, sc.constant.c, sc.constant.g, sc.constant.t))?;
    Ok(d.unbind())
}

#[pymodule]
fn snpick(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(stats, m)?)?;
    Ok(())
}
