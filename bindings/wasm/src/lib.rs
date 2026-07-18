//! WebAssembly bindings for snpick. The `snpick` library is compiled with
//! `default-features = false`, so it drops rayon (no threads) and memmap (no mmap) and runs the
//! sequential in-memory scan — suitable for the browser / Node.
//!
//! Build with wasm-pack:
//!   cd bindings/wasm && wasm-pack build --target web

use wasm_bindgen::prelude::*;

use snpick_core::fasta::{get_ref_seq, index_fasta};
use snpick_core::scan::{analyze, pass1_scan};
use snpick_core::types::build_lookup;

/// Count the variable (SNP) sites in an in-memory FASTA alignment. Returns 0 on parse error.
#[wasm_bindgen]
pub fn variable_site_count(fasta: &[u8], include_gaps: bool) -> usize {
    let lookup = build_lookup(include_gaps);
    let (records, seq_len, layout) = match index_fasta(fasta) {
        Ok(x) => x,
        Err(_) => return 0,
    };
    let bitmask = pass1_scan(fasta, &records, seq_len, layout, &lookup);
    let ref_seq = get_ref_seq(fasta, &records[0], seq_len, layout);
    let (variable, _) = analyze(&bitmask, &ref_seq, &lookup, include_gaps);
    variable.len()
}

/// The four ASC `fconst` constant-site counts (A, C, G, T) for an in-memory FASTA alignment.
#[wasm_bindgen]
pub fn fconst(fasta: &[u8], include_gaps: bool) -> Vec<u32> {
    let lookup = build_lookup(include_gaps);
    let (records, seq_len, layout) = match index_fasta(fasta) {
        Ok(x) => x,
        Err(_) => return vec![0, 0, 0, 0],
    };
    let bitmask = pass1_scan(fasta, &records, seq_len, layout, &lookup);
    let ref_seq = get_ref_seq(fasta, &records[0], seq_len, layout);
    let (_, counts) = analyze(&bitmask, &ref_seq, &lookup, include_gaps);
    let c = &counts.constant;
    vec![c.a as u32, c.c as u32, c.g as u32, c.t as u32]
}
