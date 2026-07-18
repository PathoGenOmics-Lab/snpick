#![no_main]
//! Fuzz the FASTA indexer: `index_fasta` does manual index arithmetic over untrusted,
//! machine-generated input, so it must never panic — only return `Ok`/`Err`.
//!
//! Run with a nightly toolchain and cargo-fuzz:
//!   cargo +nightly fuzz run index_fasta

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok((records, seq_len, layout)) = snpick::fasta::index_fasta(data) {
        // Exercise the reference reader on the first record too.
        let _ = snpick::fasta::get_ref_seq(data, &records[0], seq_len, layout);
    }
});
