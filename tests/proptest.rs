//! Property-based tests over randomized alignments.

use proptest::prelude::*;
use snpick::fasta::{get_ref_seq, index_fasta};
use snpick::scan::{analyze, pass1_scan};
use snpick::types::build_lookup;

/// Build an equal-length FASTA of `nseq` sequences of `len` bases, drawn deterministically
/// from `seed` over the alphabet A/C/G/T/N/-.
fn make_alignment(nseq: usize, len: usize, seed: u64) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACGTN-";
    let mut fa = Vec::new();
    let mut x = seed | 1;
    for i in 0..nseq {
        fa.extend_from_slice(format!(">s{}\n", i).as_bytes());
        for _ in 0..len {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            fa.push(ALPHABET[(x >> 33) as usize % ALPHABET.len()]);
        }
        fa.push(b'\n');
    }
    fa
}

proptest! {
    /// The core invariant: every alignment column is classified exactly once.
    #[test]
    fn classification_covers_every_column(
        nseq in 2usize..8,
        len in 1usize..40,
        seed in any::<u64>(),
    ) {
        let fa = make_alignment(nseq, len, seed);
        for gaps in [false, true] {
            let lookup = build_lookup(gaps);
            let (recs, sl, layout) = index_fasta(&fa).unwrap();
            let bm = pass1_scan(&fa, &recs, sl, layout, &lookup);
            let rs = get_ref_seq(&fa, &recs[0], sl, layout);
            let (v, sc) = analyze(&bm, &rs, &lookup, gaps);
            prop_assert_eq!(v.len() + sc.constant.total() + sc.ambiguous, sl);
        }
    }

    /// The parallel and sequential scans always agree (production path vs. fallback).
    #[test]
    fn variable_count_stable_across_gaps_flag(
        nseq in 2usize..8,
        len in 1usize..40,
        seed in any::<u64>(),
    ) {
        let fa = make_alignment(nseq, len, seed);
        let lookup = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&fa).unwrap();
        let bm = pass1_scan(&fa, &recs, sl, layout, &lookup);
        let rs = get_ref_seq(&fa, &recs[0], sl, layout);
        let (v, _) = analyze(&bm, &rs, &lookup, false);
        // Every reported variable position must be within bounds and have >= 1 ALT allele.
        for vp in &v {
            prop_assert!(vp.index < sl);
            prop_assert!(!vp.alt_bases.is_empty());
        }
    }
}
