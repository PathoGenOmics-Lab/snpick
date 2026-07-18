//! Pass 1: bitmask scan and site classification.
//!
//! Builds a per-position bitmask by OR-ing each sequence's nucleotide flags,
//! then classifies positions as variable (>1 allele), constant, or ambiguous.

use rayon::prelude::*;

use crate::fasta::FastaRecord;
use crate::types::*;

/// Minimum total scan work (records × positions) before the parallel path is
/// worth its overhead — roughly 50 seqs × 4M bp.
const PARALLEL_MIN_WORK: usize = 200_000_000;

/// Prefault mmap pages by touching one byte per OS page.
/// Eliminates soft page faults during the scan loop (~0.5s on 1 GB files).
#[inline(never)]
fn prefault(data: &[u8]) {
    const PAGE: usize = 4096;
    let mut sum = 0u8;
    let mut off = 0;
    while off < data.len() {
        sum = sum.wrapping_add(data[off]);
        off += PAGE;
    }
    std::hint::black_box(sum);
}

/// Pass 1: build bitmask of observed nucleotides at each position.
///
/// Iterates all sequences, OR-ing each base's lookup value into the bitmask.
/// Prefaults mmap pages first, then scans sequentially per sequence.
/// For multi-line FASTA, scans byte-by-byte skipping newlines.
pub fn pass1_scan(
    data: &[u8], records: &[FastaRecord], seq_length: usize,
    layout: SeqLayout, lookup: &[u8; 256],
) -> Vec<u8> {
    // Prefault all pages into RAM before the hot loop
    prefault(data);

    // Parallelism only pays off when there's enough work per thread
    // (~200M bases, e.g. 50 seqs × 4M bp). Small inputs scan sequentially.
    let num_threads = rayon::current_num_threads().min(records.len());
    let total_work = records.len() * seq_length;
    if num_threads <= 1 || total_work < PARALLEL_MIN_WORK {
        let mut bitmask = vec![0u8; seq_length];
        scan_sequential(data, records, seq_length, layout, lookup, &mut bitmask);
        bitmask
    } else {
        scan_parallel(data, records, seq_length, layout, lookup, num_threads)
    }
}

/// Parallel scan: each thread scans a disjoint chunk of records into its own
/// bitmask, then all partials are merged with OR. OR is commutative and
/// associative over the disjoint chunks, so the result is byte-for-byte
/// identical to a sequential scan.
fn scan_parallel(
    data: &[u8], records: &[FastaRecord], seq_length: usize,
    layout: SeqLayout, lookup: &[u8; 256], num_threads: usize,
) -> Vec<u8> {
    let chunk_size = records.len().div_ceil(num_threads);
    let partial_bitmasks: Vec<Vec<u8>> = records
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_bm = vec![0u8; seq_length];
            scan_sequential(data, chunk, seq_length, layout, lookup, &mut local_bm);
            local_bm
        })
        .collect();

    let mut bitmask = vec![0u8; seq_length];
    for partial in &partial_bitmasks {
        for (bm, &p) in bitmask.iter_mut().zip(partial.iter()) {
            *bm |= p;
        }
    }
    bitmask
}

/// Sequential scan of a set of records into a bitmask.
fn scan_sequential(
    data: &[u8], records: &[FastaRecord], seq_length: usize,
    layout: SeqLayout, lookup: &[u8; 256], bitmask: &mut [u8],
) {
    if layout.single_line {
        for rec in records {
            let seq = &data[rec.seq_offset..rec.seq_offset + seq_length];
            for (bm_byte, &seq_byte) in bitmask.iter_mut().zip(seq.iter()) {
                *bm_byte |= lookup[seq_byte as usize];
            }
        }
    } else {
        for rec in records {
            let mut pos = rec.seq_offset;
            let end = data.len();
            let mut i = 0;
            while i < seq_length && pos < end {
                let b = data[pos];
                pos += 1;
                if b == b'\n' || b == b'\r' { continue; }
                bitmask[i] |= lookup[b as usize];
                i += 1;
            }
        }
    }
}

/// Classify positions into variable, constant, or ambiguous-only.
pub fn analyze(
    bitmask: &[u8], ref_seq: &[u8], lookup: &[u8; 256], include_gaps: bool,
) -> (Vec<VariablePosition>, SiteCounts) {
    let mut vars = Vec::new();
    let mut cs = ConstantSiteCounts { a: 0, c: 0, g: 0, t: 0 };
    let mut ambiguous = 0usize;

    for (pos, &bits) in bitmask.iter().enumerate() {
        let ones = bits.count_ones();
        if ones > 1 {
            let rb = ref_seq[pos].to_ascii_uppercase();
            let ref_base = if lookup[rb as usize] != 0 { rb } else { bits_to_bases(bits, include_gaps)[0] };
            let alt_bases: Vec<u8> = bits_to_bases(bits, include_gaps)
                .into_iter().filter(|&b| b != ref_base).collect();
            vars.push(VariablePosition { index: pos, ref_base, alt_bases, ns: 0 });
        } else if ones == 1 {
            if bits & BIT_A != 0 { cs.a += 1; }
            else if bits & BIT_C != 0 { cs.c += 1; }
            else if bits & BIT_G != 0 { cs.g += 1; }
            else if bits & BIT_T != 0 { cs.t += 1; }
            else { ambiguous += 1; } // all-gap column (BIT_GAP only) under --include-gaps
        } else {
            ambiguous += 1;
        }
    }

    let num_variable = vars.len();
    (vars, SiteCounts { constant: cs, variable: num_variable, ambiguous })
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::index_fasta;

    #[test]
    fn parallel_matches_sequential() {
        // The parallel chunk/merge path is the production path for real
        // (whole-genome) inputs but is gated behind a size threshold that no
        // small test reaches. Drive it directly and assert it equals a
        // sequential scan byte-for-byte.
        let mut fa = Vec::new();
        for i in 0..97usize {
            fa.extend_from_slice(format!(">s{}\n", i).as_bytes());
            let mut seq = vec![b'A'; 40];
            seq[i % 40] = b'C';
            seq[(i * 7) % 40] = b'G';
            seq[(i * 13) % 40] = b'T';
            fa.extend_from_slice(&seq);
            fa.push(b'\n');
        }
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&fa).unwrap();

        let mut seq_bm = vec![0u8; sl];
        scan_sequential(&fa, &recs, sl, layout, &lk, &mut seq_bm);

        let threads = rayon::current_num_threads().min(recs.len()).max(2);
        let par_bm = scan_parallel(&fa, &recs, sl, layout, &lk, threads);

        assert_eq!(seq_bm, par_bm);
    }
}
