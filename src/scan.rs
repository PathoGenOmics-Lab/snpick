//! Pass 1: bitmask scan and site classification.
//!
//! Builds a per-position bitmask by OR-ing each sequence's nucleotide flags,
//! then classifies positions as variable (>1 allele), constant, or ambiguous.

use crate::fasta::FastaRecord;
use crate::types::*;

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
    let mut bitmask = vec![0u8; seq_length];

    // Prefault all pages into RAM before the hot loop
    prefault(data);

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

    bitmask
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
        } else {
            ambiguous += 1;
        }
    }

    let num_variable = vars.len();
    (vars, SiteCounts { constant: cs, variable: num_variable, ambiguous })
}
