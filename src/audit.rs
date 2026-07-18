//! Alignment composition audit: a per-byte histogram over the sequence data, used by
//! `--on-invalid` and `--check`. Only run on demand (it is a full extra pass).

use crate::fasta::FastaRecord;
use crate::types::SeqLayout;

/// Count every sequence byte (newlines skipped) into a 256-entry histogram.
pub fn composition(
    data: &[u8], records: &[FastaRecord], seq_length: usize, layout: SeqLayout,
) -> [u64; 256] {
    let mut h = [0u64; 256];
    for rec in records {
        if layout.single_line {
            for &b in &data[rec.seq_offset..rec.seq_offset + seq_length] {
                h[b as usize] += 1;
            }
        } else {
            let mut pos = rec.seq_offset;
            let end = data.len();
            let mut i = 0;
            while i < seq_length && pos < end {
                let b = data[pos];
                pos += 1;
                if b == b'\n' || b == b'\r' {
                    continue;
                }
                h[b as usize] += 1;
                i += 1;
            }
        }
    }
    h
}

/// Whether a byte is part of the accepted nucleotide alphabet (ACGTU, N, IUPAC ambiguity
/// codes, and gap characters). Anything else (e.g. protein residues) is "invalid".
pub fn is_valid_nucleotide(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N'
            | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
            | b'-' | b'.' | b'?'
    )
}

/// Total count of out-of-alphabet bytes.
pub fn invalid_count(h: &[u64; 256]) -> u64 {
    (0..256u16).filter(|&i| !is_valid_nucleotide(i as u8)).map(|i| h[i as usize]).sum()
}

/// A short human-readable composition breakdown (for `--check`).
pub fn summary(h: &[u64; 256]) -> String {
    let sum = |bytes: &[u8]| -> u64 { bytes.iter().map(|&b| h[b as usize]).sum() };
    let total: u64 = h.iter().sum();
    let acgt = sum(b"ACGTacgt");
    let n = sum(b"Nn");
    let gap = sum(b"-.");
    let iupac = sum(b"RYSWKMBDHVryswkmbdhvUu");
    let invalid = invalid_count(h);
    let pct = |c: u64| if total == 0 { 0.0 } else { 100.0 * c as f64 / total as f64 };
    format!(
        "total={} A/C/G/T={} ({:.2}%) N={} ({:.2}%) gap={} ({:.2}%) IUPAC={} ({:.2}%) invalid={} ({:.2}%)",
        total, acgt, pct(acgt), n, pct(n), gap, pct(gap), iupac, pct(iupac), invalid, pct(invalid),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_detection() {
        let mut h = [0u64; 256];
        h[b'A' as usize] = 10;
        h[b'N' as usize] = 2;
        h[b'-' as usize] = 1;
        h[b'E' as usize] = 3; // protein residue -> invalid
        assert_eq!(invalid_count(&h), 3);
        assert!(is_valid_nucleotide(b'R'));
        assert!(!is_valid_nucleotide(b'E'));
        assert!(summary(&h).contains("invalid=3"));
    }
}
