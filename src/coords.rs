//! Reference-anchored coordinates and BED region masking.

use std::io;

/// Ungapped reference position (1-based) for each alignment column. A column where the
/// reference base is a gap inherits the position of the preceding non-gap base.
pub fn ref_positions(ref_seq: &[u8]) -> Vec<u32> {
    let mut out = Vec::with_capacity(ref_seq.len());
    let mut rp: u32 = 0;
    for &b in ref_seq {
        if b != b'-' {
            rp += 1;
        }
        out.push(rp);
    }
    out
}

/// Parse a BED file into (start, end) 0-based half-open intervals (the chrom column is ignored;
/// snpick works over a single alignment).
pub fn parse_bed(path: &str) -> io::Result<Vec<(usize, usize)>> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        io::Error::new(e.kind(), format!("Cannot read BED '{}': {}", path, e))
    })?;
    let mut ivs = Vec::new();
    for (n, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip BED header lines, but only when the FIRST token is exactly `track`/`browser`
        // (not a chrom that merely starts with those letters, e.g. `track_scaffold_1`).
        let first = line.split_whitespace().next().unwrap_or("");
        if first == "track" || first == "browser" {
            continue;
        }
        let mut f = line.split('\t');
        let _chrom = f.next();
        let start = f.next().and_then(|s| s.trim().parse::<usize>().ok());
        let end = f.next().and_then(|s| s.trim().parse::<usize>().ok());
        match (start, end) {
            (Some(s), Some(e)) if s <= e => ivs.push((s, e)),
            _ => {
                return Err(io::Error::new(io::ErrorKind::InvalidData, format!(
                    "Malformed BED at line {}: expected 'chrom<TAB>start<TAB>end' with start <= end.",
                    n + 1
                )))
            }
        }
    }
    Ok(ivs)
}

/// Build a per-column boolean mask (true = masked). When `ref_coords` is `Some`, the intervals
/// are 0-based half-open reference positions; otherwise they are alignment columns.
pub fn build_mask(
    intervals: &[(usize, usize)], seq_length: usize, ref_coords: Option<&[u32]>,
) -> Vec<bool> {
    let mut mask = vec![false; seq_length];
    match ref_coords {
        None => {
            for &(s, e) in intervals {
                let (s, e) = (s.min(seq_length), e.min(seq_length));
                for m in &mut mask[s..e] {
                    *m = true;
                }
            }
        }
        Some(rp) => {
            for &(s, e) in intervals {
                // 0-based half-open [s, e) over the reference == 1-based positions s+1 ..= e.
                // Compare in usize (saturating on the +1) so an out-of-range or degenerate
                // interval never wraps or truncates to u32 and masks the wrong columns — the
                // None branch above is likewise robust via `.min(seq_length)`. Reference
                // positions are always <= MAX_SEQ_LENGTH, so widening `p` loses nothing.
                let lo = s.saturating_add(1);
                for (col, &p) in rp.iter().enumerate() {
                    let p = p as usize;
                    if p >= lo && p <= e {
                        mask[col] = true;
                    }
                }
            }
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ref_positions_gaps_inherit() {
        // A - C G  ->  positions 1 1 2 3
        assert_eq!(ref_positions(b"A-CG"), vec![1, 1, 2, 3]);
        assert_eq!(ref_positions(b"ACGT"), vec![1, 2, 3, 4]);
    }

    #[test]
    fn mask_alignment_columns() {
        // columns [1,3) masked
        let m = build_mask(&[(1, 3)], 5, None);
        assert_eq!(m, vec![false, true, true, false, false]);
    }

    #[test]
    fn mask_reference_coordinates() {
        // ref A-CG -> positions [1,1,2,3]; mask ref [1,2) = ref position 2 -> column 2
        let rp = ref_positions(b"A-CG");
        let m = build_mask(&[(1, 2)], 4, Some(&rp));
        assert_eq!(m, vec![false, false, true, false]);
    }

    #[test]
    fn mask_reference_out_of_range_is_noop() {
        // Reference positions are small; an interval far beyond them (or a degenerate empty one
        // at usize::MAX) must mask nothing. The old `as u32` cast wrapped/truncated these into
        // the valid range and masked the wrong columns (and panicked on the +1 in debug).
        let rp = ref_positions(b"ACGT"); // positions [1,2,3,4]
        let huge = 1usize << 40; // ≡ 0 mod 2^32
        assert_eq!(build_mask(&[(huge, huge + 4)], 4, Some(&rp)), vec![false; 4]);
        assert_eq!(build_mask(&[(usize::MAX, usize::MAX)], 4, Some(&rp)), vec![false; 4]);
        // An in-range interval still masks correctly.
        assert_eq!(build_mask(&[(1, 3)], 4, Some(&rp)), vec![false, true, true, false]);
    }
}
