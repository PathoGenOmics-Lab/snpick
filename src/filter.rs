//! Per-site filtering: recount alleles and missingness over the variable columns.
//!
//! Filtering removes *variable* sites from the output only — it never reclassifies a site as
//! constant, so the ASC `fconst` counts (constant sites) are untouched and stay valid.

use crate::fasta::FastaRecord;
use crate::types::{SeqLayout, BIT_A, BIT_C, BIT_G, BIT_T, BIT_GAP};

/// Per-variable-site allele and missingness counts across all samples.
#[derive(Clone, Copy, Default)]
pub struct SiteStat {
    /// Counts of A, C, G, T.
    pub counts: [u32; 4],
    /// Count of gaps ('-').
    pub gap: u32,
    /// Count of ambiguous / uncalled bases.
    pub missing: u32,
}

/// Thresholds for per-site filtering. `None` = no constraint.
#[derive(Clone, Copy, Default)]
pub struct SiteFilters {
    /// Max allowed fraction of missing genotypes (0.0–1.0).
    pub max_missing: Option<f64>,
    /// Min minor-allele count.
    pub mac: Option<u32>,
    /// Min minor-allele frequency (0.0–1.0).
    pub maf: Option<f64>,
    /// Min number of samples with data.
    pub min_samples: Option<u32>,
    /// Max number of distinct alleles (e.g. 2 = biallelic only).
    pub max_alleles: Option<u32>,
}

impl SiteFilters {
    /// Whether any threshold is set (skip the counting pass entirely if not).
    pub fn active(&self) -> bool {
        self.max_missing.is_some()
            || self.mac.is_some()
            || self.maf.is_some()
            || self.min_samples.is_some()
            || self.max_alleles.is_some()
    }
}

impl SiteStat {
    /// Number of samples with a called allele (gaps count only under `include_gaps`).
    pub fn ns(&self, include_gaps: bool) -> u32 {
        let acgt: u32 = self.counts.iter().sum();
        acgt + if include_gaps { self.gap } else { 0 }
    }

    /// Counts of the alleles actually present (ACGT with count > 0, plus gap if `include_gaps`).
    fn allele_counts(&self, include_gaps: bool) -> Vec<u32> {
        let mut v: Vec<u32> = self.counts.iter().copied().filter(|&c| c > 0).collect();
        if include_gaps && self.gap > 0 {
            v.push(self.gap);
        }
        v
    }

    /// Whether this site passes every set threshold.
    pub fn passes(&self, f: &SiteFilters, num_samples: u32, include_gaps: bool) -> bool {
        let ns = self.ns(include_gaps);
        let missing = num_samples - ns; // each sample lands in exactly one bucket, no underflow
        if let Some(ms) = f.min_samples {
            if ns < ms {
                return false;
            }
        }
        if let Some(mm) = f.max_missing {
            if num_samples > 0 && (missing as f64 / num_samples as f64) > mm {
                return false;
            }
        }
        let alleles = self.allele_counts(include_gaps);
        if let Some(ma) = f.max_alleles {
            if alleles.len() as u32 > ma {
                return false;
            }
        }
        let mac = alleles.iter().copied().min().unwrap_or(0);
        if let Some(m) = f.mac {
            if mac < m {
                return false;
            }
        }
        if let Some(maf) = f.maf {
            if ns == 0 || (mac as f64 / ns as f64) < maf {
                return false;
            }
        }
        true
    }
}

#[inline]
fn tally(s: &mut SiteStat, bits: u8) {
    // `bits` comes from the lookup table, so it is a single flag (or 0 for ambiguous).
    match bits {
        BIT_A => s.counts[0] += 1,
        BIT_C => s.counts[1] += 1,
        BIT_G => s.counts[2] += 1,
        BIT_T => s.counts[3] += 1,
        BIT_GAP => s.gap += 1,
        _ => s.missing += 1,
    }
}

/// Count alleles and missing calls at each variable position (sparse pass over just the
/// variable columns). Gap handling follows `lookup` (built with the same `include_gaps`).
pub fn count_sites(
    data: &[u8], records: &[FastaRecord], pos_indices: &[usize],
    layout: SeqLayout, lookup: &[u8; 256],
) -> Vec<SiteStat> {
    let num_var = pos_indices.len();
    let mut stats = vec![SiteStat::default(); num_var];
    for rec in records {
        if layout.single_line {
            let base = rec.seq_offset;
            for (vi, &p) in pos_indices.iter().enumerate() {
                tally(&mut stats[vi], lookup[data[base + p] as usize]);
            }
        } else {
            let mut pos = rec.seq_offset;
            let end = data.len();
            let mut base_idx = 0usize;
            let mut var_idx = 0usize;
            while var_idx < num_var && pos < end {
                let b = data[pos];
                pos += 1;
                if b == b'\n' || b == b'\r' {
                    continue;
                }
                if base_idx == pos_indices[var_idx] {
                    tally(&mut stats[var_idx], lookup[b as usize]);
                    var_idx += 1;
                }
                base_idx += 1;
            }
        }
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stat(a: u32, c: u32, g: u32, t: u32, gap: u32, missing: u32) -> SiteStat {
        SiteStat { counts: [a, c, g, t], gap, missing }
    }

    #[test]
    fn mac_drops_singletons() {
        let f = SiteFilters { mac: Some(2), ..Default::default() };
        // A:9 T:1 -> minor-allele count 1, dropped by --mac 2
        assert!(!stat(9, 0, 0, 1, 0, 0).passes(&f, 10, false));
        // A:8 T:2 -> minor-allele count 2, kept
        assert!(stat(8, 0, 0, 2, 0, 0).passes(&f, 10, false));
    }

    #[test]
    fn max_missing_and_min_samples() {
        // 3 of 10 missing = 0.3
        let s = stat(5, 2, 0, 0, 0, 3);
        assert!(s.passes(&SiteFilters { max_missing: Some(0.3), ..Default::default() }, 10, false));
        assert!(!s.passes(&SiteFilters { max_missing: Some(0.2), ..Default::default() }, 10, false));
        assert!(s.passes(&SiteFilters { min_samples: Some(7), ..Default::default() }, 10, false));
        assert!(!s.passes(&SiteFilters { min_samples: Some(8), ..Default::default() }, 10, false));
    }

    #[test]
    fn max_alleles_biallelic() {
        let f = SiteFilters { max_alleles: Some(2), ..Default::default() };
        assert!(stat(5, 5, 0, 0, 0, 0).passes(&f, 10, false)); // 2 alleles
        assert!(!stat(4, 3, 3, 0, 0, 0).passes(&f, 10, false)); // 3 alleles
    }

    #[test]
    fn maf_frequency() {
        let f = SiteFilters { maf: Some(0.15), ..Default::default() };
        assert!(!stat(9, 0, 0, 1, 0, 0).passes(&f, 10, false)); // maf 0.1
        assert!(stat(8, 0, 0, 2, 0, 0).passes(&f, 10, false)); // maf 0.2
    }
}
