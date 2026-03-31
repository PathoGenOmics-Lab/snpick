/// Bitmask flags for nucleotide classification.
pub const BIT_A: u8 = 0b00001;
pub const BIT_C: u8 = 0b00010;
pub const BIT_G: u8 = 0b00100;
pub const BIT_T: u8 = 0b01000;
pub const BIT_GAP: u8 = 0b10000;

/// Maximum alignment length (prevents OOM on malicious input).
pub const MAX_SEQ_LENGTH: usize = 10_000_000_000;

/// Maximum VCF genotype matrix size in bytes.
pub const MAX_VCF_GENO_BYTES: usize = 4_000_000_000;

/// I/O buffer size for BufWriter (16 MB).
pub const IO_BUF: usize = 16 * 1024 * 1024;

/// Whether all sequences are single-line (no embedded newlines).
/// When true, `data[seq_offset + pos]` gives the base at `pos` directly.
/// When false, we must scan skipping newlines for each record.
#[derive(Clone, Copy)]
pub struct SeqLayout {
    pub single_line: bool,
}

/// A variable position detected in the alignment.
pub struct VariablePosition {
    pub index: usize,
    pub ref_base: u8,
    pub alt_bases: Vec<u8>,
    pub ns: usize,
}

/// Counts of constant sites by nucleotide.
pub struct ConstantSiteCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
}

impl std::fmt::Display for ConstantSiteCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A:{} C:{} G:{} T:{}", self.a, self.c, self.g, self.t)
    }
}

impl ConstantSiteCounts {
    pub fn total(&self) -> usize { self.a + self.c + self.g + self.t }
    pub fn fconst(&self) -> String { format!("{},{},{},{}", self.a, self.c, self.g, self.t) }
}

/// Summary of site classification across the alignment.
pub struct SiteCounts {
    pub constant: ConstantSiteCounts,
    pub variable: usize,
    pub ambiguous: usize,
}

/// Build nucleotide → bitmask lookup table.
pub fn build_lookup(include_gaps: bool) -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = BIT_A; t[b'a' as usize] = BIT_A;
    t[b'C' as usize] = BIT_C; t[b'c' as usize] = BIT_C;
    t[b'G' as usize] = BIT_G; t[b'g' as usize] = BIT_G;
    t[b'T' as usize] = BIT_T; t[b't' as usize] = BIT_T;
    if include_gaps { t[b'-' as usize] = BIT_GAP; }
    t
}

/// Build lowercase → uppercase lookup table.
pub fn build_upper() -> [u8; 256] {
    let mut t = [0u8; 256];
    for (i, v) in t.iter_mut().enumerate() { *v = i as u8; }
    for c in b'a'..=b'z' { t[c as usize] = c - 32; }
    t
}

/// Convert bitmask to sorted list of bases it represents.
pub fn bits_to_bases(bits: u8, include_gaps: bool) -> Vec<u8> {
    let mut v = Vec::with_capacity(5);
    if bits & BIT_A != 0 { v.push(b'A'); }
    if bits & BIT_C != 0 { v.push(b'C'); }
    if bits & BIT_G != 0 { v.push(b'G'); }
    if bits & BIT_T != 0 { v.push(b'T'); }
    if include_gaps && bits & BIT_GAP != 0 { v.push(b'-'); }
    v
}
