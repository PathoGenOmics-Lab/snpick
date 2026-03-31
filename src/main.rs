use bio::io::fasta::{self, FastaRead};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

// =============================================================================
// CLI
// =============================================================================

#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = env!("CARGO_PKG_VERSION"),
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>",
    about = "A fast, memory-efficient tool for extracting variable sites from FASTA alignments and generating VCF output compatible with RAxML/IQ-TREE ASC."
)]
struct Args {
    /// Input FASTA alignment file
    #[arg(short, long)]
    fasta: String,

    /// Output FASTA file with variable sites
    #[arg(short, long)]
    output: String,

    /// Consider '-' (gap) as a distinct character in variant detection
    #[arg(short = 'g', long)]
    include_gaps: bool,

    /// Generate VCF file with standard GT genotypes
    #[arg(long)]
    vcf: bool,

    /// Output VCF file path (default: output.vcf)
    #[arg(long)]
    vcf_output: Option<String>,
}

// =============================================================================
// Lookup table: nucleotide byte → bitmask (branchless hot loop)
// =============================================================================

const BIT_A: u8 = 0b00001;
const BIT_C: u8 = 0b00010;
const BIT_G: u8 = 0b00100;
const BIT_T: u8 = 0b01000;
const BIT_GAP: u8 = 0b10000;

/// Maximum alignment length we'll accept (10 Gbp — larger than any known genome)
const MAX_SEQ_LENGTH: usize = 10_000_000_000;

fn build_lookup(include_gaps: bool) -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = BIT_A;
    t[b'a' as usize] = BIT_A;
    t[b'C' as usize] = BIT_C;
    t[b'c' as usize] = BIT_C;
    t[b'G' as usize] = BIT_G;
    t[b'g' as usize] = BIT_G;
    t[b'T' as usize] = BIT_T;
    t[b't' as usize] = BIT_T;
    if include_gaps {
        t[b'-' as usize] = BIT_GAP;
    }
    t
}

/// Decode bitmask bits into sorted nucleotide bytes
fn bits_to_bases(bits: u8, include_gaps: bool) -> Vec<u8> {
    let mut v = Vec::with_capacity(5);
    if bits & BIT_A != 0 { v.push(b'A'); }
    if bits & BIT_C != 0 { v.push(b'C'); }
    if bits & BIT_G != 0 { v.push(b'G'); }
    if bits & BIT_T != 0 { v.push(b'T'); }
    if include_gaps && bits & BIT_GAP != 0 { v.push(b'-'); }
    v
}

/// Check that two paths don't resolve to the same file (S3: prevent data loss)
fn check_paths_differ(input: &str, output: &str) -> io::Result<()> {
    let inp = std::fs::canonicalize(input)?;
    // Output may not exist yet — compare parent + filename
    let out_path = Path::new(output);
    let out_abs = if out_path.exists() {
        std::fs::canonicalize(output)?
    } else {
        let parent = out_path.parent().unwrap_or(Path::new("."));
        let parent_canon = std::fs::canonicalize(parent).unwrap_or(parent.to_path_buf());
        parent_canon.join(out_path.file_name().unwrap_or_default())
    };

    if inp == out_abs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Input and output paths resolve to the same file: {}. This would destroy your data.",
                inp.display()
            ),
        ));
    }
    Ok(())
}

/// Format elapsed time as human-readable string
fn format_time(secs: f64) -> String {
    if secs < 60.0 {
        format!("{:.2}s", secs)
    } else if secs < 3600.0 {
        format!("{}m {:.0}s", (secs / 60.0) as u64, secs % 60.0)
    } else {
        let h = (secs / 3600.0) as u64;
        let m = ((secs % 3600.0) / 60.0) as u64;
        format!("{}h {}m", h, m)
    }
}

// =============================================================================
// Data structures
// =============================================================================

struct VariablePosition {
    index: usize,
    ref_base: u8,
    alt_bases: Vec<u8>,
}

struct ConstantSiteCounts {
    a: usize,
    c: usize,
    g: usize,
    t: usize,
}

impl std::fmt::Display for ConstantSiteCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A:{} C:{} G:{} T:{}", self.a, self.c, self.g, self.t)
    }
}

impl ConstantSiteCounts {
    fn total(&self) -> usize {
        self.a + self.c + self.g + self.t
    }

    /// Format for IQ-TREE --fconst / RAxML --asc-corr
    fn fconst(&self) -> String {
        format!("{},{},{},{}", self.a, self.c, self.g, self.t)
    }
}

// =============================================================================
// PASS 1 — Streaming bitmask scan
// Memory: O(L) bitmask + O(L) ref_seq + O(N) names
// =============================================================================

#[derive(Debug)]
struct ScanResult {
    ref_seq: Vec<u8>,
    bitmask: Vec<u8>,
    sample_names: Vec<String>,
    seq_length: usize,
}

fn pass1_scan(input: &str, lookup: &[u8; 256]) -> io::Result<ScanResult> {
    let file = File::open(input)?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut fasta_reader = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let mut ref_seq: Vec<u8> = Vec::new();
    let mut bitmask: Vec<u8> = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();
    let mut seq_length = 0usize;
    let mut count = 0usize;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Scanned {pos} sequences")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    loop {
        fasta_reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }

        let seq = record.seq();

        if count == 0 {
            seq_length = seq.len();

            // S1: Sanity check — reject absurdly large sequences
            if seq_length > MAX_SEQ_LENGTH {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Sequence length {} exceeds maximum allowed ({} = 10 Gbp). \
                         Is this a valid alignment?",
                        seq_length, MAX_SEQ_LENGTH
                    ),
                ));
            }
            if seq_length == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "First sequence has length 0.",
                ));
            }

            bitmask = vec![0u8; seq_length];
            ref_seq = seq.to_vec();
        } else if seq.len() != seq_length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Sequence '{}' (#{}) has length {} but expected {} (same as first sequence '{}').",
                    record.id(),
                    count + 1,
                    seq.len(),
                    seq_length,
                    sample_names[0]
                ),
            ));
        }

        sample_names.push(record.id().to_string());

        // Hot loop: update bitmask with lookup table — no branches
        let bm = &mut bitmask;
        for (pos, &nuc) in seq.iter().enumerate() {
            bm[pos] |= lookup[nuc as usize];
        }

        count += 1;
        if count.is_multiple_of(1000) {
            pb.set_position(count as u64);
        }
    }

    pb.finish_and_clear();

    if count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Input FASTA file is empty.",
        ));
    }

    eprintln!(
        "[snpick] Pass 1: Scanned {} sequences of length {}.",
        count, seq_length
    );

    Ok(ScanResult {
        ref_seq,
        bitmask,
        sample_names,
        seq_length,
    })
}

// =============================================================================
// Analyze bitmask → variable positions + constant site counts
// =============================================================================

fn analyze_bitmask(
    bitmask: &[u8],
    ref_seq: &[u8],
    lookup: &[u8; 256],
    include_gaps: bool,
) -> (Vec<VariablePosition>, ConstantSiteCounts) {
    let mut vars = Vec::new();
    let mut cs = ConstantSiteCounts { a: 0, c: 0, g: 0, t: 0 };

    for (pos, &bits) in bitmask.iter().enumerate() {
        let ones = bits.count_ones();

        if ones > 1 {
            let rb_raw = ref_seq[pos].to_ascii_uppercase();
            let rb_has_bit = lookup[rb_raw as usize] != 0;

            // If ref is a standard base, use it; otherwise pick first from bitmask
            let ref_base = if rb_has_bit {
                rb_raw
            } else {
                bits_to_bases(bits, include_gaps)[0]
            };

            let alt_bases: Vec<u8> = bits_to_bases(bits, include_gaps)
                .into_iter()
                .filter(|&b| b != ref_base)
                .collect();

            vars.push(VariablePosition {
                index: pos,
                ref_base,
                alt_bases,
            });
        } else if ones == 1 {
            if bits & BIT_A != 0 { cs.a += 1; }
            else if bits & BIT_C != 0 { cs.c += 1; }
            else if bits & BIT_G != 0 { cs.g += 1; }
            else if bits & BIT_T != 0 { cs.t += 1; }
        }
    }

    (vars, cs)
}

// =============================================================================
// PASS 2 — Stream FASTA → write variable-sites FASTA + collect VCF genotypes
// Memory: O(V) positions + O(V × N) flat array if VCF
// =============================================================================

fn pass2_extract(
    input: &str,
    output: &str,
    pos_indices: &[usize],
    num_var: usize,
    total_sequences: usize,
    seq_length: usize,
    collect_vcf: bool,
) -> io::Result<Option<Vec<u8>>> {
    let file = File::open(input)?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut fasta_reader = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let out_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, out_file);

    // P1: Flat Vec<u8> of V×N instead of Vec<Vec<u8>> — better cache locality
    let mut vcf_geno: Vec<u8> = if collect_vcf {
        vec![0u8; num_var * total_sequences]
    } else {
        Vec::new()
    };

    let pb = ProgressBar::new(total_sequences as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} sequences")
            .unwrap()
            .progress_chars("█▉▊▋▌▍▎▏ "),
    );

    let mut sample_idx = 0usize;
    let mut var_seq = vec![0u8; num_var]; // Reuse buffer

    loop {
        fasta_reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }

        let seq = record.seq();

        // R1: Validate sequence length in pass 2
        if seq.len() != seq_length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Pass 2: Sequence '{}' has length {} but expected {}. \
                     Was the input file modified between passes?",
                    record.id(),
                    seq.len(),
                    seq_length
                ),
            ));
        }

        // Extract variable-site bases into reusable buffer
        for (vi, &p) in pos_indices.iter().enumerate() {
            var_seq[vi] = seq[p];
        }

        // Write FASTA record (raw bytes)
        writer.write_all(b">")?;
        writer.write_all(record.id().as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&var_seq)?;
        writer.write_all(b"\n")?;

        // Collect VCF genotypes into flat array
        if collect_vcf {
            let base = sample_idx;
            for (vi, &nuc) in var_seq.iter().enumerate() {
                vcf_geno[vi * total_sequences + base] = nuc.to_ascii_uppercase();
            }
        }

        sample_idx += 1;
        if sample_idx.is_multiple_of(1000) {
            pb.set_position(sample_idx as u64);
        }
    }

    writer.flush()?;
    pb.finish_and_clear();

    // S4: Verify we wrote the expected number of sequences
    if sample_idx != total_sequences {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Pass 2: Expected {} sequences but read {}. \
                 Was the input file modified between passes?",
                total_sequences, sample_idx
            ),
        ));
    }

    eprintln!(
        "[snpick] Pass 2: Wrote {} sequences to {}.",
        sample_idx, output
    );

    if collect_vcf {
        Ok(Some(vcf_geno))
    } else {
        Ok(None)
    }
}

// =============================================================================
// Write VCF 4.2 with standard GT genotypes
// =============================================================================

fn write_vcf(
    vcf_geno: &[u8],
    num_samples: usize,
    variable_positions: &[VariablePosition],
    vcf_path: &str,
    sample_names: &[String],
    seq_length: usize,
) -> io::Result<()> {
    let out = File::create(vcf_path)?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, out);

    // Header
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference=first_sequence")?;
    writeln!(w, "##contig=<ID=1,length={}>", seq_length)?;
    writeln!(
        w,
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;

    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for name in sample_names {
        write!(w, "\t{}", name)?;
    }
    writeln!(w)?;

    let ns = sample_names.len();

    // P2: Fixed-size allele index array instead of HashMap per position
    let mut allele_lut = [255u8; 256]; // 255 = missing

    for (vi, vp) in variable_positions.iter().enumerate() {
        let alt: String = vp
            .alt_bases
            .iter()
            .map(|&b| if b == b'-' { "*".to_string() } else { (b as char).to_string() })
            .collect::<Vec<_>>()
            .join(",");

        // Build fast allele lookup (reset only used entries)
        allele_lut[vp.ref_base as usize] = 0;
        for (i, &ab) in vp.alt_bases.iter().enumerate() {
            allele_lut[ab as usize] = (i + 1) as u8;
        }

        write!(
            w,
            "1\t{}\t.\t{}\t{}\t.\tPASS\tNS={}\tGT",
            vp.index + 1,
            vp.ref_base as char,
            alt,
            ns
        )?;

        // Read genotypes from flat array
        let row_start = vi * num_samples;
        for si in 0..num_samples {
            let gt_base = vcf_geno[row_start + si];
            let idx = allele_lut[gt_base as usize];
            if idx == 255 {
                write!(w, "\t.")?;
            } else {
                write!(w, "\t{}", idx)?;
            }
        }
        writeln!(w)?;

        // Reset LUT entries
        allele_lut[vp.ref_base as usize] = 255;
        for &ab in &vp.alt_bases {
            allele_lut[ab as usize] = 255;
        }
    }

    w.flush()?;
    Ok(())
}

// =============================================================================
// main
// =============================================================================

fn main() -> io::Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    let lookup = build_lookup(args.include_gaps);

    // S3: Verify input ≠ output before doing anything
    check_paths_differ(&args.fasta, &args.output)?;
    if let Some(ref vcf_out) = args.vcf_output {
        check_paths_differ(&args.fasta, vcf_out)?;
    }

    // ── Pass 1: stream FASTA → bitmask (O(L) memory) ────────────────────
    let scan = pass1_scan(&args.fasta, &lookup)?;
    let seq_length = scan.seq_length;
    let num_samples = scan.sample_names.len();

    // ── Analyze bitmask → variable positions + constant site counts ──────
    let (variable_positions, constant) =
        analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lookup, args.include_gaps);

    let num_var = variable_positions.len();

    eprintln!(
        "[snpick] {} variable sites, {} constant sites ({}).",
        num_var,
        constant.total(),
        constant
    );
    eprintln!("[snpick] ASC fconst: {}", constant.fconst());

    if num_var == 0 {
        eprintln!("[snpick] No variable positions found.");
        return Ok(());
    }

    // Precompute position indices (bitmask + ref_seq drop naturally here)
    let pos_indices: Vec<usize> = variable_positions.iter().map(|v| v.index).collect();

    // ── Pass 2: stream FASTA → write output FASTA (+ VCF genotypes) ─────
    let vcf_geno = pass2_extract(
        &args.fasta,
        &args.output,
        &pos_indices,
        num_var,
        num_samples,
        seq_length,
        args.vcf,
    )?;

    // ── Write VCF if requested ───────────────────────────────────────────
    if let Some(ref geno) = vcf_geno {
        let vcf_path = args.vcf_output.unwrap_or_else(|| "output.vcf".to_string());
        write_vcf(geno, num_samples, &variable_positions, &vcf_path, &scan.sample_names, seq_length)?;
        eprintln!("[snpick] VCF written to {}.", vcf_path);
    }

    let elapsed = start.elapsed();
    eprintln!(
        "[snpick] Done in {}. {} variable sites from {} sequences × {} positions.",
        format_time(elapsed.as_secs_f64()),
        num_var,
        num_samples,
        seq_length
    );

    Ok(())
}

// =============================================================================
// Tests
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn write_temp(name: &str, content: &str) -> String {
        let path = format!("/tmp/snpick_test_{}.fa", name);
        std::fs::write(&path, content).unwrap();
        path
    }

    // --- lookup table ---

    #[test]
    fn test_lookup_standard_bases() {
        let lk = build_lookup(false);
        assert_eq!(lk[b'A' as usize], BIT_A);
        assert_eq!(lk[b'a' as usize], BIT_A);
        assert_eq!(lk[b'C' as usize], BIT_C);
        assert_eq!(lk[b'G' as usize], BIT_G);
        assert_eq!(lk[b'T' as usize], BIT_T);
        assert_eq!(lk[b't' as usize], BIT_T);
    }

    #[test]
    fn test_lookup_gaps() {
        let no_gap = build_lookup(false);
        assert_eq!(no_gap[b'-' as usize], 0);
        let with_gap = build_lookup(true);
        assert_eq!(with_gap[b'-' as usize], BIT_GAP);
    }

    #[test]
    fn test_lookup_ambiguous_is_zero() {
        let lk = build_lookup(false);
        assert_eq!(lk[b'N' as usize], 0);
        assert_eq!(lk[b'R' as usize], 0);
        assert_eq!(lk[b'Y' as usize], 0);
        assert_eq!(lk[b'?' as usize], 0);
    }

    // --- bits_to_bases ---

    #[test]
    fn test_bits_to_bases() {
        assert_eq!(bits_to_bases(BIT_A | BIT_T, false), vec![b'A', b'T']);
        assert_eq!(
            bits_to_bases(BIT_A | BIT_C | BIT_G | BIT_T, false),
            vec![b'A', b'C', b'G', b'T']
        );
        assert_eq!(bits_to_bases(BIT_A | BIT_GAP, true), vec![b'A', b'-']);
        assert_eq!(bits_to_bases(BIT_A | BIT_GAP, false), vec![b'A']);
    }

    // --- check_paths_differ ---

    #[test]
    fn test_paths_differ_ok() {
        let a = write_temp("path_a", ">x\nA\n");
        assert!(check_paths_differ(&a, "/tmp/snpick_test_path_b.fa").is_ok());
        std::fs::remove_file(&a).ok();
    }

    #[test]
    fn test_paths_same_fails() {
        let a = write_temp("path_same", ">x\nA\n");
        assert!(check_paths_differ(&a, &a).is_err());
        std::fs::remove_file(&a).ok();
    }

    // --- pass1_scan ---

    #[test]
    fn test_scan_basic() {
        let path = write_temp("scan_basic2", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&path, &lk).unwrap();
        assert_eq!(r.sample_names, vec!["s1", "s2"]);
        assert_eq!(r.seq_length, 4);
        assert_eq!(r.ref_seq, b"ATGC");
        assert_eq!(r.bitmask[2].count_ones(), 2); // G + C
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_scan_empty() {
        let path = write_temp("scan_empty2", "");
        let lk = build_lookup(false);
        assert!(pass1_scan(&path, &lk).is_err());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_scan_unequal_lengths() {
        let path = write_temp("scan_uneq2", ">s1\nATGC\n>s2\nAT\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&path, &lk);
        assert!(r.is_err());
        assert!(r.unwrap_err().to_string().contains("length"));
        std::fs::remove_file(&path).ok();
    }

    // --- analyze_bitmask ---

    #[test]
    fn test_analyze_no_variants() {
        let bitmask = vec![BIT_A, BIT_T, BIT_G, BIT_C];
        let ref_seq = b"ATGC";
        let lk = build_lookup(false);
        let (vars, cs) = analyze_bitmask(&bitmask, ref_seq, &lk, false);
        assert!(vars.is_empty());
        assert_eq!(cs.total(), 4);
    }

    #[test]
    fn test_analyze_single_variant() {
        let bitmask = vec![BIT_A, BIT_T | BIT_C, BIT_G, BIT_C];
        let ref_seq = b"ATGC";
        let lk = build_lookup(false);
        let (vars, cs) = analyze_bitmask(&bitmask, ref_seq, &lk, false);
        assert_eq!(vars.len(), 1);
        assert_eq!(vars[0].index, 1);
        assert_eq!(vars[0].ref_base, b'T');
        assert_eq!(vars[0].alt_bases, vec![b'C']);
        assert_eq!(cs.total(), 3);
    }

    #[test]
    fn test_analyze_multiallelic() {
        let bitmask = vec![BIT_A | BIT_C | BIT_G | BIT_T];
        let ref_seq = b"A";
        let lk = build_lookup(false);
        let (vars, _) = analyze_bitmask(&bitmask, ref_seq, &lk, false);
        assert_eq!(vars[0].alt_bases, vec![b'C', b'G', b'T']);
    }

    #[test]
    fn test_analyze_ambiguous_ref() {
        let bitmask = vec![BIT_A | BIT_C];
        let ref_seq = b"N";
        let lk = build_lookup(false);
        let (vars, _) = analyze_bitmask(&bitmask, ref_seq, &lk, false);
        assert_eq!(vars[0].ref_base, b'A');
    }

    #[test]
    fn test_constant_site_counts_fconst() {
        let cs = ConstantSiteCounts { a: 100, c: 200, g: 150, t: 50 };
        assert_eq!(cs.fconst(), "100,200,150,50");
        assert_eq!(cs.total(), 500);
    }

    // --- end-to-end ---

    #[test]
    fn test_e2e_fasta_output() {
        let input = write_temp("e2e_fa2", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let output = "/tmp/snpick_test_e2e_fa2_out.fa";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();

        pass2_extract(&input, output, &pos_idx, vars.len(), scan.sample_names.len(), scan.seq_length, false).unwrap();

        let content = std::fs::read_to_string(output).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], ">ref");
        assert_eq!(lines[1], "TC");
        assert_eq!(lines[2], ">s1");
        assert_eq!(lines[3], "TT");
        assert_eq!(lines[4], ">s2");
        assert_eq!(lines[5], "CC");

        std::fs::remove_file(&input).ok();
        std::fs::remove_file(output).ok();
    }

    #[test]
    fn test_e2e_vcf_output() {
        let input = write_temp("e2e_vcf2_in", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fasta_out = "/tmp/snpick_test_e2e_vcf2_out.fa";
        let vcf_out = "/tmp/snpick_test_e2e_vcf2_out.vcf";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();
        let ns = scan.sample_names.len();

        let geno = pass2_extract(&input, fasta_out, &pos_idx, vars.len(), ns, scan.seq_length, true)
            .unwrap().unwrap();

        write_vcf(&geno, ns, &vars, vcf_out, &scan.sample_names, scan.seq_length).unwrap();

        let content = std::fs::read_to_string(vcf_out).unwrap();
        let data: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 1);

        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "C,G,T");
        assert_eq!(fields[9], "0");
        assert_eq!(fields[10], "1");
        assert_eq!(fields[11], "2");
        assert_eq!(fields[12], "3");

        std::fs::remove_file(&input).ok();
        std::fs::remove_file(fasta_out).ok();
        std::fs::remove_file(vcf_out).ok();
    }

    #[test]
    fn test_e2e_vcf_ambiguous_gt() {
        let input = write_temp("e2e_ambig2", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fasta_out = "/tmp/snpick_test_ambig2.fa";
        let vcf_out = "/tmp/snpick_test_ambig2.vcf";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();
        let ns = scan.sample_names.len();

        let geno = pass2_extract(&input, fasta_out, &pos_idx, vars.len(), ns, scan.seq_length, true)
            .unwrap().unwrap();
        write_vcf(&geno, ns, &vars, vcf_out, &scan.sample_names, scan.seq_length).unwrap();

        let content = std::fs::read_to_string(vcf_out).unwrap();
        let data: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(fields[11], "."); // N → missing

        std::fs::remove_file(&input).ok();
        std::fs::remove_file(fasta_out).ok();
        std::fs::remove_file(vcf_out).ok();
    }

    #[test]
    fn test_e2e_with_gaps() {
        let input = write_temp("e2e_gaps2", ">ref\nATGC\n>s1\nA-GC\n");
        let lk_no = build_lookup(false);
        let lk_yes = build_lookup(true);

        let scan = pass1_scan(&input, &lk_no).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk_no, false);
        assert!(vars.is_empty());

        let scan2 = pass1_scan(&input, &lk_yes).unwrap();
        let (vars2, _) = analyze_bitmask(&scan2.bitmask, &scan2.ref_seq, &lk_yes, true);
        assert_eq!(vars2.len(), 1);
        assert_eq!(vars2[0].index, 1);

        std::fs::remove_file(&input).ok();
    }

    // --- integrity check ---

    #[test]
    fn test_pass2_detects_modified_file() {
        let input = write_temp("integrity", ">s1\nATGC\n>s2\nATCC\n");
        let output = "/tmp/snpick_test_integrity_out.fa";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();

        // Overwrite input with fewer sequences
        std::fs::write(&input, ">s1\nATGC\n").unwrap();

        let result = pass2_extract(&input, output, &pos_idx, vars.len(), 2, scan.seq_length, false);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Expected 2"));

        std::fs::remove_file(&input).ok();
        std::fs::remove_file(output).ok();
    }
}
