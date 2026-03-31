use bio::io::fasta::{self, FastaRead};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::time::{Duration, Instant};

use chrono::Local;

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

    /// Number of threads (reserved for future use)
    #[arg(short, long, default_value_t = 1, hide = true)]
    threads: usize,

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
// Lookup table: nucleotide byte → bitmask (avoids branching in hot loop)
// =============================================================================

const BIT_A: u8 = 0b00001;
const BIT_C: u8 = 0b00010;
const BIT_G: u8 = 0b00100;
const BIT_T: u8 = 0b01000;
const BIT_GAP: u8 = 0b10000;

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
// Memory: O(L) for bitmask + O(L) for reference sequence + O(N) for names
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
            bitmask = vec![0u8; seq_length];
            ref_seq = seq.to_vec();
        } else if seq.len() != seq_length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Sequence '{}' has length {} but expected {}.",
                    record.id(),
                    seq.len(),
                    seq_length
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
        if count % 1000 == 0 {
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
        "[{}] Pass 1: Scanned {} sequences of length {}.",
        Local::now().format("%H:%M:%S"),
        count,
        seq_length
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
    let mut cs = ConstantSiteCounts {
        a: 0,
        c: 0,
        g: 0,
        t: 0,
    };

    for (pos, &bits) in bitmask.iter().enumerate() {
        let ones = bits.count_ones();

        if ones > 1 {
            // Variable position
            let rb_raw = ref_seq[pos].to_ascii_uppercase();
            let rb_has_bit = lookup[rb_raw as usize] != 0;

            // If reference is a standard base, use it; otherwise pick first from bitmask
            let ref_base = if rb_has_bit {
                rb_raw
            } else {
                // Ref is ambiguous (N, R, etc.) — pick first standard base
                let all = bits_to_bases(bits, include_gaps);
                all[0]
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
            // Constant site
            if bits & BIT_A != 0 {
                cs.a += 1;
            } else if bits & BIT_C != 0 {
                cs.c += 1;
            } else if bits & BIT_G != 0 {
                cs.g += 1;
            } else if bits & BIT_T != 0 {
                cs.t += 1;
            }
            // gap-only positions are not counted as constant for ASC
        }
        // ones == 0: all ambiguous → skip
    }

    (vars, cs)
}

// =============================================================================
// PASS 2 — Stream FASTA → write variable-sites FASTA + collect VCF genotypes
// Memory: O(V) for positions + O(V × N) if VCF requested
// =============================================================================

fn pass2_extract(
    input: &str,
    output: &str,
    pos_indices: &[usize],
    num_var: usize,
    total_sequences: usize,
    collect_vcf: bool,
) -> io::Result<Option<Vec<Vec<u8>>>> {
    let file = File::open(input)?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut fasta_reader = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let out_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, out_file);

    // VCF genotype matrix: genotypes[var_idx] = Vec of bases per sample
    let mut vcf_geno: Vec<Vec<u8>> = if collect_vcf {
        (0..num_var)
            .map(|_| Vec::with_capacity(total_sequences))
            .collect()
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

    let mut count = 0u64;

    loop {
        fasta_reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }

        let seq = record.seq();

        // Extract variable-site bases
        let var_seq: Vec<u8> = pos_indices.iter().map(|&p| seq[p]).collect();

        // Write FASTA record (raw bytes — no UTF-8 conversion needed)
        writer.write_all(b">")?;
        writer.write_all(record.id().as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&var_seq)?;
        writer.write_all(b"\n")?;

        // Collect VCF genotypes
        if collect_vcf {
            for (vi, &nuc) in var_seq.iter().enumerate() {
                vcf_geno[vi].push(nuc.to_ascii_uppercase());
            }
        }

        count += 1;
        if count % 1000 == 0 {
            pb.set_position(count);
        }
    }

    writer.flush()?;
    pb.finish_and_clear();

    eprintln!(
        "[{}] Pass 2: Wrote {} sequences to {}.",
        Local::now().format("%H:%M:%S"),
        count,
        output
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
    vcf_geno: &[Vec<u8>],
    variable_positions: &[VariablePosition],
    vcf_path: &str,
    sample_names: &[String],
    seq_length: usize,
) -> io::Result<()> {
    let out = File::create(vcf_path)?;
    let mut w = BufWriter::new(out);

    // Header
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick")?;
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

    for (vi, vp) in variable_positions.iter().enumerate() {
        let alt: String = vp
            .alt_bases
            .iter()
            .map(|&b| {
                if b == b'-' {
                    "*".to_string()
                } else {
                    (b as char).to_string()
                }
            })
            .collect::<Vec<_>>()
            .join(",");

        // Build allele index: REF=0, ALT1=1, ALT2=2, ...
        let mut aidx: HashMap<u8, usize> = HashMap::with_capacity(5);
        aidx.insert(vp.ref_base, 0);
        for (i, &ab) in vp.alt_bases.iter().enumerate() {
            aidx.insert(ab, i + 1);
        }

        write!(
            w,
            "1\t{}\t.\t{}\t{}\t.\tPASS\tNS={}\tGT",
            vp.index + 1,
            vp.ref_base as char,
            alt,
            ns
        )?;

        for &gt_base in &vcf_geno[vi] {
            match aidx.get(&gt_base) {
                Some(idx) => write!(w, "\t{}", idx)?,
                None => write!(w, "\t.")?,
            }
        }
        writeln!(w)?;
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

    // ── Pass 1: stream FASTA → bitmask (O(L) memory) ────────────────────
    let scan = pass1_scan(&args.fasta, &lookup)?;
    let seq_length = scan.seq_length;
    let num_samples = scan.sample_names.len();

    // ── Analyze bitmask → variable positions + constant site counts ──────
    let (variable_positions, constant) =
        analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lookup, args.include_gaps);

    // Free bitmask + ref_seq — no longer needed
    drop(scan.bitmask);
    drop(scan.ref_seq);

    let num_var = variable_positions.len();

    eprintln!(
        "[{}] {} variable sites, {} constant sites ({}).",
        Local::now().format("%H:%M:%S"),
        num_var,
        constant.total(),
        constant
    );
    eprintln!(
        "[{}] ASC fconst: {}",
        Local::now().format("%H:%M:%S"),
        constant.fconst()
    );

    if num_var == 0 {
        eprintln!("No variable positions found in the alignment.");
        return Ok(());
    }

    // Precompute position indices
    let pos_indices: Vec<usize> = variable_positions.iter().map(|v| v.index).collect();

    // ── Pass 2: stream FASTA → write output FASTA (+ VCF genotypes) ─────
    let vcf_geno = pass2_extract(
        &args.fasta,
        &args.output,
        &pos_indices,
        num_var,
        num_samples,
        args.vcf,
    )?;

    // ── Write VCF if requested ───────────────────────────────────────────
    if let Some(geno) = &vcf_geno {
        let vcf_path = args
            .vcf_output
            .unwrap_or_else(|| "output.vcf".to_string());
        write_vcf(
            geno,
            &variable_positions,
            &vcf_path,
            &scan.sample_names,
            seq_length,
        )?;
        eprintln!(
            "[{}] VCF written to {}.",
            Local::now().format("%H:%M:%S"),
            vcf_path
        );
    }

    let elapsed = start.elapsed();
    eprintln!(
        "[{}] Done in {:.2}s. {} variable sites from {} sequences × {} positions.",
        Local::now().format("%H:%M:%S"),
        elapsed.as_secs_f64(),
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
        assert_eq!(
            bits_to_bases(BIT_A | BIT_GAP, true),
            vec![b'A', b'-']
        );
        assert_eq!(bits_to_bases(BIT_A | BIT_GAP, false), vec![b'A']);
    }

    // --- pass1_scan ---

    #[test]
    fn test_scan_basic() {
        let path = write_temp("scan_basic", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&path, &lk).unwrap();
        assert_eq!(r.sample_names, vec!["s1", "s2"]);
        assert_eq!(r.seq_length, 4);
        assert_eq!(r.ref_seq, b"ATGC");
        // pos 0: A only, pos 1: T only, pos 2: G|C, pos 3: C only
        assert_eq!(r.bitmask[0].count_ones(), 1); // A
        assert_eq!(r.bitmask[1].count_ones(), 1); // T
        assert_eq!(r.bitmask[2].count_ones(), 2); // G + C
        assert_eq!(r.bitmask[3].count_ones(), 1); // C
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_scan_empty() {
        let path = write_temp("scan_empty", "");
        let lk = build_lookup(false);
        assert!(pass1_scan(&path, &lk).is_err());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_scan_unequal_lengths() {
        let path = write_temp("scan_uneq", ">s1\nATGC\n>s2\nAT\n");
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
        assert_eq!(cs.a, 1);
        assert_eq!(cs.t, 1);
        assert_eq!(cs.g, 1);
        assert_eq!(cs.c, 1);
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
        assert_eq!(vars.len(), 1);
        assert_eq!(vars[0].ref_base, b'A');
        assert_eq!(vars[0].alt_bases, vec![b'C', b'G', b'T']);
    }

    #[test]
    fn test_analyze_ambiguous_ref() {
        // Ref has N at a position where samples have A and C
        let bitmask = vec![BIT_A | BIT_C];
        let ref_seq = b"N";
        let lk = build_lookup(false);
        let (vars, _) = analyze_bitmask(&bitmask, ref_seq, &lk, false);
        assert_eq!(vars.len(), 1);
        // Should pick first standard base (A) as ref
        assert_eq!(vars[0].ref_base, b'A');
        assert_eq!(vars[0].alt_bases, vec![b'C']);
    }

    #[test]
    fn test_constant_site_counts_fconst() {
        let cs = ConstantSiteCounts {
            a: 100,
            c: 200,
            g: 150,
            t: 50,
        };
        assert_eq!(cs.fconst(), "100,200,150,50");
        assert_eq!(cs.total(), 500);
    }

    // --- end-to-end: FASTA output ---

    #[test]
    fn test_e2e_fasta_output() {
        let input = write_temp("e2e_fa", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let output = "/tmp/snpick_test_e2e_fa_out.fa";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();

        pass2_extract(&input, output, &pos_idx, vars.len(), scan.sample_names.len(), false)
            .unwrap();

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

    // --- end-to-end: VCF output ---

    #[test]
    fn test_e2e_vcf_output() {
        let input = write_temp("e2e_vcf_in", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fasta_out = "/tmp/snpick_test_e2e_vcf_out.fa";
        let vcf_out = "/tmp/snpick_test_e2e_vcf_out.vcf";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();

        let geno = pass2_extract(
            &input,
            fasta_out,
            &pos_idx,
            vars.len(),
            scan.sample_names.len(),
            true,
        )
        .unwrap()
        .unwrap();

        write_vcf(&geno, &vars, vcf_out, &scan.sample_names, scan.seq_length).unwrap();

        let content = std::fs::read_to_string(vcf_out).unwrap();
        let data: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 1);

        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(fields[1], "1"); // POS
        assert_eq!(fields[3], "A"); // REF
        assert_eq!(fields[4], "C,G,T"); // ALT
        assert_eq!(fields[8], "GT");
        assert_eq!(fields[9], "0"); // ref → 0
        assert_eq!(fields[10], "1"); // C → 1
        assert_eq!(fields[11], "2"); // G → 2
        assert_eq!(fields[12], "3"); // T → 3

        std::fs::remove_file(&input).ok();
        std::fs::remove_file(fasta_out).ok();
        std::fs::remove_file(vcf_out).ok();
    }

    #[test]
    fn test_e2e_vcf_ambiguous_gt() {
        let input = write_temp("e2e_ambig", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fasta_out = "/tmp/snpick_test_ambig.fa";
        let vcf_out = "/tmp/snpick_test_ambig.vcf";
        let lk = build_lookup(false);

        let scan = pass1_scan(&input, &lk).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk, false);
        let pos_idx: Vec<usize> = vars.iter().map(|v| v.index).collect();

        let geno = pass2_extract(
            &input,
            fasta_out,
            &pos_idx,
            vars.len(),
            scan.sample_names.len(),
            true,
        )
        .unwrap()
        .unwrap();

        write_vcf(&geno, &vars, vcf_out, &scan.sample_names, scan.seq_length).unwrap();

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
        let input = write_temp("e2e_gaps", ">ref\nATGC\n>s1\nA-GC\n");
        let lk_no = build_lookup(false);
        let lk_yes = build_lookup(true);

        // Without gaps
        let scan = pass1_scan(&input, &lk_no).unwrap();
        let (vars, _) = analyze_bitmask(&scan.bitmask, &scan.ref_seq, &lk_no, false);
        assert!(vars.is_empty());

        // With gaps
        let scan2 = pass1_scan(&input, &lk_yes).unwrap();
        let (vars2, _) = analyze_bitmask(&scan2.bitmask, &scan2.ref_seq, &lk_yes, true);
        assert_eq!(vars2.len(), 1);
        assert_eq!(vars2[0].index, 1);

        std::fs::remove_file(&input).ok();
    }
}
