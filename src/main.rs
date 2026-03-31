// Necessary imports
use bio::io::fasta;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;
use chrono::Local;
use sysinfo::{Pid, System};
use std::collections::{HashSet, HashMap};

/// snpick: A tool to extract variable sites from a FASTA alignment and generate a VCF with actual bases, including ambiguous bases and codons.
#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = env!("CARGO_PKG_VERSION"),
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>",
    about = "A fast and efficient tool for extracting variable sites and generating a VCF with actual bases, including ambiguous bases and codons."
)]
struct Args {
    /// Input FASTA alignment file
    #[arg(short, long, help = "Input FASTA alignment file")]
    fasta: String,

    /// Output FASTA file with variable sites
    #[arg(short, long, help = "Output FASTA file with variable sites")]
    output: String,

    /// Number of threads to use (optional)
    #[arg(short, long, default_value_t = 4, help = "Number of threads to use (optional)")]
    threads: usize,

    /// Consider the '-' (gap) symbol in variable site detection
    #[arg(short = 'g', long, help = "Consider the '-' (gap) symbol in variable site detection")]
    include_gaps: bool,

    /// Generate VCF file with variable sites
    #[arg(long, help = "Generate VCF file with variable sites")]
    vcf: bool,

    /// Output VCF file (optional)
    #[arg(long, help = "Output VCF file (optional)")]
    vcf_output: Option<String>,
}

/// Converts a nucleotide to a bitmask for variable site detection.
/// Returns None for ambiguous/unknown bases (they are excluded from variant calling).
fn nucleotide_to_bit(nuc: u8, include_gaps: bool) -> Option<u8> {
    match nuc.to_ascii_uppercase() {
        b'A' => Some(0b000001),
        b'C' => Some(0b000010),
        b'G' => Some(0b000100),
        b'T' => Some(0b001000),
        b'-' if include_gaps => Some(0b010000),
        _ => None,
    }
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let input_filename = args.fasta;
    let output_filename = args.output;
    let num_threads = args.threads;
    let include_gaps = args.include_gaps;
    let generate_vcf = args.vcf;
    let vcf_output_filename = args.vcf_output.unwrap_or_else(|| "output.vcf".to_string());

    // Configure a local thread pool for Rayon
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build Rayon thread pool");

    // Initialize system for process-level RAM usage reporting (#15)
    let system = Arc::new(Mutex::new(System::new_all()));

    pool.install(|| {
        // Step 1: Read sequences and identify variable positions
        println!("Starting Step 1: Reading sequences and identifying variable positions...");
        let (sequences, sample_names, seq_length) =
            read_sequences(&input_filename)
                .expect("Failed to read input FASTA");

        let variable_positions =
            identify_variable_positions(&sequences, seq_length, include_gaps);
        println!(
            "Step 1 Completed: Found {} variable positions.",
            variable_positions.len()
        );

        if variable_positions.is_empty() {
            eprintln!("No variable positions found in the alignment.");
            std::process::exit(0);
        }

        // Step 2: Extract variable positions and write to output file
        println!("Starting Step 2: Extracting variable positions and writing to output...");
        extract_and_write_variables(
            &sequences,
            &sample_names,
            &output_filename,
            &variable_positions,
            &system,
        )
        .expect("Failed to extract and write variable positions");
        println!("Variable positions alignment written to {}", output_filename);

        // Step 3: Generate VCF if requested
        if generate_vcf {
            println!("Generating VCF file...");
            generate_vcf_file(
                &sequences,
                &variable_positions,
                &vcf_output_filename,
                &sample_names,
                seq_length,
            )
            .expect("Failed to generate VCF file");
            println!("VCF file generated: {}", vcf_output_filename);
        }
    });

    Ok(())
}

/// Information about a variable position (no genotype storage — extracted on demand from sequences)
struct VariablePositionInfo {
    position: usize,
    reference_base: u8,
    alternate_bases: Vec<u8>,
}

/// Read all sequences from a FASTA file into memory
fn read_sequences(input_filename: &str) -> io::Result<(Vec<Vec<u8>>, Vec<String>, usize)> {
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file);
    let fasta_reader = fasta::Reader::new(reader);

    let mut sequences = Vec::new();
    let mut sample_names = Vec::new();

    for result in fasta_reader.records() {
        let record = result?;
        sample_names.push(record.id().to_string());
        sequences.push(record.seq().to_owned());
    }

    if sequences.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "The input FASTA file is empty."));
    }

    let seq_length = sequences[0].len();
    for (i, seq) in sequences.iter().enumerate() {
        if seq.len() != seq_length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Sequence '{}' has length {} but expected {} (same as first sequence).",
                    sample_names[i], seq.len(), seq_length
                ),
            ));
        }
    }

    println!(
        "[{}] Read {} sequences of length {}.",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        sequences.len(),
        seq_length
    );

    Ok((sequences, sample_names, seq_length))
}

/// Identify variable positions using the first sequence as reference.
/// Genotypes are NOT stored here — they are extracted on demand from the
/// original sequences, saving N × M bytes of redundant memory.
fn identify_variable_positions(
    sequences: &[Vec<u8>],
    seq_length: usize,
    include_gaps: bool,
) -> Vec<VariablePositionInfo> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing positions: {pos}/{len}")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    let pos_counter = AtomicUsize::new(0);
    let ref_seq = &sequences[0];

    let mut variable_positions: Vec<VariablePositionInfo> = (0..seq_length)
        .into_par_iter()
        .filter_map(|pos| {
            let mut seen = HashSet::new();

            for seq in sequences {
                let nuc = seq[pos];
                if nucleotide_to_bit(nuc, include_gaps).is_some() {
                    seen.insert(nuc.to_ascii_uppercase());
                }
            }

            // Update progress every 10K positions to reduce contention
            let current = pos_counter.fetch_add(1, Ordering::Relaxed) + 1;
            if current % 10_000 == 0 || current == seq_length {
                pb.set_position(current as u64);
                pb.set_length(seq_length as u64);
            }

            if seen.len() > 1 {
                let reference_base = ref_seq[pos].to_ascii_uppercase();

                let mut alternate_bases: Vec<u8> = seen
                    .into_iter()
                    .filter(|&nuc| nuc != reference_base)
                    .collect();
                alternate_bases.sort();

                Some(VariablePositionInfo {
                    position: pos,
                    reference_base,
                    alternate_bases,
                })
            } else {
                None
            }
        })
        .collect();

    // Sort by position for deterministic output
    variable_positions.sort_by_key(|info| info.position);

    pb.finish_with_message("Position processing completed.");

    println!(
        "[{}] Variable position identification completed.",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    println!("Total sequences processed: {}", sequences.len());

    variable_positions
}

/// Report process-level RAM usage (not system-wide) (#15)
fn report_process_memory(system: &Arc<Mutex<System>>, current: usize) {
    if let Ok(mut sys) = system.lock() {
        let pid = Pid::from_u32(std::process::id());
        sys.refresh_processes(sysinfo::ProcessesToUpdate::Some(&[pid]), true);
        if let Some(process) = sys.process(pid) {
            println!(
                "[{}] Written {} sequences. Process RSS: {} MB.",
                Local::now().format("%Y-%m-%d %H:%M:%S"),
                current,
                process.memory() / (1024 * 1024)
            );
        }
    }
}

/// Extract variable positions and write to the output FASTA file (deterministic order)
fn extract_and_write_variables(
    sequences: &[Vec<u8>],
    sample_names: &[String],
    output_filename: &str,
    variable_positions: &[VariablePositionInfo],
    system: &Arc<Mutex<System>>,
) -> io::Result<()> {
    let output_file = File::create(output_filename)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file);

    let positions: Vec<usize> = variable_positions.iter().map(|info| info.position).collect();
    let total_sequences = sequences.len();

    let pb_write = ProgressBar::new_spinner();
    pb_write.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Writing sequences: {pos}/{len}")
            .unwrap(),
    );
    pb_write.enable_steady_tick(Duration::from_millis(100));
    pb_write.set_length(total_sequences as u64);

    for (i, (seq, name)) in sequences.iter().zip(sample_names.iter()).enumerate() {
        let var_seq: Vec<u8> = positions.iter().map(|&pos| seq[pos]).collect();
        let var_seq_str = String::from_utf8_lossy(&var_seq);

        writeln!(writer, ">{}", name)?;
        writeln!(writer, "{}", var_seq_str)?;

        let current = i + 1;
        if current % 10_000 == 0 || current == total_sequences {
            pb_write.set_position(current as u64);
        }

        if current % 100_000 == 0 {
            report_process_memory(system, current);
        }
    }

    writer.flush()?;

    pb_write.finish_with_message(format!(
        "Completed: {} sequences written.",
        total_sequences
    ));

    println!(
        "[{}] Total sequences written: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        total_sequences
    );

    Ok(())
}

/// Generate a VCF 4.2 file with standard GT genotypes.
/// Genotypes are extracted on demand from the original sequences (no redundant copy).
fn generate_vcf_file(
    sequences: &[Vec<u8>],
    variable_positions: &[VariablePositionInfo],
    vcf_output_filename: &str,
    sample_names: &[String],
    seq_length: usize,
) -> io::Result<()> {
    let output_file = File::create(vcf_output_filename)?;
    let mut writer = BufWriter::new(output_file);

    // VCF header
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=snpick")?;
    writeln!(writer, "##reference=first_sequence")?;
    writeln!(writer, "##contig=<ID=1,length={}>", seq_length)?;
    writeln!(
        writer,
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;

    write!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample_name in sample_names {
        write!(writer, "\t{}", sample_name)?;
    }
    writeln!(writer)?;

    for info in variable_positions {
        let chrom = "1";
        let pos = info.position + 1; // VCF is 1-based
        let ref_base = info.reference_base as char;

        // ALT alleles separated by commas; gaps represented as * (VCF 4.2)
        let alt_strings: Vec<String> = info
            .alternate_bases
            .iter()
            .map(|&b| {
                let c = b as char;
                if c == '-' { "*".to_string() } else { c.to_string() }
            })
            .collect();
        let alt = alt_strings.join(",");

        let info_field = format!("NS={}", sample_names.len());

        // Build allele index map: REF=0, ALT1=1, ALT2=2, ...
        let mut allele_index: HashMap<u8, usize> = HashMap::new();
        allele_index.insert(info.reference_base, 0);
        for (i, &alt_base) in info.alternate_bases.iter().enumerate() {
            allele_index.insert(alt_base, i + 1);
        }

        write!(
            writer,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT",
            chrom, pos, ref_base, alt, info_field
        )?;

        // Extract genotypes on demand from original sequences
        for seq in sequences {
            let genotype_base = seq[info.position].to_ascii_uppercase();
            let gt = match allele_index.get(&genotype_base) {
                Some(idx) => idx.to_string(),
                None => ".".to_string(), // Ambiguous or unknown bases → missing (#13)
            };
            write!(writer, "\t{}", gt)?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;

    Ok(())
}

// =============================================================================
// Tests (#12)
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Helper: write a temp FASTA and return its path
    fn write_temp_fasta(name: &str, content: &str) -> String {
        let path = format!("/tmp/snpick_test_{}.fa", name);
        let mut f = File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    // --- nucleotide_to_bit ---

    #[test]
    fn test_nucleotide_to_bit_standard() {
        assert_eq!(nucleotide_to_bit(b'A', false), Some(0b000001));
        assert_eq!(nucleotide_to_bit(b'a', false), Some(0b000001));
        assert_eq!(nucleotide_to_bit(b'C', false), Some(0b000010));
        assert_eq!(nucleotide_to_bit(b'G', false), Some(0b000100));
        assert_eq!(nucleotide_to_bit(b'T', false), Some(0b001000));
        assert_eq!(nucleotide_to_bit(b't', false), Some(0b001000));
    }

    #[test]
    fn test_nucleotide_to_bit_gap() {
        assert_eq!(nucleotide_to_bit(b'-', false), None);
        assert_eq!(nucleotide_to_bit(b'-', true), Some(0b010000));
    }

    #[test]
    fn test_nucleotide_to_bit_ambiguous() {
        assert_eq!(nucleotide_to_bit(b'N', false), None);
        assert_eq!(nucleotide_to_bit(b'R', false), None);
        assert_eq!(nucleotide_to_bit(b'Y', false), None);
        assert_eq!(nucleotide_to_bit(b'?', false), None);
    }

    // --- read_sequences ---

    #[test]
    fn test_read_sequences_basic() {
        let path = write_temp_fasta("basic", ">s1\nATGC\n>s2\nATCC\n");
        let (seqs, names, len) = read_sequences(&path).unwrap();
        assert_eq!(names, vec!["s1", "s2"]);
        assert_eq!(seqs.len(), 2);
        assert_eq!(len, 4);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_read_sequences_empty() {
        let path = write_temp_fasta("empty", "");
        let result = read_sequences(&path);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("empty"));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_read_sequences_unequal_lengths() {
        let path = write_temp_fasta("unequal", ">s1\nATGC\n>s2\nAT\n");
        let result = read_sequences(&path);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("length"));
        std::fs::remove_file(&path).ok();
    }

    // --- identify_variable_positions ---

    #[test]
    fn test_identify_no_variants() {
        let seqs = vec![b"AAAA".to_vec(), b"AAAA".to_vec()];
        let result = identify_variable_positions(&seqs, 4, false);
        assert!(result.is_empty());
    }

    #[test]
    fn test_identify_single_variant() {
        let seqs = vec![b"ATGC".to_vec(), b"ACGC".to_vec()];
        let result = identify_variable_positions(&seqs, 4, false);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].position, 1);
        assert_eq!(result[0].reference_base, b'T');
        assert_eq!(result[0].alternate_bases, vec![b'C']);
    }

    #[test]
    fn test_identify_multiple_variants() {
        let seqs = vec![
            b"ATGC".to_vec(),
            b"CTGC".to_vec(),
            b"ATGA".to_vec(),
        ];
        let result = identify_variable_positions(&seqs, 4, false);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].position, 0); // A vs C
        assert_eq!(result[1].position, 3); // C vs A
    }

    #[test]
    fn test_identify_with_gaps() {
        let seqs = vec![b"ATGC".to_vec(), b"A-GC".to_vec()];
        // Without gaps: no variant at pos 1 (- is ignored)
        let result_no_gaps = identify_variable_positions(&seqs, 4, false);
        assert!(result_no_gaps.is_empty());
        // With gaps: variant at pos 1
        let result_gaps = identify_variable_positions(&seqs, 4, true);
        assert_eq!(result_gaps.len(), 1);
        assert_eq!(result_gaps[0].position, 1);
    }

    #[test]
    fn test_identify_ambiguous_ignored() {
        // N should be ignored — if only A and N at a position, it's NOT variable
        let seqs = vec![b"ATGC".to_vec(), b"ANGC".to_vec()];
        let result = identify_variable_positions(&seqs, 4, false);
        assert!(result.is_empty());
    }

    #[test]
    fn test_identify_multiallelic() {
        let seqs = vec![
            b"ATGC".to_vec(),
            b"CTGC".to_vec(),
            b"GTGC".to_vec(),
            b"TTGC".to_vec(),
        ];
        let result = identify_variable_positions(&seqs, 4, false);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].position, 0);
        assert_eq!(result[0].reference_base, b'A');
        assert_eq!(result[0].alternate_bases, vec![b'C', b'G', b'T']);
    }

    #[test]
    fn test_positions_sorted() {
        // Even with par_iter, positions must come out sorted
        let mut seqs = vec![vec![b'A'; 1000], vec![b'A'; 1000]];
        // Make positions 999, 500, 100, 0 variable
        for &pos in &[0usize, 100, 500, 999] {
            seqs[1][pos] = b'T';
        }
        let result = identify_variable_positions(&seqs, 1000, false);
        let positions: Vec<usize> = result.iter().map(|v| v.position).collect();
        assert_eq!(positions, vec![0, 100, 500, 999]);
    }

    // --- VCF format ---

    #[test]
    fn test_vcf_output_format() {
        let seqs = vec![
            b"ATGC".to_vec(),
            b"CTGC".to_vec(),
            b"GTGC".to_vec(),
        ];
        let sample_names = vec!["ref".to_string(), "s1".to_string(), "s2".to_string()];
        let var_pos = identify_variable_positions(&seqs, 4, false);

        let vcf_path = "/tmp/snpick_test_vcf.vcf";
        generate_vcf_file(&seqs, &var_pos, vcf_path, &sample_names, 4).unwrap();

        let content = std::fs::read_to_string(vcf_path).unwrap();

        // Check header
        assert!(content.starts_with("##fileformat=VCFv4.2"));
        assert!(content.contains("##FORMAT=<ID=GT,"));
        assert!(content.contains("##contig=<ID=1,length=4>"));

        // Check data line
        let data_lines: Vec<&str> = content.lines()
            .filter(|l| !l.starts_with('#'))
            .collect();
        assert_eq!(data_lines.len(), 1);

        let fields: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(fields[0], "1");        // CHROM
        assert_eq!(fields[1], "1");        // POS (1-based)
        assert_eq!(fields[3], "A");        // REF (first sequence)
        assert_eq!(fields[4], "C,G");      // ALT (comma-separated, sorted)
        assert_eq!(fields[8], "GT");       // FORMAT
        assert_eq!(fields[9], "0");        // ref → 0
        assert_eq!(fields[10], "1");       // s1 → C = 1
        assert_eq!(fields[11], "2");       // s2 → G = 2

        std::fs::remove_file(vcf_path).ok();
    }

    #[test]
    fn test_vcf_ambiguous_genotype() {
        let seqs = vec![
            b"ATGC".to_vec(),
            b"CTGC".to_vec(),
            b"NTGC".to_vec(), // N at variable position
        ];
        let sample_names = vec!["ref".to_string(), "s1".to_string(), "s2".to_string()];
        let var_pos = identify_variable_positions(&seqs, 4, false);

        let vcf_path = "/tmp/snpick_test_vcf_ambig.vcf";
        generate_vcf_file(&seqs, &var_pos, vcf_path, &sample_names, 4).unwrap();

        let content = std::fs::read_to_string(vcf_path).unwrap();
        let data_lines: Vec<&str> = content.lines()
            .filter(|l| !l.starts_with('#'))
            .collect();
        let fields: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(fields[11], "."); // N → missing

        std::fs::remove_file(vcf_path).ok();
    }

    // --- FASTA output ---

    #[test]
    fn test_fasta_output_deterministic() {
        let seqs = vec![
            b"ATGCATGC".to_vec(),
            b"ATGTATGC".to_vec(),
            b"ACGCATGC".to_vec(),
        ];
        let names = vec!["ref".to_string(), "s1".to_string(), "s2".to_string()];
        let var_pos = identify_variable_positions(&seqs, 8, false);
        let system = Arc::new(Mutex::new(System::new_all()));

        let fasta_path = "/tmp/snpick_test_fasta.fa";
        extract_and_write_variables(&seqs, &names, fasta_path, &var_pos, &system).unwrap();

        let content = std::fs::read_to_string(fasta_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        // Order must match input order
        assert_eq!(lines[0], ">ref");
        assert_eq!(lines[1], "TC");
        assert_eq!(lines[2], ">s1");
        assert_eq!(lines[3], "TT");
        assert_eq!(lines[4], ">s2");
        assert_eq!(lines[5], "CC");

        std::fs::remove_file(fasta_path).ok();
    }
}
