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
use sysinfo::{System};
use std::collections::{HashSet, HashMap};

/// snpick: A tool to extract variable sites from a FASTA alignment and generate a VCF with actual bases, including ambiguous bases and codons.
#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = "1.0.0",
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

/// Converts a nucleotide to a bitmask
fn nucleotide_to_bit(nuc: u8, include_gaps: bool) -> Option<u8> {
    match nuc.to_ascii_uppercase() {
        b'A' => Some(0b000001),
        b'C' => Some(0b000010),
        b'G' => Some(0b000100),
        b'T' => Some(0b001000),
        b'-' if include_gaps => Some(0b010000),
        _ => None, // Exclude ambiguous bases
    }
}

fn main() -> io::Result<()> {
    // Parse command-line arguments
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

    // Initialize system for RAM usage reporting
    let system = Arc::new(Mutex::new(System::new_all()));

    // Execute processing within the thread pool
    pool.install(|| {
        // Step 1: Read sequences and identify variable positions
        println!("Starting Step 1: Reading sequences and identifying variable positions...");
        let (sequences, sample_names, seq_length) =
            read_sequences(&input_filename)
                .expect("Failed to read input FASTA");
        let variable_positions_info =
            identify_variable_positions(&sequences, seq_length, include_gaps);
        println!(
            "Step 1 Completed: Found {} variable positions.",
            variable_positions_info.len()
        );

        if variable_positions_info.is_empty() {
            eprintln!("No variable positions found in the alignment.");
            std::process::exit(0);
        }

        // Step 2: Extract variable positions and write to output file (deterministic order)
        println!("Starting Step 2: Extracting variable positions and writing to output...");
        extract_and_write_variables(
            &sequences,
            &sample_names,
            &output_filename,
            &variable_positions_info,
            &system,
        )
        .expect("Failed to extract and write variable positions");
        println!("Variable positions alignment written to {}", output_filename);

        // Step 3: Generate VCF if requested
        if generate_vcf {
            println!("Generating VCF file...");
            generate_vcf_file(
                &variable_positions_info,
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

/// Structure to store information about variable positions
struct VariablePositionInfo {
    position: usize,
    reference_base: u8,
    alternate_bases: Vec<u8>,
    genotypes: Vec<u8>, // Bases at this position for each sample
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

/// Identify variable positions using the first sequence as reference (FIX #2)
fn identify_variable_positions(
    sequences: &[Vec<u8>],
    seq_length: usize,
    include_gaps: bool,
) -> Vec<VariablePositionInfo> {
    let total_sequences = sequences.len();

    // Initialize a progress spinner
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing positions: {pos}/{len}")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    let pos_counter = AtomicUsize::new(0);

    // Use first sequence as reference (standard VCF convention)
    let ref_seq = &sequences[0];

    let mut variable_positions_info: Vec<VariablePositionInfo> = (0..seq_length)
        .into_par_iter()
        .filter_map(|pos| {
            let mut seen = HashSet::new();
            let mut genotypes = Vec::with_capacity(total_sequences);

            for seq in sequences {
                let nuc = seq[pos];
                genotypes.push(nuc);
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

                // Alternate bases: everything that's not the reference, sorted for determinism
                let mut alternate_bases: Vec<u8> = seen
                    .into_iter()
                    .filter(|&nuc| nuc != reference_base)
                    .collect();
                alternate_bases.sort();

                Some(VariablePositionInfo {
                    position: pos,
                    reference_base,
                    alternate_bases,
                    genotypes,
                })
            } else {
                None
            }
        })
        .collect();

    // Sort by position for deterministic output (FIX #1: par_iter may reorder)
    variable_positions_info.sort_by_key(|info| info.position);

    pb.finish_with_message("Position processing completed.");

    println!(
        "[{}] Variable position identification completed.",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    println!("Total sequences processed: {}", total_sequences);

    variable_positions_info
}

/// Step 2: Extract variable positions and write to the output FASTA file
/// FIX #1: Write sequentially from in-memory sequences to guarantee deterministic order
fn extract_and_write_variables(
    sequences: &[Vec<u8>],
    sample_names: &[String],
    output_filename: &str,
    variable_positions_info: &[VariablePositionInfo],
    system: &Arc<Mutex<System>>,
) -> io::Result<()> {
    let output_file = File::create(output_filename)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file);

    let variable_positions: Vec<usize> = variable_positions_info.iter().map(|info| info.position).collect();
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
        // Extract variable nucleotides
        let var_seq: Vec<u8> = variable_positions.iter().map(|&pos| seq[pos]).collect();
        let var_seq_str = String::from_utf8_lossy(&var_seq);

        writeln!(writer, ">{}", name)?;
        writeln!(writer, "{}", var_seq_str)?;

        let current = i + 1;
        if current % 10_000 == 0 || current == total_sequences {
            pb_write.set_position(current as u64);
        }

        // Report RAM usage every 100,000 sequences
        if current % 100_000 == 0 {
            if let Ok(mut sys) = system.lock() {
                sys.refresh_memory();
                let total_memory = sys.total_memory();
                let used_memory = sys.used_memory();
                println!(
                    "[{}] Written {} sequences. RAM usage: {} KB used / {} KB total.",
                    Local::now().format("%Y-%m-%d %H:%M:%S"),
                    current,
                    used_memory,
                    total_memory
                );
            }
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

/// Step 3: Generate a VCF 4.2 file with standard GT genotypes
/// FIX #3: ALT alleles separated by commas (VCF standard)
/// FIX #4: Standard GT format with allele indices instead of custom BASE format
fn generate_vcf_file(
    variable_positions_info: &[VariablePositionInfo],
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

    // Generate VCF entries
    for info in variable_positions_info {
        let chrom = "1";
        let pos = info.position + 1; // VCF is 1-based
        let id = ".";
        let ref_base = info.reference_base as char;

        // ALT alleles separated by commas; gaps represented as * (VCF 4.2 standard)
        let alt_strings: Vec<String> = info
            .alternate_bases
            .iter()
            .map(|&b| {
                let c = b as char;
                if c == '-' { "*".to_string() } else { c.to_string() }
            })
            .collect();
        let alt = alt_strings.join(",");

        let qual = ".";
        let filter = "PASS";
        let info_field = format!("NS={}", sample_names.len());

        // Build allele index map: REF=0, ALT1=1, ALT2=2, ...
        let mut allele_index: HashMap<u8, usize> = HashMap::new();
        allele_index.insert(info.reference_base, 0);
        for (i, &alt_base) in info.alternate_bases.iter().enumerate() {
            allele_index.insert(alt_base, i + 1);
        }

        // Standard GT genotypes (index-based)
        let genotypes: Vec<String> = info
            .genotypes
            .iter()
            .map(|&genotype_base| {
                let upper = genotype_base.to_ascii_uppercase();
                match allele_index.get(&upper) {
                    Some(idx) => idx.to_string(),
                    None => ".".to_string(), // Missing or ambiguous base
                }
            })
            .collect();

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT",
            chrom, pos, id, ref_base, alt, qual, filter, info_field
        )?;
        for gt in genotypes {
            write!(writer, "\t{}", gt)?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;

    Ok(())
}
