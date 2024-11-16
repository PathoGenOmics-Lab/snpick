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
use sysinfo::{System, SystemExt};
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
        // Step 1: Identify variable positions and extract individual genotypes
        println!("Starting Step 1: Identifying variable positions...");
        let (variable_positions_info, total_sequences, sample_names) =
            identify_variable_positions(&input_filename, &system, include_gaps)
                .expect("Failed to identify variable positions");
        println!(
            "Step 1 Completed: Found {} variable positions.",
            variable_positions_info.len()
        );

        if variable_positions_info.is_empty() {
            eprintln!("No variable positions found in the alignment.");
            std::process::exit(0);
        }

        // Step 2: Extract variable positions and write to output file
        println!("Starting Step 2: Extracting variable positions and writing to output...");
        extract_and_write_variables(
            &input_filename,
            &output_filename,
            &variable_positions_info,
            total_sequences,
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
    alternate_bases: HashSet<u8>,
    genotypes: Vec<u8>, // Bases at this position for each sample
}

/// Step 1: Identify variable positions and extract individual genotypes
fn identify_variable_positions(
    input_filename: &str,
    system: &Arc<Mutex<System>>,
    include_gaps: bool,
) -> io::Result<(Vec<VariablePositionInfo>, usize, Vec<String>)> {
    // Open the input FASTA file
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file);
    let mut fasta_reader = fasta::Reader::new(reader);

    // Read all sequences and store names and sequences
    let mut sequences = Vec::new();
    let mut sample_names = Vec::new();

    for result in fasta_reader.records() {
        let record = result?;
        let id = record.id().to_string();
        let seq = record.seq().to_owned();
        sample_names.push(id);
        sequences.push(seq);
    }

    let total_sequences = sequences.len();

    if total_sequences == 0 {
        eprintln!("The input FASTA file is empty.");
        return Ok((Vec::new(), 0, Vec::new()));
    }

    let seq_length = sequences[0].len();

    // Verify that all sequences have the same length
    for seq in &sequences {
        if seq.len() != seq_length {
            eprintln!("All sequences must have the same length.");
            return Ok((Vec::new(), 0, Vec::new()));
        }
    }

    println!(
        "[{}] Processing {} sequences of length {}.",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        total_sequences,
        seq_length
    );

    // Initialize a progress spinner
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing positions: {pos}/{len}")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    // Atomic counter for processed positions
    let pos_counter = AtomicUsize::new(0);

    // Identify variable positions and extract genotypes
    let variable_positions_info: Vec<VariablePositionInfo> = (0..seq_length)
        .into_par_iter()
        .filter_map(|pos| {
            let mut counts = HashMap::new();
            let mut genotypes = Vec::with_capacity(total_sequences);

            // Collect bases at this position for all samples
            for seq in &sequences {
                let nuc = seq[pos];
                genotypes.push(nuc);

                if nucleotide_to_bit(nuc, include_gaps).is_some() {
                    counts.entry(nuc).and_modify(|c| *c += 1).or_insert(1);
                }
            }

            // Count types of nucleotides A, C, G, T (and gap if included)
            let nucleotide_types = counts.len();

            // Update the progress spinner
            let current = pos_counter.fetch_add(1, Ordering::SeqCst) + 1;
            pb.set_position(current as u64);
            pb.set_length(seq_length as u64);

            if nucleotide_types > 1 {
                // Determine the reference base (most frequent)
                let mut max_count = 0;
                let mut reference_base = b'N';
                for (&nuc, &count) in &counts {
                    if count > max_count {
                        max_count = count;
                        reference_base = nuc;
                    }
                }

                // Alternate bases
                let alternate_bases: HashSet<u8> = counts
                    .keys()
                    .cloned()
                    .filter(|&nuc| nuc != reference_base)
                    .collect();

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

    pb.finish_with_message("Position processing completed.");

    println!(
        "[{}] Variable position identification completed.",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    println!("Total sequences processed: {}", total_sequences);

    Ok((variable_positions_info, total_sequences, sample_names))
}

/// Step 2: Extract variable positions and write to the output FASTA file
fn extract_and_write_variables(
    input_filename: &str,
    output_filename: &str,
    variable_positions_info: &[VariablePositionInfo],
    total_sequences: usize,
    system: &Arc<Mutex<System>>,
) -> io::Result<()> {
    // Open the input FASTA file again for the second step
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file);
    let fasta_reader = fasta::Reader::new(reader);

    // Open the output FASTA file for writing
    let output_file = File::create(output_filename)?;
    let writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file);
    let writer = Arc::new(Mutex::new(writer));

    // Get the variable positions
    let variable_positions: Vec<usize> = variable_positions_info.iter().map(|info| info.position).collect();

    // Initialize a progress spinner for Step 2
    let pb_write = ProgressBar::new_spinner();
    pb_write.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Writing sequences: {pos}/{len}")
            .unwrap(),
    );
    pb_write.enable_steady_tick(Duration::from_millis(100));

    // Atomic counter for written sequences
    let seq_counter = AtomicUsize::new(0);

    // Process and write sequences in parallel
    fasta_reader
        .records()
        .par_bridge()
        .try_for_each(|result| -> io::Result<()> {
            let record = result?;
            let seq = record.seq();

            // Extract variable nucleotides based on variable_positions
            let var_seq: Vec<u8> = variable_positions.iter().map(|&pos| seq[pos]).collect();

            // Convert to String
            let var_seq_str = String::from_utf8_lossy(&var_seq).to_string();

            // Safely write to the output file
            {
                let mut writer = writer.lock().unwrap();
                writeln!(writer, ">{}", record.id())?;
                writeln!(writer, "{}", var_seq_str)?;
            }

            // Update the progress spinner
            let current = seq_counter.fetch_add(1, Ordering::SeqCst) + 1;
            pb_write.set_position(current as u64);
            pb_write.set_length(total_sequences as u64);

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

            Ok(())
        })?;

    // Finish the progress spinner for Step 2
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

/// Step 3: Generate the VCF file with actual bases, including ambiguous bases and codons
fn generate_vcf_file(
    variable_positions_info: &[VariablePositionInfo],
    vcf_output_filename: &str,
    sample_names: &[String],
) -> io::Result<()> {
    // Open the output VCF file for writing
    let output_file = File::create(vcf_output_filename)?;
    let mut writer = BufWriter::new(output_file);

    // Write the VCF header
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=snpick")?;
    writeln!(writer, "##reference=.")?;
    writeln!(writer, "##contig=<ID=1,length={}>", variable_positions_info.len())?;
    writeln!(
        writer,
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=BASE,Number=1,Type=String,Description=\"Observed base at this position\">"
    )?;

    // Write the header line with sample names
    write!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample_name in sample_names {
        write!(writer, "\t{}", sample_name)?;
    }
    writeln!(writer)?;

    // Generate VCF entries
    for info in variable_positions_info {
        let chrom = "1"; // Change this if you have chromosome information
        let pos = info.position + 1; // Positions in VCF are 1-based
        let id = ".";
        let ref_base = info.reference_base as char;
        let mut alt_bases: Vec<char> = info.alternate_bases.iter().map(|&b| b as char).collect();
        alt_bases.sort();
        let alt = alt_bases.iter().collect::<String>().replace('-', ".");
        let qual = ".";
        let filter = "PASS";
        let info_field = format!("NS={}", sample_names.len());
        let format_field = "BASE"; // Using a custom field

        // Generate genotypes for each sample, using the actual bases
        let genotypes: Vec<String> = info
            .genotypes
            .iter()
            .map(|&genotype_base| {
                let base_char = genotype_base as char;
                base_char.to_string()
            })
            .collect();

        // Write the VCF line
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chrom, pos, id, ref_base, alt, qual, filter, info_field, format_field
        )?;
        for gt in genotypes {
            write!(writer, "\t{}", gt)?;
        }
        writeln!(writer)?;
    }

    // Ensure all data is written correctly
    writer.flush()?;

    Ok(())
}
