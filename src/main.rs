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

/// snpick: A tool to extract SNP sites from a FASTA alignment, considering IUPAC ambiguous nucleotides.
#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = "0.1.0",
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>",
    about = "A fast and memory-efficient SNP extraction tool."
)]
struct Args {
    /// Input FASTA alignment file
    #[arg(short, long, help = "Input FASTA alignment file")]
    fasta: String,

    /// Output FASTA file with SNPs
    #[arg(short, long, help = "Output FASTA file with SNPs")]
    output: String,

    /// Number of threads to use (optional)
    #[arg(short, long, default_value_t = 4, help = "Number of threads to use (optional)")]
    threads: usize,
}

static NUCLEOTIDE_BITS: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 0b0001;
    table[b'C' as usize] = 0b0010;
    table[b'G' as usize] = 0b0100;
    table[b'T' as usize] = 0b1000;
    table[b'R' as usize] = 0b0101; // A or G
    table[b'Y' as usize] = 0b1010; // C or T
    table[b'S' as usize] = 0b0110; // G or C
    table[b'W' as usize] = 0b1001; // A or T
    table[b'K' as usize] = 0b1100; // G or T
    table[b'M' as usize] = 0b0011; // A or C
    table[b'B' as usize] = 0b1110; // C or G or T
    table[b'D' as usize] = 0b1101; // A or G or T
    table[b'H' as usize] = 0b1011; // A or C or T
    table[b'V' as usize] = 0b0111; // A or C or G
    table[b'N' as usize] = 0b1111; // Any nucleotide (A, C, G, T)
    table[b'X' as usize] = 0b1111; // Any nucleotide (A, C, G, T)
    table
};

fn nucleotide_to_bits(nuc: u8) -> u8 {
    NUCLEOTIDE_BITS[nuc.to_ascii_uppercase() as usize]
}

fn main() -> io::Result<()> {
    // Parse command-line arguments
    let args = Args::parse();

    let input_filename = args.fasta;
    let output_filename = args.output;
    let num_threads = args.threads;

    // Configure the number of threads for Rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Initialize system for RAM usage reporting
    let system = Arc::new(Mutex::new(System::new_all()));

    // Step 1: Identify SNP positions
    println!("Starting Step 1: Identifying SNP positions...");
    let (snp_positions, total_sequences) = identify_snp_positions(&input_filename, &system)?;
    println!(
        "Step 1 Complete: Found {} SNP positions.",
        snp_positions.len()
    );

    if snp_positions.is_empty() {
        eprintln!("No SNPs found in the alignment.");
        std::process::exit(0);
    }

    // Step 2: Extract SNPs and write to output
    println!("Starting Step 2: Extracting SNPs and writing to output...");
    extract_and_write_snps(
        &input_filename,
        &output_filename,
        &snp_positions,
        total_sequences,
        &system,
    )?;
    println!("SNP alignment written to {}", output_filename);

    Ok(())
}

/// Step 1: Identify SNP positions by iterating through all sequences
fn identify_snp_positions(
    input_filename: &str,
    system: &Arc<Mutex<System>>,
) -> io::Result<(Vec<usize>, usize)> {
    // Open the input FASTA file
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file); // 16 MB buffer
    let fasta_reader = fasta::Reader::new(reader);

    // Read the first sequence to initialize
    let mut records = fasta_reader.records();
    let first_record = match records.next() {
        Some(Ok(rec)) => rec,
        Some(Err(e)) => {
            eprintln!("Error reading the first sequence: {}", e);
            return Ok((Vec::new(), 0));
        }
        None => {
            eprintln!("Input FASTA is empty.");
            return Ok((Vec::new(), 0));
        }
    };
    let first_seq = first_record.seq();
    let seq_length = first_seq.len();

    println!(
        "[{}] nucleotide_presence vector initialized with length {}.",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        seq_length
    );

    // Initialize a progress spinner
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing sequences: {msg}")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    // Atomic counter for sequences processed (already processed first)
    let seq_counter = Arc::new(AtomicUsize::new(1));

    // Reopen the FASTA file for processing all sequences
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file);
    let fasta_reader = fasta::Reader::new(reader);

    // Process sequences in parallel using Rayon
    let nucleotide_presence = fasta_reader
        .records()
        .par_bridge()
        .fold(
            || vec![0u8; seq_length],
            |mut acc, result| {
                let record = match result {
                    Ok(rec) => rec,
                    Err(e) => {
                        eprintln!("Error reading a sequence: {}", e);
                        return acc;
                    }
                };
                let seq = record.seq();

                for (i, &nuc) in seq.iter().enumerate() {
                    let bits = nucleotide_to_bits(nuc);
                    acc[i] |= bits;
                }

                // Increment the sequence counter
                let current = seq_counter.fetch_add(1, Ordering::SeqCst) + 1;
                if current % 1000 == 0 {
                    pb.set_message(format!("{}", current));
                }

                // Report RAM usage every 100,000 sequences
                if current % 100_000 == 0 {
                    if let Ok(mut sys) = system.lock() {
                        sys.refresh_memory();
                        let total_memory = sys.total_memory(); // in KB
                        let used_memory = sys.used_memory(); // in KB
                        println!(
                            "[{}] Processed {} sequences. RAM Usage: {} KB used / {} KB total.",
                            Local::now().format("%Y-%m-%d %H:%M:%S"),
                            current,
                            used_memory,
                            total_memory
                        );
                    }
                }

                acc
            },
        )
        .reduce(
            || vec![0u8; seq_length],
            |mut acc1, acc2| {
                for (i, &bit) in acc2.iter().enumerate() {
                    acc1[i] |= bit;
                }
                acc1
            },
        );

    // Finalize the progress spinner
    let total_sequences = seq_counter.load(Ordering::SeqCst);
    pb.finish_with_message(format!("{}", total_sequences));
    println!(
        "[{}] Processing completed: {} sequences.",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        total_sequences
    );

    // Identify SNP positions
    println!(
        "[{}] Identifying SNP positions...",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    let snp_positions: Vec<usize> = nucleotide_presence
        .par_iter()
        .enumerate()
        .filter_map(|(i, &bits)| {
            if bits.count_ones() > 1 {
                Some(i)
            } else {
                None
            }
        })
        .collect();

    println!(
        "[{}] SNP identification completed.",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    println!("Total sequences processed: {}", total_sequences);

    Ok((snp_positions, total_sequences))
}

/// Step 2: Extract SNPs based on identified positions and write to output FASTA
fn extract_and_write_snps(
    input_filename: &str,
    output_filename: &str,
    snp_positions: &[usize],
    total_sequences: usize,
    system: &Arc<Mutex<System>>,
) -> io::Result<()> {
    // Open the input FASTA file again for the second step
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file); // 16 MB buffer
    let fasta_reader = fasta::Reader::new(reader);

    // Open the output FASTA file for writing
    let output_file = File::create(output_filename)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file); // 16 MB buffer

    // Sequence counter
    let mut seq_count = 0usize;

    // Initialize a progress bar for Step 2
    let pb_write = ProgressBar::new(total_sequences as u64);
    pb_write.set_style(
        ProgressStyle::default_bar()
            .template(
                "[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} SNP sequences written",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Iterate through each record (sequence) in the FASTA file sequentially
    for result in fasta_reader.records() {
        let record = match result {
            Ok(rec) => rec,
            Err(e) => {
                eprintln!("Error reading a sequence: {}", e);
                continue;
            }
        };
        let seq = record.seq();

        // Extract SNP nucleotides based on snp_positions
        let snp_seq: Vec<u8> = snp_positions.iter().map(|&pos| seq[pos]).collect();

        // Write the header
        if let Err(e) = writeln!(writer, ">{}", record.id()) {
            eprintln!("Error writing header: {}", e);
            continue;
        }

        // Write the SNP sequence
        if let Err(e) = writeln!(writer, "{}", String::from_utf8_lossy(&snp_seq)) {
            eprintln!("Error writing SNP sequence: {}", e);
            continue;
        }

        seq_count += 1;
        pb_write.inc(1);

        // Report RAM usage every 100,000 sequences
        if seq_count % 100_000 == 0 {
            if let Ok(mut sys) = system.lock() {
                sys.refresh_memory();
                let total_memory = sys.total_memory(); // in KB
                let used_memory = sys.used_memory(); // in KB
                println!(
                    "[{}] Written {} SNP sequences. RAM Usage: {} KB used / {} KB total.",
                    Local::now().format("%Y-%m-%d %H:%M:%S"),
                    seq_count,
                    used_memory,
                    total_memory
                );
            }
        }
    }

    // Finalize the progress bar for Step 2
    pb_write.finish_with_message(format!(
        "Completed writing {} SNP sequences.",
        seq_count
    ));

    println!(
        "[{}] Total SNP sequences written: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        seq_count
    );

    writer.flush()?;
    Ok(())
}
