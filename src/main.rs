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

    /// Consider the '-' symbol (gap) in SNP detection
    #[arg(short = 'g', long, help = "Consider the '-' symbol (gap) in SNP detection")]
    include_gaps: bool,

    /// Only consider A, T, G, C nucleotides in SNP detection
    #[arg(short = 'a', long, help = "Only consider A, T, G, C nucleotides in SNP detection")]
    only_atgc: bool,
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
    let include_gaps = args.include_gaps;
    let only_atgc = args.only_atgc;

    // Configurar un pool de hilos local para Rayon
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build Rayon thread pool");

    // Inicializar sistema para reportar uso de RAM
    let system = Arc::new(Mutex::new(System::new_all()));

    // Ejecutar el procesamiento dentro del pool de hilos
    pool.install(|| {
        // Paso 1: Identificar posiciones SNP
        println!("Starting Step 1: Identifying SNP positions...");
        let (snp_positions, total_sequences) =
            identify_snp_positions(&input_filename, &system, include_gaps, only_atgc)
                .expect("Failed to identify SNP positions");
        println!(
            "Step 1 Complete: Found {} SNP positions.",
            snp_positions.len()
        );

        if snp_positions.is_empty() {
            eprintln!("No SNPs found in the alignment.");
            std::process::exit(0);
        }

        // Paso 2: Extraer SNPs y escribir en el archivo de salida
        println!("Starting Step 2: Extracting SNPs and writing to output...");
        extract_and_write_snps(
            &input_filename,
            &output_filename,
            &snp_positions,
            total_sequences,
            &system,
        )
        .expect("Failed to extract and write SNPs");
        println!("SNP alignment written to {}", output_filename);
    });

    Ok(())
}

/// Paso 1: Identificar posiciones SNP iterando a través de todas las secuencias
fn identify_snp_positions(
    input_filename: &str,
    system: &Arc<Mutex<System>>,
    include_gaps: bool,
    only_atgc: bool,
) -> io::Result<(Vec<usize>, usize)> {
    // Abrir el archivo FASTA de entrada
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file); // Buffer de 16 MB
    let fasta_reader = fasta::Reader::new(reader);

    // Leer la primera secuencia para inicializar
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

    // Inicializar un spinner de progreso
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing sequences: {msg}")
            .unwrap(),
    );
    pb.enable_steady_tick(Duration::from_millis(100));

    // Contador atómico para secuencias procesadas (ya se procesó la primera)
    let seq_counter = Arc::new(AtomicUsize::new(1));

    // Reabrir el archivo FASTA para procesar todas las secuencias
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file);
    let fasta_reader = fasta::Reader::new(reader);

    // Procesar secuencias en paralelo usando Rayon
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

                    // Si solo se consideran ATGC, ignorar nucleótidos no estándar
                    if only_atgc {
                        // Verificar si el nucleótido es A, T, G o C
                        if bits == 0b0001 || bits == 0b0010 || bits == 0b0100 || bits == 0b1000 {
                            acc[i] |= bits;
                        } else {
                            // Ignorar nucleótidos no ATGC
                            continue;
                        }
                    } else {
                        acc[i] |= bits;
                    }
                }

                // Incrementar el contador de secuencias
                let current = seq_counter.fetch_add(1, Ordering::SeqCst) + 1;
                if current % 1000 == 0 {
                    pb.set_message(format!("{}", current));
                }

                // Reportar uso de RAM cada 100,000 secuencias
                if current % 100_000 == 0 {
                    if let Ok(mut sys) = system.lock() {
                        sys.refresh_memory();
                        let total_memory = sys.total_memory(); // en KB
                        let used_memory = sys.used_memory(); // en KB
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

    // Finalizar el spinner de progreso
    let total_sequences = seq_counter.load(Ordering::SeqCst);
    pb.finish_with_message(format!("{}", total_sequences));
    println!(
        "[{}] Processing completed: {} sequences.",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        total_sequences
    );

    // Identificar posiciones SNP
    println!(
        "[{}] Identifying SNP positions...",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    );
    let snp_positions: Vec<usize> = nucleotide_presence
        .par_iter()
        .enumerate()
        .filter_map(|(i, &bits)| {
            if include_gaps {
                if bits.count_ones() > 1 {
                    Some(i)
                } else {
                    None
                }
            } else {
                // Excluir el bit de gap (bit 5) de bits
                let bits_without_gap = bits & 0b1111;
                if bits_without_gap.count_ones() > 1 {
                    Some(i)
                } else {
                    None
                }
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

/// Paso 2: Extraer SNPs basados en las posiciones identificadas y escribir en el archivo FASTA de salida
fn extract_and_write_snps(
    input_filename: &str,
    output_filename: &str,
    snp_positions: &[usize],
    total_sequences: usize,
    system: &Arc<Mutex<System>>,
) -> io::Result<()> {
    // Abrir el archivo FASTA de entrada nuevamente para el segundo paso
    let input_file = File::open(input_filename)?;
    let reader = BufReader::with_capacity(16 * 1024 * 1024, input_file); // Buffer de 16 MB
    let fasta_reader = fasta::Reader::new(reader);

    // Abrir el archivo FASTA de salida para escribir
    let output_file = File::create(output_filename)?;
    let writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file); // Buffer de 16 MB
    let writer = Arc::new(Mutex::new(writer)); // Proteger el escritor con un Mutex para acceso concurrente

    // Inicializar una barra de progreso para el Paso 2
    let pb_write = ProgressBar::new(total_sequences as u64);
    pb_write.set_style(
        ProgressStyle::default_bar()
            .template(
                "[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} SNP sequences written",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Contador atómico para secuencias escritas
    let seq_count = Arc::new(AtomicUsize::new(0));

    // Procesar y escribir secuencias en paralelo
    fasta_reader
        .records()
        .par_bridge()
        .for_each(|result| {
            let record = match result {
                Ok(rec) => rec,
                Err(e) => {
                    eprintln!("Error reading a sequence: {}", e);
                    return;
                }
            };
            let seq = record.seq();

            // Extraer nucleótidos SNP basados en snp_positions
            let snp_seq: Vec<u8> = snp_positions.iter().map(|&pos| seq[pos]).collect();

            // Convertir a String
            let snp_seq_str = String::from_utf8_lossy(&snp_seq).to_string();

            // Escribir de forma segura en el archivo de salida
            {
                let mut writer = writer.lock().unwrap();
                if let Err(e) = writeln!(writer, ">{}", record.id()) {
                    eprintln!("Error writing header: {}", e);
                    return;
                }
                if let Err(e) = writeln!(writer, "{}", snp_seq_str) {
                    eprintln!("Error writing SNP sequence: {}", e);
                    return;
                }
            }

            // Incrementar el contador y actualizar la barra de progreso
            let current = seq_count.fetch_add(1, Ordering::SeqCst) + 1;
            pb_write.inc(1);

            // Reportar uso de RAM cada 100,000 secuencias
            if current % 100_000 == 0 {
                if let Ok(mut sys) = system.lock() {
                    sys.refresh_memory();
                    let total_memory = sys.total_memory(); // en KB
                    let used_memory = sys.used_memory(); // en KB
                    println!(
                        "[{}] Written {} SNP sequences. RAM Usage: {} KB used / {} KB total.",
                        Local::now().format("%Y-%m-%d %H:%M:%S"),
                        current,
                        used_memory,
                        total_memory
                    );
                }
            }
        });

    // Finalizar la barra de progreso para el Paso 2
    let final_count = seq_count.load(Ordering::SeqCst);
    pb_write.finish_with_message(format!(
        "Completed writing {} SNP sequences.",
        final_count
    ));

    println!(
        "[{}] Total SNP sequences written: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S"),
        final_count
    );

    // Asegurar que todos los datos se hayan escrito correctamente
    {
        let mut writer = writer.lock().unwrap();
        writer.flush()?;
    }

    Ok(())
}
