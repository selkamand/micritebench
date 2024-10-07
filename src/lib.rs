use bio::io::fastq;
use std::{collections::HashSet, process::Command, usize};

pub struct FastqPath {
    r1: String,
    r2: String,
}
/// Simulate paired-end illumina reads by combining simulated reads from 2 different reference genomes
///
/// Generates simulated datasets using Mason.
///
/// Mason must be found on path
///
/// ref1: path to reference genome 1
/// ref2: path to reference genome 2
/// ref1_key: key to encode reference genome 1 in spliced fastq filename
/// ref2_key: key to encode reference genome 2 in spliced fastq filename
/// ref1_reads: number of reads to produce from refgenome 1
/// ref2_reads: number of reads to produce from refgenome 2
pub fn simulate_spliced(
    ref1: &str,
    ref1_key: &str,
    ref2: &str,
    ref2_key: &str,
    ref1_reads: u32,
    ref2_reads: u32,
) {
}

/// Simulate paired-end illumina reads by combining simulated reads from 2 different reference genomes
///
/// Generates simulated datasets using Mason.
///
/// Mason must be found on path
///
/// * ref1: path to reference genome
///
/// * ref1_key: key to encode reference genome name fastq filename
///
/// * ref1_reads: number of reads to produce from refgenome 1
///
/// # Returns
/// Path the the fastq file produced
pub fn simulate(ref_path: &str, ref_key: &str, reads: u32, workdir: &str) -> FastqPath {
    // Create Working Directory
    std::fs::create_dir_all(workdir).expect("Failed to create working directory");

    // Setup outfile prefixes
    let outfile_prefix = format!("{workdir}/{ref_key}_{reads}");
    let outfile_prefix_pre_subsample = format!("{outfile_prefix}.pre_subsample.R");

    // Check ART is Installed
    let art_path = which::which("art_illumina").unwrap_or_else(|err|{
        panic!("Failed to find art_illumina binary. Please ensure ART is installed and added to path. Error: {err}")
    }
    );
    eprintln!("Found ART Path: {:#?}", art_path);

    // Run ART simulation
    eprintln!("Running ART read simulation");
    std::process::Command::new(art_path)
        .arg("--noALN")
        .arg("--paired")
        .arg("--quiet")
        .arg("--len")
        .arg("150")
        .arg("--id")
        .arg(ref_key)
        .arg("--rcount")
        .arg(reads.to_string())
        .arg("--mflen")
        .arg("400")
        .arg("--sdev")
        .arg("10")
        .arg("--in")
        .arg(ref_path)
        .arg("--out")
        .arg(&outfile_prefix_pre_subsample)
        .output()
        .unwrap_or_else(|err| {
            panic!("Mason simulation failed for reference {ref_key}. Error: {err}")
        });

    // Since if 'reads' = 100 ART will produce 100 reads per contig, we then need to subsample down
    // to the absolute value of 'reads'
    eprintln!("Downsampling simulated reads to {reads} reads");
    let outfile_r1 = format!("{outfile_prefix_pre_subsample}1.fq");
    let outfile_r2 = format!("{outfile_prefix_pre_subsample}2.fq");

    let outfile_r1_subsampled = format!("{outfile_prefix}.R1.fq");
    let outfile_r2_subsampled = format!("{outfile_prefix}.R2.fq");

    subsample_fastq(&outfile_r1, &outfile_r1_subsampled, reads);
    subsample_fastq(&outfile_r2, &outfile_r2_subsampled, reads);

    // Return the path to the produced fastq files
    FastqPath {
        r1: outfile_r1,
        r2: outfile_r2,
    }
}

/// Subsample a fastq file
///
/// Creates a subsampled fastq file with `n_reads`
pub fn subsample_fastq(path_in: &str, path_out: &str, n_reads: u32) {
    // First pass to count how many reads we have
    let n_reads_observed = {
        let reader = fastq::Reader::from_file(path_in).expect("Failed to read fastq file");
        reader.records().count()
    };

    // Check if we have enough reads to downsample to n_reads
    if (n_reads_observed as u32) < n_reads {
        panic!("Cannot subsample fastq file [{path_in}] down to {n_reads} reads since it only contains {n_reads_observed} reads to begin with");
    }

    // Randomly generate indices for subsampling
    let mut rng = rand::thread_rng();
    let indices = rand::seq::index::sample(&mut rng, n_reads_observed, n_reads as usize);

    // Convert IndexVec into a HashSet for O(1) lookup
    let indices_set: HashSet<usize> = indices.into_iter().collect();

    // Second pass to read the records and write the sampled ones
    let reader = fastq::Reader::from_file(path_in).expect("Failed to read fastq file");
    let mut writer = fastq::Writer::to_file(path_out).expect("Failed to create fastq writer");

    // Subsample records
    for (i, record_result) in reader.records().enumerate() {
        if indices_set.contains(&i) {
            let record = record_result.expect("Error reading record");
            writer
                .write_record(&record)
                .expect("Error writing subsampled fastq read");
        }
    }
}
