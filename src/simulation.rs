use crate::datastructures::FastqPath;
use crate::utils::prefix;

use bio::io::fastq;
use std::{collections::HashSet, fs::File, io, path::Path, usize};

/// Simulates and align multiple combinations of 2 reference genomes to a third 'alignment' genome.
///
/// This function uses ART to simulate paired-end reads from two different reference genomes and then combines them together.
/// Supplying multiple sets of read counts for each reference genome will generate a distinct set of paired-end fastq files for each combination.
///
/// Each set of read pairs
/// # Arguments
///
/// * `ref1_path` - Path to the first reference genome file.
/// * `ref1_key` - Identifier for the first reference genome (used in naming the output files).
/// * `ref1_reads` - Vector containing the number of reads to generate from the first reference genome for each paired combination.
/// * `ref2_path` - Path to the second reference genome file.
/// * `ref2_key` - Identifier for the second reference genome (used in naming the output files).
/// * `ref2_reads` - Vector containing the number of reads to generate from the second reference genome for each paired combination.
/// * `workdir` - Directory where intermediate files will be stored (automatically cleaned).
/// * `outdir` - Directory where the final fastq files will be stored.
/// * `alignment_refgenome` - Path to reference genome used for alignment
/// * `refgenome_key` - Identifier for the alignment reference genome (used in naming of output files)
/// * `threads` - Number of threads to use for alignment & bam sorting
///
/// # Returns
///
/// A vector of `FastqPath` structs, each representing the paths to the paired fastq files for each combination.
///
/// # Panics
///
/// The function will panic if `ref1_reads` and `ref2_reads` vectors have different lengths or if Mason is not installed.
///
/// # Examples
///
/// Generate 3 combinations EBV and HPV genomes containing
///
/// 1. 10 EBV reads + 20 HPV reads
/// 2. 100 EBV reads + 30 HPV reads
/// 3. 1000 EBV reads + 40 HPV reads
///
/// Plus coord-sorted BAMs created by alignment each set of fastqs to the 'metagenome.fna' reference using bwa-mem2
///
/// ```
/// micritebench::simulation::simulate_and_align(
///     "genomes/ebv.fna",
///     "ebv",
///     vec![10, 100, 1000],
///     "genomes/hpv16.fna",
///     "hpv16",
///     vec![20, 30, 40],
///     "workdir",
///     "outdir",
///     "genomes/metagenome.fna",
///     "metagenome",
///     2
/// );
/// ```
/// ///
///
/// Combines simulated reads from reference genome `ref1_path` (n=`ref1_reads`) with
/// those from `ref2_path` (n=`ref2_reads`) to create a spliced fastq file with a filename that encodes
/// its origin (genomes are indicated from `key` arguments).
///
/// The paired fastq files are then aligned to `alignment_refgenome` with bwa-mem2 and coordinate sorted
/// using samtools.
///  
///
#[allow(clippy::too_many_arguments)]
pub fn simulate_and_align(
    ref1_path: &str,
    ref1_key: &str,
    ref1_reads: Vec<u32>,
    ref2_path: &str,
    ref2_key: &str,
    ref2_reads: Vec<u32>,
    workdir: &str,
    outdir: &str,
    alignment_refgenome: &str,
    refgenome_key: &str,
    threads: u8,
) {
    //simulate multiple fastqs

    //align all of them
    let fastqs = crate::simulation::simulate_spliced_many(
        ref1_path, ref1_key, ref1_reads, ref2_path, ref2_key, ref2_reads, workdir, outdir,
    );

    for fastq in fastqs {
        crate::alignment::align(&fastq, alignment_refgenome, refgenome_key, outdir, threads)
    }
}

/// Simulates multiple combinations of 2 reference genomes.
///
/// This function uses ART to simulate paired-end reads from two different reference genomes and then combines them together.
/// Supplying multiple sets of read counts for each reference genome will generate a distinct set of paired-end fastq files for each combination.
///
/// # Arguments
///
/// * `ref1_path` - Path to the first reference genome file.
/// * `ref1_key` - Identifier for the first reference genome (used in naming the output files).
/// * `ref1_reads` - Vector containing the number of reads to generate from the first reference genome for each paired combination.
/// * `ref2_path` - Path to the second reference genome file.
/// * `ref2_key` - Identifier for the second reference genome (used in naming the output files).
/// * `ref2_reads` - Vector containing the number of reads to generate from the second reference genome for each paired combination.
/// * `workdir` - Directory where intermediate files will be stored (automatically cleaned).
/// * `outdir` - Directory where the final fastq files will be stored.
///
/// # Returns
///
/// A vector of `FastqPath` structs, each representing the paths to the paired fastq files for each combination.
///
/// # Panics
///
/// The function will panic if `ref1_reads` and `ref2_reads` vectors have different lengths or if Mason is not installed.
///
/// # Examples
///
/// Generate 3 combinations EBV and HPV genomes containing
///
/// 1. 10 EBV reads + 20 HPV reads
/// 2. 100 EBV reads + 30 HPV reads
/// 3. 1000 EBV reads + 40 HPV reads
///
/// ```
/// micritebench::simulation::simulate_spliced_many(
///     "genomes/ebv.fna",
///     "ebv",
///     vec![10, 100, 1000],
///     "genomes/hpv16.fna",
///     "hpv16",
///     vec![20, 30, 40],
///     "workdir",
///     "outdir",
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn simulate_spliced_many(
    ref1_path: &str,
    ref1_key: &str,
    ref1_reads: Vec<u32>,
    ref2_path: &str,
    ref2_key: &str,
    mut ref2_reads: Vec<u32>,
    workdir: &str,
    outdir: &str,
) -> Vec<FastqPath> {
    // Assert length of ref1 and ref2 reads are the same
    assert!(
        ref1_reads.len() == ref2_reads.len(),
        "ref1_reads must be the same length as ref2_reads. {} != {}",
        ref1_reads.len(),
        ref2_reads.len()
    );

    assert!(
        !ref1_reads.is_empty(),
        "ref1_reads must be at least 1 element long."
    );

    // Create a vector to hold FastqPath objects
    let mut fastqpaths: Vec<FastqPath> = Vec::new();
    for it in ref1_reads.iter().zip(ref2_reads.iter_mut()) {
        let (nreads_ref1, nreads_ref2) = it;

        eprintln!(
            "Simulating a paired-end fastq ({} reads from {} + {} reads from {})",
            nreads_ref1, ref1_key, nreads_ref2, ref2_key,
        );

        let spliced = simulate_spliced(
            ref1_path,
            ref1_key,
            *nreads_ref1,
            ref2_path,
            ref2_key,
            *nreads_ref2,
            workdir,
            outdir,
        );
        fastqpaths.push(spliced);
    }

    fastqpaths
}

/// Simulate paired-end illumina reads by combining simulated reads from 2 different reference genomes
///
/// Generates simulated datasets using Mason.
///
/// Mason must be found on path
///
/// * ref1_path: path to reference genome 1
/// * ref2_path: path to reference genome 2
/// * ref1_key: key to encode reference genome 1 in spliced fastq filename
/// * ref2_key: key to encode reference genome 2 in spliced fastq filename
/// * ref1_reads: number of reads to produce from refgenome 1
/// * ref2_reads: number of reads to produce from refgenome 2
/// * workdir: working directory where intermediary files are created
/// * cleanup: should all but the final fastqs be deleted from the workdir
#[allow(clippy::too_many_arguments)]
pub fn simulate_spliced(
    ref1_path: &str,
    ref1_key: &str,
    ref1_reads: u32,
    ref2_path: &str,
    ref2_key: &str,
    ref2_reads: u32,
    workdir: &str,
    outdir: &str,
) -> FastqPath {
    let ref1_fastqs = simulate(ref1_path, ref1_key, ref1_reads, workdir);
    let ref2_fastqs = simulate(ref2_path, ref2_key, ref2_reads, workdir);

    let spliced_fastqs = cat_fastqs(ref1_fastqs, ref2_fastqs, workdir);

    // Move files from workdir to outdir
    let moved_spliced_fastq = mv_fastq_to_dir(spliced_fastqs, outdir)
        .expect("Failed to move spliced fastq to output directory");

    // Rest of this function is just cleanup (deleting unused directories)
    let workdir_files = std::fs::read_dir(workdir).expect("Failed to list workdir contents");
    for file_result in workdir_files {
        let file = file_result.expect("Failed to read workdir contents");
        std::fs::remove_file(file.path()).expect("Failed to delete file");
    }

    moved_spliced_fastq
}

fn mv_fastq_to_dir(fastq: FastqPath, dirpath: &str) -> Result<FastqPath, std::io::Error> {
    let r1_newpath = mv_file_to_dir(std::path::Path::new(&fastq.r1), dirpath)?;
    let r2_newpath = mv_file_to_dir(std::path::Path::new(&fastq.r2), dirpath)?;

    Ok(FastqPath {
        r1: r1_newpath,
        r2: r2_newpath,
    })
}

/// Move a file to a new directory
fn mv_file_to_dir(filepath: &Path, dirpath: &str) -> Result<String, std::io::Error> {
    std::fs::create_dir_all(dirpath).expect("failed to create output directory");
    let filename = filepath.file_name().unwrap_or_default();
    let newfilename = format!("{}/{}", dirpath, filename.to_string_lossy());

    std::fs::rename(filepath, &newfilename)?;

    Ok(newfilename)
}

/// Concatenates fastq files
///
/// Concatenates two fastq files and names output file
/// by pasting the two together
///
/// Only works on uncompressed fastq files
///
/// Returns a FastqPath of the concatenated fastqs
fn cat_fastqs(fastq1: FastqPath, fastq2: FastqPath, workdir: &str) -> FastqPath {
    let f1r1_prefix = prefix(&fastq1.r1);
    let f1r2_prefix = prefix(&fastq1.r2);

    let f2r1_prefix = prefix(&fastq2.r1);
    let f2r2_prefix = prefix(&fastq2.r2);

    let output_filename_r1 = format!("{workdir}/{f1r1_prefix}_{f2r1_prefix}.R1.fq");
    let output_filename_r2 = format!("{workdir}/{f1r2_prefix}_{f2r2_prefix}.R2.fq");

    // Concatenate Read 1 FASTQs
    concatenate_files(&fastq1.r1, &fastq2.r1, &output_filename_r1)
        .expect("Failed to concatenate Read 1 FASTQs");

    // Concatenate Read 2 FASTQs
    concatenate_files(&fastq1.r2, &fastq2.r2, &output_filename_r2)
        .expect("Failed to concatenate Read 2 FASTQs");

    FastqPath {
        r1: output_filename_r1,
        r2: output_filename_r2,
    }
}

fn concatenate_files(source1: &str, source2: &str, destination: &str) -> io::Result<()> {
    //Open file1
    let mut file1 = std::fs::File::open(source1)?;

    let mut dest = File::options()
        .create(true)
        .write(true)
        .truncate(true)
        .open(destination)?;

    // Copy the contents of the first file into the destination file
    std::io::copy(&mut file1, &mut dest)?;

    // Open the second source file for reading
    let mut file2 = File::open(source2)?;

    // Reopen file in append mode
    let mut dest = File::options()
        .append(true)
        .truncate(false)
        .open(destination)?;

    // Copy the contents of the second file into the destination file
    std::io::copy(&mut file2, &mut dest)?;

    Ok(())

    // Open
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
            panic!("ART simulation failed for reference {ref_key}. Error: {err}")
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
fn subsample_fastq(path_in: &str, path_out: &str, n_reads: u32) {
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
