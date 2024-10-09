use std::fs;

use std::path::Path;

use crate::datastructures::FastqPath;
use crate::utils::prefix;

pub fn align(fastq: &FastqPath, refgenome: &str, refgenome_key: &str, outdir: &str, threads: u8) {
    // Set up output filepaths
    let fastq_prefix = prefix(&fastq.r1);
    let outfile_sam = format!("{outdir}/{fastq_prefix}.{refgenome_key}.bwamem2.sam");

    // Check BWA MEM 2 is Installed
    let bwa_path = which::which("bwa-mem2").unwrap_or_else(|err|{
           panic!("Failed to find bwa-mem2 binary. Please ensure bwa-mem2 is installed and added to path. Error: {err}")
       });

    // Check reference genomes exist

    // Check the genome is indexed
    assert_file_exists(refgenome);
    assert_indexed(refgenome);

    //    bwa-mem2 mem [options] <idxbase> <in1.fq> [in2.fq]
    let bwamem = std::process::Command::new(bwa_path)
        .arg("mem")
        .arg("-t")
        .arg(threads.to_string())
        .arg("-o")
        .arg(&outfile_sam)
        .arg(refgenome)
        .arg(&fastq.r1)
        .arg(&fastq.r2)
        .output()
        .expect("Failed to run bwa-mem2 command to align fastq files");

    if !bwamem.status.success() {
        let stderr = String::from_utf8_lossy(&bwamem.stderr);
        panic!(
            "BWA mem2 alignment ran with exit code {}. \nStderr:\n{}",
            bwamem.status.code().unwrap(),
            stderr
        )
    } else {
        eprintln!("Alignment successful")
    }

    // Now we need to sort and
    let bampath = sam_to_sorted_bam(&outfile_sam, true);

    eprintln!("bampath: {bampath}")
    // let taken = bwamem.stdout.take();
}

fn assert_indexed(idxname: &str) {
    let expected_index_file = format!("{idxname}.bwt.2bit.64");
    let exists = std::path::Path::new(&expected_index_file).exists();

    assert!(
        exists,
        "Could not find reference index {expected_index_file}. Please run `bwa-mem2 index {idxname}` to create the index"
    )
}

fn assert_file_exists(file: &str) {
    assert!(
        std::path::Path::new(file).exists(),
        "File does not exist: [{}]",
        file
    )
}

fn sam_to_sorted_bam(path_sam: &str, delete_sam: bool) -> String {
    // Filenames
    let inpath = Path::new(path_sam);
    assert!(inpath.exists(), "Could not find SAM file {}", path_sam);

    // Get the parent directory (or current directory if none)
    let parent_dir = inpath.parent().unwrap_or(Path::new("."));
    // Get the file stem (filename without extension)
    let stem = inpath
        .file_stem()
        .expect("Failed to get the stem of the SAM file");

    // Construct the output BAM file path in the same directory
    let outbam = parent_dir.join(format!("{}.bam", stem.to_string_lossy()));

    eprintln!("Will save file to: {}", outbam.display());

    // Check samtools is installed
    let samtools_path = which::which("samtools").unwrap_or_else(|err| {
        panic!(
            "Failed to find samtools binary. Please ensure samtools is installed and added to PATH. Error: {}",
            err
        )
    });

    // Run samtools
    let output = std::process::Command::new(samtools_path)
        .arg("sort")
        .args(["-O", "BAM"])
        .args(["-o", outbam.to_str().unwrap()])
        .arg(path_sam)
        .output()
        .unwrap_or_else(|err| {
            panic!(
                "Samtools sorting failed for SAM file {}. Error: {}",
                path_sam, err
            )
        });

    // Check Samtools run was successful
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "Samtools sort ran with exit code {}. \nStderr:\n{}",
            output.status.code().unwrap_or(-1),
            stderr
        )
    }

    // Delete Sam Input File
    if delete_sam {
        fs::remove_file(path_sam).expect("Failed to delete sam file after use");
    }

    // Return the output BAM file path as a String
    outbam.to_str().unwrap().to_string()
}
