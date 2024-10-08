/// Struct describing paths to paired fastq data
///
/// * r1: path to forward reads
/// * r2: path to reverse reads
///
#[derive(Debug)]
pub struct FastqPath {
    pub r1: String,
    pub r2: String,
}
