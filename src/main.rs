use clap::Parser;

#[derive(Parser)]
#[command(name = "micritebenchmark")]
#[command(version = "0.0.1")]
#[command(about = "Create simulated datasets for benchmarking microbial detection toolkits", long_about = None)]
struct Cli {
    /// Path to the first reference genome file.
    #[arg(long)]
    ref1: String,

    /// Identifier for the first reference genome (used when naming outputs)
    #[arg(long)]
    ref1_key: String,

    /// Comma delimited list describing number of reads to generate from the first reference genome for each paired combination.
    #[arg(long, num_args(1..=10), value_delimiter=',')]
    ref1_reads: Vec<u32>,

    /// Path to the second reference genome file.
    #[arg(long)]
    ref2: String,

    /// Identifier for the second reference genome (used when naming outputs)
    #[arg(long)]
    ref2_key: String,

    /// Comma delimited list describing number of reads to generate from the second reference genome for each paired combination (e.g. 10,20,40)
    #[arg(long, num_args(1..=10), value_delimiter=',')]
    ref2_reads: Vec<u32>,

    /// Path to reference genome used for alignment file.
    #[arg(long)]
    alignment_genome: String,

    /// Identifier for the alignment reference genome (used when naming outputs)
    #[arg(long)]
    alignment_key: String,

    // Number of cores
    #[arg(long)]
    threads: Option<u8>,

    /// Outut directory
    #[arg(short, long)]
    outdir: String,
}

fn main() {
    let cli = Cli::parse();

    //Preprocess Args
    let threads = cli.threads.unwrap_or(1);
    micritebench::simulation::simulate_and_align(
        cli.ref1.as_str(),
        cli.ref1_key.as_str(),
        cli.ref1_reads,
        cli.ref2.as_str(),
        cli.ref2_key.as_str(),
        cli.ref2_reads,
        "workdir",
        cli.outdir.as_str(),
        cli.alignment_genome.as_str(),
        cli.alignment_key.as_str(),
        threads,
    );

    eprintln!("======================",);
    eprintln!("See [{}] for results", cli.outdir);
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert();
}
