fn main() {
    micritebench::simulation::simulate_and_align(
        "genomes/ebv.fna",
        "ebv",
        vec![10, 100, 1000],
        "genomes/hpv16.fna",
        "hpv16",
        vec![20, 30, 40],
        "workdir",
        "outdir",
        "genomes/test.fna",
        "test",
        2,
    );
}
