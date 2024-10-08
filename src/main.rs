fn main() {
    // micritebench::simulate("genomes/test.fna", "test", 100, "workdir");
    // let spliced_fastqs = micritebench::simulation::simulate_spliced(
    //     "genomes/ebv.fna",
    //     "ebv",
    //     100,
    //     "genomes/hpv16.fna",
    //     "hpv16",
    //     10,
    //     "workdir",
    //     true,
    // );

    // eprintln!("{:#?}", spliced_fastqs);

    let spliced_fastqs = micritebench::simulation::simulate_spliced_many(
        "genomes/ebv.fna",
        "ebv",
        vec![10, 100, 1000],
        "genomes/hpv16.fna",
        "hpv16",
        vec![20, 30, 40],
        "workdir",
        "outdir",
    );

    eprintln!("{:#?}", spliced_fastqs)
}
