# micritebench

micritebench is a collection of simulated datasets used for benchmarking tools that detect microbial DNA/RNA against a human genome background.


## FASTQs

For each commonly pathogenic microbes we generate an *in silico* dilution series by combining reads simulated from human & microbial reference genomes then combining them in different proportions. This produces a series of fastqs whose contents can be identified by their names:

* **genome1**\_nreads\_**genome2**\_nreads\_**genome3**\_nreads\_etc.*R1*\_or\_*R2*.fastq


For example, FASTQ files containing 900 reads simulated from the telomere to telomere human reference genome + 100 from an EBV reference genome would look like:

* **humanT2T**\_900\_**ebv**\_100.R1.fastq

* **humanT2T**\_900\_**ebv**\_100.R2.fastq

## BAMs

We then align each FASTQ to different human reference genome analysis sets using common aligners. BAMs names encode the alignment strategy and reference genomes

**fastq_name**.**refgenome**.**alignment_strategy**.bam

For example, if we used bwa-mem2 to align the fastq above against the GDC GRCH38 reference genome which includes viral decoy sequences that can mess up some microbial identification pipelines ([Description](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files); [download](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834)), the resulting bam would be named:


**humanT2T_900_ebv_100.gdc_grch38_d1_vd1.bwamem2.bam**


Full descriptions of all [reference genomes](genomes.csv) and [alignment strategies](alignment_strategies.csv) used are available as CSV files.


## Read simulation approach

Art for read simulation. Due to ease of installation accross many operating systems. All simulations produce paired-reads. Because ART simulates reads per contig we also add a downsampling step to get the exact number of reads we're after.

This repo includes a rust CLI tool that automates the read simulation, downsampling, inter-species splicing, and alignment.

To function correctly it requires ART, bwa-mem2 & samtools to be installed and available on PATH.

## CLI example

Goal: Generate 3 combinations EBV and HPV genomes containing

1. 10 EBV reads + 20 HPV reads
2. 100 EBV reads + 30 HPV reads
3. 1000 EBV reads + 40 HPV reads

Plus coord-sorted BAMs created by alignment each set of fastqs to the 'metagenome.fna' reference (contains both an EBV and HPV reference genome) using bwa-mem2.

You can run this 
```
micritebenchmark \
    --ref1 genomes/ebv.fna \
    --ref1-key ebv \
    --ref1-reads 10,100,1000 \
    --ref2 genomes/hpv16.fna \
    --ref2-key hpv16 \
    --ref2-reads 20,30,40 \
    --alignment-genome genomes/metagenome.fna \
    --alignment-key metagenome \
    --threads 2 \
    --outdir outdir

```