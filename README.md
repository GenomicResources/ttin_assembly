Repository name: Illumina Assembly
Maintainer: Jorge Langa - jorge.langa.arranz at gmail dot com


This repository contains tools and scritps for:
- Illumina adaptor finding (custom script)
- Read trimming (trimmomatic and seqclean)
- Transcriptome assembly (trinity)

The entire pipeline is packed into a Snakefile (a mix between of GNU make and Python3)


Requirements
------------
- GNU Parallel (for simple parallelization)
- Seqclean and associated binaries (psx, cdbyank, cdbfasta, etc) present in the path (please test the missing libraries of each binary with `ldd`; I think it was libc6-dev the one missing)
- Trinity associated software (bowtie1|2, samtools, etc.)
- Snakemake (pip3 install snakemake; python3 required)


Usage
----------
1. Create a folder called `data/fastq_raw`
2. Use the `scripts/find_adaptors.sh` script to determine the set of Illumina adaptors used. Output is on `data/adaptors` folder
3. Modify the Snakefile in a fashion that:
    - The `SAMPLES` variable contains the names of the samples to be analyzed, i.e., if your files are named sample_1.fastq.gz and sample_2.fastq.gz, then `SAMPLES="sample1"`
    - The `ADAPTORS` variable points towards the desired adaptors file in `./src/trimmomatic-0.33/adapters/file.fa`
4. Type `snakemake` to run the pipeline

