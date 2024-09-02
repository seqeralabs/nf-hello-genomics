# nf-hello-gatk

Welcome to the **nf-hello-gatk** repository! This is a demonstration pipeline built using [Nextflow](https://www.nextflow.io/), designed to showcase a basic genomic analysis workflow incorporating the [GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/). The pipeline is ideal for demonstrating how to handle input and output files via channels and pass them between processes effectively.

## Overview

The **nf-hello-gatk** pipeline performs a variant calling analysis using GATK HaplotypeCaller on a set of BAM files. It comes with test data located in [./data](./data/) that should run in seconds, allowing you to demonstrate the pipeline quickly. Furthermore, it has built-in support for Docker, which simplifies dependency management and ensures consistent execution environments.

## Getting Started

### Prerequisites

To run this pipeline, you need to have the following software installed:

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) (optional but recommended for containerized execution)
- [GATK](https://gatk.broadinstitute.org/) (handled internally by the pipeline)

### Pre-requisites

1. Ensure that you have Nextflow installed and accessible in your `PATH`.
1. (Preferred) Ensure Docker is installed and available
    a. If you do not wish to use Docker you must use an alternative method of handling software dependencies

### Running the Pipeline

To run the pipeline, use the following command:

```bash
nextflow run seqeralabs/nf-hello-gatk
```

If you wish you can manually supply your own parameters using command line options. These are the defaults specified from the root of the repository:

```bash
nextflow run seqeralabs/nf-hello-gatk \
    --bams "./data/bam/*.bam" \
    --reference ./data/ref/ref.fasta \
    --reference_index ./data/ref/ref.fasta.fai \
    --reference_dict ./data/ref/ref.dict \
    --calling_intervals data/ref/intervals.bed \
    --cohort_name my_cohort
```

This will run the pipeline using the supplied BAM files and reference data.

#### Parameters

The pipeline allows for the following input parameters:

- `--bams`: A glob pattern to specify the input BAM files.
- `--reference`: Path to the reference genome FASTA file.
- `--reference_index`: Path to the index file (`.fai`) of the reference genome.
- `--reference_dict`: Path to the dictionary file (`.dict`) of the reference genome.
- `--calling_intervals`: Path to the intervals file for variant calling.
- `--cohort_name`: A name for the cohort being analyzed (used in naming output files).

Example of running the pipeline:

```bash
nextflow run seqeralabs/nf-hello-gatk \
    --bams "./data/bams/*.bam" \
    --reference ./data/ref/hg38.fasta \
    --reference_index ./data/ref/hg38.fasta.fai \
    --reference_dict ./data/ref/hg38.dict \
    --calling_intervals ./data/ref/intervals.bed \
    --cohort_name sample_cohort
```

### Running with Docker

By default, the pipeline will use Docker for each process. This is enabled via the configuration option in the [nextflow.config](./nextflow.config). Nextflow will handle downloading the necessary Docker images and running the pipeline within containers.

If you wish to disable this, you can use the following configuration option:

```config
docker.enabled = false
```

Note: You will need to provide the software dependencies yourself or use an alternative method to manage them.

### Advanced Usage

For more advanced usage, such as customizing the workflow, modifying the process definitions, or integrating additional tools, you can edit the `main.nf` file or create custom configurations.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue if you have suggestions for improvements or find any bugs.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This pipeline is developed and maintained by [Seqera Labs](https://seqera.io/). It is inspired by the broader community of bioinformatics and computational biology professionals who contribute to open-source genomics tools.
