# Lexogen pipeline

> Developed by: Álvaro Herrero Reiriz  
> GitHub: [a-hr](https://github.com/a-hr)

- [Lexogen pipeline](#lexogen-pipeline)
  - [Introduction](#introduction)
  - [How to use the pipeline](#how-to-use-the-pipeline)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Installing the pipeline (for command line use)](#installing-the-pipeline-for-command-line-use)
    - [Installing the pipeline (for EPI2ME use)](#installing-the-pipeline-for-epi2me-use)
    - [Dependencies](#dependencies)
  - [Running the pipeline](#running-the-pipeline)
    - [Input files](#input-files)
    - [Usage for Nextflow/nf-core](#usage-for-nextflownf-core)
    - [Usage for EPI2ME](#usage-for-epi2me)
  - [Pipeline output](#pipeline-output)
  - [Available arguments](#available-arguments)
    - [Input Options](#input-options)
    - [Output Options](#output-options)
    - [Sample Options](#sample-options)
    - [Advanced Options](#advanced-options)

## Introduction

This pipeline is designed to process the sequencing results of targeted RNA experiments. It includes the following steps:

- Quality control of the samples using FastQC
- Demultiplexing of the samples by barcode using cutadapt
- UMI extraction using umi_tools
- Read alignment using STAR (optionally, the index can be created)
- UMI demultiplication using umi_tools
- Optional indexing of the BAM files to be viewed in IGV (enabled with the `get_bams: true` parameter)
- Target gene quantification using featureCounts. A SAF (tab-separated table) containing the target genes and their coordinates is required.

## How to use the pipeline

The pipeline is written in Nextflow, a workflow manager that allows to run the pipeline in a wide variety of systems. It is configured to be run either on a SLURM-managed HPC cluster or a local machine, though it can be run on a cloud instance or using other workload managers by editing the configuration file according to ![Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes).

There are two main ways to run the pipeline:

- Cloning the repository and running the pipeline manually through the command line using Nextflow or nf-core. (straightforward, but requires bioinformatics knowledge).
- Installing EPI2ME and running the pipeline through the EPI2ME interface. (user-friendly, but requires to install EPI2ME).

## Requirements

In order to run, at least 32GB of RAM are required (human reference index is about 27GB). It can run on Windows, MacOS or Linux. If using Windows, we recommend installing the Windows Subsystem for Linux (WSL).

## Installation

If using EPI2ME, [install it](https://labs.epi2me.io/installation/) on your system and follow the instructions on the [EPI2ME documentation](https://community.nanoporetech.com/protocols/epi2me-labs) to run the pipeline. To add it to your saved workflows simply copy this repository's URL and paste it on the "Add workflow" section of the EPI2ME interface.
> WARNING: as of July 2023, there is a bug in Docker Desktop that might cause pipelines to fail. For this reason, as stated by the development team, it is recommended to [install docker manually](https://docs.docker.com/engine/install/ubuntu/) through the command line (just a couple of commands are needed to have it up and running).

If using Nextflow/nf-core, clone the repository and install the dependencies. The easiest way to do so is using conda. The pipeline can be run on any system that supports Docker or Singularity. If using Windows, we recommend using the Windows Subsystem for Linux (WSL).

The pipeline is especially tailored to be run on a HPC cluster, though it can seamlessly be run on a local machine and, with some configuration, on a cloud instance.

### Installing the pipeline (for command line use)

The pipeline can be installed on any directory. Simply *cd* to the desired directory and clone the repository:

```bash
git clone https://github.com/a-hr/lexogen_pipeline.git
```

### Installing the pipeline (for EPI2ME use)

1. Open EPI2ME and go to the "View Workflows" tab.
2. On the top right corner, click on "Import workflow".
3. Paste the repository's URL (https://github.com/a-hr/lexogen_pipeline) on the pop-up window and click on "Download".
4. The workflow will be added to your saved workflows.

### Dependencies

If you are using EPI2ME, the program will guide you through the dependencies installation.
> WARNING: as of July 2023, there is a bug in Docker Desktop that might cause pipelines to fail. For this reason, as stated by the development team, it is recommended to [install docker manually](https://docs.docker.com/engine/install/ubuntu/) through the command line (just a couple of commands are needed to have it up and running).

If you are using Nextflow/nf-core, you will need to install them manually. The pipeline requires [Docker](https://docs.docker.com/engine/install/ubuntu/) or [Singularity](https://github.com/sylabs/singularity/releases) (mainly for HPC clusters) to run.

In order to install Nextflow, you can use the following command:

```bash
conda create -n Nextflow -c bioconda nextflow
```

Process specific dependencies are automatically downloaded as Docker or Singularity containers. Docker images will be automatically pulled while running. On the other hand, Singularity images are cached locally (default location: `{pipeline_dir}/containers`), and can be manually installed using the provided Makefile:

```bash
make pull
```

## Running the pipeline

### Input files

- FASTQ files: Illumina pair-end FASTQ files (either gzipped or not).
- fiveprime.csv: maps file names (without pair-end suffix "_R" and without extensions) to sample groups, if a single file pair is present, add it as well. This will rename files to classify samples for the output table (e.g. mHTT_WT-1a).
    > fiveprime.csv

    ```csv
    File;Name
    group_1;mHTT
    group_2;control
    ```

    Note that *File* is the name of the fastq file without the suffix (e.g. ``lib1_R1.fastq.gz`` -> ``lib1``) and *Name* is the arbitrary name of the group or experiment.

- threeprime.csv: maps barcodes to sample names for the demultiplexing step.
    > threeprime.csv

    ```csv
    Barcode;Sequence;Sample
    BC1;ATCG;WT-1a
    BC2;TGCA;WT-1b
    BC3;ACGT;KO-1a
    BC4;TGAC;KO-1b
    ```

    Note that *Barcode* is the arbitrary name of the barcode, *Sequence* that will be searched on the R2 reads (3' barcodes have higher quality in the reverse strand) and *Sample* is the name of the specific sample.

- SAF file: tab-separated table containing the target genes and their coordinates. The pipeline will use featureCounts to quantify the transcripts mapping to these genes.
    > target_genes.saf

    ```csv
    GeneID Chr Start End Strand
    gene1 chr1 1000 2000 +
    gene2 chr1 3000 4000 +
    gene3 chr1 5000 6000 +
    gene4 chr1 7000 8000 +
    ```

    Note that *GeneID* is the arbitrary name of the gene, *Chr* is the chromosome, *Start* and *End* are the coordinates of the gene and *Strand* is the orientation of the strand.

### Usage for Nextflow/nf-core

Running the pipeline is really straightforward:

1. First, add the input files to the input directory (default: `./inputs`):
   - The CSV files (./inputs/csvs) will contain the barcodes (``threeprime.csv``) and the libraries (``fiveprime.csv``). Additionally, the ``annotations.saf``file will also be placed here.
   - The FASTQ files (./inputs/fastqs) will contain the FASTQ files to be processed.
   - Locate an existing STAR reference index (verions 2.7.10b) or create a new one using the pipeline (requires reference genome `.fa` and annotation `.gtf` files).
 2. Then, fill in the `input_params.yaml` file with the desired parameters. The available parameters are described in the [Available arguments](#available-arguments) section.
 3. Next, load the conda environment. If you have Nextflow installed system-wide, you can skip this step.

    ```bash
    conda activate Nextflow
    ```

    Now, decide what profile to use:
      - If running on a local machine, use the `local_docker` profile. You can optionally use the `local_singularity` profile if you have Singularity installed.
      - If running on a HPC cluster, use the `cluster` profile. It uses the SLURM executor by default. If your cluster happens to use a different scheduler, you can change the executor on the `./confs/cluster.config` file. Note that it requires to have Singularity available on the main thread's environment.
 4. Finally, run the pipeline:

    ```bash
    nextflow run main.nf -resume -profile cluster -params-file input_params.yaml
    ```

    The `-resume` flag will resume the pipeline if it was interrupted.  
    If you are working on a cluster, it is good practice to launch the main script as a job itself (make sure that computing nodes can launch jobs in your cluster, contact your system admins if you are not sure). A sample job launching script is provided in launch_cluster.sh. You can use it as follows:

    ```bash
    sbatch launch_cluster.sh
    ```

### Usage for EPI2ME

The pipeline can be run on EPI2ME without any bioinformatics kwnoledge:

1. Open the EPI2ME interface and select "Add a new workflow". Add this custom workflow by pasting the repo link (https://github.com/a-hr/lexogen_pipeline.git).
2. Select the workflow and click on "Run workflow".
3. Fill in the parameters on the interface. The available parameters are described in the [Available arguments](#available-arguments) section.
4. Additionally, on the *Nextflow Configuration* tab, enter a specific profile:
    - `local_docker` for running on a local machine with Docker.
    - `local_singularity` for running on a local machine with Singularity.
    - `cluster` for running on a SLURM cluster with Singularity.
    - additional profiles can be added by creating a copy of any of the files in `./confs` and modifying it. Do no forget to add the new profile name and link its `.conf` in the `nextflow.config` file.
5. Run the workflow.

## Pipeline output

The pipeline will generate the following:

- a directory containing the expression tables both before and after deduplication.
- a directory containing `seqkit stat` runs after every step of the pipeline.
- a ``multiqc`` report of the pipeline processes.
- a performance report of the pipeline processes.
- optionally, a directory containing the final deduplicated BAM files, to be able to check the alignment on IGV.

## Available arguments

Arguments can be passed to the pipeline through the input_params.yaml file.

### Input Options

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `csv_dir` | .csv and .saf directory path <details><summary>Help</summary><small>Path to the directory where fiveprime.csv, three prime.csv and annotations.saf are located.<br></small></details>| `string` | lexogen_pipeline/inputs/csvs | True |
| `fastq_dir` | .fastq(.gz) files' directory path <details><summary>Help</summary><small>Path to the directory containing the FASTQ files (either gzipped or not) to be used.<br></small></details>| `string` | lexogen_pipeline/inputs/fastqs | True |
| `index_dir` | reference genome directory path <details><summary>Help</summary><small>Directory containing the reference STAR genome if already created. If not, the index created by the pipeline could be saved here.</small></details>| `string` | lexogen_pipeline/index | True |
| `suffix` | pair-end pattern <details><summary>Help</summary><small>Common pattern to detect both pair-end files, without the 1 and 2. For example: file_R01; file_R02 --> suffix = _R0</small></details>| `string` | _R |  |
| `extension` | FASTQ file extension <details><summary>Help</summary><small>Either .fastq or .fq. + gzip compression .gz</small></details>| `string` | .fastq.gz |  |

### Output Options

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `output_dir` | output directory path <details><summary>Help</summary><small>Path to the directory that will contain the output tables, BAMs, QC reports and logs.</small></details>| `string` | lexogen_pipeline/output |  |
| `get_bams` | whether to output aligned BAMs or not <details><summary>Help</summary><small>Useful to visualize individual samples with IGV to explore the reads.</small></details>| `boolean` | True |  |
| `save_index` | whether to save the generated index to the index_dir path | `boolean` |  |  |

### Sample Options

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `umi_length` | length of the UMI olives | `integer` | 6 | True |

### Advanced Options

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `create_index` | whether to create the STAR index or not <details><summary>Help</summary><small>Not necessary to create one if an index is already available (link it with index_dir parameter).</small></details>| `boolean` |  |  |
| `ref_gen` | directory containing index creation files <details><summary>Help</summary><small>A directory containing the reference FASTA and GTF annotations to create the STAR index.</small></details>| `string` | lexogen_pipeline/ref_gen |  |
| `map_size` | STAR parameter to create index | `integer` | 100 |  |
| `sharedMemory` | False when running on a cluster, True else. <details><summary>Help</summary><small>Enables the use of shared memory to optimize RAM usage of STAR by sharing the reference index across parallel processes. HPC clusters usually do not allow its use.</small></details>| `boolean` |  |  |
