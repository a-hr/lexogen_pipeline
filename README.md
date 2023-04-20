# Lexogen pipeline

> Developed by: Álvaro Herrero Reiriz  
> GitHub: [a-hr](https://github.com/a-hr)

- [Lexogen pipeline](#lexogen-pipeline)
  - [Introduction](#introduction)
  - [Installation](#installation)
    - [Installing dependencies](#installing-dependencies)
    - [Installing the pipeline](#installing-the-pipeline)
  - [Usage](#usage)
  - [Available arguments](#available-arguments)

## Introduction

This pipeline is designed to process the sequencing results of target RNA experiments. It includes the following steps:

* Quality control of the samples using FastQC
* Demultiplexing of the samples by barcode using cutadapt
* UMI extraction using umi_tools
* Read alignment using STAR (optionally, the index can be created)
* UMI demultiplication using umi_tools
* Optional indexing of the bam files to be viewed in IGV (enabled with the `get_bams: true` parameter)
* Target gene quantification using featureCounts. A SAF (tab-separated table) containing the target genes and their coordinates is required.

## Installation

To install the pipeline, clone the repository and install the dependencies. The easiest way to do so is using conda. The pipeline can be run on any system that supports Docker or Singularity. If using Windows, we recommend using the Windows Subsystem for Linux (WSL).

The pipeline is especially tailored to be run on a HPC cluster, though it can seamlessly be run on a local machine and, with some configuration, on a cloud instance.

### Installing dependencies

In order to install the dependencies, you can use the following commands:

```bash
conda create -n Nextflow -c bioconda nextflow
```

Process specific dependencies are automatically downloaded as Docker or Singularity containers. The latter are cached locally (default location: `pipeline_dir/containers`), and can be manually installed using the provided Makefile:

```bash
make pull
```

### Installing the pipeline

The pipeline can be installed on any directory. Simply *cd* to the desired directory and clone the repository:

```bash
git clone https://github.com/a-hr/lexogen_pipeline.git
```


## Usage

Running the pipeline is really straightforward:

1. First, add the input files to the input directory (default: `./inputs`):
   * The CSV files (./inputs/csvs) will contain the barcodes (``threeprime.csv``) and the libraries (``fiveprime.csv``). The CSV files must be semi-colon separated and have the following format:
        > fiveprime.csv

        ```csv
         File;Name
         lib1;WT
         lib2;KO
        ```

        Note that *File* is the name of the fastq file without the suffix (e.g. ``lib1_R1.fastq.gz`` -> ``lib1``) and *Name* is the arbitrary name of the group or experiment.

        > threeprime.csv

        ```csv
        Barcode;Sequence;Sample
        BC1;ATCG;WT-1a
        BC2;TGCA;WT-1b
        BC3;ACGT;KO-1a
        BC4;TGAC;KO-1b
        ```

        Note that *Barcode* is the arbitrary name of the barcode, *Sequence* that will be searched on the R2 reads (3' barcodes have higher quality in the reverse strand) and *Sample* is the name of the specific sample.
        
    * The SAF file (./inputs/csvs) will contain the target genes and their coordinates. The SAF file must be tab-separated and have the following format:

        ```saf
        Geneid	Chr	Start	End Strand
        gene1	5	35069758	35069874	+
        gene2	6	125138815	125138970	-
        ```
        Note that *Geneid* is the arbitrary name of the gene, *Chr* the chromosome, *Start* and *End* the coordinates and *Strand* the strand sense.

    * The fastq files (./inputs/fastqs) must be paired-end. Can be either gzipped or not. The pipeline will automatically detect the suffix and the extension of the files. The default suffix is ``_R`` and the default extension is ``.fastq.gz``. These can be changed in the ``input_params.yaml`` file.
 
 2. Then, fill in the `input_params.yaml` file with the desired parameters. The available parameters are described in the [Available arguments](#available-arguments) section.
 
 3. Next, load the conda environment. If you have Nextflow installed system-wide, you can skip this step.
 
    ```bash
    conda activate Nextflow
    ```

4. Now, decide what profile to use:
   * If running on a local machine, use the `local_docker` profile. You can optionally use the `local_singularity` profile if you have Singularity installed.
   * If running on a HPC cluster, use the `cluster` profile. It uses the SLURM executor by default. If your cluster happens to use a different scheduler, you can change the executor on the `./confs/cluster.config` file. Note that it requires to have Singularity available on the main thread environment.

 5. Finally, run the pipeline:
 
    ```bash
    nextflow run main.nf -resume -profile cluster -params-file input_params.yaml
    ```

    The `-resume` flag will resume the pipeline if it was interrupted.  
    If you are working on a cluster, it is good practice to launch the main script as a job itself (make sure that computing nodes can launch jobs in your cluster, contact your systems admins if you are not sure). A sample job launching script is provided in launch_cluster.sh. You can use it as follows:

    ```bash
    sbatch launch_cluster.sh
    ```


## Available arguments

Arguments can be passed to the pipeline through the input_params.yaml file.

The mandatory arguments are:
* csv_dir: Directory containing the csv barcode files and the SAF file (default: "./inputs/csvs")
* fasq_dir: Directory containing the fastq files (default: "./inputs/fastqs")
* output_dir: Directory where the output files will be stored (default: "./output")
* index_dir: Directory containing the STAR index files (or where to store them if created)
* suffix: Suffix of the paired fastq files (default: "_R")
* extension: Extension of the fastq files (default: ".fastq.gz")

The optional arguments are:
* get_bams: If set to true, the bam and bai files will be stored in the `output/bams/` (default: false)
* create_index: If set to true, the STAR index will be created. Note that it will require for the reference fasta and annotations to be provided through *ref_gen* (default: false)
* ref_gen: Path to the directory containing the reference .fa file and the .gtf file (default: "./ref_gen")
* sharedMemory: Whether to run STAR with the genome loaded in shared memory accross processes or not. It should e turned off if running on a cluster (default: false)

