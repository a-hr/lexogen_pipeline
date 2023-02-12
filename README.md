# Lexogen pipeline

## Introduction

## Installation

## Usage

## Available arguments

Arguments can be passed to the pipeline trough the input_params.yaml file.

The mandatory arguments are:
* csv_dir: Directory containing the csv barcode files and the SAF file
* fasq_dir: Directory containing the fastq files
* output_dir: Directory where the output files will be stored
* index_dir: Directory containing the STAR index files (or where to store them if created)
* suffix: Suffic of the paired fastq files (default: "_R")
* extension: Suffix of the fastq files (default: ".fastq.gz")

The optional arguments are:
* get_bams: If set to true, the bam and bai files will be stored in the output dir (default: false)
* create_index: If set to true, the STAR index will be created (default: false)
* ref_gen: Path to the directory containing the reference .fa file and the .gtf file (default: "ref_gen")

