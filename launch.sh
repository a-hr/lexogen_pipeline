#!/usr/bin/env bash

#SBATCH --time=20:00:00
#SBATCH --mem=8G

# load dependencies
module load Singularity Nextflow

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx8G"

# disable mmap for leveldb to avoid cluster-specific error
export NXF_OPTS="-Dleveldb.mmap=false"

# launch the main process
nextflow run main.nf -profile cluster -params-file input_params.yaml