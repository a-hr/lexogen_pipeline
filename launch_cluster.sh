#!/usr/bin/env bash

#SBATCH --time=20:00:00
#SBATCH --mem=8G

# load dependencies (if using singularity on a cluster, the host machine should have
#  Singularity and GoLang installed to download and run the containers)
module load Singularity Nextflow Go

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx8G"

# launch the main process
nextflow run main.nf -resume -profile cluster -params-file input_params.yaml