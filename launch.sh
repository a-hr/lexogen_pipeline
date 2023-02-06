#!/usr/bin/env bash

#SBATCH --time=20:00:00
#SBATCH --mem=2G

# load dependencies
module load Singularity Python Nextflow
# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx2G"
# launch the main process
nextflow run main.nf -resume -with-trace -profile cluster -params-file input_params.yaml