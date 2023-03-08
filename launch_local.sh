#!/usr/bin/env bash

nextflow run main.nf -resume -profile local_docker -params-file input_params.yaml