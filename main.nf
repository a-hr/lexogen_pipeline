#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { stats } from './modules/seqkit'
include { demultiplex; trim_adapters } from './modules/cutadapt'
include { fastqc } from './modules/fastqc'
include { multiqc } from './modules/multiqc'
include { STAR_INDEX; STAR_ALIGN } from './modules/STAR'
include { extract_UMI; dedup_UMI } from './modules/UMI'
include { BAM_INDEX } from './modules/samtools'
include { featureCounts } from './modules/subread'
include { group_results; plot_results } from './modules/visualization'


// VARIABLE DECLARATIONS
Channel
    .fromFilePairs("$params.fastq_dir/*$params.suffix{1,2}*$params.extension", checkIfExists: true)
    .set { paired_fastq_files }
Channel
    .value(params.csv_dir)
    .set { csv_dir }
Channel
    .value(params.index_dir)
    .set { index_dir }

if (params.create_index) {
    Channel
        .fromPath("$params.ref_gen/*.fa", checkIfExists: true)
        .set { index_fa }
    Channel
        .fromPath("$params.ref_gen/*.gtf", checkIfExists: true) 
        .set { index_annot }
}
else {
    Channel
        .empty()
        .set { index_fa }
    Channel
        .empty()
        .set { index_annot }
}

process input_logger {
    publishDir "${params.output_dir}/input_logs", mode: 'copy'

    input:
        path lib_csv
        path bc_csv
        path saf
        path input_params
        path fastq_dir

    output:
        path "*.txt"
        path "*.csv", includeInputs: true
        path "*.yaml", includeInputs: true
    
    """
    ls $fastq_dir > input_files.txt
    """
}


workflow {
    /* INPUT LOGGING */
    input_logger(
        Channel.fromPath("$params.csv_dir/fiveprime.csv"),
        Channel.fromPath("$params.csv_dir/threeprime.csv"),
        Channel.fromPath("$params.csv_dir/*.saf"),
        Channel.fromPath("$baseDir/*.yaml"),
        params.fastq_dir
    )

    // star index creation
    STAR_INDEX(index_fa, index_annot)

    // TODO trim adapters
    paired_fastq_files \
    | trim_adapters \
    | set { trimmed_fastq_files }

    // raw file qc
    fastqc(trimmed_fastq_files)
    fastqc_multiqc = fastqc.out.collect()
    
    // fastq demultiplexing
    demultiplex(trimmed_fastq_files, csv_dir)
    demultiplex_multiqc = demultiplex.out.logs.collect()
    
    stats(
        demultiplex.out.fastqs.collect(),
        "bc_demultiplexed"
    )

    // UMI sequence extraction
    extract_UMI(demultiplex.out.fastqs.flatten())

    // read alignment
    index_dir = params.create_index ? STAR_INDEX.out.index_dir : params.index_dir

    STAR_ALIGN(extract_UMI.out, index_dir)
    alignment_multiqc = STAR_ALIGN.out.logs.collect()

    // deduplication of aligned bams
    dedup_UMI(STAR_ALIGN.out.bams)
    dedup_multiqc = dedup_UMI.out.logs.collect()
    
    // reindexing of bams for IGV
    BAM_INDEX(dedup_UMI.out.dedup_bams)

    // expression quantification
    featureCounts(
        dedup_UMI.out.dedup_bams.collect(),
        STAR_ALIGN.out.bams.collect(),
        csv_dir
    ) 
    featureCounts_multiqc = featureCounts.out.logs.collect()

    // result visualization
    featureCounts.out.counts.collect() \
    | group_results

    featureCounts.out.counts.collect() \
    | plot_results

    // qc report
    multiqc(
        fastqc_multiqc,
        demultiplex_multiqc,
        alignment_multiqc,
        dedup_multiqc,
        featureCounts_multiqc
    )
}
