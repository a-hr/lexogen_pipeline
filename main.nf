#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { stats } from './modules/seqkit.nf'
include { demultiplex } from './modules/cutadapt.nf'
include { fastqc } from './modules/fastqc.nf'
include { STAR_INDEX; STAR_ALIGN } from './modules/STAR.nf'
include { extract_UMI; dedup_UMI } from './modules/umi_tools.nf'
include { BAM_INDEX } from './modules/samtools.nf'
include { featureCounts } from './modules/subread.nf'


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


workflow {
    // star index creation
    STAR_INDEX(index_fa, index_annot)

    // raw file qc
    fastqc(paired_fastq_files)
    fastqc_multiqc = fastqc.out.collect()

    // fastq demultiplexing
    demultiplex(paired_fastq_files, csv_dir)
    demultiplex_multiqc = demultiplex.out.logs.collect()
    
    stats(
        demultiplex.out.collect(),
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

    // qc report
    multiqc(
        fastqc_multiqc,
        demultiplex_multiqc,
        alignment_multiqc,
        dedup_multiqc,
        featureCounts_multiqc
    )
}
