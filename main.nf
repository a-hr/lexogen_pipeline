#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqc {
    label 'process_low'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
    output:
        path "*.zip"

    """
    fastqc $fastqs -o .
    """
}

process multiqc {
    label 'process_low'
    publishDir "${params.output}", mode: 'copy'

    input:
        path fastqc_logs
        path alignment_logs

    output:
        path "*.html"
    
    """
    multiqc . -n Lexogen_pipeline_multiqc_report.html
    """
}

process demultiplex {
    label 'process_high'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
        each csv_dir
    output:
        path "*.gz"

    """
    demultiplex.py -i $csv_dir -p . -s ${params.suffix}1 ${params.suffix}2
    """
}

process seqkit_demultiplex {
    label 'process_high'
    publishDir "$params.output", pattern: "*.txt", mode: "copy"

    input:
        path dmplx_fastqs
    output:
        path "*.gz", includeInputs: true, emit: fastqs 
        path "*.txt", emit: logs

    """
    seqkit stat *R1.fastq.gz > demultiplex_qc.txt
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}

process STAR_INDEX {
    label 'process_max'
    publishDir "${params.index_dir}", mode: 'move'

    input:
        path fa_file
        path annot_file
    output:
        path "star_index/", emit: index_path
        val true, emit: finished
    when:
        params.create_index == true
    
    """
    mkdir star_index
    STAR \\
        --runThreadN 12 \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles $fa_file \\
        --sjdbGTFfile  $annot_file \\
        --sjdbOverhang $params.map_size
    """
}

process extract_UMI {
    label 'process_medium'
    tag "$fastq"
    maxForks 25

    input:
        path fastq
    output:
        path "*_UMI.*"
    
    """
    umi.py -m extract -I $fastq
    """
}

process STAR_ALIGN {
    label 'process_max'
    tag "$fastq"
    maxForks 25

    input:
        path fastq
        path genome_index
        val index_created
    output:
        path "*.bam", emit: bams
        path "*Log.final.out", emit: logs

    """
    # star_utils.py -t align -I $genome_index -f $fastq
    STAR \\
        --runThreadN $task.cpus \\
        --genomeDir $genome_index \\
        --readFilesIn $fastq \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outReadsUnmapped Fastx \\
        --limitBAMsortRAM 10000000000 \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outFilterMatchNminOverLread 0.33 \\
        --alignEndsType Local \\
        --outSAMattributes Standard 
    """
}

process dedup {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output}/bams", mode: 'copy'
    maxForks 25

    input:
        path bam
    output:
        path "dedup_*"
    
    """
    samtools index $bam
    umi_tools dedup -I $bam -S dedup_$bam
    """
}

process BAM_INDEX {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output}/bams", mode: 'move'
    maxForks 25

    input:
        path bam
    output:
        path "*.bai"

    """
    samtools index $bam
    """
}

process featureCounts {
    label 'process_low'
    publishDir "${params.output}", mode: 'move'
    input:
        path bams
        path csv_dir
    output:
        path "counts.txt"

    """
    featureCounts -a $csv_dir/*.saf -o counts.txt $bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/_R1_UMIAligned.sortedByCoord.out.bam//g' counts.txt
    """
}


// VARIABLE DECLARATIONS
paired_fastq_files = channel.fromFilePairs("$params.fastq_files/*$params.suffix{1,2}*$params.extension", checkIfExists: true)
csv_dir            = channel.fromPath(params.csv_dir, type: 'dir', checkIfExists: true)
index_dir          = channel.fromPath(params.index_dir, type: 'dir', checkIfExists: true)

if (params.create_index) {
    index_fa    = fromPath("$params.ref_gen/*.fa", checkIfExists: true)
    index_annot = fromPath("$params.ref_gen/*.gtf", checkIfExists: true) 
}
else {
    index_fa    = channel.empty()
    index_annot = channel.empty()
}


workflow {
    /* //TODO
        - random bus error or input/output error on any process:
            * as a patch, retry allowed as a strategy (the error occurs randomly and doesnt usually repeat on the same process fork)
            * limited the maxForks in case the error comes from memory, little effect
            * removed the shared memory strategy on STAR genome, little effect
        - ask the user whether to save the created index or not and modify its publishDir 
        - add dry-run for dependency download 
        - add demultiplex log files to multiqc report
    */

    // star index creation
    STAR_INDEX(index_fa, index_annot)

    // raw file qc
    fastqc(paired_fastq_files)
    fastqc_multiqc_inputs = fastqc.out.collect()
    
    // raw fastq demultiplexing
    demultiplex(paired_fastq_files, csv_dir) | collect \
    | seqkit_demultiplex

    // UMI extraction
    seqkit_demultiplex.out.fastqs | flatten \
    | extract_UMI

    // read alignment
    align_trigger = params.create_index ? STAR_INDEX.out.finished : true
    index_path    = params.create_index ? STAR_INDEX.out.index_path : params.index_dir

    STAR_ALIGN(extract_UMI.out, index_path, align_trigger)
    alignment_multiqc_inputs = STAR_ALIGN.out.logs.collect()

    // qc report
    multiqc(
        fastqc_multiqc_inputs,
        alignment_multiqc_inputs
    )

    // deduplication of aligned bams
    dedup(STAR_ALIGN.out.bams)
    
    // reindexing of bams for IGV
    BAM_INDEX(dedup.out)

    // expression quantification
    featureCounts(dedup.out.collect(), csv_dir) 
}
