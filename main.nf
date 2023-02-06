#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastqc {
    label 'process_low'
    input:
        tuple val(sample_id), path(fastqs)
    output:
        path "*.zip"

    """
    fastqc $params.input/$sample_id* -o .
    """
}

process multiqc {
    label 'process_low'
    publishDir "${params.output}", mode: 'copy'

    input:
        path report_infiles
        val outname

    output:
        path "*.html"
    
    """
    multiqc . -n $outname
    """
}

process demultiplex {
    // todo: pass fromfilepairs.collect() as argument and replace -p inpuuts with "."
    label 'process_high'
    output:
        path "*.gz"

    """
    demultiplex.py -i $params.input -p $params.input -s ${params.suffix}1 ${params.suffix}2
    touch $params.output/demultiplex_qc.txt
    seqkit stat *R1.fastq.gz > $params.output/demultiplex_qc.txt
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}

process STAR_INDEX {
    label 'process_max'
    publishDir "${params.index}", mode: 'move'
    output:
        path "*" 
    when:
        params.create_index == true
    
    """
    mkdir -p index
    star_utils.py -t index -I . -g ${params.ref_gen} -L ${params.map_size}
    """
}

process extract_UMI {
    label 'process_medium'
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
    input:
        path fastq
        val index_created // to wait for the index to be created if needed
    output:
        path "*.bam", emit: bams
        path "*[!.bam][!.gz]", emit: logs

    """
    star_utils.py -t align -I $params.index -f $fastq
    """
}

process dedup {
    label 'process_low'
    publishDir "${params.output}/bams", mode: 'copy'
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
    publishDir "${params.output}/bams", mode: 'move'
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
    output:
        path "counts.txt"

    """
    featureCounts -a $params.input/*.saf -o counts.txt $bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/_R1_UMIAligned.sortedByCoord.out.bam//g' counts.txt
    """
}


workflow {
    /* //TODO
        - add dry-run for dependency download 
    */

    // star index creation
    STAR_INDEX()

    // raw file qc
    channel.fromFilePairs("$params.input/*$params.suffix{1,2}*$params.extension") \
    | fastqc

    multiqc_inputs = fastqc.out.collect() // holds data to feed multiqc
    
    // raw file demultiplexing and aligning of resulting .fastqs
    demultiplex | flatten \
    | extract_UMI
    STAR_ALIGN(extract_UMI.out, params.create_index ? STAR_INDEX.out : true)

    // qc report
    multiqc_inputs.concat(STAR_ALIGN.out.logs.collect())
    multiqc(multiqc_inputs, "qc_report")

    // deduplication of aligned bams
    dedup(STAR_ALIGN.out.bams)
    
    // reindexing of bams to be viewable through IGV
    BAM_INDEX(dedup.out)

    // expression quantification
    featureCounts(dedup.out.collect()) 
}
