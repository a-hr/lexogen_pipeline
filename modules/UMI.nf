
process extract_UMI {
    label 'process_medium'
    tag "$fastq"

    input:
        path fastq
    output:
        path "UMI_$fastq"
    
    """
    umi_tools extract --bc-pattern=NNNNNN -I $fastq -S UMI_${fastq}
    """
}

process dedup_UMI {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams",
        mode: 'copy',
        pattern: "dedup_$bam",
        enabled: params.publish_mapped

    input:
        path bam
    output:
        path "dedup_$bam", emit: dedup_bams
        path "${bam.simpleName}.log", emit: logs
    
    """
    samtools index $bam
    umi_tools dedup -I $bam -S dedup_$bam > ${bam.simpleName}.log
    """
}