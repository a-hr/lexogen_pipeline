process BAM_INDEX {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams", mode: 'copy', enabled: params.get_bams

    input:
        path bam
    output:
        path "*.bai"
    when:
        params.get_bams == true

    """
    samtools index $bam
    """
}
