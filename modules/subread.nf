
process featureCounts {
    label 'process_low'
    publishDir "${params.output_dir}",
        mode: 'copy',
        pattern: "*.tsv"

    input:
        path dedup_bams
        path raw_bams
        path csv_dir
    output:
        path "*.tsv", emit: counts
        path "*.summary", emit: logs

    """
    featureCounts -a $csv_dir/*.saf -o dedup_counts.tsv $dedup_bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_UMI_//g ; 2 s/_R1Aligned.sortedByCoord.out.bam//g' dedup_counts.tsv

    featureCounts -a $csv_dir/*.saf -o raw_counts.tsv $raw_bams -F 'SAF'
    sed -i '1d ; 2 s/UMI_//g ; 2 s/_R1Aligned.sortedByCoord.out.bam//g' raw_counts.tsv
    """
}
