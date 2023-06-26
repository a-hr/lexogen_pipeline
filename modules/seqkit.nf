process stats {
    tag "demultiplex_stats"

    publishDir "${params.out_dir}/seqkit", mode: 'copy'
    
    input:
        path fastq_files
        val name
    output:
        path "*stats.tsv"
    
    """
    seqkit stats -a -T *fastq.gz > ${name}_demultiplex_stats.tsv
    """
}
