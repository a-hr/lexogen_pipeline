process multiqc {
    label 'process_low'
    publishDir "${params.output_dir}", mode: 'copy'
    containerOptions '--user $(id -u):$(id -g) --group-add 100'

    input:
        path fastqc_logs
        path demultiplex_logs
        path alignment_logs
        path dedup_logs
        path featureCounts_logs

    output:
        path "*.html"
    
    """
    multiqc . -n Lexogen_pipeline_multiqc_report.html
    """
}
