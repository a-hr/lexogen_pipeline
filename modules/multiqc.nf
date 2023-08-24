process multiqc {
    label 'process_low'
    publishDir "${params.out_dir}", mode: 'copy'
    containerOptions ${workflow.containerEngine == "docker" ? '--user $(id -u):$(id -g) --group-add 100' : ''}

    input:
        path fastqc_logs
        path trim_fastqc_logs
        path demultiplex_logs
        path alignment_logs
        path dedup_logs
        path featureCounts_logs

    output:
        path "*.html"
    
    """
    multiqc . -n lexogen_pipeline_multiqc_report.html
    """
}
