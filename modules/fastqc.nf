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
