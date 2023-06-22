process fastqc {
    label 'process_medium'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
    output:
        path "*.zip"

    """
    fastqc $fastqs -o .
    """
}
