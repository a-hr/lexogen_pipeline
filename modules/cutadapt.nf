process demultiplex {
    label 'process_high'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
        path csv_dir
    output:
        path "*.gz", emit: fastqs
        path "demultiplex${sample_id}.log", emit: logs

    """
    demultiplex_PE.py -i $csv_dir -p . -s ${params.suffix}1 ${params.suffix}2 > demultiplex${sample_id}.log
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}
