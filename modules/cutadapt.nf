process demultiplex {
    label 'process_high'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
        path csv_dir
    output:
        path "*.gz", emit: fastqs
        path "*.cutadapt.json", emit: logs

    """
    demultiplex_PE.py -b $csv_dir -f . -s _R1 _R2
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}

process adapter_trimming {
    label 'process_high'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
        val fw_adapter
        val rv_adapter
    output:
        tuple val(sample_id), path("*.gz"), emit: fastqs
        path "trim_${sample_id}.cutadapt.json", emit: logs

    """
    cutadapt \\
        -g $fw_adapter -G $rv_adapter \\
        -e 0.1 \\
        -j 8 \\
        -o trim_${sample_id}_R1.fastq.gz -p trim_${sample_id}_R2.fastq.gz \\
        --json=trim_${sample_id}.cutadapt.json \\
        ${fastqs}
    """
}
