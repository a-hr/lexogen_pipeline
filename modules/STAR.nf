process STAR_INDEX {
    label 'process_max'
    publishDir "${params.index_dir}", mode: 'copy', enabled: params.save_index

    input:
        path fa_file
        path annot_file
    output:
        path "star_index/", emit: index_dir
    when:
        params.create_index == true
    
    """
    mkdir star_index
    STAR \\
        --runThreadN 12 \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles $fa_file \\
        --sjdbGTFfile $annot_file \\
        --sjdbOverhang $params.map_size
    """
}

process STAR_ALIGN {
    label 'process_max'
    tag "$fastq"

    input:
        path fastq
        path genome_index

    output:
        path "*.bam", emit: bams
        path "*Log.final.out", emit: logs

    script:
    def outName = fastq.simpleName
    """
    STAR \\
        --runThreadN $task.cpus \\
        --genomeDir $genome_index \\
        --readFilesIn $fastq \\
        --readFilesCommand zcat \\
        --outFileNamePrefix $outName \\
        --outSAMtype BAM SortedByCoordinate \\
        --outReadsUnmapped Fastx \\
        --limitBAMsortRAM 10000000000 \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outFilterMatchNminOverLread 0.33 \\
        --alignEndsType Local \\
        --outSAMattributes Standard \\
        ${executor.name == 'local' ? '--genomeLoad LoadAndRemove' : ''} 
    """
}