#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

process multiqc {
    label 'process_low'
    publishDir "${params.output_dir}", mode: 'copy'

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

process demultiplex {
    label 'process_high'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(fastqs)
        each csv_dir
    output:
        path "*.gz", emit: fastqs
        path "demultiplex${sample_id}.log", emit: logs

    """
    demultiplex.py -i $csv_dir -p . -s ${params.suffix}1 ${params.suffix}2 > demultiplex${sample_id}.log
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}

process STAR_INDEX {
    label 'process_max'
    publishDir "${params.index_dir}", mode: 'copy', enabled: params.save_index

    input:
        path fa_file
        path annot_file
    output:
        path "star_index/", emit: index_dir
        val true, emit: finished
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

process extract_UMI {
    label 'process_medium'
    tag "$fastq"

    input:
        path fastq
    output:
        path "*_UMI.*"
    
    """
    umi.py -m extract -I $fastq
    """
}

process STAR_ALIGN {
    label 'process_max'
    tag "$fastq"

    input:
        path fastq
        each genome_index

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
            --outSAMattributes Standard 
        """
}

process dedup {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams", mode: 'copy', enabled: params.get_bams

    input:
        path bam
    output:
        path "dedup_*", emit: dedup_bams
        path "${logname}.log", emit: logs
    
    script:
        logname = bam.simpleName
        """
        samtools index $bam
        umi_tools dedup -I $bam -S dedup_$bam > ${logname}.log
        """
}

process BAM_INDEX {
    label 'process_low'
    tag "$bam"
    publishDir "${params.output_dir}/bams", mode: 'move', enabled: params.get_bams

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

process featureCounts {
    label 'process_low'
    publishDir "${params.output_dir}", mode: 'copy', pattern: "*.txt"

    input:
        path bams
        path csv_dir
    output:
        path "counts.txt", emit: counts
        path "*.summary", emit: logs

    """
    featureCounts -a $csv_dir/*.saf -o counts.txt $bams -F 'SAF'
    sed -i '1d ; 2 s/dedup_//g ; 2 s/_R1_UMIAligned.sortedByCoord.out.bam//g' counts.txt
    """
}


// VARIABLE DECLARATIONS
paired_fastq_files = channel.fromFilePairs("$params.fastq_dir/*$params.suffix{1,2}*$params.extension", checkIfExists: true)
csv_dir            = channel.fromPath(params.csv_dir, type: 'dir', checkIfExists: true)
index_dir          = channel.fromPath(params.index_dir, type: 'dir', checkIfExists: true)

if (params.create_index) {
    index_fa    = channel.fromPath("$params.ref_gen/*.fa", checkIfExists: true)
    index_annot = channel.fromPath("$params.ref_gen/*.gtf", checkIfExists: true) 
}
else {
    index_fa    = channel.empty()
    index_annot = channel.empty()
}


workflow {
    // star index creation
    STAR_INDEX(index_fa, index_annot)

    // raw file qc
    fastqc(paired_fastq_files)
    fastqc_multiqc = fastqc.out.collect()

    // fastq demultiplexing
    demultiplex(paired_fastq_files, csv_dir)
    demultiplex_multiqc = demultiplex.out.logs.collect()
    
    // UMI sequence extraction
    extract_UMI(demultiplex.out.fastqs.flatten())

    // read alignment
    index_dir = params.create_index ? STAR_INDEX.out.index_dir : params.index_dir

    STAR_ALIGN(extract_UMI.out, index_dir)
    alignment_multiqc = STAR_ALIGN.out.logs.collect()

    // deduplication of aligned bams
    dedup(STAR_ALIGN.out.bams)
    dedup_multiqc = dedup.out.logs.collect()
    
    // reindexing of bams for IGV
    BAM_INDEX(dedup.out.dedup_bams)

    // expression quantification
    featureCounts(dedup.out.dedup_bams.collect(), csv_dir) 
    featureCounts_multiqc = featureCounts.out.logs.collect()

    // qc report
    multiqc(
        fastqc_multiqc,
        demultiplex_multiqc,
        alignment_multiqc,
        dedup_multiqc,
        featureCounts_multiqc
    )
}
