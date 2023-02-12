#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process pullContainers {
    publishDir "${baseDir}/containers", mode: 'move'

    output:
        path "containers"
    
    exec:
    // parse the config file and get the container names
    def config = new File("nextflow.config")
    def containers = []
    config.eachLine { line ->
        if (line.contains("container =")) {
            containers.add(line.split("=")[1].trim())
        }
    }

    // remove the duplicates
    containers = containers.unique()

    // create the containers directory
    def containersDir = new File("${baseDir}containers")
    containersDir.mkdirs()

    // pull the containers
    containers.each { container ->
        def containerName = container.replace("/", "_")
        def command = "singularity pull --name " + $containerName + ".img --dir " + $baseDir + "/containers docker://" + $container
        println command
        def process = command.execute()
        process.waitFor()
    }
}


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
        path ".command.log", emit: logs

    """
    demultiplex.py -i $csv_dir -p . -s ${params.suffix}1 ${params.suffix}2
    rm -f *R2.fastq.gz *unknown_R1.fastq.gz
    """
}

process STAR_INDEX {
    label 'process_max'
    publishDir "${params.index_dir}", mode: 'move'

    input:
        path fa_file
        path annot_file
    output:
        path "star_index/", emit: index_path
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
        --sjdbGTFfile  $annot_file \\
        --sjdbOverhang $params.map_size
    """
}

process extract_UMI {
    label 'process_medium'
    tag "$fastq"
    maxForks 25

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
    maxForks 25

    input:
        path fastq
        path genome_index
        val index_created
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
    maxForks 25

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
    maxForks 25

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
    publishDir "${params.output_dir}", mode: 'move'
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
    index_fa    = fromPath("$params.ref_gen/*.fa", checkIfExists: true)
    index_annot = fromPath("$params.ref_gen/*.gtf", checkIfExists: true) 
}
else {
    index_fa    = channel.empty()
    index_annot = channel.empty()
}


workflow {
    /* BUS ERROR
    - random bus error or input/output error on any process:
            * as a patch, 'retry' allowed as a strategy (the error occurs randomly and doesnt usually repeat on the same process fork)
            * limited the maxForks in case the error comes from RAM memory, little effect
            * removed the shared memory strategy on STAR genome, little effect
            * might be related to disk space, according to different sources
    */

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
    align_trigger = params.create_index ? STAR_INDEX.out.finished : true
    index_path    = params.create_index ? STAR_INDEX.out.index_path : params.index_dir

    STAR_ALIGN(extract_UMI.out, index_path, align_trigger)
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
