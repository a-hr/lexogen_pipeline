process group_results {
    publishDir "${params.output_dir}/counts", mode: 'copy'

    input:
        path result_tables
    output:
        path "counts/*.tsv"
    when:
        params.enable_grouping

    script:
    def mouses = params.mouses.collect { "-m ${it}" }.join(" ")
    def tissues = params.tissues.collect { "-t ${it}" }.join(" ")
    """
    group_results.py ${result_tables} -o counts ${mouses} ${tissues}
    """
}

process plot_results {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    input:
        path result_tables
    output:
        path "plots/*"
    when:
        params.enable_plotting

    script:
    def libs = params.libraries.collect { "-l ${it}" }.join(" ")
    def samps = params.samples.collect { "-s ${it}" }.join(" ")
    def tars = params.targets.collect { "-t ${it}" }.join(" ")
    def ctrls = params.controls.collect { "-c ${it}" }.join(" ")
    """
    plot_results.py ${result_tables} -o plots -p ${params.groups} ${libs} ${samps} ${tars} ${ctrls}
    """
}