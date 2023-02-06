#!/usr/bin/env nextflow
nextflow.enable.dsl=2
clean = true
multiqc_inputs = []

multiqc_inputs << "notes.md"
multiqc_inputs << "launch.sh"

process a {
    input:
        val x

    output:
        path "*.txt"
    
    """
    touch a.txt
    touch b.txt
    """
}

process b {
    input:
        val y

    output:
        path "*.txt"
    
    """
    touch c${y}.txt
    touch blabla${y}.txt
    """
}

process c {
    input:
        path infiles
    output:
        stdout
    
    """
    echo *.txt
    """
}

workflow {
    a("hola")
    b(Channel.of([1, 2, 3]).flatten())
    input_channel = a.out.collect().concat(b.out.collect())
    c(input_channel.collect()).view {it}
}