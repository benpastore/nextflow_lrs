

process NANOPLOT_RAW {

    label 'nanoplot'

    tag "${sampleID}_nanoplot_raw"

    publishDir "${params.results}/nanoplot_raw", mode : 'copy', pattern : "*"

    input : 
        tuple val(sampleID), val(bam), val(fastq)
    
    output : 
        path("*")
    
    script : 
    """
    #!/bin/bash

    name=\$(basename ${fastq} .fastq.gz)

    NanoPlot \\
        -t ${task.cpus} \\
        -o \${PWD}/\${name} \\
        --fastq ${fastq} \\
        --loglength \\
        --N50 \\
        --tsv_stats

    """
}