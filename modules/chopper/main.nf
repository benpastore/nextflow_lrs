

process CHOPPER {

    label 'chopper'

    tag "${sampleID}_chopper"

    publishDir "${params.results}/chopper", mode : 'copy'


    input : 
        tuple val(sampleID), val(fastq)
    
    output : 
        tuple val(sampleID), path("*.chopper.fastq.gz"), emit : fq
        path("*")
    
    script : 
    """
    #!/bin/bash

    name=\$(basename ${fastq} .fastq.gz)

    zcat ${fastq} | ${params.bin}/chopper-linux \
        --trim-approach trim-by-quality \
        --cutoff 10 \
        --quality 10 \
        --minlength 1000 \
        --headcrop 25 \
        --tailcrop 50 | gzip > \${name}.chopper.fastq.gz

    """

}