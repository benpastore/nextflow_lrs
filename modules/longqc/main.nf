process LONGQC_RAW {

    label 'qc'

    tag "${sampleID}_longqc_raw"

    publishDir "${params.results}/longqc_raw", mode : 'copy', pattern : "*"

    input : 
        tuple val(sampleID), val(fastq)
    
    output : 
        path("*")
    
    script : 
    """
    #!/bin/bash

    name=\$(basename ${fastq} .fastq.gz)

    zcat ${fastq} | \\
        sampleqc \\
        -x ont-ligation \\
        -p ${task.cpus} \\
        -o \${PWD}/\${name} \\
        -
    """

}