process PORECHOP {

    label 'porechop'

    publishDir "${params.results}/porechop", mode: params.publish_mode

    input:
        tuple val(sampleID), val(fastq)

    output:
        tuple val(sampleID), path("*.fastq.gz"), emit : porechop_ch

    script:
    """
    #!/bin/bash
    
    name=\$(basename ${fastq} .fastq.gz)

    porechop -t ${task.cpus} -i ${fastq} -o \$name.porechop.fastq.gz

    """
}