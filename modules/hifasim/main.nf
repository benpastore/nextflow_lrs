process HIFASIM { 

    label 'hifasim'

    publishDir "${params.results}/hifasim", mode: params.publish_mode

    input :
        tuple val(sampleID), val(reads) 

    output : 
        tuple val(sampleID), val(reads), path("*hap1.p_ctg.gfa"), path("*hap2.p_ctg.gfa"), emit : hifasim_asm
        path("*")

    script : 
    """
    #!/bin/bash
    
    name=\$(basename ${reads} .fastq.gz)
    hifiasm \\
        -t64 \\
        --ont \\
        -o \$name \\
        ${reads}
    """
}