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

    singularity run \\
        docker://benpasto/hifiasm:latest \\
        -t64 \\
        --ont \\
        -o ${sampleID} \\
        ${reads}

    #hifiasm -t64 -o ${sampleID}.asm ${reads}
    #hap1.p_ctg.gfa

    """
}