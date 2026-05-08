process HIFASIM { 

    label 'hifasim'

    publishDir "${params.results}/hifasim", mode: params.publish_mode

    input :
        tuple val(sampleID), val(reads) 

    output : 
        tuple val(sampleID), val(reads), path("*.asm"), emit : hifasim_asm

    script : 
    """
    #!/bin/bash

    hifiasm -t64 -o ${sampleID}.asm ${reads}

    """
}