process GFA_CONVERT { 

    label 'gfatools'

    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input : 
        tuple val(sampleID), val(asm)

    output : 
        tuple val(sampleID), path("*.asm.fa"), path("*.asm.fa.fai"), emit : fasta_asm

    script : 
    """
    #!/bin/bash

    gfatools gfa2fa ${asm} > ${sampleID}.asm.fa

    samtools faidx ${sampleID}.asm.fa

    """
}