process GFA_CONVERT { 

    label 'gfatools'

    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input : 
        tuple val(sampleID), val(hap1), val(hap2)

    output : 
        tuple val(sampleID), path("*.hap1.fa"), path("*.hap2.fa"), path("*.hap1.fa.fai"), path("*.hap2.fa.fai"), emit : fasta_asm

    script : 
    """
    #!/bin/bash

    gfatools gfa2fa ${hap1} > ${sampleID}.hap1.fa
    gfatools gfa2fa ${hap2} > ${sampleID}.hap2.fa

    samtools faidx ${sampleID}.hap1.fa
    samtools faidx ${sampleID}.hap2.fa
    """
}