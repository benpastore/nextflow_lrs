process GFA_CONVERT { 
    label 'gfatools'


    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input: 
        tuple val(sampleID), val(reads), path(hap1), path(hap2)

    output: 
        tuple val(sampleID), val(reads), path("${sampleID}.haps.combined.fa"), emit: combined_fa

    script:
    """
    #!/bin/bash
    gfatools gfa2fa ${hap1} > ${sampleID}.hap1.fa
    gfatools gfa2fa ${hap2} > ${sampleID}.hap2.fa

    sed 's/^>/>hap1_/' ${sampleID}.hap1.fa > ${sampleID}.haps.combined.fa
    sed 's/^>/>hap2_/' ${sampleID}.hap2.fa >> ${sampleID}.haps.combined.fa
    """
}

process GFA_FAIDX { 
    label 'samtools'

    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input:
        tuple val(sampleID), val(reads), path(combined_fa)

    output:
        tuple val(sampleID), val(reads), path(combined_fa), path("${combined_fa}.fai"), emit: fasta_asm
        tuple val(sampleID), val(reads), path("*hap1.fa"), path("*hap2.fa"), emit: haplotype_fasta_asm


    script:
    """
    #!/bin/bash
    samtools faidx ${combined_fa}

    awk '/^>hap1_/ {p=1} /^>hap2_/ {p=0} p' ${combined_fa} > ${sampleID}.hap1.fa
    awk '/^>hap2_/ {p=1} /^>hap1_/ {p=0} p' ${combined_fa} > ${sampleID}.hap2.fa

    samtools faidx ${sampleID}.hap1.fa
    samtools faidx ${sampleID}.hap2.fa


    """
}