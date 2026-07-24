process GFA_CONVERT { 
    label 'gfatools'
    errorStrategy 'ignore'


    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input: 
        tuple val(sampleID), val(bam), val(reads), path(hap1), path(hap2)

    output: 
        tuple val(sampleID), val(bam), val(reads), path("*.haps.combined.fa"), emit: combined_fa

    script:
    """
    #!/bin/bash
    name=\$(basename ${reads} .fastq.gz)

    gfatools gfa2fa ${hap1} > \${name}.hap1.fa
    gfatools gfa2fa ${hap2} > \${name}.hap2.fa

    sed 's/^>/>hap1_/' \${name}.hap1.fa > \${name}.haps.combined.fa
    sed 's/^>/>hap2_/' \${name}.hap2.fa >> \${name}.haps.combined.fa
    """
}

process GFA_FAIDX { 

    label 'samtools'
    errorStrategy 'ignore'

    publishDir "${params.results}/gfatools", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), val(reads), path(combined_fa)

    output:
        tuple val(sampleID), val(bam), val(reads), path(combined_fa), path("${combined_fa}.fai"), emit: fasta_asm
        tuple val(sampleID), val(bam), val(reads), path("*hap1.fa"), path("*hap2.fa"), emit: haplotype_fasta_asm


    script:
    """
    #!/bin/bash
    samtools faidx ${combined_fa}

    name=\$(basename ${combined_fa} .fa)
    
    awk '/^>hap1_/ {p=1} /^>hap2_/ {p=0} p' ${combined_fa} > \$name.hap1.fa
    awk '/^>hap2_/ {p=1} /^>hap1_/ {p=0} p' ${combined_fa} > \$name.hap2.fa

    samtools faidx \$name.hap1.fa
    samtools faidx \$name.hap2.fa
    """
}