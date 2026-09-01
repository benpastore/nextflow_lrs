process HIFASIM { 

    label 'hifasim'

    publishDir "${params.results}/hifasim", mode: params.publish_mode

    input :
        tuple val(sampleID), val(bam), val(reads) 

    output : 
        tuple val(sampleID), val(bam), val(reads), path("*hap1.p_ctg.gfa"), path("*hap2.p_ctg.gfa"), emit : hifasim_asm
        path("*")

    script :
    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${reads} .fastq.gz)

    # No --ont: herro-corrected reads are FASTA (no quality scores, see
    # HERRO_INFERENCE) under a misleading .fastq.gz name. --ont mode hard-
    # requires real FASTQ and refuses FASTA regardless of extension/naming
    # ("is in fasta format rather than fastq format"); hifiasm's default
    # HiFi-mode path has always accepted plain FASTA, so treat the
    # corrected reads as HiFi-equivalent input instead.
    hifiasm \\
        -t64 \\
        -o \$name \\
        ${reads}
    """
}

process HIFASIM_TRIO {

    tag "$sampleID"
    label 'hifasim'

    publishDir "${params.results}/hifasim_trio", mode: params.publish_mode

    input :
        tuple val(sampleID), val(bam), val(reads), val(pat_yak), val(mat_yak)

    output :
        tuple val(sampleID), val(bam), val(reads), path("*hap1.p_ctg.gfa"), path("*hap2.p_ctg.gfa"), emit : hifasim_trio_asm
        path("*")

    script :
    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${reads} .fastq.gz)

    # No --ont: herro-corrected reads are FASTA (no quality scores, see
    # HERRO_INFERENCE) under a misleading .fastq.gz name. --ont mode hard-
    # requires real FASTQ and refuses FASTA regardless of extension/naming
    # ("is in fasta format rather than fastq format"); hifiasm's default
    # HiFi-mode path has always accepted plain FASTA, so treat the
    # corrected reads as HiFi-equivalent input instead.
    hifiasm \\
        -t64 \\
        -1 ${pat_yak} \\
        -2 ${mat_yak} \\
        -o \$name \\
        ${reads}
    """
}