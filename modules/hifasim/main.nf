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

    # herro-corrected reads are FASTA content under a .fastq.gz name (see
    # HERRO_INFERENCE) -- minimap2 tolerates this via content auto-detect,
    # but hifiasm validates extension against content and refuses
    # ("is in fasta format rather than fastq format"). Symlink to a
    # correctly-named copy rather than touching the shared herro output.
    ln -s ${reads} \${name}.fasta.gz

    hifiasm \\
        -t64 \\
        --ont \\
        -o \$name \\
        \${name}.fasta.gz
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

    # herro-corrected reads are FASTA content under a .fastq.gz name (see
    # HERRO_INFERENCE) -- minimap2 tolerates this via content auto-detect,
    # but hifiasm validates extension against content and refuses
    # ("is in fasta format rather than fastq format"). Symlink to a
    # correctly-named copy rather than touching the shared herro output.
    ln -s ${reads} \${name}.fasta.gz

    hifiasm \\
        -t64 \\
        --ont \\
        -1 ${pat_yak} \\
        -2 ${mat_yak} \\
        -o \$name \\
        \${name}.fasta.gz
    """
}