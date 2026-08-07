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

    name=\$(basename ${reads} .fastq.gz)
    hifiasm \\
        -t64 \\
        --ont \\
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

    name=\$(basename ${reads} .fastq.gz)
    hifiasm \\
        -t64 \\
        --ont \\
        -1 ${pat_yak} \\
        -2 ${mat_yak} \\
        -o \$name \\
        ${reads}
    """
}