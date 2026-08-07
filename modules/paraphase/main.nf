process PARAPHASE {

    tag "$sampleID"
    label 'paraphase'

    publishDir "${params.results}/SMA_ASSEMBLY", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref_fa), val(ref_fai)

    output:
        tuple val(sampleID),
            path("*.paraphase.json"),
            path("*.paraphase.bam"),
            path("*.paraphase.bam.bai"),
            path("${sampleID}_paraphase_vcfs"),
            emit: paraphase_ch

    script:
    def gene_arg = params.paraphase_genes ? "-g ${params.paraphase_genes}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    paraphase \\
        -b ${bam} \\
        -r ${ref_fa} \\
        -o . \\
        -p ${sampleID} \\
        --genome ${params.paraphase_genome_build ?: 'chm13'} \\
        -t ${task.cpus} \\
        ${gene_arg} \\
        ${params.paraphase_args ?: ''}
    """
}
