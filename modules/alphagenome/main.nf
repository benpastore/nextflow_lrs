process ALPHAGENOME {

    tag "$sampleID"
    label 'alphagenome'
    errorStrategy 'ignore'

    publishDir "${params.results}/alphagenome", mode: params.publish_mode

    input:
        tuple val(sampleID), val(vcf)
        tuple path(ref), path(ref_fai)

    output:
        tuple val(sampleID), path("*.alphagenome.tsv"), emit: alphagenome_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    python3 ${params.bin}/run_alphagenome.py \\
        --sample ${sampleID} \\
        --vcf ${vcf} \\
        --reference ${ref} \\
        --weights ${params.alphagenome_weights} \\
        --output ${sampleID}.alphagenome.tsv \\
        ${params.alphagenome_args ?: ''}
    """
}
