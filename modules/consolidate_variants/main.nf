process CONSOLIDATE_VARIANTS {

    tag "$sampleID"
    label 'low'
    errorStrategy 'ignore'

    publishDir "${params.results}/variants/consolidated", mode: params.publish_mode

    input:
        tuple val(sampleID),
            val(clair3_longphase_vcf),
            val(sniffles_vcf),
            val(dipcall_vcf),
            val(hapdiff_vcf)

    output:
        tuple val(sampleID), path("*.master.vcf"), emit: consolidated_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    source activate rnaseq

    python3 ${params.bin}/consolidate_variants.py \\
        --sample ${sampleID} \\
        --caller clair3_longphase:alignment:${clair3_longphase_vcf} \\
        --caller sniffles:alignment:${sniffles_vcf} \\
        --caller dipcall:assembly:${dipcall_vcf} \\
        --caller hapdiff:assembly:${hapdiff_vcf} \\
        --output ${sampleID}.master.vcf
    """
}
