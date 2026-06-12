
process SNIFFLES {

    label 'sniffles'
    tag "${sampleID}"
    publishDir "${params.results}/variants/sniffles", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.sniffles.vcf.gz"), path("*.sniffles.vcf.gz.tbi"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    sniffles \
      --input ${bam} \
      --vcf \$name.sniffles.vcf \
      --threads ${task.cpus} \
      --allow-overwrite
    
    bgzip -f \$name.sniffles.vcf
    tabix -f -p vcf \$name.sniffles.vcf.gz
    
    """
}

process INDEX_SNIFFLES_VCF {

    label "bcftools"

    input:
        tuple val(sampleID), val(vcf)

    output:
        tuple val(sampleID),
            path("*.sniffles.vcf.gz"),
            path("*.sniffles.vcf.gz.tbi"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    set -euo pipefail

    if [[ "${vcf}" == *.vcf.gz ]]; then
        tabix -f -p vcf ${vcf}
    else
        bgzip -f ${vcf}
        tabix -f -p vcf ${vcf}.gz
    fi
    """
}