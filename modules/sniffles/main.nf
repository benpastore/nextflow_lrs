
process SNIFFLES {

    label 'sniffles'
    tag "${sampleID}"
    publishDir "${params.results}/variants/sniffles", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.sniffles.vcf"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    sniffles \
      --input ${bam} \
      --vcf \$name.sniffles.vcf \
      --threads ${task.cpus} \
      --allow-overwrite
    """
}

process INDEX_SNIFFLES_VCF {

    label "bcftools"

    input:
        tuple val(sampleID), val(vcf)

    output:
        tuple val(sampleID),
            path("${sampleID}.sniffles.vcf.gz"),
            path("${sampleID}.sniffles.vcf.gz.tbi"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    bgzip -f ${vcf}
    tabix -f -p vcf ${vcf}.sniffles.vcf.gz
    """
}