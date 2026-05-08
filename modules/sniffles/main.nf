
process SNIFFLES {

    label 'sniffles'
    tag "${sampleID}"
    publishDir "${params.results}/variants/sniffles2", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.sniffles.vcf.gz"), path("*.sniffles.vcf.gz.tbi"), emit : sniffles_ch

    script:
    """
    #!/bin/bash

    name=\$(${bam} .bam)
    sniffles \
      --input ${bam} \
      --vcf \$name.sniffles.vcf \
      --threads ${task.cpus} \
      --allow-overwrite

    bgzip -f \$name.sniffles.vcf
    tabix -p vcf \$name.sniffles.vcf
    """
}