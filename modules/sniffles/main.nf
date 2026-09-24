
process SNIFFLES {

    label 'sniffles'
    tag "${sampleID}"

    input:
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.sniffles.vcf"), emit : sniffles_raw_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    # --phase uses the HP/PS tags already on this BAM (haplotagged by
    # longphase upstream) to emit a per-SV PHASE= field in INFO.
    #
    # bgzip/tabix aren't in this container (the biocontainers sniffles
    # image only ships sniffles itself, no htslib CLI tools) -- compression
    # + indexing is done by INDEX_SNIFFLES_VCF (bcftools container) instead.
    sniffles \
      --input ${bam} \
      --vcf \$name.sniffles.vcf \
      --threads ${task.cpus} \
      --phase \
      --allow-overwrite
    """
}

process INDEX_SNIFFLES_VCF {

    tag "${sampleID}"
    label "bcftools"
    publishDir "${params.results}/06_variants/sniffles", mode: params.publish_mode

    input:
        tuple val(sampleID), path(vcf)

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