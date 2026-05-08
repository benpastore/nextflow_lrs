
process WHATSHAP_PHASE {

    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/variants/phased", mode: params.publish_mode

    input:
        val ref 
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID), val(bam), val(bai), path("*whatshapphase.vcf.gz"), path("whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_ch

    script:
    """
    #!/bin/bash

    name=\$(${vcf} .vcf.gz)

    whatshap phase \
      --reference ${ref} \
      --output phased.vcf.gz \
      ${vcf} \
      ${bam}

    tabix -p vcf phased.vcf.gz

    mv phased.vcf.gz \$name.whatshapphase.vcf.gz
    mv phased.vcf.gz.tbi \$name.whatshapphase.vcf.gz.tbi
    """
}

process WHATSHAP_HAPLOTAG {
    
    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/alignment/haplotagged", mode: params.publish_mode

    input:
        val ref
        tuple val(sampleID), val(bam), val(bai), val(phased_vcf), val(phased_tbi)

    output:
        tuple val(sampleID), path("*.happlotagged.bam"), path("*.happlotagged.bam.bai"), emit : whathap_haplotag_ch

    script:
    """
    #!/bin/bash

    name=\$(${bam} .bam)

    whatshap haplotag \
      --reference ${ref} \
      --output \$name.happlotagged.bam \
      --output-haplotag-list \$name.happlotaggs.tsv \
      ${phased_vcf} \
      ${bam}

    samtools index \$name.happlotagged.bam
    """
}