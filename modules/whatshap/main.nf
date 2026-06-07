
process WHATSHAP_PHASE {

    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/variants/whatshap_phased", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID), val(bam), val(bai), path("*whatshapphase.vcf.gz"), path("*whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${vcf} .vcf.gz)

    whatshap phase \\
      --reference ${ref} \\
      --output \$name.whatshapphase.vcf.gz \\
      --ignore-read-groups \\
      ${vcf} \\
      ${bam}

    bcftools index -f -t \$name.whatshapphase.vcf.gz

    #mv phased.vcf.gz \$name.whatshapphase.vcf.gz
    #mv phased.vcf.gz.tbi \$name.whatshapphase.vcf.gz.tbi
    """
}

process WHATSHAP_HAPLOTAG {
    
    label 'whatshap'
    tag "${sampleID}"
    publishDir "${params.results}/alignment/whatshap_haplotagged", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(bam), val(bai), val(phased_vcf), val(phased_tbi)

    output:
        tuple val(sampleID), path("*.happlotagged.bam"), path("*.happlotagged.bam.bai"), emit : whathap_haplotag_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    whatshap haplotag \\
      --reference ${ref} \\
      --ignore-read-groups \\
      --output \$name.happlotagged.bam \\
      --output-haplotag-list \$name.happlotagged.tsv \\
      ${phased_vcf} \\
      ${bam}

    samtools index \$name.happlotagged.bam
    """
}