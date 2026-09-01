
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

process WHATSHAP_PHASE_TRIO {

    label 'whatshap'
    tag "${child}"
    publishDir "${params.results}/variants/whatshap_trio_phased", mode: params.publish_mode

    // Additive, trio-aware refinement of the child's small-variant phasing,
    // run alongside (not instead of) the existing longphase-based
    // WHATSHAP_PHASE/longphase_sv path. Parents contribute genotype-only
    // Mendelian constraints from their own Illumina VCF -- no parental BAM
    // is required for whatshap's pedigree mode to improve child phasing.
    input:
        tuple val(ref), val(ref_fai)
        tuple val(child), val(father), val(mother), val(child_vcf), val(child_bam), val(child_bai), val(father_vcf), val(mother_vcf)

    output:
        tuple val(child), path("*.trio.whatshapphase.vcf.gz"), path("*.trio.whatshapphase.vcf.gz.tbi"), emit : whatshap_phase_trio_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    printf "%s\\t%s\\t%s\\t%s\\t0\\t0\\n" "${child}_family" "${child}" "${father}" "${mother}" > ${child}.ped
    printf "%s\\t%s\\t0\\t0\\t1\\t0\\n" "${child}_family" "${father}" >> ${child}.ped
    printf "%s\\t%s\\t0\\t0\\t2\\t0\\n" "${child}_family" "${mother}" >> ${child}.ped

    # No --ignore-read-groups: whatshap rejects it together with --ped
    # (pedigree mode identifies samples by name, which is the opposite of
    # what --ignore-read-groups is for). Requires child_vcf/father_vcf/
    # mother_vcf's sample columns to actually match ${child}/${father}/
    # ${mother} -- see CLAIR3's --sample_name.
    whatshap phase \\
      --ped ${child}.ped \\
      --reference ${ref} \\
      --output ${child}.trio.whatshapphase.vcf.gz \\
      ${child_vcf} \\
      ${father_vcf} \\
      ${mother_vcf} \\
      ${child_bam}

    bcftools index -f -t ${child}.trio.whatshapphase.vcf.gz
    """
}