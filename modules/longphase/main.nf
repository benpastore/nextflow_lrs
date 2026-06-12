process LONGPHASE {

    label 'longphase'
    tag "$sampleID"

    publishDir "${params.results}/longphase", mode: 'copy'

    input:
        tuple val(ref_fa), val(ref_fai)
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID),
            path("*.longphase.vcf.gz"),
            path("*.longphase.vcf.gz.tbi"),
            path("*.longphase.haplotagged.bam"),
            path("*.longphase.haplotagged.bam.bai")
        
        tuple val(sampleID), path("*.longphase.haplotagged.bam"), path("*.longphase.haplotagged.bam.bai"), emit : longphase_ch

    script:
    """
    set -euo pipefail

    # LongPhase can read VCF, but keeping bgzipped/indexed output is nicer downstream.
    # Use --indels if the Clair3 VCF contains SNPs + small indels.
    name=\$(basename ${vcf})
    name=\${name%.vcf.gz}
    name=\${name%.vcf}

    longphase phase \\
        -s ${vcf} \\
        -b ${bam} \\
        -r ${ref_fa} \\
        -t ${task.cpus} \\
        -o \$name.longphase \\
        --indels \\
        --ont

    # LongPhase usually writes: <prefix>_SNP.vcf
    PHASED_VCF="\$name.longphase.vcf"

    bgzip -f "\$PHASED_VCF"
    tabix -f -p vcf "\${PHASED_VCF}.gz"

    longphase haplotag \\
        -s "\${PHASED_VCF}.gz" \\
        -b ${bam} \\
        -r ${ref_fa} \\
        -t ${task.cpus} \\
        -o \$name.longphase.haplotagged \\
        --log

    # LongPhase usually writes: <prefix>.bam
    samtools index -@ ${task.cpus} \$name.longphase.haplotagged.bam
    """
}

process LONGPHASE_SV {

    label 'longphase'
    tag "$sampleID"

    publishDir "${params.results}/longphase_sv", mode: 'copy'

    input:
        tuple val(ref_fa), val(ref_fai)
        tuple val(sampleID), val(snp_vcf), val(snp_tbi), val(sv_vcf), val(sv_tbi), val(bam), val(bai)

    output:
        tuple val(sampleID),
            path("*.longphase_sv.vcf.gz"),
            path("*.longphase_sv.vcf.gz.tbi"),
            path("*.longphase_sv.haplotagged.bam"),
            path("*.longphase_sv.haplotagged.bam.bai")
        
        tuple val(sampleID), path("*.longphase_sv.haplotagged.bam"), path("*.longphase_sv.haplotagged.bam.bai"), emit : longphase_sv_ch

    script:
    """
    set -euo pipefail

    # LongPhase can read VCF, but keeping bgzipped/indexed output is nicer downstream.
    # Use --indels if the Clair3 VCF contains SNPs + small indels.
    name=\$(basename ${sv_vcf})
    name=\${name%.vcf.gz}
    name=\${name%.vcf}

    longphase phase \\
        -s ${snp_vcf} \\
        --sv-file ${sv_vcf} \\
        -b ${bam} \\
        -r ${ref_fa} \\
        -t ${task.cpus} \\
        -o \$name.longphase_sv \\
        --indels \\
        --ont

    # LongPhase usually writes: <prefix>_SNP.vcf
    PHASED_VCF="\$name.longphase_sv.vcf"

    bgzip -f "\$PHASED_VCF"
    tabix -f -p vcf "\${PHASED_VCF}.gz"

    longphase haplotag \\
        -s "\${PHASED_VCF}.gz" \\
        -b ${bam} \\
        -r ${ref_fa} \\
        -t ${task.cpus} \\
        -o \$name.longphase_sv.haplotagged \\
        --log

    # LongPhase usually writes: <prefix>.bam
    samtools index -@ ${task.cpus} \$name.longphase_sv.haplotagged.bam
    """
}