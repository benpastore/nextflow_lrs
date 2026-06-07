process LONGPHASE {

    label 'longphase'
    tag "$sampleID"

    cpus 8
    memory '32 GB'
    time '12h'

    publishDir "${params.outdir}/longphase", mode: 'copy'

    input:
        tuple val(ref_fa), val(ref_fai)
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID),
            path("*.longphase.phased.vcf.gz"),
            path("*.longphase.phased.vcf.gz.tbi"),
            path("*.longphase.haplotagged.bam"),
            path("*.longphase.haplotagged.bam.bai")
        
        tuple val(sampleID), path("*.happlotagged.bam"), path("*.happlotagged.bam.bai"), emit : longphase_ch


    script:
    """
    set -euo pipefail

    # LongPhase can read VCF, but keeping bgzipped/indexed output is nicer downstream.
    # Use --indels if the Clair3 VCF contains SNPs + small indels.
    name=\$(basename ${vcf} .vcf)

    longphase phase \
        -s ${vcf} \
        -b ${bam} \
        -r ${ref_fa} \
        -t ${task.cpus} \
        -o \$name.longphase \
        --indels \
        --ont

    # LongPhase usually writes: <prefix>_SNP.vcf
    PHASED_VCF="\$name.longphase_SNP.vcf"

    bgzip -f "\$PHASED_VCF"
    tabix -f -p vcf "\$name.longphase_SNP.vcf.gz"

    longphase haplotag \
        -s "\$name.longphase_SNP.vcf.gz" \
        -b ${bam} \
        -r ${ref_fa} \
        -t ${task.cpus} \
        -o \$name.longphase.haplotagged \
        --log

    # LongPhase usually writes: <prefix>.bam
    mv \$name.longphase_SNP.vcf.gz \$name.longphase.phased.vcf.gz
    mv \$name.longphase_SNP.vcf.gz.tbi \$name.longphase.phased.vcf.gz.tbi

    if [[ -f \$name.longphase.haplotagged.bam ]]; then
        :
    elif [[ -f \$name.longphase.haplotagged_HAPLOTAG.bam ]]; then
        mv \$name.longphase.haplotagged_HAPLOTAG.bam \$name.longphase.haplotagged.bam
    else
        echo "ERROR: Could not find LongPhase haplotagged BAM output" >&2
        ls -lh >&2
        exit 1
    fi

    samtools index -@ ${task.cpus} \$name.longphase.haplotagged.bam
    """
}