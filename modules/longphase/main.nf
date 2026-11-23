process LONGPHASE {

    label 'longphase'
    tag "$sampleID"

    publishDir "${params.results}/05_alignment/longphase", mode: 'copy'

    input:
        tuple val(ref_fa), val(ref_fai)
        tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)

    output:
        tuple val(sampleID),
            path("*.longphase.vcf.gz"),
            path("*.longphase.vcf.gz.tbi"),
            path("*.longphase.haplotagged.bam"),
            path("*.longphase.haplotagged.bam.bai"), emit : longphase_full_ch

        tuple val(sampleID), path("*.longphase.haplotagged.bam"), path("*.longphase.haplotagged.bam.bai"), emit : longphase_ch

        tuple val(sampleID), path("*.longphase.vcf.gz"), path("*.longphase.vcf.gz.tbi"), emit : longphase_vcf_ch

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

    publishDir "${params.results}/05_alignment/longphase_sv", mode: 'copy'

    input:
        tuple val(ref_fa), val(ref_fai)
        tuple val(sampleID), val(snp_vcf), val(snp_tbi), val(sv_vcf), val(sv_tbi), val(bam), val(bai)

    output:
        tuple val(sampleID),
            path("*.longphase_sv.vcf.gz"),
            path("*.longphase_sv.vcf.gz.tbi"),
            path("*.longphase_sv.haplotagged.bam"),
            path("*.longphase_sv.haplotagged.bam.bai"), emit : longphase_sv_full_ch

        tuple val(sampleID), path("*.longphase_sv.haplotagged.bam"), path("*.longphase_sv.haplotagged.bam.bai"), emit : longphase_sv_ch

        tuple val(sampleID), path("*.longphase_sv.vcf.gz"), path("*.longphase_sv.vcf.gz.tbi"), emit : longphase_sv_vcf_ch

        // longphase's phased SVs (same phase blocks as the SNV VCF + BAM above);
        // the only SV VCF downstream uses (SV_TRIO_ORIGIN, consolidation).
        // Named *.sniffles.vcf.gz so consolidate_variants.py still sees it as Sniffles.
        tuple val(sampleID), path("*.longphase_sv.SV.sniffles.vcf.gz"), path("*.longphase_sv.SV.sniffles.vcf.gz.tbi"), emit : longphase_sv_sv_vcf_ch

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

    # LongPhase writes the phased SNVs/indels to <prefix>.vcf and, with
    # --sv-file, the phased SVs to <prefix>_SV.vcf
    PHASED_VCF="\$name.longphase_sv.vcf"
    SV_VCF="\$name.longphase_sv_SV.vcf"
    if [ ! -s "\$SV_VCF" ]; then
        SV_VCF=\$(ls *_SV.vcf 2>/dev/null | head -n 1 || true)
    fi
    if [ -z "\$SV_VCF" ] || [ ! -s "\$SV_VCF" ]; then
        echo "[LONGPHASE_SV] longphase wrote no phased SV VCF (<prefix>_SV.vcf); files here:" >&2
        ls -l >&2
        exit 1
    fi
    n_in=\$(zcat -f ${sv_vcf} | grep -vc '^#' || true)
    n_out=\$(grep -vc '^#' "\$SV_VCF" || true)
    echo "[LONGPHASE_SV] \$SV_VCF: \$n_out SV records (Sniffles input: \$n_in)" >&2

    bgzip -f "\$PHASED_VCF"
    tabix -f -p vcf "\${PHASED_VCF}.gz"
    bgzip -c "\$SV_VCF" > \$name.longphase_sv.SV.sniffles.vcf.gz
    tabix -f -p vcf \$name.longphase_sv.SV.sniffles.vcf.gz

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