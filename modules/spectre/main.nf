process SPECTRE {

    tag "$sampleID"
    label 'spectre'

    publishDir "${params.results}/06_variants/spectre", mode: 'copy'

    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref_fa), val(ref_fai)

    output:
        tuple val(sampleID), path("*.spectre.vcf.gz"), path("*.spectre.vcf.gz.tbi"), emit: spectre_vcf_ch
        tuple val(sampleID), path("*.spectre.bed.gz"), path("*.spectre.bed.gz.tbi"), emit: spectre_bed_ch
        tuple val(sampleID), path("*.spectre.spc"), emit: spectre_spc_ch
        path("*mosdepth.regions.bed.gz")

    script:
    """
    set -euo pipefail

    name=\$(basename ${bam} .bam)

    mosdepth \
        -t ${task.cpus} \
        -x \
        -b ${params.spectre_bin_size ?: 1000} \
        -Q ${params.spectre_mapq ?: 20} \
        \${name}.mosdepth \
        ${bam}

    tabix -f -0 -s 1 -b 2 -e 3 \${name}.mosdepth.regions.bed.gz

    spectre CNVCaller \
        --coverage \${name}.mosdepth.regions.bed.gz \
        --sample-id \${name}.spectre \
        --output-dir . \
        --reference ${ref_fa} \
        --threads ${task.cpus} \
        ${params.spectre_args ?: ''}

    # spectre names outputs after --sample-id, but doesn't guarantee the
    # extension/compression state (may emit .vcf or .vcf.gz, plain .bed,
    # and may or may not index) -- find whatever it actually wrote and
    # normalize to *.spectre.vcf.gz(+.tbi)/*.spectre.bed.gz(+.tbi)/*.spectre.spc
    # rather than assuming a fixed shape.
    # strip the "./" find prefixes off -- otherwise "./foo" vs "foo" never
    # compares equal below, and mv refuses a same-file no-op rename (exit 1)
    SPECTRE_VCF=\$(find . -maxdepth 1 \\( -name "*.vcf" -o -name "*.vcf.gz" \\) | head -n 1)
    SPECTRE_VCF=\${SPECTRE_VCF#./}
    SPECTRE_BED=\$(find . -maxdepth 1 \\( -name "*.bed" -o -name "*.bed.gz" \\) -not -name "*mosdepth*" | head -n 1)
    SPECTRE_BED=\${SPECTRE_BED#./}
    SPECTRE_SPC=\$(find . -maxdepth 1 -name "*.spc" | head -n 1)
    SPECTRE_SPC=\${SPECTRE_SPC#./}

    if [[ "\$SPECTRE_VCF" == *.gz ]]; then
        [ "\$SPECTRE_VCF" = "\${name}.spectre.vcf.gz" ] || mv "\$SPECTRE_VCF" \${name}.spectre.vcf.gz
    else
        bgzip -f -c "\$SPECTRE_VCF" > \${name}.spectre.vcf.gz
    fi
    tabix -f -p vcf \${name}.spectre.vcf.gz

    if [[ "\$SPECTRE_BED" == *.gz ]]; then
        [ "\$SPECTRE_BED" = "\${name}.spectre.bed.gz" ] || mv "\$SPECTRE_BED" \${name}.spectre.bed.gz
    else
        bgzip -f -c "\$SPECTRE_BED" > \${name}.spectre.bed.gz
    fi
    tabix -f -0 -s 1 -b 2 -e 3 \${name}.spectre.bed.gz

    [ "\$SPECTRE_SPC" = "\${name}.spectre.spc" ] || mv "\$SPECTRE_SPC" \${name}.spectre.spc
    """
}