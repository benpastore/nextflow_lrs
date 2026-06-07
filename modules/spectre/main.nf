process SPECTRE {

    tag "$sampleID"
    label 'spectre'

    publishDir "${params.results}/spectre", mode: 'copy'

    input:
        tuple val(sampleID), path(bam), path(bai), path(ref_fa), path(ref_fai)
    output:
        tuple val(sampleID), path("*.spectre.*"), emit: spectre_cnv_ch
            //path("*.spectre.vcf.gz.tbi"),
            //path("*.spectre.bed.gz"),
            //path("*.spectre.spc"),
            //emit: spectre_cnv_ch

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
        --sample-id \${name} \
        --output-dir . \
        --reference ${ref_fa} \
        --threads ${task.cpus} \
        ${params.spectre_args ?: ''}

    #SPECTRE_VCF=\$(find . -maxdepth 1 -name "*.vcf.gz" | head -n 1)
    #SPECTRE_BED=\$(find . -maxdepth 1 -name "*.bed" | head -n 1)
    #SPECTRE_SPC=\$(find . -maxdepth 1 -name "*.spc" | head -n 1)

    #mv "\$SPECTRE_VCF" \${name}.spectre.vcf.gz

    #if [ -f "\${SPECTRE_VCF}.tbi" ]; then
    #    mv "\${SPECTRE_VCF}.tbi" \${name}.spectre.vcf.gz.tbi
    #else
    #    tabix -f -p vcf \${name}.spectre.vcf.gz
    #fi

    #bgzip -f -c "\$SPECTRE_BED" > \${name}.spectre.bed.gz
    #tabix -f -0 -s 1 -b 2 -e 3 \${name}.spectre.bed.gz

    #mv "\$SPECTRE_SPC" \${name}.spectre.spc
    """
}