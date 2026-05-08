process SPECTRE {

    tag "$sampleID"
    label 'spectre'

    publishDir "${params.results}/spectre", mode: 'copy'

    input:
        tuple val(sampleID), path(bam), path(bai), path(ref_fa), path(ref_fai), emit : spectre_ch

    output:
        tuple val(sampleID),
            path("${sampleID}.spectre.vcf.gz"),
            path("${sampleID}.spectre.vcf.gz.tbi"),
            path("${sampleID}.spectre.bed.gz"),
            path("${sampleID}.spectre.spc"),
            emit: spectre_cnv_ch

    script:
    """
    set -euo pipefail

    mosdepth \\
        -t ${task.cpus} \\
        -x \\
        -b ${params.spectre_bin_size ?: 1000} \\
        -Q ${params.spectre_mapq ?: 20} \\
        ${sampleID}.mosdepth \\
        ${bam}

    spectre CNVCaller \\
        --coverage ${sampleID}.mosdepth.regions.bed.gz \\
        --sample-id ${sampleID} \\
        --output-dir . \\
        --reference ${ref_fa} \\
        --threads ${task.cpus} \\
        ${params.spectre_args ?: ''}

    # Normalize expected filenames
    SPECTRE_VCF=\$(ls *.vcf.gz | head -n 1)
    SPECTRE_BED=\$(ls *.bed | head -n 1)
    SPECTRE_SPC=\$(ls *.spc | head -n 1)

    mv "\$SPECTRE_VCF" ${sampleID}.spectre.vcf.gz

    if [ -f "\${SPECTRE_VCF}.tbi" ]; then
        mv "\${SPECTRE_VCF}.tbi" ${sampleID}.spectre.vcf.gz.tbi
    else
        tabix -f -p vcf ${sampleID}.spectre.vcf.gz
    fi

    bgzip -f -c "\$SPECTRE_BED" > ${sampleID}.spectre.bed.gz
    mv "\$SPECTRE_SPC" ${sampleID}.spectre.spc
    """
}