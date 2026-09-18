process SNPEFF_BUILD {

    tag "${params.snpeff_db}"
    label 'snpeff'

    // published under params.snpeff_data_dir so the next run finds
    // snpEffectPredictor.bin already there and skips rebuilding (see
    // the snpeff_db_exists check in the snpeff subworkflow, subworkflows/ont.nf)
    publishDir "${params.snpeff_data_dir}", mode: params.publish_mode

    input:
        path(genome)
        path(gtf)

    output:
        path("${params.snpeff_db}"), emit: snpeff_db_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    mkdir -p ${params.snpeff_db}
    cp ${genome} ${params.snpeff_db}/sequences.fa
    cp ${gtf} ${params.snpeff_db}/genes.gtf

    snpEff build \\
        -gtf22 \\
        -dataDir \$(pwd) \\
        -configOption ${params.snpeff_db}.genome=${params.snpeff_db} \\
        -v ${params.snpeff_db} \\
        -noCheckCds \\
        -noCheckProtein
    """
}

process SNPEFF {

    tag "$sampleID"
    label 'snpeff'
    errorStrategy 'ignore'

    publishDir "${params.results}/variants/snpeff", mode: params.publish_mode

    input:
        tuple val(sampleID), val(vcf)
        path(snpeff_db_dir)

    output:
        tuple val(sampleID), path("${sampleID}.snpeff.vcf"), emit: snpeff_vcf_ch
        tuple val(sampleID), path("${sampleID}.snpeff.csv"), path("${sampleID}.snpeff.html"), emit: snpeff_stats_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    snpEff \\
        -dataDir \$(pwd) \\
        -configOption ${params.snpeff_db}.genome=${params.snpeff_db} \\
        -csvStats ${sampleID}.snpeff.csv \\
        -htmlStats ${sampleID}.snpeff.html \\
        ${params.snpeff_args ?: ''} \\
        ${params.snpeff_db} \\
        ${vcf} \\
        > ${sampleID}.snpeff.vcf
    """
}
