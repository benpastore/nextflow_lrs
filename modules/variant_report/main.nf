// Final per-variant report for the recessive/dominant-modifier PASS
// candidates -- joins modifier-search info, cohort carriers (mild/severe),
// gnomAD AF, snpEff, AlphaMissense and AlphaGenome into one table. See
// bin/build_variant_report.py. The hg38 VCF / AlphaGenome inputs are
// optional: pass file("NO_FILE") and their columns are left empty (each is
// staged into its own subdir so multiple NO_FILE placeholders don't collide).

process VARIANT_REPORT {

    label 'recessive_modifier'
    publishDir "${params.results}/09_variant_report", mode: params.publish_mode

    input:
        path(recessive_tsv)
        path(dominant_tsv)
        path(family_json)
        path(vcf_manifest)
        path(vcfs)
        path(hg38_vcf, stageAs: 'hg38/*')
        path(alphagenome_tsv, stageAs: 'alphagenome_local/*')
        path(alphagenome_api_tsv, stageAs: 'alphagenome_api/*')
        path(common_script)
        path(script)

    output:
        path("cohort.variant_report.tsv"), emit: variant_report_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -recessive_tsv ${recessive_tsv} \\
        -dominant_tsv ${dominant_tsv} \\
        -family_json ${family_json} \\
        -vcf_manifest ${vcf_manifest} \\
        -hg38_vcf ${hg38_vcf} \\
        -alphagenome_tsv ${alphagenome_tsv} \\
        -alphagenome_api_tsv ${alphagenome_api_tsv} \\
        ${params.alphamissense_tsv ? "-alphamissense_tsv ${params.alphamissense_tsv}" : ''} \\
        -output cohort.variant_report.tsv
    """
}
