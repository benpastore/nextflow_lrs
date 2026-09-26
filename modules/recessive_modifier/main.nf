// Recessive-modifier candidate search over the cohort's consolidated +
// snpEff VCFs (see bin/recessive_modifier_consolidated.py) -- derives
// mild/severe sibling groupings from family.json and requires compound-het
// candidates to be confirmed in trans via each variant's own calling
// caller's phase evidence. Single process for the whole cohort, not
// per-sample, since the mild-vs-severe comparison is inherently cross-sample.

process RECESSIVE_MODIFIER {

    label 'recessive_modifier'
    publishDir "${params.results}/06_variants/recessive_modifier", mode: params.publish_mode

    input:
        path(family_json)
        path(vcf_manifest)
        path(vcfs)
        path(script)

    output:
        path("*.recessive_modifier.tsv"), emit: recessive_modifier_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -family_json ${family_json} \\
        -vcf_manifest ${vcf_manifest} \\
        -sv_merge_dist ${params.recessive_modifier_sv_merge_dist ?: 500} \\
        ${params.gnomad_chm13_vcf ? "-gnomad_chm13_vcf ${params.gnomad_chm13_vcf}" : ""} \\
        -output cohort.recessive_modifier.tsv
    """
}
