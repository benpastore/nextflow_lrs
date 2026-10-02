// Recessive- and dominant-modifier candidate search over the cohort's
// consolidated + snpEff VCFs (see bin/recessive_modifier_consolidated.py,
// bin/dominant_modifier_consolidated.py, and their shared parsing/family-
// grouping code in bin/_modifier_common.py) -- derives mild/severe sibling
// groupings from family.json. The recessive search requires compound-het
// candidates to be confirmed in trans via each variant's own calling
// caller's phase evidence; the dominant search needs only a single
// qualifying variant. Both run here, back to back, against the same
// staged inputs, with distinctly-named outputs. Single process for the
// whole cohort, not per-sample, since the mild-vs-severe comparison is
// inherently cross-sample.

process RECESSIVE_MODIFIER {

    label 'recessive_modifier'
    publishDir "${params.results}/06_variants/recessive_modifier", mode: params.publish_mode

    input:
        path(family_json)
        path(vcf_manifest)
        path(vcfs)
        path(common_script)
        path(recessive_script)
        path(dominant_script)

    output:
        path("*.recessive_modifier.tsv"), emit: recessive_modifier_ch
        path("*.recessive_modifier.gene_summary.tsv"), emit: recessive_modifier_gene_summary_ch
        path("*.dominant_modifier.tsv"), emit: dominant_modifier_ch
        path("*.dominant_modifier.gene_summary.tsv"), emit: dominant_modifier_gene_summary_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${recessive_script} \\
        -family_json ${family_json} \\
        -vcf_manifest ${vcf_manifest} \\
        -sv_merge_dist ${params.recessive_modifier_sv_merge_dist ?: 500} \\
        ${params.gnomad_chm13_vcf ? "-gnomad_chm13_vcf ${params.gnomad_chm13_vcf}" : ""} \\
        -output cohort.recessive_modifier.tsv

    python3 ${dominant_script} \\
        -family_json ${family_json} \\
        -vcf_manifest ${vcf_manifest} \\
        -sv_merge_dist ${params.recessive_modifier_sv_merge_dist ?: 500} \\
        ${params.gnomad_chm13_vcf ? "-gnomad_chm13_vcf ${params.gnomad_chm13_vcf}" : ""} \\
        -output cohort.dominant_modifier.tsv
    """
}
