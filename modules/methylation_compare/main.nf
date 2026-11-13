// Mild vs severe promoter methylation per parental haplotype (trio families,
// --run_methylation). See bin/methylation_compare.py.
//
// Children's modkit HP-partitioned bedMethyl + the longphase_sv VCF those HP
// tags came from, the parents' Clair3 VCFs (to assign each phase block's HP1
// to a parent), and the sibling-sharing segments (shared / not_shared
// haplotype per side). Whole cohort in one process, files staged next to
// two manifests.

process METHYLATION_COMPARE {

    label 'recessive_modifier'
    publishDir "${params.results}/09_methylation_compare", mode: params.publish_mode

    input:
        path(family_json)
        path(manifest)          // child<TAB>hp1<TAB>hp2<TAB>ungrouped<TAB>phased vcf
        path(child_files)
        path(parent_manifest)   // parent<TAB>clair3 vcf
        path(parent_vcfs)
        path(parent_segments)
        path(gtf, stageAs: 'gtf/*')
        path(modifier_common)
        path(cnv_repeat_common)
        path(script)

    output:
        path("cohort.promoter_methylation.tsv"), emit: promoter_methylation_ch
        path("cohort.promoter_methylation.flagged.tsv"), emit: flagged_ch
        path("cohort.methylation_phase_blocks.tsv"), emit: phase_blocks_ch
        path("cohort.methylation_imprinting_check.tsv"), emit: imprinting_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -family_json ${family_json} \\
        -manifest ${manifest} \\
        -parent_manifest ${parent_manifest} \\
        -sib_parent_segments ${parent_segments} \\
        -gtf ${gtf} \\
        -upstream ${params.meth_promoter_upstream} \\
        -downstream ${params.meth_promoter_downstream} \\
        -min_cov ${params.meth_min_cov} \\
        -min_cpgs ${params.meth_min_cpgs} \\
        -min_delta ${params.meth_min_delta} \\
        -min_block_votes ${params.meth_min_block_votes} \\
        -exclude_chroms '${params.cnv_repeat_exclude_chroms ?: 'chrM'}' \\
        -output_prefix cohort
    """
}
