// Mild-vs-severe sibling comparisons of copy-number variants (Spectre) and
// repeat expansions (Straglr), cohort-wide, using family.json's mild/severe
// groupings (same as the modifier search, see load_family_groups in
// bin/_modifier_common.py). Gene association comes from the same GTF
// snpEff is built from (params.cnv_repeat_gtf, default params.snpeff_gtf).
//   CNV_COMPARE     bin/cnv_compare.py    -> 10_cnv_repeat_compare/cnv
//   REPEAT_COMPARE  bin/repeat_compare.py -> 10_cnv_repeat_compare/repeats
// Single process each for the whole cohort (the comparison is inherently
// cross-sample); every sample's caller output is staged next to a
// sampleID<TAB>filename manifest.

process CNV_COMPARE {

    label 'low'
    publishDir "${params.results}/10_cnv_repeat_compare/cnv", mode: params.publish_mode

    input:
        path(family_json)
        path(manifest)
        path(vcfs)
        path(gtf, stageAs: 'gtf/*')
        path(modifier_common)
        path(common_script)
        path(script)

    output:
        path("cohort.cnv_gene_summary.tsv"), emit: gene_summary_ch
        path("cohort.cnv_mild_vs_severe.tsv"), emit: pairs_ch
        path("cohort.cnv_gene_cn.tsv"), emit: gene_cn_ch

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
        -gtf ${gtf} \\
        -min_size ${params.cnv_compare_min_size ?: 0} \\
        -min_gene_overlap ${params.cnv_min_gene_overlap ?: 0} \\
        ${params.cnv_include_noncoding ? '-include_noncoding' : ''} \\
        -exclude_chroms '${params.cnv_repeat_exclude_chroms ?: ''}' \\
        -output_prefix cohort
    """
}

process REPEAT_COMPARE {

    label 'low'
    publishDir "${params.results}/10_cnv_repeat_compare/repeats", mode: params.publish_mode

    input:
        path(family_json)
        path(manifest)
        path(beds)
        path(gtf, stageAs: 'gtf/*')
        path(known_loci_bed, stageAs: 'known/*')
        path(modifier_common)
        path(common_script)
        path(script)

    output:
        path("cohort.repeat_locus_summary.tsv"), emit: locus_summary_ch
        path("cohort.repeat_gene_summary.tsv"), emit: gene_summary_ch
        path("cohort.repeat_pattern_summary.tsv"), emit: pattern_summary_ch
        path("cohort.repeat_mild_vs_severe.tsv"), emit: pairs_ch
        path("cohort.repeat_sample_burden.tsv"), emit: sample_burden_ch
        path("cohort.repeat_loci.tsv"), emit: loci_ch

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
        -gtf ${gtf} \\
        -known_loci_bed ${known_loci_bed} \\
        -merge_dist ${params.repeat_merge_dist ?: 50} \\
        -min_copy_diff ${params.repeat_min_copy_diff} \\
        -min_frac_diff ${params.repeat_min_frac_diff} \\
        -promoter_bp ${params.repeat_promoter_bp ?: 2000} \\
        ${params.repeat_absent_as_missing ? '-absent_as_missing' : ''} \\
        -exclude_chroms '${params.cnv_repeat_exclude_chroms ?: ''}' \\
        -output_prefix cohort
    """
}
