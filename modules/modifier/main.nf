// Modifier-gene search (variant-based), see bin/modifier_variants.py and
// bin/gene_recurrence.py. Plain genotype rules -- no statistical model and no
// allele-frequency filter.
//
// MODIFIER_VARIANTS  one task for the whole cohort: every discordant
//                    family's mild child vs its severe sib(s), with the
//                    concordant (severe-only) families' severe children as
//                    extra exclusion controls -> per family, passing genes (homozygous / compound_het /
//                    regulatory), the mild + severe VCFs of the supporting
//                    variants with hap-ibd IBD state + AF, a variant table,
//                    every mild AND severe variant in each passing gene
//                    (gene_variants.tsv.gz) and a per-step count summary
// GENE_RECURRENCE    whole cohort: how many mild children each gene passes
//                    in; for genes passing in >= modifier_min_mild of them,
//                    all mild / severe variants in the gene as two deduplicated
//                    AlphaGenome / AlphaMissense input VCFs (mild, severe)

process MODIFIER_VARIANTS {

    label 'cohort_python'
    publishDir "${params.results}/10_modifier_analysis/variants/families", mode: params.publish_mode

    input:
        val(sample_vcf_args)
        path(vcfs)
        path(family_json)
        path(ibd_segments, stageAs: 'ibd/*')    // cohort.sib_ibd_segments.tsv, or NO_FILE
        path(gtf)
        path(modifier_common)
        path(cnv_repeat_common)
        path(script)

    output:
        // <family>/<file>, published as families/<family>/<file>
        path("*/*.passing_genes.tsv"), emit: genes_ch
        path("*/*.modifier_variants.tsv.gz"), emit: variants_ch
        path("*/*.gene_variants.tsv.gz"), emit: gene_variants_ch
        path("*/*.mild.vcf.gz*"), emit: mild_vcf_ch
        path("*/*.severe.vcf.gz*"), emit: severe_vcf_ch
        path("*/*.modifier_summary.tsv"), emit: summary_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -family_json ${family_json} \\
        ${sample_vcf_args} \\
        -ibd_segments ${ibd_segments} \\
        -gtf ${gtf} \\
        -sv_match_dist ${params.modifier_sv_match_dist}

    # bgzip + tabix when the container has WORKING ones (the rnaseq image's
    # bgzip is on PATH but missing libcrypto.so.1.0.0), plain gzip otherwise
    use_bgzip=0
    if command -v bgzip >/dev/null 2>&1 && printf '' | bgzip -c >/dev/null 2>&1; then use_bgzip=1; fi
    for f in */*.mild.vcf */*.severe.vcf; do
        if [ \$use_bgzip -eq 1 ]; then
            bgzip -f "\$f"
            if command -v tabix >/dev/null 2>&1; then tabix -f -p vcf "\$f.gz" || true; fi
        else
            gzip -f "\$f"
        fi
    done
    """
}

process GENE_RECURRENCE {

    label 'cohort_python'
    publishDir "${params.results}/10_modifier_analysis/variants", mode: params.publish_mode

    input:
        path(passing_genes, stageAs: 'genes/*')
        path(gene_variants, stageAs: 'gene_variants/*')
        path(script)

    output:
        path("cohort.passing_genes.tsv"), emit: passing_genes_ch
        path("cohort.gene_recurrence.tsv"), emit: recurrence_ch
        path("cohort.candidate_gene_variants.tsv.gz"), emit: candidate_table_ch
        path("candidate_variants.mild.chm13.vcf"), emit: mild_candidate_vcf_ch
        path("candidate_variants.severe.chm13.vcf"), emit: severe_candidate_vcf_ch

    script:
    def max_af = params.modifier_candidate_max_af != null && params.modifier_candidate_max_af != '' ? "-max_af ${params.modifier_candidate_max_af}" : ''
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -passing_genes ${passing_genes} \\
        -gene_variants ${gene_variants} \\
        -min_mild ${params.modifier_min_mild} \\
        ${max_af}
    """
}
