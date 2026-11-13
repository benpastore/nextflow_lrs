// Mild vs severe haplotype comparison + genetic models (trio/quad families only).
//
// Each variant's parental haplotype (INFO/ORIGIN) is already in the
// consolidated+snpEff VCF (bin/consolidate_variants.py). These two steps add
// the sibling comparison:
//
//   SIB_HAPLOTYPE_SHARING  which parental haplotypes each mild/severe sibling
//                          pair shares along the genome (IBD0/1/2), from
//                          genotypes alone (bin/sib_haplotype_sharing.py)
//   HAPLOTYPE_MODIFIER     labels every mild-child variant shared /
//                          mild_specific, then runs the recessive / dominant
//                          models on those haplotypes (bin/haplotype_modifier.py).
//                          Writes the same cohort.{recessive,dominant}_modifier*.tsv
//                          files as RECESSIVE_MODIFIER, so everything downstream
//                          (AlphaGenome, AlphaMissense, report) is unchanged.

process SIB_HAPLOTYPE_SHARING {

    label 'recessive_modifier'
    publishDir "${params.results}/07_haplotypes", mode: params.publish_mode

    input:
        path(family_json)
        path(vcf_manifest)   // sampleID<TAB>vcf name: parents' Clair3 + children's WhatsHap-trio VCFs
        path(vcfs)
        path(script)

    output:
        path("cohort.sib_parent_segments.tsv"), emit: parent_segments_ch
        path("cohort.sib_ibd_segments.tsv"), emit: ibd_segments_ch
        path("cohort.sib_ibd_summary.tsv"), emit: ibd_summary_ch

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
        -error_rate ${params.hap_ibd_error_rate} \\
        -switch_rate ${params.hap_ibd_switch_rate} \\
        -min_segment_sites ${params.hap_ibd_min_segment_sites} \\
        -exclude_chroms ${params.hap_exclude_chroms} \\
        -output_prefix cohort
    """
}

process HAPLOTYPE_MODIFIER {

    label 'recessive_modifier'
    publishDir "${params.results}/07_haplotypes", mode: params.publish_mode

    input:
        path(family_json)
        path(vcf_manifest)                              // childID<TAB>snpEff VCF name
        path(vcfs)
        path(parent_segments)
        path(common_script)
        path(script)

    output:
        path("cohort.recessive_modifier.tsv"), emit: recessive_modifier_ch
        path("cohort.recessive_modifier.gene_summary.tsv"), emit: recessive_modifier_gene_summary_ch
        path("cohort.dominant_modifier.tsv"), emit: dominant_modifier_ch
        path("cohort.dominant_modifier.gene_summary.tsv"), emit: dominant_modifier_gene_summary_ch
        path("cohort.haplotype_variants.tsv"), emit: haplotype_variants_ch
        path("cohort.haplotype_funnel.tsv"), emit: haplotype_funnel_ch

    // gnomAD / 1000 Genomes AF are read from the VCFs' INFO (ANNOTATE_AF)
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
        -sib_parent_segments ${parent_segments} \\
        -sv_merge_dist ${params.recessive_modifier_sv_merge_dist ?: 500} \\
        -output_prefix cohort
    """
}
