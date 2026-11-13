// Population allele frequencies written into every child's consolidated +
// snpEff VCF, once, right after snpEff (bin/annotate_af.py):
//   INFO/gnomAD_AF_NFE  gnomAD AF_nfe, CHM13 projection (params.gnomad_chm13_vcf)
//   INFO/KG_AF_EUR      1000 Genomes EUR, re-called on CHM13 (kg_af subworkflow)
// The two cover each other's gaps (the gnomAD projection misses sites where
// CHM13 carries the minor allele). Annotation only -- nothing is filtered.
// Whole cohort in one task so each population file is streamed once.

process ANNOTATE_AF {

    label 'recessive_modifier'
    publishDir "${params.results}/06_variants/snpeff_af", mode: params.publish_mode

    input:
        path(manifest)                          // sampleID<TAB>snpEff VCF name
        path(vcfs)
        path(kg_af_vcf, stageAs: 'kg_af/*')     // or NO_FILE
        path(modifier_common)
        path(script)

    output:
        path("*.snpeff.af.vcf.gz"), emit: af_vcfs_ch

    script:
    def kg_arg = kg_af_vcf.toString().tokenize('/')[-1] != 'NO_FILE' ? "-kg_af_vcf ${kg_af_vcf}" : ''
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -manifest ${manifest} \\
        ${params.gnomad_chm13_vcf ? "-gnomad_vcf ${params.gnomad_chm13_vcf}" : ""} \\
        ${kg_arg}
    """
}
