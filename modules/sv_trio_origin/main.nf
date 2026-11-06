// Parent of origin for each child's Sniffles SVs (bin/sv_trio_origin.py):
// transmission (is the SV in the father's / mother's own Sniffles calls?)
// and haplotype (the SV's longphase haplotype, mapped to a parent from the
// child's phased SNVs where the parents' Clair3 genotypes differ), with a
// final origin + which evidence gave it. The SNV VCF must be the phased
// VCF of the BAM Sniffles ran on -- LONGPHASE round 1 in this pipeline --
// so the phase-block IDs match.

process SV_TRIO_ORIGIN {

    tag "${child}"
    label 'low'
    publishDir "${params.results}/06_variants/sv_trio_origin", mode: params.publish_mode

    input:
        tuple val(child), val(child_sv), val(child_snv), val(father), val(mother),
              val(father_snv), val(mother_snv), val(father_sv), val(mother_sv)
        path(script)

    output:
        tuple val(child), path("${child}.sv_trio_origin.tsv"), path("${child}.sv_trio_blocks.tsv"), emit: sv_trio_origin_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -child ${child} \\
        -child_sv ${child_sv} \\
        -father_sv ${father_sv} \\
        -mother_sv ${mother_sv} \\
        -child_snv ${child_snv} \\
        -father_snv ${father_snv} \\
        -mother_snv ${mother_snv} \\
        -merge_dist ${params.recessive_modifier_sv_merge_dist ?: 500} \\
        -output ${child}.sv_trio_origin.tsv \\
        -block_output ${child}.sv_trio_blocks.tsv
    """
}
