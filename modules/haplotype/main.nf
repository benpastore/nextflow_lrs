// Mild vs severe sibling haplotype sharing with hap-ibd, then the mild
// child's variants on copies the severe sibling didn't inherit
// (trio/quad families only).
//
//   ORIENT_PHASE             per family: WhatsHap family VCF -> children only,
//                            phase oriented paternal|maternal chromosome-wide
//                            (bin/orient_family_phase.py)
//   HAP_IBD                  per family: hap-ibd on that VCF (uniform 1 cM/Mb map)
//   IBD_SEGMENTS             cohort: hap-ibd segments -> per pair, per parent,
//                            shared / not_shared + IBD0/1/2 (bin/hapibd_segments.py)
//   MILD_HAPLOTYPE_VARIANTS  cohort: mild variants on a not_shared copy
//                            (bin/mild_haplotype_variants.py)

process ORIENT_PHASE {

    tag "${family}"
    label 'cohort_python'
    publishDir "${params.results}/07_haplotypes/oriented", mode: params.publish_mode

    input:
        tuple val(family), path(family_vcf), path(family_tbi)
        path(family_json)
        path(script)

    output:
        tuple val(family), path("${family}.oriented.vcf.gz"), emit: oriented_ch
        path("${family}.orient_blocks.tsv"), emit: blocks_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -vcf ${family_vcf} \\
        -family_json ${family_json} \\
        -min_gq ${params.hapibd_min_gq} \\
        -exclude_chroms '${params.hap_exclude_chroms}' \\
        -output ${family}.oriented.vcf.gz \\
        -blocks ${family}.orient_blocks.tsv
    """
}

process HAP_IBD {

    tag "${family}"
    label 'hapibd'
    publishDir "${params.results}/07_haplotypes/hapibd", mode: params.publish_mode

    // hap-ibd (Zhou, Browning & Browning 2020) needs a genetic map in PLINK
    // format on the VCF's coordinates. A uniform 1 cM/Mb map is generated
    // from the reference .fai: sibling IBD segments are tens of Mb, so local
    // recombination-rate detail barely matters. min-mac 1: a variant only one
    // child carries is still evidence their haplotypes differ there.
    input:
        tuple val(family), path(oriented_vcf)
        tuple val(ref), val(ref_fai)

    output:
        tuple val(family), path("${family}.hapibd.ibd.gz"), emit: ibd_ch
        path("${family}.hapibd.log"), emit: log_ch

    script:
    if (!params.hapibd_jar) { error "params.hapibd_jar not set -- download hap-ibd.jar (https://faculty.washington.edu/browning/hap-ibd.jar) and point hapibd_jar at it" }
    """
    #!/bin/bash
    set -euo pipefail

    awk 'BEGIN{OFS="\\t"} { printf "%s\\t.\\t0\\t1\\n%s\\t.\\t%.6f\\t%d\\n", \$1, \$1, \$2/1e6, \$2 }' ${ref_fai} > uniform_1cM_per_Mb.map

    java -Xmx${params.hapibd_java_mem} -jar ${params.hapibd_jar} \\
        gt=${oriented_vcf} \\
        map=uniform_1cM_per_Mb.map \\
        out=${family}.hapibd \\
        min-seed=${params.hapibd_min_seed} \\
        min-extend=${params.hapibd_min_extend} \\
        min-output=${params.hapibd_min_output} \\
        max-gap=${params.hapibd_max_gap} \\
        min-mac=1 \\
        nthreads=${task.cpus}
    """
}

process IBD_SEGMENTS {

    label 'cohort_python'
    publishDir "${params.results}/07_haplotypes", mode: params.publish_mode

    input:
        path(family_json)
        val(families)
        path(oriented_vcfs)
        path(ibd_files)
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
        -families ${families.join(' ')} \\
        -bridge_cm ${params.hapibd_bridge_cm} \\
        -output_prefix cohort
    """
}

process MILD_HAPLOTYPE_VARIANTS {

    label 'cohort_python'
    publishDir "${params.results}/08_mild_haplotype_variants", mode: params.publish_mode

    input:
        path(family_json)
        path(vcf_manifest)       // childID<TAB>annotated master VCF name
        path(vcfs)
        path(parent_segments)
        path(script)

    output:
        path("cohort.mild_haplotype_variants.tsv.gz"), emit: variants_ch
        path("cohort.mild_haplotype_genes.tsv"), emit: genes_ch
        path("cohort.mild_haplotype_funnel.tsv"), emit: funnel_ch

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
        -sv_merge_dist ${params.sv_merge_dist ?: 500} \\
        -rare_af ${params.mild_variants_rare_af} \\
        -output_prefix cohort
    """
}
