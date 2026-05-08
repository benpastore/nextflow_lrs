process STRAGLR {

    tag "$sampleID"
    label 'straglr'

    publishDir "${params.results}/straglr", mode: 'copy'

    input:
    tuple val(sampleID),
          path(bam),
          path(bai),
          path(ref_fa),
          path(ref_fai),
          path(loci_bed), emit : straglr_ch

    output:
    tuple val(sampleID),
          path("${sampleID}.straglr.tsv"),
          path("${sampleID}.straglr.bed"),
          emit: straglr_ch

    script:
    def loci_arg = loci_bed ? "--loci ${loci_bed}" : ""

    """
    set -euo pipefail

    straglr.py \\
        ${bam} \\
        ${ref_fa} \\
        ${sampleID}.straglr \\
        ${loci_arg} \\
        --genotype_in_size \\
        --min_support ${params.straglr_min_support ?: 2} \\
        --max_str_len ${params.straglr_max_str_len ?: 100} \\
        --max_num_clusters ${params.straglr_max_num_clusters ?: 2} \\
        --nprocs ${task.cpus} \\
        ${params.straglr_args ?: ''}

    mv ${sampleID}.straglr.tsv ${sampleID}.straglr.tsv
    mv ${sampleID}.straglr.bed ${sampleID}.straglr.bed
    """
}