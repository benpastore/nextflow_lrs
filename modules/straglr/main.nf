process STRAGLR {

    tag "$sampleID"
    label 'straglr'

    publishDir "${params.results}/variants/straglr", mode: 'copy'

    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref_fa), val(ref_fai)
        
    output:
        tuple val(sampleID),
            path("*.straglr.tsv"),
            path("*.straglr.bed"),
            emit: straglr_ch

    script:
    def loci_arg = params.straglr_loci_bed ? "--loci ${params.straglr_loci_bed}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${bam} .bam)
    straglr.py \\
        ${bam} \\
        ${ref_fa} \\
        \${name}.straglr \\
        ${loci_arg} \\
        --genotype_in_size \\
        --min_support ${params.straglr_min_support ?: 2} \\
        --max_str_len ${params.straglr_max_str_len ?: 100} \\
        --max_num_clusters ${params.straglr_max_num_clusters ?: 2} \\
        --nprocs ${task.cpus} \\
        ${params.straglr_args ?: ''}
    """
}