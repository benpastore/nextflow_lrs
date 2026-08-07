process MODKIT_PILEUP {

    tag "$sampleID"
    label 'modkit'

    publishDir "${params.results}/methylation/modkit", mode: params.publish_mode

    // Reference-coordinate methylation pileup, split by haplotype. Runs on
    // longphase_sv's haplotagged bam, so --partition-tag HP gives per-haplotype
    // bedMethyl output for free (haplotype 1, haplotype 2, and an ungrouped
    // file for reads that couldn't be phased) -- no separate phasing step
    // needed here. Depends on MM/ML tags surviving the bam -> fastq -> bam
    // round trip (see SAMTOOLS_CONVERT_BAM_TO_FASTQ's -T flag and
    // MINIMAP2_ALIGN's -y flag); if those were ever removed, this would
    // start silently producing empty/no-modification pileups.
    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref), val(ref_fai)

    output:
        tuple val(sampleID), path("${sampleID}*.bed.gz"), emit: modkit_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    modkit pileup \\
        ${bam} \\
        ${sampleID}.modkit.bed \\
        --ref ${ref} \\
        --partition-tag HP \\
        --threads ${task.cpus} \\
        --log-filepath ${sampleID}.modkit.log

    for f in ${sampleID}.modkit*.bed; do
        bgzip "\$f"
    done
    """
}
