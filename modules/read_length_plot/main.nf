// Publication-quality read-length histogram (reads vs length, N50 marked)
// per sample -- see bin/plot_read_length.py. Reads the same trimmed FASTQ
// NANOPLOT_RAW does, but is fully independent of it: nothing downstream
// consumes its output, so adding/changing it never invalidates any other
// step's -resume cache.

process READ_LENGTH_PLOT {

    tag "${sampleID}"
    label 'read_length_plot'
    publishDir "${params.results}/00_qc/read_length", mode: params.publish_mode

    input:
        tuple val(sampleID), val(bam), path(fastq)
        path(script)

    output:
        tuple val(sampleID), path("${sampleID}.read_length.png"), path("${sampleID}.read_length.pdf"), emit: plot_ch
        tuple val(sampleID), path("${sampleID}.read_length.stats.tsv"), emit: stats_ch
        path("${sampleID}.read_length.counts.tsv.gz"), emit: counts_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    # \$HOME is read-only inside the container
    export MPLCONFIGDIR=\$PWD/.matplotlib

    python3 ${script} \\
        -fastq ${fastq} \\
        -sample ${sampleID} \\
        -prefix ${sampleID}.read_length \\
        ${params.read_length_plot_args ?: ''}
    """
}
