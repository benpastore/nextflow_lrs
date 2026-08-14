process SORT_READS_FOR_ROPEBWT2 {

    tag "$sampleID"
    label 'sort_reads'

    input:
        tuple val(sampleID), val(r1), val(r2)

    output:
        tuple val(sampleID), path("${sampleID}.sorted_seqs.txt"), emit: sorted_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # runs on the general-purpose container (real GNU coreutils, unlike the
    # minimal ropebwt2 biocontainer -- see ROPEBWT2 below) so -S/--parallel/-T
    # actually work and sort can spill to disk instead of crashing in-memory
    # on WGS-scale read counts.
    zcat -f ${r1} ${r2} \\
        | awk 'NR % 4 == 2' \\
        | sort -S 224G --parallel=${task.cpus} -T . \\
        | tr NT TN > ${sampleID}.sorted_seqs.txt
    """
}

process ROPEBWT2 {

    tag "$sampleID"
    label 'ropebwt2'

    input:
        tuple val(sampleID), path(sorted_seqs)

    output:
        tuple val(sampleID), path("*.ropebwt2.txt"), emit: ropebwt2_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    ropebwt2 -LR < ${sorted_seqs} > ${sampleID}.ropebwt2.txt
    """
}
