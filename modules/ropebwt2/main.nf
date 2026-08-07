process ROPEBWT2 {

    tag "$sampleID"
    label 'ropebwt2'

    input:
        tuple val(sampleID), val(r1), val(r2)

    output:
        tuple val(sampleID), path("*.ropebwt2.txt"), emit: ropebwt2_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    zcat -f ${r1} ${r2} \\
        | awk 'NR % 4 == 2' \\
        | sort \\
        | tr NT TN \\
        | ropebwt2 -LR > ${sampleID}.ropebwt2.txt
    """
}
