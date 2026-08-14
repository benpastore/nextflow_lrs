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

    # -S caps sort's own memory well below the SLURM --mem allocation (see
    # nextflow.config, label 'ropebwt2') -- GNU sort otherwise grabs RAM
    # unboundedly and ignores the cgroup limit, which is what was OOM-killing
    # this step; -T . spills to the task work dir instead of a possibly-small
    # shared /tmp once past that cap.
    zcat -f ${r1} ${r2} \\
        | awk 'NR % 4 == 2' \\
        | sort -S 48G --parallel=${task.cpus} -T . \\
        | tr NT TN \\
        | ropebwt2 -LR > ${sampleID}.ropebwt2.txt
    """
}
