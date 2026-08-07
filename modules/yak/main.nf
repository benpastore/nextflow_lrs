process YAK_COUNT {

    tag "$sampleID"
    label 'yak'

    publishDir "${params.results}/yak", mode: params.publish_mode

    input:
        tuple val(sampleID), val(r1), val(r2)

    output:
        tuple val(sampleID), path("*.yak"), emit: yak_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    yak count -k31 -b37 -t ${task.cpus} -o ${sampleID}.yak ${r1} ${r2}
    """
}
