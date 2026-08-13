process YAK_COUNT {

    tag "$sampleID"
    label 'yak'

    publishDir "${params.results}/yak", mode: params.publish_mode

    // fastqs is a list (one ONT fastq.gz, or a pair of Illumina R1/R2) --
    // yak count just wants a whitespace-separated list of read files.
    input:
        tuple val(sampleID), val(fastqs)

    output:
        tuple val(sampleID), path("*.yak"), emit: yak_ch

    script:
    def fastq_list = (fastqs instanceof List ? fastqs : [fastqs]).join(' ')
    """
    #!/bin/bash
    set -euo pipefail

    yak count -k31 -b37 -t ${task.cpus} -o ${sampleID}.yak ${fastq_list}
    """
}
