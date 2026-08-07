process FMLRC2_CONVERT {

    tag "$sampleID"
    label 'fmlrc2'

    input:
        tuple val(sampleID), val(ropebwt2_txt)

    output:
        tuple val(sampleID), path("*.msbwt.npy"), emit: fmlrc2_convert_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    cat ${ropebwt2_txt} | tr NT TN | fmlrc2-convert ${sampleID}.msbwt.npy
    """
}

process FMLRC2_CORRECT {

    tag "$sampleID"
    label 'fmlrc2'

    input:
        tuple val(sampleID), val(msbwt), val(fastq)

    output:
        tuple val(sampleID), path("*.fmlrc2.corrected.fasta"), emit: fmlrc2_correct_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    zcat -f ${fastq} > long_reads.fastq

    fmlrc2 -t ${task.cpus} ${msbwt} long_reads.fastq ${sampleID}.fmlrc2.corrected.fasta
    """
}
