p
rocess HAPDUP_PHASE {

    label 'hapdup'
    tag "${sampleID}_hapdup_phase"
    publishDir "${params.outdir}/hapdup/${sampleID}", mode: params.publish_mode

    input:
        tuple val(sampleID), val(asm_fa), val(asm_fai), val(bam), val(bai)

    output:
        tuple val(sampleID), path("hapdup/hap1.fa"), path("hapdup/hap2.fa"), emit : hapdup_ch
        path("*")

    script:
    """
    #!/bin/bash

    mkdir -p hapdup

    hapdup \
      --assembly ${asm_fa} \
      --bam ${bam} \
      --out-dir hapdup \
      --threads ${task.cpus}
    


    """
}