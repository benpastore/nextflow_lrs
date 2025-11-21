
process BAM_TO_BW {

    tag "${sampleID}_BAM_to_BW"

    label 'deeptools'

    publishDir "$params.results/BigWig", mode: 'copy', pattern : "*.bw"

    input :
        tuple val(sampleID), val(bam), val(bai)

    output :
        tuple val(sampleID), path("*bw"), emit : bws_ch
    
    script :
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    singularity run \\
        docker://cabimerbioinfo/deeptools:v1 \\
        bamCoverage -b ${bam} \\
            -o \$name.bw \\
            -p ${task.cpus} \\
            --outFileFormat bigwig \\
            --binSize 10 \\
            --normalizeUsing CPM 
        
    """
}