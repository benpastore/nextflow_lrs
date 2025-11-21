process MINIMAP2_INDEX { 

    label 'minimap2'

    publishDir "${index_dir}", mode : 'copy'

    input : 
        val genome 
        val index_dir
        val index_path
    
    output : 
        path("*mmi"), emit : minimap_index
        val(index_path), emit : index_ch


    script : 
    """
    #!/bin/bash

    genome_name=\$(basename ${genome} .fa)
    cp ${genome} \${genome_name}

    singularity run \\
        docker://staphb/minimap2:latest \\
        minimap2 \\
        -d \${genome_name}.mmi \\
        \${genome_name}

    """
}

process MINIMAP2_ALIGN { 

    label 'minimap2'

    publishDir "${params.results}/minimap2", mode : 'copy'

    input : 
        val index 
        tuple val(sampleID), val(fastq) 
    
    output : 
        tuple val(sampleID), path("*.bam"), path("*.bai"), emit : bam_ch

    script : 
    """
    #!/bin/bash
    
    name=\$(basename ${fastq} .fastq.gz)

    singularity run \\
        docker://staphb/minimap2:latest \\
        minimap2 \\
        -t ${task.cpus} \\
        -x map-ont \\
        -a \\
        -Y \\
        --MD \\
        ${index} \\
        ${fastq} > alignment.sam

    singularity run \\
        docker://biocontainers/samtools:v1.9-4-deb_cv1 \\
        samtools sort -@ ${task.cpus} -o \${name}.minimap2.sorted.bam alignment.sam

    singularity run \\
        docker://biocontainers/samtools:v1.9-4-deb_cv1 \\
        samtools index -@ ${task.cpus} \${name}.minimap2.sorted.bam

    """

}