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

    genome_name=\$(basename "${genome}")
    genome_name=\${genome_name%.fa}
    genome_name=\${genome_name%.fna}
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
        tuple val(sampleID), path("*.10pct.bam"), path("*.10pct.bam.bai"), emit : bam_ch

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

    singularity run \\
        docker://biocontainers/samtools:v1.9-4-deb_cv1 \\
        samtools view -@ ${task.cpus} -s 42.10 -b \${name}.minimap2.sorted.bam > \${name}.minimap2.sorted.10pct.bam
    
    singularity run \\
        docker://biocontainers/samtools:v1.9-4-deb_cv1 \\
        samtools index -@ ${task.cpus} \${name}.minimap2.sorted.10pct.bam

    """

}

process MAP_READS_TO_ASSEMBLY {

    tag "${params.sample_id}"
    label 'minimap2'
    publishDir "${params.results}/assembly/read_to_assembly", mode: params.publish_mode

    input:
        tuple val(sampleID), val(fastq), val(asm_fa), val(asm_fai)

    output:
        tuple val(sampleID), path(asm_fa), path(asm_fai), path("*.bam"), path("*.bai"), emit : alignment

    script:
    """
    minimap2 -t ${task.cpus} \
        -ax map-ont ${asm} \
        ${reads} | samtools sort -@ ${task.cpus} -o ${sampleID}_align_reads_to_assembly.sorted.bam

    samtools index ${sampleID}_align_reads_to_assembly.sorted.bam
    """
}

process ALIGN_HAPLOTYPES_TO_REFERENCE {
    
    tag "${sampleID}"
    publishDir "${params.results}/assembly/haplotype_alignments", mode: params.publish_mode

    input:
        val(ref)
        tuple val(sampleID), val(hap1), val(hap2)

    output:
        tuple val(sampleID), path("*hap1*"), path("*hap2*"), emit : aligned_haplotypes
    path "*hap1_to_ref.paf"
    path "*hap2_to_ref.paf"

    script:
    """
    minimap2 -t ${task.cpus} -x asm5 ${ref} ${hap1} > ${sampleID}_hap1_to_ref.paf
    minimap2 -t ${task.cpus} -x asm5 ${ref} ${hap2} > ${sampleID}_hap2_to_ref.paf
    """
}