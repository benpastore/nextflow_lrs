process BWA_INDEX {

    tag "${genome}_BWA_INDEX"

    label 'low'

    publishDir "${params.index_path}", mode : 'copy'
    
    input : 
        val genome
        val index

    output : 
        path("*")
        val("${params.index}/bwa/${params.genome_name}"), emit : bwa_index_ch

    script :
    """
    #!/bin/bash

    source activate rnaseq

    prefix=\$(basename ${genome} .fa)
    bwa index -p \$prefix -a ${params.bwa_index_alg} ${genome}
    
    samtools faidx ${genome}
    cut -f1,2 ${genome}.fai > chrom_sizes.txt
    """
}

process BWA_MEM {

    time { 5.hour * task.attempt } 
    errorStrategy 'retry'
    maxRetries 3 

    tag "${sampleID}_BWA_MEM"

    label 'high'

    publishDir "$params.results/bwa/bam", mode : 'copy', pattern : '*sorted.bam'
    publishDir "$params.results/bwa/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/bwa/logs", mode : 'copy', pattern : '*.summary.txt'

    input : 
        val bwa_idx
        tuple val(sampleID), val(fastq)

    output : 
        tuple val(sampleID), path("*sorted.bam"), emit : bwa_bam_ch
        tuple val(sampleID), path("*sorted.bam"), path("*bai"), emit : bwa_bam_bai_ch
        path("*.summary.txt")
    
    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}" 
    RG = ["@RG",
        "ID:${sampleID}",
        "SM:${sampleID}",
        "LB:1}",
        "PL:illumina"].join("\\t")
    """
    #!/bin/bash

    source activate rnaseq

    bwa mem ${bwa_idx} \\
        -t ${task.cpus} \\
        ${fastq_command} | samtools view -@ ${task.cpus} -b -h -O BAM -o ${sampleID}.bam -

    samtools sort  -@ ${task.cpus} ${sampleID}.bam -T ${sampleID} -o ${sampleID}.sorted.bam

    samtools index -@ ${task.cpus} ${sampleID}.sorted.bam

    samtools flagstat ${sampleID}.sorted.bam > ${sampleID}.sorted.bam.summary.txt
    """
}