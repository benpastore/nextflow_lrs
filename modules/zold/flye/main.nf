process FLYE_ASSEMBLE {

    label 'flye'
    publishDir "${params.outdir}/assembly/flye", mode: params.publish_mode

    input:
        tuple val(sampleID), val(fastq) 
        val genome

    output:
        tuple val(sampleID), path("${sampleID}_assembly.fasta"), path("${sampleID}_assembly.fasta.fai"), emit : flye_assembly

    script:
    """
    flye \
      --${params.flye_model} ${fastq} \
      --genome-size ${params.genome_size} \
      --threads ${task.cpus} \
      --out-dir flye_out

    cp flye_out/assembly.fasta ${sampleID}_assembly.fasta

    samtools faidx ${sampleID}_assembly.fasta
    """
}