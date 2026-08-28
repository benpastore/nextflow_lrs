/*
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

    minimap2 \\
        -d \${genome_name}.mmi \\
        \${genome_name}

    """
}
*/

process MINIMAP2_INDEX { 
    label 'minimap2'
    publishDir "${index_dir}", mode: 'copy'

    input: 
        path genome
        val index_dir

    output: 
        path("*.mmi"), emit: minimap_index

    script: 
    """
    set -euo pipefail

    idx="${index_dir}/${genome.simpleName}.mmi"

    if [[ -s "\$idx" ]]; then
        ln -s "\$idx" ${genome.simpleName}.mmi
    else
        minimap2 -d ${genome.simpleName}.mmi ${genome}
    fi
    """
}

process MINIMAP2_ALIGN { 

    label 'minimap2'

    publishDir "${params.results}/minimap2", mode : 'copy'

    input : 
        path index 
        tuple val(sampleID), val(fastq) 
    
    output : 
        tuple val(sampleID), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit : bam_ch

    script :
    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${fastq} .fastq.gz)

    # No -y here: this only ever runs on herro-corrected reads (see
    # main.nf, ont_ch_corrected), whose FASTA headers carry a literal
    # "None" comment (herro's own writer, not a real MM/ML tag) rather
    # than the methylation tags -y expects -- copying it into the SAM aux
    # field produces unparseable output ("Incomplete aux field"). herro
    # also rewrites/consensus-corrects the sequence, so original per-base
    # modification calls wouldn't have valid coordinates on it anyway.
    # Piped straight into sort rather than written to an intermediate SAM
    # file -- with --MD on a full WGS alignment the uncompressed SAM can run
    # to 100s of GB, and a full/quota-exceeded disk truncates it silently.
    minimap2 \\
        -t ${task.cpus} \\
        -x map-ont \\
        -a \\
        -Y \\
        --MD \\
        ${index} \\
        ${fastq} \\
      | samtools sort -@ ${task.cpus} -o \${name}.minimap2.sorted.bam -

    samtools index -@ ${task.cpus} \${name}.minimap2.sorted.bam

    samtools view -@ ${task.cpus} -s 42.10 -b \${name}.minimap2.sorted.bam > \${name}.minimap2.sorted.10pct.bam
    
    samtools index -@ ${task.cpus} \${name}.minimap2.sorted.10pct.bam

    """

}

process MINIMAP2_ALIGN_METHYLATION {

    label 'minimap2'

    publishDir "${params.results}/minimap2_methylation", mode : 'copy'

    input :
        path index
        tuple val(sampleID), val(fastq)

    output :
        tuple val(sampleID), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit : bam_ch

    script :
    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${fastq} .fastq.gz)

    # Same as MINIMAP2_ALIGN but on the RAW (pre-herro) reads, which is
    # the only place real MM/ML methylation tags survive -- -y copies them
    # (carried as fastq comments by SAMTOOLS_CONVERT_BAM_TO_FASTQ's -T
    # flag) onto the aligned reads. This bam has no HP (haplotype) tags of
    # its own; see TRANSFER_HP_TAGS, which projects them over from the
    # herro-corrected/phased bam before modkit runs.
    minimap2 \\
        -t ${task.cpus} \\
        -x map-ont \\
        -a \\
        -Y \\
        -y \\
        --MD \\
        ${index} \\
        ${fastq} \\
      | samtools sort -@ ${task.cpus} -o \${name}.minimap2.sorted.bam -

    samtools index -@ ${task.cpus} \${name}.minimap2.sorted.bam
    """

}

process MAP_ASM_TO_REF {

    tag "$sampleID"
    label 'minimap2'

    publishDir "${params.results}/minimap2_map_asm_to_ref", mode : 'copy'

    input:
        path ref_mmi 
        tuple val(sampleID), val(reads), val(hap1_fa), val(hap2_fa)
    

    output:
    tuple val(sampleID),
          path("${sampleID}.hap1.ref.bam"),
          path("${sampleID}.hap1.ref.bam.bai"),
          path("${sampleID}.hap2.ref.bam"),
          path("${sampleID}.hap2.ref.bam.bai"),
          emit: hap_bams

    script:
    """
    #!/bin/bash
    set -euo pipefail

    minimap2 -ax asm20 -t ${task.cpus} ${ref_mmi} ${hap1_fa} \
      | samtools sort -@ ${task.cpus} -o ${sampleID}.hap1.ref.bam

    samtools index -@ ${task.cpus} ${sampleID}.hap1.ref.bam

    minimap2 -ax asm20 -t ${task.cpus} ${ref_mmi} ${hap2_fa} \
      | samtools sort -@ ${task.cpus} -o ${sampleID}.hap2.ref.bam

    samtools index -@ ${task.cpus} ${sampleID}.hap2.ref.bam
    """
}