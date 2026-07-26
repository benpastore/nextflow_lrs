process DORADO { 

    label 'dorado'

    publishDir "${params.results}/dorado", mode: params.publish_mode

    input : 
        tuple val(sampleID), val(unal_bam), val(reads), val(hapfasta), val(hapfai)
    
    output : 
        tuple val(sampleID), val(reads), path("*hap1.doradopolish.fa"), path("*hap2.doradopolish.fa"), emit : dorado_output_ch
    
    script : 
    """
    #!/bin/bash
    
    export TF_FORCE_UNIFIED_MEMORY='1'

    name=\$(basename ${reads} .fastq.gz)
    # Align unmapped reads to a reference using dorado aligner, sort and index
    
    dorado aligner ${hapfasta} ${unal_bam} | samtools sort --threads ${task.cpus} > aligned_reads.raw.bam

    # dorado polish requires a single read group, but the merged input bam can carry
    # multiple RGs from its source files -- overwrite them all with one uniform RG.
    # Preserve the DS:basecall_model=... field from the original @RG line, since
    # dorado polish needs it to auto-resolve its model and addreplacerg -w would
    # otherwise drop it when synthesizing a bare ID/SM-only RG.
    ds_field=\$(samtools view -H ${unal_bam} | awk -F'\\t' '\$1=="@RG"' | head -n1 | tr '\\t' '\\n' | grep '^DS:' || true)
    if [ -n "\$ds_field" ]; then
        rg_line=\$(printf '@RG\\tID:%s\\tSM:%s\\t%s' "${sampleID}" "${sampleID}" "\$ds_field")
    else
        rg_line=\$(printf '@RG\\tID:%s\\tSM:%s' "${sampleID}" "${sampleID}")
    fi
    samtools addreplacerg -r "\$rg_line" -w -@ ${task.cpus} -o aligned_reads.bam aligned_reads.raw.bam
    samtools index aligned_reads.bam

    # Call consensus
    dorado polish aligned_reads.bam ${hapfasta} > polished_assembly.fasta

    awk '/^>hap1_/ {p=1} /^>hap2_/ {p=0} p' polished_assembly.fasta > \${name}.hap1.doradopolish.fa
    awk '/^>hap2_/ {p=1} /^>hap1_/ {p=0} p' polished_assembly.fasta > \${name}.hap2.doradopolish.fa

    """
}

process DORADO_TRIM { 

    label 'dorado'

    publishDir "${params.results}/dorado", mode: params.publish_mode

    input : 
        tuple val(sampleID), val(fastq), val(unal_bam)
    
    output : 
        tuple val(sampleID), path("*dorado.trim.bam"), emit : dorado_trim_output_ch
    
    script : 
    """
    #!/bin/bash

    name=\$(basename ${unal_bam} .bam)

    dorado trim --sequencing-kit ${params.default_dorado_seq_kit} ${unal_bam} > \$name.dorado.trim.bam

    
    """
}