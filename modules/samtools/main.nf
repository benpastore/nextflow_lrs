process INDEX_BAM {

    tag "${condition}_index_bam"

    label 'high'

    errorStrategy 'ignore'

    publishDir "$params.results/bams", mode : 'copy', pattern : '*.bam*'

    input : 
        tuple val(condition), val(bam)

    output : 
        tuple val(condition), path("*.bam"), path("*bai")
    
    script :
    """
    #!/bin/bash

    source activate rnaseq

    cp ${bam} .
    name=\$(basename ${bam})
    samtools index \$name 
    
    """
}

process FILTER_BAM {

    errorStrategy 'retry'
    maxRetries 3
    time { 5.hour * task.attempt } 
    
    tag "${condition}_filter_merge_bam"

    label 'high'

    publishDir "$params.results/filtered/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/filtered/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/filtered/log", mode : 'copy', pattern : '*stat'
    publishDir "$params.results/filtered/log", mode : 'copy', pattern : '*stat'

    input : 
        tuple val(condition), val(bam), val(bai)

    output : 
        tuple val(condition), path("*filtered.sorted.bam"), emit : filtered_bams_ch
        tuple val(condition), path("*filtered.sorted.bam"), path("*filtered.sorted.bam.bai"), emit : bam_to_bw_input
        path("*stats")
    
    script:
    filter_params = params.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 20'
    blacklist_params = params.blacklist ? "-L $bed" : ''
    remove_orphan = params.remove_orphan && !params.single_end ? "-f 0x2" : ''
    select_primary = params.select_primary ? "-F 0x0100" : ''
    """
    #!/bin/bash

    source activate rnaseq

    name=\$(basename $bam .bam)

    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
        $remove_orphan \\
        $select_primary \\
        -h \\
        -b ${bam} > \$name.filtered.bam

    samtools sort  -@ ${task.cpus} \$name.filtered.bam -T \$name -o \$name.filtered.sorted.bam
    samtools index \$name.filtered.sorted.bam
    samtools flagstat \$name.filtered.sorted.bam > \$name.filtered.sorted.flagstat
    samtools idxstats \$name.filtered.sorted.bam > \$name.filtered.sorted.idxstats
    samtools stats \$name.filtered.sorted.bam > \$name.filtered.sorted.stats
    """   
}

process INDEX_REFERENCE {

    label 'high'
    tag "index_reference"

    publishDir "${params.results}/reference", mode: params.publish_mode

    input:
        path(ref)

    output:
        tuple path(ref), path("${ref}.fai"), emit : ref_indexed_ch

    script:
    """
    #!/bin/bash

    source activate rnaseq

    samtools faidx ${ref}
    """
}

process INDEX_FLYE_ASSEMBLY {

    tag "index_flye_assembly"
    label 'samtools'

    publishDir "${params.results}/assembly/flye", mode: params.publish_mode

    input:
        path assembly

    output:
        tuple path(assembly), path("${assembly}.fai"), emit : flye_assembly_index_ch

    script:
    """
    #!/bin/bash
    
    samtools faidx ${assembly}
    """
}


process SAMTOOLS_CONVERT_BAM_TO_FASTQ {

    tag "convert_bam2fastq"
    label 'samtools_high'

    publishDir "${params.results}/samtools", mode: params.publish_mode

    input:
        tuple val(sampleID), val(unal_bam)

    output:
        tuple val(sampleID), val(unal_bam), path("*.fastq.gz"), emit : samtool_convert_unalbam2fastq

    script:
    """
    #!/bin/bash

    set -euo pipefail
    name=\$(basename ${unal_bam} .bam)
    # -T MM,ML,MN carries dorado's methylation-call tags through as fastq
    # comments; minimap2 -y (see MINIMAP2_ALIGN) then copies them back onto
    # the aligned reads. Without this here, they're gone before minimap2
    # ever sees the reads, regardless of that flag.
    samtools fastq -@ ${task.cpus} -T MM,ML,MN ${unal_bam} | gzip > \$name.fastq.gz

    """
}

process SAMTOOLS_FASTQ_TO_BAM {

    tag "$sampleID"
    label 'samtools_high'

    publishDir "${params.results}/samtools", mode: params.publish_mode

    input:
        tuple val(sampleID), path(fastq), path(orig_bam)

    output:
        tuple val(sampleID), path("*.unal.bam"), path(fastq), emit : fastq_to_bam_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    name=\$(basename ${fastq} .fastq.gz)

    # reuse the real @RG line (incl. basecall_model=... that dorado polish needs to
    # auto-resolve its model) from the original bam -- samtools import alone can only
    # synthesize a bare ID/SM line, which loses that field entirely
    rg_line=\$(samtools view -H ${orig_bam} | awk -F'\\t' '\$1=="@RG"' | head -n1)
    if [ -z "\$rg_line" ]; then
        rg_line=\$(printf '@RG\\tID:%s\\tSM:%s' "${sampleID}" "${sampleID}")
    fi

    samtools import -r "\$rg_line" ${fastq} -o \${name}.unal.bam
    """
}

process MERGE_INPUTS {

    tag "merge_inputs"
    label 'samtools_high'

    input:
    tuple val(sample),
          path(fastqs),
          path(bams)

    output:
    tuple val(sample),
          path("${sample}.fastq.gz"),
          path("${sample}.bam"), emit : merged_inputs_ch

    script:
    """
    # Merge FASTQs
    cat ${fastqs.join(' ')} > ${sample}.fastq.gz

    # Merge BAMs
    samtools merge -@ ${task.cpus} -o ${sample}.bam ${bams.join(' ')}

    samtools index -@ ${task.cpus} ${sample}.bam
    """
}