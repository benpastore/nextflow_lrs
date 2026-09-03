process BAM_COMPARE {

    tag "${ip}_vs_${control}_BAM_Compare"

    label 'high'

    publishDir "$params.results/BamCompare", mode: 'copy', pattern : "*.bw"

    input : 
        tuple val(ip), path(ip_bam), val(control), path(control_bam)

    output :
        path("*bw"), emit : bws_ch
    
    script :
    normalize = params.normalize_bw ? "--normalizeUsing ${params.normalize_bw} --scaleFactorsMethod None" : '--normalizeUsing CPM --scaleFactorsMethod None'
    smooth = params.smooth_length ? "--smoothLength ${params.smooth_length}" : '--smoothLength 10' 
    exact = params.exact_scaling ? "--exactScaling" : ''
    filter_strand = params.filter_strand ? '--filterRNAstrand forward' : ''
    operation = params.bw_compare_op ? "--operation ${params.bw_compare_op}" : "--operation subtract"
    """
    #!/bin/bash

    source activate rnaseq

    samtools index -@ ${task.cpus} ${ip_bam}
    samtools index -@ ${task.cpus} ${control_bam}

    source activate DeepTools

    name=\$(basename ${ip_bam} .bam)

    bamCompare -b1 ${ip_bam} \\
        -b2 ${control_bam} \\
        -o ${ip}_vs_${control}_compare.bw \\
        ${normalize} \\
        ${smooth} \\
        ${exact} \\
        ${operation} \\
        ${filter_strand}
    """
}

process BAM_TO_BW {

    tag "${sampleID}_BAM_to_BW"

    label 'high'

    publishDir "$params.results/BigWig", mode: 'copy', pattern : "*.bw"

    input :
        tuple val(sampleID), path(bam), path(bai)

    output :
        tuple val(sampleID), path("*bw"), emit : bws_ch
    
    script :
    normalize = params.normalize_bw ? "--normalizeUsing ${params.normalize_bw}" : '--normalizeUsing CPM'
    smooth = params.smooth_length ? "--smoothLength ${params.smooth_length}" : '--smoothLength 10' 
    exact = params.exact_scaling ? "--exactScaling" : ''
    filter_strand = params.filter_strand ? '--filterRNAstrand forward' : ''
    bin_size = params.bin_size ? "--binSize ${params.bin_size}" : '--binSize 5'
    """
    #!/bin/bash

    source activate DeepTools

    name=\$(basename ${bam} .bam)

    bamCoverage -b ${bam} \\
        -o \$name.bw \\
        -p ${task.cpus} \\
        ${normalize} \\
        ${smooth} \\
        ${exact} \\
        ${filter_strand} \\
        ${bin_size}
        
    """
}

process BW_COMPARE {

    tag "${ip}_vs_${control}_MACS2"

    label 'high'

    publishDir "$params.results/subtract_BigWig", mode: 'copy', pattern : "*.bw"

    input : 
        tuple val(control), val(ip), path(ip_bw), path(control_bw)

    output :
        path("*bw")
    
    script :
    pseudocount = params.bw_compare_pseudocount ? "--pseudocount ${params.bw_compare_pseudocount}" : ''
    operation = params.bw_compare_operation ? "--operation ${params.bw_compare_operation}" : '--operation ratio' 
    bin_size = params.bin_size ? "--binSize ${params.bin_size}" : '--binSize 5'
    """
    #!/bin/bash

    source activate DeepTools

    name=${ip}.subtract_${control}.bw

    bigwigCompare -b1 ${ip_bw} \\
        -b2 ${control_bw} \\
        -o \$name \\
        -p ${task.cpus} \\
        ${operation} \\
        ${pseudocount} \\
        ${bin_size}
        
    """
}

process MERGE_BW {

    tag "${condition}_merge_bw"

    label 'high'

    publishDir "$params.results/merge_BigWig", mode: 'copy', pattern : "*.bw"

    input : 
        tuple val(sampleIDs), path(bws), val(condition)
        val chrom_sizes

    output : 
        tuple val(condition), path("*.bw"), emit : merge_bw_ch
    
    script : 
    """
    #!/bin/bash

    source activate rnaseq

    # make array with bigwigs 
    bws=(${bws.join(' ')})

    # get length of array
    N_bw=\${#bws[@]}

    echo \$N_bw

    if [ "\$N_bw" -gt 1 ]; then 

        # merge bigwig(s) --> bedgraph
        bigWigMerge ${bws.join(' ')} ${condition}.tmp

        # divide counts column by N samples to average
        cat ${condition}.tmp | awk -F'\\t' -v OFS='\\t' -v nsamp=\$N_bw '{ print \$1,\$2,\$3,\$4/nsamp }' > ${condition}.bedGraph

        # sort the bed file 
        bedSort ${condition}.bedGraph ${condition}.bedGraph.sorted

        # convert bedGraph --> bigwig
        bedGraphToBigWig ${condition}.bedGraph.sorted ${chrom_sizes} ${condition}.bw
    else 
        cp ${bws.join(' ')} ${condition}.bw
    fi
    """

}