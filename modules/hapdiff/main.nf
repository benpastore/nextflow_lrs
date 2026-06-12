process HAPDIFF { 

    label 'hapdiff'

    publishDir "${params.results}/variants/hapdiff", mode: params.publish_mode

    input :
        tuple val(sampleID), val(reads), val(hap1), val(hap2)
        tuple path(ref), path(ref_fai)

    output : 
        tuple val(sampleID), path("*phased.vcf.gz"), emit : hapdiff_ch

    script : 
    """
    #!/bin/bash
    
    name=\$(basename ${reads} .fastq.gz)
    
    hapdiff.py --reference ${ref} \\
        --pat ${hap1} \\
        --mat ${hap2} \\
        --out-dir \$(pwd) \\
        -t ${task.cpus}

    mv hapdiff_unphased.vcf.gz \${name}.hapdiff.unphased.vcf.gz
    mv hapdiff_phased.vcf.gz \${name}.hapdiff.phased.vcf.gz
   """
}