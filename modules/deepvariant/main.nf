process DEEPVARIANT_CALL_VARIANTS {

    time { 5.hour * task.attempt } 
    errorStrategy 'retry'
    maxRetries 3 

    tag "${condition}_deepvariant"

    label 'DeepVariant'

    publishDir "$params.results/deepvariant", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai)
        val genome

    output :
        // narrowed to the small-variant vcf specifically -- "*.vcf.gz" also
        // matches "${condition}.dv.gvcf.gz" since gzip's glob upstream
        // catches both, which would otherwise make this a 2-file list
        tuple val(condition), path("${condition}.dv.vcf.gz"), path("${condition}.dv.vcf.gz.tbi"), emit : deepvariant_vcf_ch
        //tuple val(condition), val("*.gvcf.gz"), emit : deepvariant_gvcf_ch
        //path("*")
        //path("*vcf.gz"), emit : vcf
    
    script:
    """
    #!/bin/bash

    export TF_FORCE_UNIFIED_MEMORY='1'

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type="WGS" \\
        --ref=${genome} \\
        --reads=${bam} \\
        --output_vcf=${condition}.dv.vcf \\
        --output_gvcf=${condition}.dv.gvcf \\
        --num_shards=${task.cpus} \\
        --logging_dir=\$PWD/${condition}_dv_logs

    gzip *vcf*

    bgzip -@ 8 *vcf
    tabix -p vcf *.vcf.gz
    """   
}

process SPLIT_DV_VCF {

    tag "${condition}_filter_dv"

    label 'bcftools_low'

    input : 
        tuple val(condition), val(vcf)

    output : 
        tuple val(condition), path("*.vcf.gz"), emit : vcfs
    
    script:
    """
    #!/bin/bash

    zcat ${vcf} | grep "^#"  > header

    zcat ${vcf} | awk '(\$7=="PASS")' > passing

    cat header passing | bgzip > vcf.gz

    tabix vcf.gz

    chroms=\$(zcat vcf.gz | grep -v "^#" | cut -f1 | sort | uniq | grep chr)

    for chrom in \$chroms; do
        echo \$chrom
        bcftools view -r \$chrom vcf.gz -Oz -o \${chrom}.vcf.gz
    done

    """  

}