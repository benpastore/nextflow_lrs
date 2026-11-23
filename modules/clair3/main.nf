process CLAIR3 {

    label 'clair3'
    tag "${sampleID}"
    publishDir "${params.results}/06_variants/clair3", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.clair3.vcf.gz"), path("*.clair3.vcf.gz.tbi"), emit : clair3_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    MODEL_NAME="r1041_e82_400bps_hac_v600"
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/\${MODEL_NAME} \\
        --output=\$(pwd) \\
        --include_all_ctgs \\
        --sample_name=${sampleID} \\
        --use_gpu
    
    mv merge_output.vcf.gz \$name.clair3.vcf.gz
    mv merge_output.vcf.gz.tbi \$name.clair3.vcf.gz.tbi
    """
}