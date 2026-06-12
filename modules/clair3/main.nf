process CLAIR3 {

    label 'clair3'
    tag "${sampleID}"
    publishDir "${params.results}/variants/clair3", mode: params.publish_mode

    input:
        tuple val(ref), val(ref_fai)
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("*.clair3.vcf.gz"), path("*.clair3.vcf.gz.tbi"), emit : clair3_ch

    script:
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    MODEL_NAME="r941_prom_sup_g5014"
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/\${MODEL_NAME} \\
        --output=\$(pwd) \\
        --include_all_ctgs \\
        --use_gpu
    
    mv merge_output.vcf.gz \$name.clair3.vcf.gz
    mv merge_output.vcf.gz.tbi \$name.clair3.vcf.gz.tbi

    #singularity run --nv --cleanenv --env TMPDIR=/tmp \\
    #docker://hkubal/clair3:v2.0.1_gpu \\
    #mkdir -p clair3
    #model="r1041_e82_400bps_hac_v520_with_mv"
    #run_clair3.sh \
    #    --bam_fn=${bam} \
    #    --ref_fn=${ref} \
    #    --threads=${task.cpus} \
    #    --platform=ont \
    #    --model_path=/opt/models/\$model \
    #    --output=clair3 \
    #    --include_all_ctgs \
    #    --use_gpu
    #test -s clair3/merge_output.vcf.gz
    #test -s clair3/merge_output.vcf.gz.tbi
    """
}