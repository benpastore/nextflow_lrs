process CLAIR3 {

    label 'clair3'
    tag "${sampleID}"
    publishDir "${params.results}/variants/clair3", mode: params.publish_mode

    input:
        val ref
        tuple val(sampleID), val(bam), val(bai)

    output:
        tuple val(sampleID), path("${sampleID}.clair3.vcf.gz"), path("${sampleID}.clair3.vcf.gz.tbi"), emit : clair3_ch

    script:
    """
    mkdir -p clair3

    model="r1041_e82_400bps_sup_v500"
    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --threads=${task.cpus} \
        --platform=ont \
        --model_path=/opt/models/\$model \
        --output=clair3 \
        --use_gpu

    test -s clair3/merge_output.vcf.gz
    test -s clair3/merge_output.vcf.gz.tbi

    mv clair3/merge_output.vcf.gz ${sampleID}.clair3.vcf.gz
    mv clair3/merge_output.vcf.gz.tbi ${sampleID}.clair3.vcf.gz.tbi
    """
}