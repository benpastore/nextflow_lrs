process DIPCALL {

    label 'dipcall'
    tag "${sampleID}"
    
    publishDir "${params.results}/dipcall/${sampleID}", mode: params.publish_mode

    input:
        tuple val(sampleID), val(reads), val(hap1), val(hap2)
        tuple path(ref), path(ref_fai)

    output:
        tuple val(sampleID), path("*.dip.vcf.gz"), path("*.dip.vcf.gz.tbi"),     emit: dipcall_ch
        tuple val(sampleID), path("*.dip.vcf.gz"),     emit: vcf
        tuple val(sampleID), path("*.dip.vcf.gz.tbi"), emit: tbi
        tuple val(sampleID), path("*.dip.bed"),         emit: bed
        tuple val(sampleID), path("*.pair.vcf.gz"),     emit: pair_vcf
        tuple val(sampleID), path("*.hap1.bam"),        emit: hap1_bam
        tuple val(sampleID), path("*.hap2.bam"),        emit: hap2_bam
        path "versions.yml",                        emit: versions

    script:
    def prefix   = sampleID
    par_flag = params.par_opt ? params.par_opt : ""
    """
    # bgzip inputs if not already compressed (dipcall expects .gz)
    HAP1=${hap1}
    HAP2=${hap2}
    if [[ "${hap1}" != *.gz ]]; then
        bgzip -c ${hap1} > hap1.fa.gz && HAP1=hap1.fa.gz
    fi
    if [[ "${hap2}" != *.gz ]]; then
        bgzip -c ${hap2} > hap2.fa.gz && HAP2=hap2.fa.gz
    fi

    # Index if missing
    [ -f "\${HAP1}.fai" ] || samtools faidx "\${HAP1}"
    [ -f "\${HAP2}.fai" ] || samtools faidx "\${HAP2}"

    # Generate the Makefile
    run-dipcall \\
        ${par_flag} \\
        ${prefix} \\
        ${ref} \\
        \${HAP1} \\
        \${HAP2} \\
        > ${prefix}.mak

    # Execute the pipeline
    make -j ${task.cpus} -f ${prefix}.mak

    # Tabix-index the output VCF if not already done
    [ -f "${prefix}.dip.vcf.gz.tbi" ] || tabix -p vcf ${prefix}.dip.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dipcall: \$(run-dipcall 2>&1 | head -1 | sed 's/.*dipcall-//;s/ .*//')
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}