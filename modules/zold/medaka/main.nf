process MEDAKA_CONSENSUS {

    tag "$sample_id"
    label 'medaka'

    publishDir "${params.outdir}/medaka", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.consensus.fasta"), emit: consensus
    tuple val(sample_id), path("${sample_id}.variants.vcf"), optional: true, emit: vcf
    tuple val(sample_id), path("${sample_id}.consensus.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.consensus.bam.bai"), emit: bai
    tuple val(sample_id), path("${sample_id}_medaka"), emit: medaka_dir
    path "versions.yml", emit: versions

    script:
    def threads = task.cpus ?: 4
    def extra_args = params.medaka_args ?: ''

    """
    set -euo pipefail

    mkdir -p ${sample_id}_medaka

    medaka_consensus \\
        -i ${reads} \\
        -d ${reference} \\
        -o ${sample_id}_medaka \\
        -t ${threads} \\
        ${extra_args}

    if [[ ! -f ${sample_id}_medaka/consensus.fasta ]]; then
        echo "ERROR: Medaka did not produce consensus.fasta" >&2
        exit 1
    fi

    cp ${sample_id}_medaka/consensus.fasta ${sample_id}.consensus.fasta

    if [[ -f ${sample_id}_medaka/round_1.vcf ]]; then
        cp ${sample_id}_medaka/round_1.vcf ${sample_id}.variants.vcf
    elif [[ -f ${sample_id}_medaka/variants.vcf ]]; then
        cp ${sample_id}_medaka/variants.vcf ${sample_id}.variants.vcf
    fi

    minimap2 -a -x asm5 -t ${threads} ${reference} ${sample_id}.consensus.fasta | \\
        samtools sort -@ ${threads} -o ${sample_id}.consensus.bam

    samtools index -@ ${threads} ${sample_id}.consensus.bam

    cat <<-EOF > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>&1 | head -n 1)
        minimap2: \$(minimap2 --version 2>&1 | head -n 1)
        samtools: \$(samtools --version 2>&1 | head -n 1)
    EOF
    """
}