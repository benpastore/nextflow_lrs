// lbcb-sci/herro (https://github.com/lbcb-sci/herro) -- GPU deep-learning
// ONT self-correction. NOT the same tool as the old hybrid Illumina+ONT
// "HERO" (kangxiongbin/HERO) this replaced -- herro never touches Illumina
// data, it corrects ONT reads via all-vs-all self-alignment + neural-net
// inference, so there's no illumina_design-style per-sample gating here.

process HERRO_PREPROCESS {

    tag "$sampleID"
    label 'herro_preprocess'

    input:
        tuple val(sampleID), val(unal_bam), val(fastq)

    output:
        tuple val(sampleID), path("${sampleID}.herro_preprocessed.fastq.gz"), emit: preprocessed_ch

    script:
    // "parts to split job into" -- herro's own knob for bounding memory on
    // large read sets (chunks the porechop/duplex_tools trimming step).
    def parts = params.herro_preprocess_parts ?: 4
    """
    #!/bin/bash
    set -euo pipefail

    /opt/herro/scripts/preprocess.sh \\
        ${fastq} \\
        ${sampleID}.herro_preprocessed \\
        ${task.cpus} \\
        ${parts}
    """
}

process HERRO_ALIGN_BATCHES {

    tag "$sampleID"
    label 'herro_align'

    input:
        tuple val(sampleID), path(preprocessed_fastq)

    output:
        tuple val(sampleID), path(preprocessed_fastq), path("${sampleID}_alignment_batches"), emit: batches_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    seqkit seq -ni ${preprocessed_fastq} > ${sampleID}.read_ids.txt

    mkdir -p ${sampleID}_alignment_batches
    /opt/herro/scripts/create_batched_alignments.sh \\
        ${preprocessed_fastq} \\
        ${sampleID}.read_ids.txt \\
        ${task.cpus} \\
        ${sampleID}_alignment_batches
    """
}

process HERRO_INFERENCE {

    tag "$sampleID"
    label 'herro_inference'

    input:
        tuple val(sampleID), path(preprocessed_fastq), path(alignment_batches)

    output:
        tuple val(sampleID), path("*.herro.corrected.fastq.gz"), emit: herro_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    herro inference \\
        --read-alns ${alignment_batches} \\
        -t ${task.cpus} \\
        -d ${params.herro_gpu_ids ?: '0'} \\
        -m ${params.herro_model} \\
        -b ${params.herro_batch_size ?: 64} \\
        ${preprocessed_fastq} \\
        ${sampleID}.herro.corrected.fasta

    # herro writes corrected sequence as FASTA (no per-base quality, since
    # inference re-derives consensus rather than keeping raw basecalls) --
    # gzip and keep the .fastq.gz name anyway so basename(...) .fastq.gz
    # calls downstream (hifiasm/minimap2, which both auto-detect FASTA vs
    # FASTQ by content, not extension) keep producing clean sample-name-based
    # output prefixes.
    gzip -c ${sampleID}.herro.corrected.fasta > ${sampleID}.herro.corrected.fastq.gz
    """
}
