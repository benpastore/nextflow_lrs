process HERO {

    tag "$sampleID"
    label 'hero'

    publishDir "${params.results}/hero_corrected", mode: params.publish_mode

    input:
        tuple val(sampleID), val(fmlrc2_fasta), val(r1), val(r2)

    output:
        tuple val(sampleID), path("*.hero.corrected.fastq.gz"), emit: hero_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    zcat -f ${r1} ${r2} > combined_short_reads.fastq

    python3 ${params.hero_bin}/HERO.py \\
        -r combined_short_reads.fastq \\
        -lc ${fmlrc2_fasta} \\
        -p \\
        -o ${sampleID}.hero.corrected.fa \\
        -s ${params.hero_chunk_size ?: 30}

    # HERO writes corrected sequence as FASTA (no per-base quality, since
    # the OLC step re-derives consensus sequence rather than keeping raw
    # basecalls) -- gzip and keep the .fastq.gz name anyway so the
    # basename(...) .fastq.gz calls downstream (hifiasm/minimap2, which
    # both auto-detect FASTA vs FASTQ by content, not extension) keep
    # producing clean sample-name-based output prefixes.
    gzip -c ${sampleID}.hero.corrected.fa > ${sampleID}.hero.corrected.fastq.gz
    """
}
