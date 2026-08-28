process TRANSFER_HP_TAGS {

    tag "$sampleID"
    label 'samtools_high'

    // herro-corrected reads carry no real methylation tags (see
    // MINIMAP2_ALIGN's comment) and are also rewritten/consensus-corrected,
    // so real MM/ML data can only come from a separate alignment of the
    // RAW reads (MINIMAP2_ALIGN_METHYLATION). But haplotype phasing runs
    // on the herro-corrected alignment for accuracy, so that raw bam has
    // no HP tags of its own. This projects HP from the phased
    // herro-corrected bam onto the raw bam by read ID, so modkit can still
    // do --partition-tag HP on real methylation data.
    //
    // herro sometimes splits one raw read into multiple corrected reads
    // (IDs suffixed :0, :1, ...) -- stripping that suffix maps each
    // fragment back to its raw parent read. If a read's fragments ever
    // disagree on HP (split exactly at a phase-switch point), whichever
    // fragment is encountered first in the bam wins; this is rare and
    // only affects that one read's haplotype label, never its
    // methylation call.
    input:
        tuple val(sampleID), path(phased_bam), path(phased_bai), path(raw_bam), path(raw_bai)

    output:
        tuple val(sampleID), path("${sampleID}.raw.haplotagged.bam"), path("${sampleID}.raw.haplotagged.bam.bai"), emit: haplotagged_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    samtools view ${phased_bam} \\
      | awk 'BEGIN{FS=OFS="\\t"} {
            n = split(\$1, a, ":"); base = a[1];
            if (base in seen) next;
            for (i = 12; i <= NF; i++) {
                if (\$i ~ /^HP:i:/) { print base, \$i; seen[base] = 1; next }
            }
        }' > read_hp_map.tsv

    samtools view -h -@ ${task.cpus} ${raw_bam} \\
      | awk -v OFS='\\t' -v map=read_hp_map.tsv '
            BEGIN {
                while ((getline line < map) > 0) {
                    split(line, a, "\\t"); hp[a[1]] = a[2]
                }
            }
            /^@/ { print; next }
            { if (\$1 in hp) print \$0, hp[\$1]; else print \$0 }
        ' \\
      | samtools view -b -@ ${task.cpus} -o ${sampleID}.raw.haplotagged.bam -

    samtools index -@ ${task.cpus} ${sampleID}.raw.haplotagged.bam
    """
}

process MODKIT_PILEUP {

    tag "$sampleID"
    label 'modkit'

    publishDir "${params.results}/methylation/modkit", mode: params.publish_mode

    // Reference-coordinate methylation pileup, split by haplotype. Runs on
    // TRANSFER_HP_TAGS's output -- the RAW (pre-herro) read alignment,
    // haplotagged with HP from the herro-corrected/phased bam -- so
    // --partition-tag HP gives per-haplotype bedMethyl output (haplotype
    // 1, haplotype 2, and an ungrouped file for reads that couldn't be
    // phased). Depends on MM/ML tags surviving the bam -> fastq -> bam
    // round trip (see SAMTOOLS_CONVERT_BAM_TO_FASTQ's -T flag and
    // MINIMAP2_ALIGN_METHYLATION's -y flag); if those were ever removed,
    // this would start silently producing empty/no-modification pileups.
    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref), val(ref_fai)

    output:
        tuple val(sampleID), path("${sampleID}*.bed.gz"), emit: modkit_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    modkit pileup \\
        ${bam} \\
        ${sampleID}.modkit.bed \\
        --ref ${ref} \\
        --partition-tag HP \\
        --threads ${task.cpus} \\
        --log-filepath ${sampleID}.modkit.log

    for f in ${sampleID}.modkit*.bed; do
        bgzip "\$f"
    done
    """
}
