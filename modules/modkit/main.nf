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
    // do --phased (HP) on real methylation data.
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

    publishDir "${params.results}/07_methylation/modkit", mode: params.publish_mode

    // Reference-coordinate methylation pileup, split by haplotype. Runs on
    // TRANSFER_HP_TAGS's output -- the RAW (pre-herro) read alignment,
    // haplotagged with HP from the herro-corrected/phased bam -- so
    // --phased (modkit >= 0.6; replaces --partition-tag HP) gives
    // per-haplotype bedMethyl output (HP=1, HP=2, and a combined file over
    // ALL reads, phased or not). Depends on MM/ML tags surviving the bam -> fastq -> bam
    // round trip (see SAMTOOLS_CONVERT_BAM_TO_FASTQ's -T flag and
    // MINIMAP2_ALIGN_METHYLATION's -y flag); if those were ever removed,
    // this would start silently producing empty/no-modification pileups.
    input:
        tuple val(sampleID), val(bam), val(bai)
        tuple val(ref), val(ref_fai)

    //
    // CpG-only, strands combined (--cpg --combine-strands): one row per CpG
    // per mod code. With --phased, modkit writes a DIRECTORY with
    // hp1.bedmethyl, hp2.bedmethyl and combined.bedmethyl (all reads -- NOT
    // just the unphased ones, unlike the old "ungrouped" file); they're
    // renamed to fixed names here so downstream (METHYLATION_COMPARE) can
    // rely on them. A missing file becomes an empty one.
    output:
        tuple val(sampleID),
            path("${sampleID}.modkit.hp1.bed.gz"),
            path("${sampleID}.modkit.hp2.bed.gz"),
            path("${sampleID}.modkit.combined.bed.gz"), emit: modkit_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    modkit pileup \\
        ${bam} \\
        modkit_out \\
        --ref ${ref} \\
        --cpg \\
        --combine-strands \\
        --phased \\
        --threads ${task.cpus} \\
        --log-filepath ${sampleID}.modkit.log

    if [ ! -d modkit_out ]; then
        echo "modkit did not write a phased output directory (modkit_out)" >&2
        exit 1
    fi
    ls -l modkit_out >&2
    for out in hp1 hp2 combined; do
        f=\$(ls modkit_out/*\${out}*.bedmethyl* 2>/dev/null | head -n 1 || true)
        if [ -n "\$f" ]; then
            case "\$f" in
                *.gz) cp "\$f" ${sampleID}.modkit.\$out.bed.gz ;;
                *)    gzip -c "\$f" > ${sampleID}.modkit.\$out.bed.gz ;;
            esac
        else
            echo "warning: no \$out bedmethyl in modkit_out, writing an empty file" >&2
            printf '' | gzip -c > ${sampleID}.modkit.\$out.bed.gz
        fi
    done
    """
}
