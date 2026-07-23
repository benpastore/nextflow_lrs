process GENOME_STATS {

    tag "$sampleID"
    label 'samtools'
    errorStrategy 'ignore'

    publishDir "${params.results}/genome_stats", mode: params.publish_mode

    input:
        tuple val(sampleID), val(reads), val(hap1_fa), val(hap2_fa)

    output:
        tuple val(sampleID), path("*.genome_stats.tsv"), emit: genome_stats_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # number of reads and bases in the read set that went into assembly
    read n_reads n_bases <<< \$(zcat -f ${reads} | awk 'NR % 4 == 2 { n++; b += length(\$0) } END { print n, b }')

    # assembled genome size = sum of contig lengths across both haplotypes
    samtools faidx ${hap1_fa}
    samtools faidx ${hap2_fa}

    hap1_size=\$(awk '{ sum += \$2 } END { print sum + 0 }' ${hap1_fa}.fai)
    hap2_size=\$(awk '{ sum += \$2 } END { print sum + 0 }' ${hap2_fa}.fai)
    asm_size=\$(( hap1_size + hap2_size ))

    coverage=\$(awk -v b="\$n_bases" -v g="\$asm_size" 'BEGIN { if (g > 0) printf "%.2f", b / g; else print "NA" }')

    {
        printf "sample\\tnum_reads\\tnum_bases\\thap1_size\\thap2_size\\tassembled_genome_size\\tcoverage_vs_assembly\\n"
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${sampleID}" "\$n_reads" "\$n_bases" "\$hap1_size" "\$hap2_size" "\$asm_size" "\$coverage"
    } > ${sampleID}.genome_stats.tsv
    """
}

process CHROM_COVERAGE {

    tag "$sampleID"
    label 'samtools'
    errorStrategy 'ignore'

    publishDir "${params.results}/genome_stats", mode: params.publish_mode

    input:
        tuple val(sampleID), path(hap1_bam), path(hap1_bai), path(hap2_bam), path(hap2_bai)

    output:
        tuple val(sampleID), path("*.chrom_coverage.tsv"), emit: chrom_coverage_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    {
        samtools coverage ${hap1_bam} | head -n1 | awk -v OFS='\\t' '{print "sample","hap",\$0}'
        samtools coverage ${hap1_bam} | tail -n +2 | awk -v OFS='\\t' -v s="${sampleID}" -v h="hap1" '{print s,h,\$0}'
        samtools coverage ${hap2_bam} | tail -n +2 | awk -v OFS='\\t' -v s="${sampleID}" -v h="hap2" '{print s,h,\$0}'
    } > ${sampleID}.chrom_coverage.tsv
    """
}
