// Post-hoc haplotype annotation for callers with no native phasing mode
// (see bin/annotate_haplotype.py for why: Straglr and Spectre don't read the
// HP:i:1/HP:i:2 tags already present on the haplotagged BAM they're called
// against). Clair3 (phased via longphase/longphase_sv) and Sniffles (--phase)
// need no equivalent step -- they already carry real phase.

process HAPLOTYPE_ANNOTATE_STRAGLR {

    tag "$sampleID"
    // needs python3 + samtools (annotate_haplotype.py shells out to
    // `samtools view -d HP:...`) -- the plain samtools biocontainer doesn't
    // reliably ship python3, so use the general-purpose env other bin/
    // python scripts in this pipeline use (see CONSOLIDATE_VARIANTS,
    // modules/consolidate_variants/main.nf), which BAM_COMPARE
    // (modules/deeptools/main.nf) already confirms also has samtools.
    label 'low'

    publishDir "${params.results}/06_variants/straglr", mode: params.publish_mode

    input:
        tuple val(sampleID), path(tsv), path(bed), path(bam), path(bai)
        path(script)

    output:
        tuple val(sampleID), path("*.straglr.haplotagged.tsv"), emit: straglr_haplotagged_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    name=\$(basename ${tsv} .straglr.tsv)

    python3 ${script} \\
        --bam ${bam} \\
        --calls ${tsv} \\
        --mode straglr \\
        --output \${name}.straglr.haplotagged.tsv \\
        --flank ${params.haplotype_flank_bp ?: 1000}
    """
}

process HAPLOTYPE_ANNOTATE_SPECTRE {

    tag "$sampleID"
    label 'low' // see HAPLOTYPE_ANNOTATE_STRAGLR above

    input:
        tuple val(sampleID), path(bed_gz), path(bed_gz_tbi), path(bam), path(bai)
        path(script)

    output:
        tuple val(sampleID), path("*.spectre.haplotagged.bed"), emit: spectre_haplotagged_raw_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    name=\$(basename ${bed_gz} .spectre.bed.gz)

    zcat ${bed_gz} > \${name}.spectre.bed

    python3 ${script} \\
        --bam ${bam} \\
        --calls \${name}.spectre.bed \\
        --mode spectre \\
        --output \${name}.spectre.haplotagged.bed \\
        --flank ${params.haplotype_flank_bp ?: 1000}
    """
}

process INDEX_SPECTRE_HAPLOTAGGED_BED {

    tag "$sampleID"
    // the rnaseq container's bgzip is broken (missing libcrypto.so.1.0.0),
    // same underlying issue INDEX_SNIFFLES_VCF works around -- compress +
    // index under the bcftools container instead.
    label "bcftools"
    publishDir "${params.results}/06_variants/spectre", mode: params.publish_mode

    input:
        tuple val(sampleID), path(bed)

    output:
        tuple val(sampleID), path("*.spectre.haplotagged.bed.gz"), path("*.spectre.haplotagged.bed.gz.tbi"), emit: spectre_haplotagged_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    bgzip -f ${bed}
    tabix -f -0 -s 1 -b 2 -e 3 ${bed}.gz
    """
}
