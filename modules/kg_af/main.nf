// 1000 Genomes allele frequencies native to T2T-CHM13: the T2T consortium's
// recall of all 3,202 high-coverage 1KGP samples against CHM13 (Aganezov et
// al. 2022), reduced once to a small sites-only VCF whose INFO/AF is
// recomputed over the unrelated EUR samples (comparable to gnomAD AF_nfe).
// Read by the modifier search (kg_AF_EUR column, see GnomadNFEAnnotator in
// bin/_modifier_common.py) alongside the Ensembl gnomAD projection, which
// misses sites where CHM13 carries the minor allele. "Build once, reuse
// forever" -- published under kg_af_dir, which the kg_af subworkflow
// (subworkflows/ont.nf) checks first.

process KG_AF_SITES {

    // one task per per-chromosome 1KGP VCF
    tag "${vcf.name}"
    label 'kg_af'

    input:
        path(vcf)
        path(panel)

    output:
        path("${stem}.EUR.sites.vcf.gz"), emit: sites_ch

    script:
    // not simpleName: per-chromosome files share everything before the
    // first dot (e.g. 1kgp.chr1.recalibrated...vcf.gz)
    stem = vcf.name.replaceAll(/\.vcf(\.gz)?$/, '')
    """
    #!/bin/bash
    set -euo pipefail

    # unrelated EUR samples: no parents listed in the 1KGP pedigree/population
    # panel (drops the trio children), columns located by header name
    awk 'NR == 1 { for (i = 1; i <= NF; i++) c[\$i] = i
                   if (!c["SampleID"] || !c["Superpopulation"] || !c["FatherID"] || !c["MotherID"]) {
                       print "ERROR: panel needs SampleID/FatherID/MotherID/Superpopulation columns" > "/dev/stderr"; exit 1 }
                   next }
         \$c["Superpopulation"] == "EUR" && \$c["FatherID"] == "0" && \$c["MotherID"] == "0" { print \$c["SampleID"] }' \\
        ${panel} > eur_samples.txt
    echo "[KG_AF_SITES] \$(wc -l < eur_samples.txt) unrelated EUR sample(s) in the panel" >&2

    # subset -> one ALT per record -> AF/AC/AN over the subset -> sites only
    bcftools view --threads ${task.cpus} -S eur_samples.txt --force-samples -Ou ${vcf} \\
        | bcftools norm -m- -Ou \\
        | bcftools +fill-tags -Ou -- -t AF,AC,AN \\
        | bcftools annotate -x ^INFO/AF,INFO/AC,INFO/AN -Ou \\
        | bcftools view -G --threads ${task.cpus} -Oz -o ${stem}.EUR.sites.vcf.gz
    """
}

process KG_AF_MERGE {

    label 'kg_af'
    publishDir "${params.kg_af_dir ?: "${params.results}/00_resources/1kg_chm13"}", mode: params.publish_mode

    input:
        path(sites)

    output:
        tuple path("1kgp.chm13.EUR.sites.vcf.gz"), path("1kgp.chm13.EUR.sites.vcf.gz.tbi"), emit: kg_af_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # natural sort keeps chromosomes contiguous (chr1, chr2, ..., chr10, ...)
    bcftools concat --threads ${task.cpus} -Oz -o 1kgp.chm13.EUR.sites.vcf.gz \$(ls -v *.EUR.sites.vcf.gz)
    bcftools index -t 1kgp.chm13.EUR.sites.vcf.gz
    """
}
