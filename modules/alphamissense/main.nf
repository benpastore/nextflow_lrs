// AlphaMissense missense-pathogenicity lookup for the recessive/dominant-
// modifier candidates that land in protein-coding sequence (see
// bin/run_alphamissense.py). Scores the same candidate set as AlphaGenome
// (same gnomAD AF_nfe pre-filter, params.alphagenome_max_af) but runs
// independently: it lifts its own copy and never touches the alphagenome
// processes, so turning it on doesn't invalidate their cache.
// Four steps: index the AlphaMissense table once (PREPARE_ALPHAMISSENSE_INDEX,
// skipped when an index already exists), union the modifier TSVs into a VCF
// (BUILD_ALPHAMISSENSE_INPUT), lift the SNVs CHM13 -> hg38
// (LIFTOVER_ALPHAMISSENSE_VARIANTS), then look them up (ALPHAMISSENSE).

process PREPARE_ALPHAMISSENSE_INDEX {

    // "Build once, reuse forever" (same pattern as SNPEFF_BUILD): published
    // under alphamissense_index_dir, which the alphamissense subworkflow
    // (subworkflows/ont.nf) checks first, so this only runs on the first run.
    label 'alphamissense_index'
    publishDir "${params.alphamissense_index_dir ?: "${params.results}/00_resources/alphamissense"}", mode: params.publish_mode

    input:
        path(alphamissense_tsv)
        path(script)

    output:
        tuple path("*.tabix.tsv.gz"), path("*.tabix.tsv.gz.tbi"), emit: index_ch

    script:
    def stem = alphamissense_tsv.name.replaceAll(/\.tsv(\.gz)?$/, '')
    """
    #!/bin/bash
    set -euo pipefail

    python3 ${script} -input ${alphamissense_tsv} -output ${stem}.tabix.tsv
    """
}

process BUILD_ALPHAMISSENSE_INPUT {

    // same union + AF filter as BUILD_ALPHAGENOME_INPUT (same script, same
    // max_af), so both tools score the same candidates
    label 'low'
    publishDir "${params.results}/08_annotation/alphamissense", mode: params.publish_mode

    input:
        path(recessive_tsv)
        path(dominant_tsv)
        path(script)

    output:
        path("candidate_variants.chm13.vcf"), emit: candidate_vcf_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    set +u
    source activate rnaseq
    set -u

    python3 ${script} \\
        -recessive_tsv ${recessive_tsv} \\
        -dominant_tsv ${dominant_tsv} \\
        ${params.alphagenome_max_af != null && params.alphagenome_max_af != '' ? "-max_af ${params.alphagenome_max_af}" : ''} \\
        -output candidate_variants.chm13.vcf
    """
}

process LIFTOVER_ALPHAMISSENSE_VARIANTS {

    // SNVs only (AlphaMissense scores nothing else). Plain CrossMap vcf:
    // the Fail(REF==ALT) sites that LIFTOVER_ALPHAGENOME_VARIANTS recovers
    // by flipping are left in the .unmap here -- hg38's reference residue
    // there IS the candidate allele, so there's no missense change to score
    // (run_alphamissense.py reports them as hg38_ref_is_candidate).
    label 'crossmap'
    publishDir "${params.results}/08_annotation/alphamissense", mode: params.publish_mode

    input:
        path(vcf)
        path(chain)
        tuple path(hg38_genome), path(hg38_genome_fai)

    output:
        path("candidate_snvs.hg38.vcf"), emit: lifted_vcf_ch
        path("candidate_snvs.hg38.vcf.unmap"), emit: unmapped_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    awk -F'\\t' '/^#/ || (length(\$4)==1 && length(\$5)==1)' ${vcf} > candidate_snvs.chm13.vcf
    CrossMap vcf ${chain} candidate_snvs.chm13.vcf ${hg38_genome} candidate_snvs.hg38.vcf
    touch candidate_snvs.hg38.vcf.unmap

    n_in=\$(grep -vc '^#' candidate_snvs.chm13.vcf || true)
    n_unmapped=\$(grep -vc '^#' candidate_snvs.hg38.vcf.unmap || true)
    echo "[LIFTOVER_ALPHAMISSENSE_VARIANTS] \${n_unmapped}/\${n_in} SNV(s) not lifted CHM13 -> hg38 (see candidate_snvs.hg38.vcf.unmap)" >&2
    """
}

process ALPHAMISSENSE {

    // crossmap container: has pysam for the tabix lookup
    label 'crossmap'
    publishDir "${params.results}/08_annotation/alphamissense", mode: params.publish_mode

    input:
        path(candidate_vcf, stageAs: 'candidate/*')
        path(hg38_vcf, stageAs: 'hg38/*')
        path(unmap, stageAs: 'unmap/*')
        tuple path(alphamissense_tsv), path(alphamissense_tbi)
        path(gene_tsv, stageAs: 'gene/*')
        path(script)

    output:
        path("cohort.alphamissense.tsv"), emit: alphamissense_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    python3 ${script} \\
        -candidate_vcf ${candidate_vcf} \\
        -hg38_vcf ${hg38_vcf} \\
        -unmap ${unmap} \\
        -alphamissense ${alphamissense_tsv} \\
        -gene_tsv ${gene_tsv} \\
        -output cohort.alphamissense.tsv
    """
}
