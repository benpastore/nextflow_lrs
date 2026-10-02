// AlphaGenome functional scoring of the recessive/dominant-modifier
// candidate variants (see modules/recessive_modifier/main.nf) -- for this
// SMA study, restricted to motor-neuron-relevant tissue (see
// bin/run_alphagenome.py). Three steps, chained: union the two modifier
// TSVs into one small CHM13-coordinate VCF (BUILD_ALPHAGENOME_INPUT),
// liftover to hg38 since AlphaGenome is trained on that build
// (LIFTOVER_ALPHAGENOME_VARIANTS, via CrossMap), then score locally against
// a pre-downloaded Kaggle AlphaGenome weights directory (ALPHAGENOME).
// Single process chain for the whole cohort, not per-sample, matching
// RECESSIVE_MODIFIER's own shape -- the modifier TSVs are already
// cohort-level.

process FETCH_CHM13_TO_HG38_CHAIN {

    // Only runs when params.chm13_to_hg38_chain is unset -- fetches
    // UCSC's own official CHM13v2(hs1)->hg38 liftOver chain rather than
    // generating one from scratch via alignment (UCSC already did that
    // carefully; re-deriving it would be slower and more error-prone for
    // no benefit). Cached under results permanently ("build once, reuse
    // forever", same pattern as SNPEFF_BUILD / the snpeff_db_exists check
    // in the snpeff subworkflow, subworkflows/ont.nf) so subsequent runs
    // don't redownload.
    //
    // This runs as a Slurm job, i.e. on a compute node -- if your cluster
    // blocks outbound internet from compute nodes (common), this will fail
    // loudly rather than hang. If it does, fetch the file once yourself
    // from a login/data-transfer node:
    //   wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/chm13v2-hg38.over.chain.gz
    // and point params.chm13_to_hg38_chain at it to skip this step entirely.
    label 'low'
    publishDir "${params.results}/00_resources/liftover", mode: params.publish_mode

    output:
        path("chm13v2-hg38.over.chain.gz"), emit: chain_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    wget -q -O chm13v2-hg38.over.chain.gz \\
        "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/chm13v2-hg38.over.chain.gz" \\
        || { echo "ERROR: failed to download the liftover chain -- if your cluster blocks internet on compute nodes, fetch it manually from a login node instead (see comment above this process in modules/alphagenome/main.nf) and set params.chm13_to_hg38_chain" >&2; exit 1; }

    # sanity-check: a truncated download or an HTML error page saved as
    # .gz would otherwise be cached and silently reused on every future run.
    size=\$(stat -c%s chm13v2-hg38.over.chain.gz 2>/dev/null || stat -f%z chm13v2-hg38.over.chain.gz)
    if [ "\$size" -lt 100000 ]; then
        echo "ERROR: downloaded chain file is suspiciously small (\${size} bytes) -- likely a failed/partial download" >&2
        exit 1
    fi
    zcat chm13v2-hg38.over.chain.gz | head -1 | grep -q '^chain ' || {
        echo "ERROR: downloaded file doesn't look like a chain file (first line isn't a 'chain' record)" >&2
        exit 1
    }
    """
}

process BUILD_ALPHAGENOME_INPUT {

    label 'low'
    publishDir "${params.results}/08_annotation/alphagenome", mode: params.publish_mode

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
        -output candidate_variants.chm13.vcf
    """
}

process LIFTOVER_ALPHAGENOME_VARIANTS {

    label 'crossmap'
    publishDir "${params.results}/08_annotation/alphagenome", mode: params.publish_mode

    input:
        path(vcf)
        path(chain)
        tuple path(hg38_genome), path(hg38_genome_fai)

    output:
        path("candidate_variants.hg38.vcf"), emit: lifted_vcf_ch
        path("candidate_variants.hg38.vcf.unmap"), optional: true, emit: unmapped_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    CrossMap vcf ${chain} ${vcf} ${hg38_genome} candidate_variants.hg38.vcf

    n_unmapped=0
    if [ -f candidate_variants.hg38.vcf.unmap ]; then
        n_unmapped=\$(grep -vc '^#' candidate_variants.hg38.vcf.unmap || true)
    fi
    echo "[LIFTOVER_ALPHAGENOME_VARIANTS] \${n_unmapped} variant(s) failed to lift CHM13 -> hg38 (see candidate_variants.hg38.vcf.unmap)" >&2
    """
}

process ALPHAGENOME {

    label 'alphagenome'
    publishDir "${params.results}/08_annotation/alphagenome", mode: params.publish_mode

    input:
        path(candidate_vcf)
        path(vcf)
        tuple path(hg38_genome), path(hg38_genome_fai)
        path(script)

    output:
        path("cohort.alphagenome.tsv"), emit: alphagenome_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    python3 ${script} \\
        -candidate_vcf ${candidate_vcf} \\
        -vcf ${vcf} \\
        -reference ${hg38_genome} \\
        -weights ${params.alphagenome_weights} \\
        -quantile_threshold ${params.alphagenome_quantile_threshold ?: 0.9} \\
        ${params.alphagenome_args ?: ''} \\
        -output cohort.alphagenome.tsv
    """
}
