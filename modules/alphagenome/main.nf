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
    # pipefail off in a subshell: head exits after one line, so zcat gets
    # SIGPIPE (141) and would otherwise fail this check on a valid file.
    if ! (set +o pipefail; zcat chm13v2-hg38.over.chain.gz | grep -v '^#' | head -n 1 | grep -q '^chain '); then
        echo "ERROR: downloaded file doesn't look like a chain file (first non-comment line isn't a 'chain' record)" >&2
        exit 1
    fi
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
        ${params.alphagenome_max_af != null && params.alphagenome_max_af != '' ? "-max_af ${params.alphagenome_max_af}" : ''} \\
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
        path(flip_script)

    output:
        path("candidate_variants.hg38.vcf"), emit: lifted_vcf_ch
        path("candidate_variants.hg38.vcf.unmap"), optional: true, emit: unmapped_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    CrossMap vcf ${chain} ${vcf} ${hg38_genome} candidate_variants.hg38.raw.vcf
    touch candidate_variants.hg38.raw.vcf.unmap

    # Fail(REF==ALT): hg38's reference base IS the candidate's ALT (CHM13 and
    # hg38 carry different alleles there). CrossMap 0.7.4 drops these before
    # --no-comp-alleles is even checked, and the .unmap keeps only the CHM13
    # record -- so lift just their positions with bed mode (coordinates only,
    # no allele logic; strand in column 6) and let flip_ref_alt_alleles.py
    # rewrite them as hg38 REF (= candidate allele) > CHM13 REF.
    awk -F'\\t' 'BEGIN{OFS="\\t"} \$NF=="Fail(REF==ALT)" && length(\$4)==1 && length(\$5)==1 {print \$1, \$2-1, \$2, \$3, 0, "+"}' \\
        candidate_variants.hg38.raw.vcf.unmap > ref_eq_alt.chm13.bed
    if [ -s ref_eq_alt.chm13.bed ]; then
        CrossMap bed ${chain} ref_eq_alt.chm13.bed ref_eq_alt.hg38.bed
    fi
    touch ref_eq_alt.hg38.bed

    python3 ${flip_script} \\
        -original_vcf ${vcf} \\
        -lifted_vcf candidate_variants.hg38.raw.vcf \\
        -unmap_in candidate_variants.hg38.raw.vcf.unmap \\
        -lifted_bed ref_eq_alt.hg38.bed \\
        -reference ${hg38_genome} \\
        -output candidate_variants.hg38.vcf \\
        -unmap_out candidate_variants.hg38.vcf.unmap
    rm candidate_variants.hg38.raw.vcf candidate_variants.hg38.raw.vcf.unmap

    n_unmapped=0
    if [ -f candidate_variants.hg38.vcf.unmap ]; then
        n_unmapped=\$(grep -vc '^#' candidate_variants.hg38.vcf.unmap || true)
    fi
    echo "[LIFTOVER_ALPHAGENOME_VARIANTS] \${n_unmapped} variant(s) failed to lift CHM13 -> hg38 (see candidate_variants.hg38.vcf.unmap)" >&2
    """
}

process ALPHAGENOME {

    // One task per chunk of candidates (params.alphagenome_chunk_size), so
    // a cohort's scoring spreads across several GPU jobs instead of one
    // long one that hits the time limit -- the model scores one variant at a
    // time (~1 s/variant on GPU), there's no batching within a process. A
    // failed chunk reruns alone on -resume. Chunk outputs aren't published;
    // MERGE_ALPHAGENOME_CHUNKS joins them into cohort.alphagenome.tsv.
    label 'alphagenome'
    tag "chunk ${chunk_index + 1}/${n_chunks}"

    input:
        tuple val(chunk_index), val(n_chunks), path(candidate_vcf), path(vcf), path(hg38_genome), path(hg38_genome_fai)
        path(script)

    output:
        path("cohort.alphagenome.chunk_*.tsv"), emit: alphagenome_chunk_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # \$HOME is read-only inside the container
    export MPLCONFIGDIR=\$PWD/.matplotlib

    python3 ${script} \\
        -candidate_vcf ${candidate_vcf} \\
        -vcf ${vcf} \\
        -reference ${hg38_genome} \\
        -weights ${params.alphagenome_weights} \\
        ${params.alphagenome_resources ? "-resources ${params.alphagenome_resources}" : ''} \\
        -quantile_threshold ${params.alphagenome_quantile_threshold ?: 0.99} \\
        -chunk_index ${chunk_index} \\
        -n_chunks ${n_chunks} \\
        ${params.alphagenome_args ?: ''} \\
        -output cohort.alphagenome.chunk_${String.format('%05d', chunk_index)}.tsv
    """
}

process MERGE_ALPHAGENOME_CHUNKS {

    // Chunks are contiguous slices of the VID-sorted candidates and their
    // filenames are zero-padded, so concatenating in filename order gives
    // the same row order a single unchunked run would.
    label 'low'
    publishDir "${params.results}/08_annotation/alphagenome", mode: params.publish_mode

    input:
        path(chunks)

    output:
        path("cohort.alphagenome.tsv"), emit: alphagenome_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    files=\$(ls cohort.alphagenome.chunk_*.tsv | sort)
    head -n 1 \$(echo "\$files" | head -n 1) > cohort.alphagenome.tsv
    for f in \$files; do
        tail -n +2 "\$f" >> cohort.alphagenome.tsv
    done
    echo "[MERGE_ALPHAGENOME_CHUNKS] \$(echo "\$files" | wc -l) chunk(s), \$(( \$(wc -l < cohort.alphagenome.tsv) - 1 )) row(s)" >&2
    """
}

process ALPHAGENOME_API {

    // Backup to ALPHAGENOME: same script and inputs, but scores via
    // DeepMind's hosted API instead of local GPU inference. CPU-only, but
    // needs internet from the compute node. Publishes to its own dir so it
    // can run alongside ALPHAGENOME without clobbering its output.
    //
    // API key via a Nextflow secret, so it never lands in .command.sh or
    // the work dir. Set it once on OSC:
    //   nextflow secrets set ALPHAGENOME_API_KEY <your key>
    label 'alphagenome_api'
    publishDir "${params.results}/08_annotation/alphagenome_api", mode: params.publish_mode
    secret 'ALPHAGENOME_API_KEY'

    input:
        path(candidate_vcf)
        path(vcf)
        tuple path(hg38_genome), path(hg38_genome_fai)
        path(script)

    output:
        path("cohort.alphagenome_api.tsv"), emit: alphagenome_api_ch

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # \$HOME is read-only inside the container
    export MPLCONFIGDIR=\$PWD/.matplotlib

    python3 ${script} \\
        -backend api \\
        -candidate_vcf ${candidate_vcf} \\
        -vcf ${vcf} \\
        -reference ${hg38_genome} \\
        -quantile_threshold ${params.alphagenome_quantile_threshold ?: 0.99} \\
        ${params.alphagenome_args ?: ''} \\
        -output cohort.alphagenome_api.tsv
    """
}
