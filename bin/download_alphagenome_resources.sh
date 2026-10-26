#!/bin/bash
# One-time download of the hg38 annotation/calibration files AlphaGenome
# otherwise fetches from Google Cloud Storage on every run (see
# RESOURCE_FILES in bin/run_alphagenome.py -- filenames must match). Run from
# a login/data-transfer node, then point params.alphagenome_resources at the
# output directory. ~390MB total. The FASTA isn't fetched: the pipeline's own
# params.hg38_genome is used instead.
#
#   bash bin/download_alphagenome_resources.sh /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/models/alphagenome_resources
set -euo pipefail

outdir=${1:?usage: $0 <output_dir>}
mkdir -p "$outdir"

base=https://storage.googleapis.com/alphagenome
for path in \
    reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather \
    reference/exon/hg38/polyadb_human_v3_exon3_contiguous_gtfv46.feather \
    reference/gencode/hg38/gencode.v46.splice_sites_starts.feather \
    reference/gencode/hg38/gencode.v46.splice_sites_ends.feather \
    data/hg38/calibration_scores.pb
do
    fname=$(basename "$path")
    echo "[download_alphagenome_resources] $fname"
    # download to .part first so an interrupted fetch never leaves a
    # truncated file that run_alphagenome.py would pick up as valid
    wget -q -O "$outdir/$fname.part" "$base/$path"
    mv "$outdir/$fname.part" "$outdir/$fname"
done

ls -lh "$outdir"
