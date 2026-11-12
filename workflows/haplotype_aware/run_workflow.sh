#!/bin/bash

# --run_trio_analysis + --family_json (mild/severe phenotypes) enable
# SIB_HAPLOTYPE_SHARING and HAPLOTYPE_MODIFIER; --run_snpeff is required by
# the modifier search.
sh run.sh \
    --design /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family_A_design.txt \
    --family_json /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family.json \
    --run_trio_analysis true \
    --run_snpeff true \
    --results /fs/scratch/PAS0631/00_HAPLOTYPE_AWARE_ANALYSIS/FAMILY_A_1022 \
    --run_alphagenome true \
    -resume
