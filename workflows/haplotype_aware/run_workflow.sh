#!/bin/bash

# --run_trio_analysis + --family_json (mild/severe phenotypes) enable family
# phasing, parent of origin, hap-ibd sibling IBD and the mild-copy variant
# table (needs --hapibd_jar or params.hapibd_jar); --run_methylation adds
# the mild vs severe promoter methylation comparison.
sh run.sh \
    --design /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family_A_design.txt \
    --family_json /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family.json \
    --run_trio_analysis true \
    --run_snpeff true \
    --run_methylation false \
    --results /fs/scratch/PAS0631/00_HAPLOTYPE_AWARE_ANALYSIS/FAMILY_A_1022_November \
    -resume
