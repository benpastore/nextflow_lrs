#!/bin/bash

sh run.sh \
    --design /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family_A_design.txt \
    --family_json /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family.json \
    --run_trio_analysis true \
    --illumina_design /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/illumina_design.csvt \
    --results /fs/scratch/PAS0631/00_HAPLOTYPE_AWARE_ANALYSIS/FAMILY_A \
    -resume