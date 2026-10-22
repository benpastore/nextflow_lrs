#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --ntasks-per-node=1 --time=24:00:00 --mem=16gb --cluster=ascend --partition=nextgen

nextflow /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/main.nf --design /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family_A_design.txt --family_json /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/workflows/haplotype_aware/family.json --run_trio_analysis true --results /fs/scratch/PAS0631/00_HAPLOTYPE_AWARE_ANALYSIS/FAMILY_A -resume


