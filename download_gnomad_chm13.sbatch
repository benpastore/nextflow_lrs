#!/bin/bash
#SBATCH --account=PCON0160
#SBATCH --job-name=gnomad_chm13_dl
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --output=%x_%j.log

set -euo pipefail

cd /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/gnomAD

module load miniconda3/24.1.2-py310

URL="https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz"

aria2c -x 16 -s 16 -k 100M --file-allocation=none --max-tries=0 --retry-wait=10 \
    -d /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/gnomAD "$URL"
