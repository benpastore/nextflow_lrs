#!/bin/bash
#SBATCH --account=PAS0631
#SBATCH --job-name=gnomad_chm13_dl
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --output=%x_%j.log

set -euo pipefail

cd /fs/ess/PAS0631/BEN_PIPELINE/nextflow_lrs/gnomAD

URL="https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz"
EXPECTED_BYTES=46053099774

wget -c "$URL"

actual_bytes=$(stat -c%s Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz)
if [ "$actual_bytes" -ne "$EXPECTED_BYTES" ]; then
    echo "ERROR: downloaded ${actual_bytes} bytes, expected ${EXPECTED_BYTES} -- download incomplete/corrupt" >&2
    exit 1
fi

gzip -t Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz
echo "Download verified: correct size and valid gzip stream."
