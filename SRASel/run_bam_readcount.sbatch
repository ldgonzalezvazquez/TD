#!/bin/bash
#SBATCH --job-name=bam_readcount
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=fast
#SBATCH --output=logs/bam_readcount_%j.out
#SBATCH --error=logs/bam_readcount_%j.err

set -euo pipefail

module purge
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/02_rc

reference="MN908947.3/MN908947.3.fna"

for bam in processing/01_rg/*.bam; do
  base=$(basename "$bam" ".bam")
  bam-readcount \
    -q 20 \
    -b 20 \
    -w 0 \
    -f "$reference" \
    "$bam" \
    > "processing/02_rc/${base}_readcount.txt"
done

