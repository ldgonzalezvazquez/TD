#!/bin/bash
#SBATCH --job-name=add_rg
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=fast
#SBATCH --output=logs/add_rg_%j.out
#SBATCH --error=logs/add_rg_%j.err

set -euo pipefail

module purge
module load cesga/2020 picard/2.25.5

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/01_rg

for bam in processing/ptrim_renombrados/*.bam; do
  base=$(basename "$bam" ".bam")
  java -jar "$EBROOTPICARD/picard.jar" AddOrReplaceReadGroups \
    I="$bam" \
    O="processing/01_rg/${base}_RG.bam" \
    RGID="$base" \
    RGLB="$base" \
    RGPL=ILLUMINA \
    RGPU="$base" \
    RGSM="$base" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
done

