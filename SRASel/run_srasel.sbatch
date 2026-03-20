#!/bin/bash
#SBATCH --job-name=srasel
#SBATCH --cpus-per-task=48
#SBATCH --mem=60G
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --output=logs/srasel_%j.out
#SBATCH --error=logs/srasel_%j.err

set -euo pipefail

module purge
module load cesga/2020
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/06_selection_analysis

# Get all _genes.tsv files
TSV_FILES=$(find processing/05_annotated_tables -name "*_genes.tsv" -type f | sort | tr '\n' ' ')

if [ -z "$TSV_FILES" ]; then
    echo "Error: No *_genes.tsv files found in processing/05_annotated_tables/" >&2
    exit 1
fi

echo "Processing $(echo $TSV_FILES | wc -w) patient files..."

# Execute SRASel.py with all _genes.tsv files
python3 SRASel.py \
    $TSV_FILES \
    --outdir processing/06_selection_analysis \
    --times dates.tsv \
    --nperm 2000 \
    --jobs 48 \
    --freq-delta 0.03 \
    --thr-clear 0.95 \
    --eps-tie 0.05 \
    --first-k-min 10

echo "✓ Selection analysis complete"
echo "  Results in: processing/06_selection_analysis/"

