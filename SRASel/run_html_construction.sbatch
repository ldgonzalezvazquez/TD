#!/bin/bash
#SBATCH --job-name=html_construction
#SBATCH --cpus-per-task=48
#SBATCH --mem=30G
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --output=logs/html_construction_%j.out
#SBATCH --error=logs/html_construction_%j.err

set -euo pipefail

module purge
module load cesga/2020
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/07_html

# Find GenBank file (check common locations)
GB_FILE=""
if [ -f "NC_045512.2.gb" ]; then
    GB_FILE="NC_045512.2.gb"
elif [ -f "NC_045512.2/NC_045512.2.gb" ]; then
    GB_FILE="NC_045512.2/NC_045512.2.gb"
elif [ -f "MN908947.3/MN908947.3.gb" ]; then
    GB_FILE="MN908947.3/MN908947.3.gb"
else
    echo "Error: GenBank file (.gb) not found. Please specify with --gb" >&2
    exit 1
fi

echo "Using GenBank file: $GB_FILE"
echo "Input directory: processing/06_selection_analysis"

echo "Copying files with gene information ..."
cp processing/05_annotated_tables/*_genes.tsv processing/06_selection_analysis/

echo "Output HTML: processing/07_html/selection_UPGMA.html"

# Execute HTML construction script
python3 Selection_HTML_construction.py \
    --dir processing/06_selection_analysis \
    --gb "$GB_FILE" \
    --out processing/07_html/selection_UPGMA.html \
    --workers 48 \
    --min 1

echo "✓ HTML construction complete"
echo "  Results in: processing/07_html/"

