#!/bin/bash
#SBATCH --job-name=annotate_tables
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=logs/annotate_tables_%j.out
#SBATCH --error=logs/annotate_tables_%j.err

set -euo pipefail

module purge
module load cesga/2020
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/05_annotated_tables

# Find GenBank file (check common locations)
GB_FILE=""
if [ -f "NC_045512.2/NC_045512.2.gb" ]; then
    GB_FILE="NC_045512.2/NC_045512.2.gb"
elif [ -f "NC_045512.2.gb" ]; then
    GB_FILE="NC_045512.2.gb"
elif [ -f "MN908947.3/MN908947.3.gb" ]; then
    GB_FILE="MN908947.3/MN908947.3.gb"
elif [ -f "MN908947.3.gb" ]; then
    GB_FILE="MN908947.3.gb"
else
    echo "Error: GenBank file (.gb) not found. Please specify with --gb" >&2
    exit 1
fi

echo "Using GenBank file: $GB_FILE"

# Execute annotation script
python3 annotate_patient_tables.py \
    --tsv-dir processing/04_position_time_matrices \
    --gb "$GB_FILE"

# Move annotated files (_genes.tsv) to output directory
python3 - << 'EOF'
from pathlib import Path
import shutil

input_dir = Path("processing/04_position_time_matrices")
output_dir = Path("processing/05_annotated_tables")

# Find all _genes.tsv files
annotated_files = sorted(list(input_dir.glob("*_genes.tsv")))

print(f"Moving {len(annotated_files)} annotated files to {output_dir}...")

moved = 0
for annotated_file in annotated_files:
    dest = output_dir / annotated_file.name
    shutil.move(str(annotated_file), str(dest))
    moved += 1

print(f"✔ Moved {moved} files to {output_dir}")
print(f"✓ Annotated tables ready in processing/05_annotated_tables/")
EOF

