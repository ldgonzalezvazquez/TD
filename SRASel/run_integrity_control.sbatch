#!/bin/bash
#SBATCH --job-name=integrity_control
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=5:00:00
#SBATCH --output=logs/integrity_control_%j.out
#SBATCH --error=logs/integrity_control_%j.err

set -euo pipefail

module purge
module load cesga/2020
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/03_integrity_checked

python3 integrity_control.py

