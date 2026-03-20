#!/bin/bash
#SBATCH --job-name=built_tables
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=logs/built_tables_%j.out
#SBATCH --error=logs/built_tables_%j.err

set -euo pipefail

module purge
module load cesga/2020
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bam-readcount

cd /mnt/lustre/scratch/nlsas/home/uvi/cm/lgv

mkdir -p logs
mkdir -p processing/04_position_time_matrices

# Script to process all patients
python3 - << 'EOF'
from pathlib import Path
import subprocess
import sys

input_dir = Path("processing/03_integrity_checked")
output_dir = Path("processing/04_position_time_matrices")
script_path = Path("built_time_mutation_tsv.py")

# Get all patient directories
patient_dirs = sorted([d for d in input_dir.iterdir() if d.is_dir()])

print(f"Processing {len(patient_dirs)} patients...")
print("="*60)

failed = []
success = []

for patient_dir in patient_dirs:
    patient = patient_dir.name
    
    # Check if .time file exists
    time_file = next(patient_dir.glob("*.time"), None)
    if not time_file:
        print(f"⚠️  {patient}: No .time file found, skipping", file=sys.stderr)
        failed.append((patient, "No .time file"))
        continue
    
    # Execute built_time_mutation_tsv.py for this patient
    # The script writes output to the patient directory
    try:
        result = subprocess.run(
            [sys.executable, str(script_path), str(patient_dir)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # The output file is created in the patient directory as {patient}.tsv
        output_file = patient_dir / f"{patient}.tsv"
        
        if output_file.exists():
            # Move to output directory (not copy, to avoid leaving files in 03)
            dest = output_dir / f"{patient}.tsv"
            import shutil
            shutil.move(str(output_file), str(dest))
            success.append(patient)
            print(f"✔  {patient}: {result.stdout.strip()} (moved to {output_dir})")
        else:
            print(f"⚠️  {patient}: Script ran but output file not found", file=sys.stderr)
            failed.append((patient, "Output file not found"))
            
    except subprocess.CalledProcessError as e:
        print(f"❌ {patient}: Error running script", file=sys.stderr)
        print(f"   stdout: {e.stdout}", file=sys.stderr)
        print(f"   stderr: {e.stderr}", file=sys.stderr)
        failed.append((patient, f"Script error: {e.returncode}"))
    except Exception as e:
        print(f"❌ {patient}: Unexpected error: {e}", file=sys.stderr)
        failed.append((patient, str(e)))

print()
print("="*60)
print(f"Summary:")
print(f"  Total patients processed: {len(patient_dirs)}")
print(f"  Successful: {len(success)}")
print(f"  Failed: {len(failed)}")
print("="*60)

if failed:
    print(f"\n⚠️  Failed patients ({len(failed)}):")
    for patient, reason in failed:
        print(f"  {patient}: {reason}")

if len(success) > 0:
    print(f"\n✔  Successfully processed {len(success)} patients")
    print(f"   Output files copied to: {output_dir}")

EOF

