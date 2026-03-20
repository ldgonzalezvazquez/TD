#!/usr/bin/env python3
"""
Integrity control and case selection - Step 3

Reorganization of outputs and safety checks:
- Number of positions check (>29,000 positions)
- Identification of time points with different tokens but sharing a date (1,2; not 1a, 1b)
- Only samples/patients that pass these checks proceed to downstream analyses.
"""
import re
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

MIN_POSITIONS = 29000
DATE_FMT = "%d/%m/%y"

# Token regex: matches "-1", "-2", "-1a", "-2b", etc.
TOKEN_RE = re.compile(r"-(\d+)([a-z]?)$")

def parse_sample_id(filename):
    """Extract patient and token from filename like '007-A-1a_RG_readcount.txt'"""
    # Remove _RG_readcount.txt suffix
    stem = filename.replace("_RG_readcount.txt", "")
    # Extract token (e.g., "1a", "2b", "3")
    parts = stem.rsplit("-", 1)
    if len(parts) == 2:
        patient = parts[0]  # e.g., "007-A"
        token = parts[1]     # e.g., "1a"
        return patient, token
    return None, None

def count_positions(filepath):
    """Count number of positions in a readcount file"""
    count = 0
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 4:
                        try:
                            pos = int(fields[1])
                            count += 1
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return 0
    return count

def check_time_point_conflicts(patient_data, dates_tsv):
    """
    Check for time points with different tokens but sharing a date (1,2; not 1a, 1b).
    Returns True if patient passes the check, False otherwise.
    """
    # Group samples by date
    date_to_tokens = defaultdict(list)
    
    for sample_id, date_str in dates_tsv.items():
        if sample_id in patient_data:
            try:
                date = datetime.strptime(date_str, DATE_FMT).date()
                token = patient_data[sample_id]['token']
                # Extract numeric part of token (1 from "1a", 2 from "2", etc.)
                match = TOKEN_RE.search(f"-{token}")
                if match:
                    numeric_token = match.group(1)
                    # Only check tokens without letter suffix (1, 2, 3, not 1a, 1b)
                    if not match.group(2):  # No letter suffix
                        date_to_tokens[date].append((numeric_token, token, sample_id))
            except (ValueError, KeyError):
                continue
    
    # Check for conflicts: same date, different numeric tokens
    conflicts = []
    for date, token_list in date_to_tokens.items():
        numeric_tokens = set([t[0] for t in token_list])
        if len(numeric_tokens) > 1:
            conflicts.append((date, token_list))
    
    return len(conflicts) == 0, conflicts

def main():
    input_dir = Path("processing/02_rc")
    output_dir = Path("processing/03_integrity_checked")
    dates_file = Path("dates.tsv")
    
    # Read dates.tsv
    dates_tsv = {}
    try:
        with open(dates_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        sample_id = parts[0]
                        date_str = parts[1]
                        dates_tsv[sample_id] = date_str
    except Exception as e:
        print(f"Error reading dates.tsv: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Group files by patient
    patient_files = defaultdict(dict)
    
    for rc_file in sorted(input_dir.glob("*_readcount.txt")):
        patient, token = parse_sample_id(rc_file.name)
        if not patient or not token:
            print(f"Warning: Could not parse sample ID from {rc_file.name}", file=sys.stderr)
            continue
        
        # Remove _RG_readcount.txt suffix to get sample_id
        # Sample IDs in dates.tsv don't have _RG suffix
        sample_id_base = rc_file.name.replace("_RG_readcount.txt", "")
        sample_id_clean = sample_id_base.replace("_RG", "") if sample_id_base.endswith("_RG") else sample_id_base
        
        patient_files[patient][sample_id_clean] = {
            'filepath': rc_file,
            'token': token,
            'sample_id_for_file': sample_id_base  # Keep original for file naming
        }
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Statistics
    all_patients = []
    patients_with_issues = []
    total_samples_checked = 0
    total_samples_passed = 0
    total_samples_failed = 0
    
    # Process each patient
    for patient, samples in sorted(patient_files.items()):
        patient_dir = output_dir / patient
        patient_dir.mkdir(exist_ok=True)
        
        patient_issues = []
        failed_samples = []
        passed_samples = []
        
        # Check each sample for minimum positions
        for sample_id, data in samples.items():
            total_samples_checked += 1
            positions = count_positions(data['filepath'])
            
            if positions < MIN_POSITIONS:
                issue_msg = f"{sample_id}: Only {positions} positions (< {MIN_POSITIONS})"
                print(f"⚠️  {patient}/{issue_msg}", file=sys.stderr)
                patient_issues.append(issue_msg)
                failed_samples.append((sample_id, f"positions={positions}"))
                total_samples_failed += 1
            else:
                passed_samples.append(sample_id)
                total_samples_passed += 1
        
        # Check time point conflicts (check all samples, not just passed ones)
        no_conflicts, conflicts = check_time_point_conflicts(
            {sid: data for sid, data in samples.items()},
            dates_tsv
        )
        
        if not no_conflicts:
            conflict_details = []
            for date, token_list in conflicts:
                tokens_str = ", ".join([f"{t[1]} ({t[0]})" for t in token_list])
                conflict_details.append(f"Date {date}: {tokens_str}")
                print(f"⚠️  {patient}: Time point conflicts detected:", file=sys.stderr)
                print(f"   Date {date}: {tokens_str}", file=sys.stderr)
            issue_msg = f"Time point conflicts: {'; '.join(conflict_details)}"
            patient_issues.append(issue_msg)
        
        # Always copy ALL samples to output directory (even if there are issues)
        all_sample_ids = sorted(samples.keys())
        
        if all_sample_ids:
            # Create .time file for this patient (with ALL samples)
            time_file = patient_dir / f"{patient}.time"
            with open(time_file, 'w') as tf:
                for sample_id in all_sample_ids:
                    if sample_id in dates_tsv:
                        tf.write(f"{sample_id}\t{dates_tsv[sample_id]}\n")
            
            # Copy and rename readcount files (ALL samples, including failed ones)
            for sample_id in all_sample_ids:
                if sample_id in samples:
                    source = samples[sample_id]['filepath']
                    # Use sample_id for file naming (without _RG if present)
                    # built_time_mutation_tsv.py expects files named like: sample_id.readcounts.tsv
                    dest = patient_dir / f"{sample_id}.readcounts.tsv"
                    shutil.copy2(source, dest)
            
            all_patients.append(patient)
            if patient_issues:
                patients_with_issues.append((patient, patient_issues, len(passed_samples), len(failed_samples)))
                print(f"⚠️  {patient}: {len(passed_samples)} samples passed, {len(failed_samples)} failed - ISSUES DETECTED (all samples copied)")
            else:
                print(f"✔ {patient}: {len(passed_samples)} samples passed (no issues, all samples copied)")
        else:
            # Patient has no samples (shouldn't happen)
            all_patients.append(patient)
            patients_with_issues.append((patient, patient_issues, 0, len(failed_samples)))
            print(f"❌ {patient}: No samples found", file=sys.stderr)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Total samples checked: {total_samples_checked}")
    print(f"  Total samples passed: {total_samples_passed}")
    print(f"  Total samples failed: {total_samples_failed}")
    print(f"  Total patients processed: {len(all_patients)}")
    print(f"  Patients with issues: {len(patients_with_issues)}")
    print(f"  Patients without issues: {len(all_patients) - len(patients_with_issues)}")
    print(f"{'='*60}")
    
    if patients_with_issues:
        print(f"\n⚠️  Patients with issues detected:")
        print(f"{'='*60}")
        for patient, issues, passed_count, failed_count in patients_with_issues:
            print(f"\n  Patient: {patient}")
            print(f"    Samples passed: {passed_count}")
            print(f"    Samples failed: {failed_count}")
            print(f"    Issues ({len(issues)}):")
            for issue in issues:
                print(f"      - {issue}")
        print(f"\n{'='*60}")

if __name__ == "__main__":
    main()

