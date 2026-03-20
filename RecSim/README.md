# RecSim: Viral Recombination Simulation and Detection Pipeline

**RecSim** is an integrated pipeline for simulating viral recombination events using a **modified version of [Santa-Sim](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407609/)** and comparing the detection capabilities of multiple phylogenetic recombination analysis tools: [**ClonalFrameML**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041), [**BEAST-Bacter**](https://academic.oup.com/genetics/article/205/2/857/6066433), and [**RDP5**](https://academic.oup.com/ve/article/7/1/veaa087/6020281?login=true) (Recombination Detection Program).

---

## Table of Contents

1. [Overview](#overview)
2. [Modifications to Santa-Sim](#modifications-to-santa-sim)
3. [System Requirements](#system-requirements)
4. [Installation](#installation)
5. [Project Structure](#project-structure)
6. [Using the Pipeline](#using-the-pipeline)
7. [Command-Line Arguments](#command-line-arguments)
8. [Usage Examples](#usage-examples)
9. [Pipeline Outputs](#pipeline-outputs)
10. [RDP Integration](#rdp-integration)
11. [SLURM Cluster Execution](#slurm-cluster-execution)
12. [Authors and Contact](#authors-and-contact)

---

## Overview

This project aims to **evaluate the ability of different phylogenetic tools to detect recombination events** in viral genomes through simulation-based benchmarking. The workflow:

1. **Simulates** viral populations with recombination using a modified Santa-Sim
2. **Reconstructs** maximum-likelihood phylogenetic trees with IQ-TREE
3. **Detects** recombination events using:
   - **ClonalFrameML**
   - **BEAST-Bacter**
   - **RDP**
4. **Compares** detected events against the simulated ground truth
5. **Generates** Excel spreadsheets with detailed metrics on precision, sensitivity, and breakpoint concordance

---

## Modifications to Santa-Sim

[Santa-Sim](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407609/) is a viral sequence evolution simulator. For this project, **some modifications were made to the Java source code** to capture comprehensive recombination metadata.

### 1. **Enhanced Recombination Event Tracking**

[REVIEW - Ask Darren for the pre-modified version] The original Santa-Sim only recorded which sequences were recombinant descendants and a single (internal) breakpoint per event. After substantial development work, the modified version now tracks:

#### **For every sequence in the simulation:**
- **Events where it descends from the recombinant** (`R:` tags)
- **Events where it descends from the recombinant's grandfather** (`P:` tags), including cases where the recombination occurred *after* that ancestral lineage had already branched

#### **For every recombination event:**
- **Both genomic segments** donated by each parent (not just the internal breakpoint)
- **Ancestral recombination history**: All upstream recombination events (`R:` tags) in which each parent itself descends from a recombinant

This comprehensive tracking enables **reconstruction of the complete recombination genealogy** and crucially allows us to determine **which breakpoints remain detectable** in the final sample.

---

### 2. **New Output Files**

#### **a) `rec_events.txt`** (Recombination Events Log)
A real-time log written during simulation containing detailed information for each recombination event:

**Format:**
```
Event*Gen*P1_frags*P1_tags*P2_frags*P2_tags
1*5234*0-7355**7356-29903*
206*536*0-13062,23325-29903*R:145@378*13062-23325*
```

**Columns:**
- `EventID` - Unique numeric identifier for the event
- `Generation` - Simulation generation when the event occurred
- `Par1_fragments` - Genomic coordinates inherited from Parent 1 (format: `start-end`)
- `Par1_tags` - Ancestral recombination events for Parent 1 (format: `R:id@gen`)
- `Par2_fragments` - Genomic coordinates inherited from Parent 2
- `Par2_tags` - Ancestral recombination events for Parent 2

This file is essential for tracing event genealogies and computing detectability.

---

#### **b) `sequence_events_map_final_sample.txt`** (Sequence-to-Event Mapping)
Maps each sampled sequence to the recombination events in which it participates:

**Format:**
```
SeqID*[R:5@59 ; P:3@5234@0-7355 ; P:7@6012@15000-29903]
23*[R:5@59 ; P:3@5234@0-7355]
```

**Tag types:**
- `R:id@gen` → Sequence **is a recombinant** descendant of event `id`
- `P:id@gen@region` → Ancestors of the sequence acted as a **parental** contributor to event `id`, donating the specified genomic region

This allows identification of:
- Direct recombinant descendants
- Parental contributors to other events
- Carriers of ancestral fragments from multiple events

---

#### **c) Enhanced `recombination_events_final_sample.txt`**
[REVIEW - Ask Darren for the pre-modified version] The existing events file was augmented with **three additional columns** containing the complete sequences of:

1. **Ancestral recombinant sequence** (`anc_rec_seq`)
2. **Ancestral Parent 1 sequence** (`anc_rep_par1`)
3. **Ancestral Parent 2 sequence** (`anc_rep_par2`)

These sequences enable calculation of **sequence similarities between true ancestral parents and sampled sequences**, which is critical for evaluating whether a tool correctly identified the parental lineages.

---

### 3. **Detectability Algorithm**

The pipeline implements a **detectability propagation algorithm** that determines which genomic fragments from each parent remain visible (detectable) in the final sample, accounting for nested and overlapping recombination events.

#### **Why is this necessary?**

Consider this example:
- **Event A** (generation 5000): A 1000-bp recombinant with Parent 1 contributing bases 1–500, Parent 2 contributing 500–1000
- **Event B** (generation 7000): A descendant of Event A acts as a parent, donating only bases 250–750
- **Sampling**: Only descendants of Event B are sampled

In this scenario:
- Bases **1–250 and 750–1000** from Event A are **invisible** (masked by Event B)
- The **detectable breakpoints** for Event A are now **250–750**, completely different from the original 1–500 / 500–1000
- Even entire parents or the whole event might disappear from detection

#### **How the algorithm works:**

1. **Parses relationships** between events using `R:` tags from `rec_events.txt`
2. **Propagates fragments** backward through the genealogy, starting from the most recent (closest to sample) events
3. **Computes intersections** between inherited fragments and original parental contributions
4. **Generates detectability values** (`det_par1`, `det_par2`) which can be:
   - Single region: `0-7355`
   - Multiple regions: `0-990,1500-2000`
   - Multiple alternatives (separated by `/`): `0-990/0-990,1500-2000`
   - Not detectable: `ND`

This ensures **fair tool comparison**: tools can only be evaluated on events that are actually detectable in the final sample.

---

### 4. **Modified Java Source Files**

The modifications were made to:
- src/santa/simulator/samplers/AlignmentSampler.java
- src/santa/simulator/replicators/RecombinantReplicatorWithHotSpots.java
- src/santa/simulator/replicators/RecombinantReplicator.java
- src/santa/simulator/replicators/RecombinantTracker.java
- src/santa/simulator/replicators/RecombinationEvent.java
- src/santa/simulator/samplers/SamplingSchedule.java
- src/santa/simulator/genomes/SimpleGenome.java
- src/santa/simulator/population/Population.java
- src/santa/simulator/Virus.java


**The modified code is included in `software/santa-sim-Recomb_and_align_Luis_MOD/`** and is **pre-compiled** (`dist/santa.jar`).

---

## System Requirements

### Required Software
- **Python 3.8+** with packages:
  - `numpy`, `pandas`, `openpyxl`, `biopython`, `argparse`, `pathlib`, `subprocess`

- **Java 8 or higher** (to run Santa-Sim and BEAST)
- **Apache Ant** (only needed if recompiling Santa-Sim)
- **gcc** and **make** (to compile ClonalFrameML)
- **git** (to download ClonalFrameML)
- **wget** (to download BEAST and IQ-TREE)

### External Tools (Auto-downloadable)
The script can automatically download and compile:
- **Santa-Sim** (modified version, already included in `software/`)
- **IQ-TREE 2**
- **ClonalFrameML**
- **BEAST 2.7.7** with **Bacter** plugin

---

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/ldgonzalezvazquez/RecSim.git
cd RecSim
```

### 2. Install Python Dependencies
```bash
pip install numpy pandas openpyxl biopython
```

### 3. Download and Install External Software

**Option A: Automatic Download (Recommended)**
```bash
# Caution: These installation is designed for a Linux environment; you may need to adjust the versions by consulting the tools download page – see Option B to do manually.
python3 RecSimulator.py --download-software
```

This command:
- Downloads IQ-TREE 2.2.0 for Linux
- Clones and compiles ClonalFrameML from GitHub
- Downloads BEAST 2.7.7 and installs Bacter, BEASTLabs, and MODEL_SELECTION plugins
- **Note:** Santa-Sim is already included and pre-compiled at `software/santa-sim-Recomb_and_align_Luis_MOD/dist/santa.jar`

**Option B: Manual Installation**

If you prefer manual installation:

```bash
# Caution: These installations are designed for a Linux environment; you may need to adjust the versions by consulting the tools download page.

# IQ-TREE
wget https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz
tar -xzf iqtree-2.2.0-Linux.tar.gz -C software/

# ClonalFrameML
cd software/
git clone --depth 1 https://github.com/xavierdidelot/ClonalFrameML
cd ClonalFrameML/src
make
cd ../../..

# BEAST 2.7.7
cd software/
wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz
tar -xzf BEAST.v2.7.7.Linux.x86.tgz
mv beast beast2.7.7
./beast2.7.7/bin/packagemanager -add bacter
./beast2.7.7/bin/packagemanager -add BEASTLabs
./beast2.7.7/bin/packagemanager -add MODEL_SELECTION
cd ..
```

### 4. Verify Installation

**Santa-Sim:**
```bash
java -jar software/santa-sim-Recomb_and_align_Luis_MOD/dist/santa.jar
```
You should see Santa-Sim's help message.

**IQ-TREE:**
```bash
software/iqtree-2.2.0-Linux/bin/iqtree2 --version
```

**ClonalFrameML:**
```bash
software/ClonalFrameML/src/ClonalFrameML
```

**BEAST:**
```bash
software/beast2.7.7/bin/beast -version
```

---

## Project Structure

```
RecSim/
├── RecSimulator.py              # Main pipeline script
├── README.md                    # This file
├── howthemainworks.md           # Technical workflow explanation
├── LICENSE                      # Project license
├── NC_045512.2.fna              # Reference genome (SARS-CoV-2)
│
├── software/                    # External tools
│   ├── santa-sim-Recomb_and_align_Luis_MOD/
│   │   ├── dist/santa.jar       # Modified Santa-Sim (pre-compiled)
│   │   └── src/                 # Modified Java source code
│   ├── iqtree-2.2.0-Linux/
│   ├── ClonalFrameML/
│   └── beast2.7.7/
│
└── results/                     # Output directory (auto-created)
    ├── RDP/                     # FASTA files for RDP analysis
    ├── 1_50_NC_045512.2/        # Each replicate folder
    │   ├── santa/               # Santa-Sim outputs
    │   ├── iqtree/              # ML trees
    │   ├── cfml/                # ClonalFrameML results
    │   └── bacter/              # BEAST-Bacter results
    ├── recombination_summary.xlsx    # Raw events by tool
    └── santa_vs_tools.xlsx           # Detailed comparisons
```

---

## Using the Pipeline

### Mode 1: Full Pipeline (Simulation + Detection)

Run a complete simulation with default parameters:

```bash
python3 RecSimulator.py --start 1 --end 10 --num_seqs 50
```

This executes 10 independent replicates, each with:
- Simulation of 50 sequences with Santa-Sim
- Tree reconstruction with IQ-TREE
- Detection with ClonalFrameML
- Detection with BEAST-Bacter

### Mode 2: Resume BEAST-Bacter Only

If you have previous simulations and only want to run/resume BEAST:

```bash
python3 RecSimulator.py --start 1 --end 10 --resume-beast-only
```

### Mode 3: Extract Results from a Replicate

To re-analyze and extract events from an existing results folder:

```bash
python3 RecSimulator.py --extract-results results/5_50_NC_045512.2
```

### Mode 4: Integrate RDP Results

After analyzing FASTA files with RDP (external GUI), integrate the results:

```bash
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 10
```

---

## Command-Line Arguments

### **General Execution**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--download-software` | flag | - | Download and install all required tools in `./software/` |
| `--start` | int | 1 | Starting index for simulation replicates |
| `--end` | int | 1 | Ending index for simulation replicates |
| `--num_seqs` | int | 50 | Number of sequences to sample per simulation |
| `--ref` | str | `NC_045512.2.fna` | Reference genome in FASTA format (added as "root") |

**Example:**
```bash
python3 RecSimulator.py --start 1 --end 100 --num_seqs 100 --ref my_genome.fasta
```

---

### **Simulation Parameters (Santa-Sim)**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--popSize` | int | 10000 | Viral population size |
| `--mutationRate` | str | `1.1E-6` | Per-site mutation rate per generation |
| `--recombProb` | float | 0.002 | Genome-wide recombination probability |
| `--dualInfectProb` | float | 0.02 | Co-infection probability (required for recombination) |
| `--genCount` | int | 10000 | Number of generations to simulate |

**Example:**
```bash
python3 RecSimulator.py --recombProb 0.01 --dualInfectProb 0.05 --genCount 15000
```

---

### **Software Paths**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--soft_dir` | str | `software` | Base directory for installed software |
| `--santa_bin` | str | `software/.../santa.jar` | Path to Santa-Sim JAR |
| `--iqtree_bin` | str | `software/.../iqtree2` | Path to IQ-TREE binary |
| `--cfml_bin` | str | `software/.../ClonalFrameML` | Path to ClonalFrameML binary |
| `--bacter_bin` | str | `software/.../beast` | Path to BEAST binary |

---

### **BEAST-Bacter Parameters**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--bacter_chain` | int | 50000 | MCMC chain length |
| `--bacter_delta` | float | 100.0 | Initial δ value (mean recombination fragment length, in bp) |
| `--bacter_rho` | float | 0.002 | Initial ρ value (recombination rate parameter) |
| `--beast_mem` | int | 3.5 | Maximum JVM heap size for BEAST (in GB) |
| `--burnin` | float | 0.1 | Burn-in fraction (0–1) for ESS computation |

**Note:** BEAST automatically resumes until ESS ≥ 200 for both `rho` and `delta`.

**Example (longer chains):**
```bash
python3 RecSimulator.py --bacter_chain 1000000 --beast_mem 16
```

---

### **Special Execution Modes**

| Argument | Type | Description |
|----------|------|-------------|
| `--extract-results DIR` | str | Extract recombination events from a completed run directory |
| `--rdp-dir DIR` | str | Parse all `*.fasta.csv` RDP results files in the specified directory |
| `--resume-beast-only` | flag | Skip Santa-Sim, IQ-TREE, and ClonalFrameML; only run/resume BEAST-Bacter |

**Example (extract results):**
```bash
python3 RecSimulator.py --extract-results results/42_100_NC_045512.2
```

**Example (integrate RDP):**
```bash
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 50
```

---

### **Parallelization**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--threads` | int | 1 | Number of threads for IQ-TREE, ClonalFrameML, and BEAST |

**Example (8 threads):**
```bash
python3 RecSimulator.py --threads 8
```

---

### **Filtering and Comparison**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--par-threshold` | float | 1e-8 | Maximum p-distance to consider a sequence as a matching parental |
| `--exclude-RDP` | str | `""` | Comma-separated list of RDP breakpoint flags to ignore (`~,*,$,^`) |

**Example (exclude recombinants marked with ^ in RDP):**
```bash
python3 RecSimulator.py --rdp-dir results/RDP --exclude-RDP "^"
```

---

### **SLURM Cluster Script Generation**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--script_cluster` | flag | - | Generate a `submit_cluster.sh` SLURM script and exit |
| `--cluster_jobname` | str | `RDP_job` | SLURM job name |
| `--cluster_time` | str | `6:00:00` | Time limit (HH:MM:SS) |
| `--cluster_mem` | str | `32G` | Memory per task |
| `--cluster_cpus` | int | 16 | CPUs per task |
| `--cluster_array` | str | `1-100` | SLURM array range |

**Example:**
```bash
python3 RecSimulator.py --script_cluster --cluster_array 1-500 --cluster_cpus 32 --cluster_mem 64G
```

Then submit:
```bash
sbatch submit_cluster.sh
```

---

## Usage Examples

### Example 1: Quick Test Simulation
```bash
python3 RecSimulator.py --start 1 --end 3 \
    --num_seqs 30 \
    --genCount 5000 \
    --bacter_chain 10000 \
    --threads 4
```

### Example 2: High Recombination Simulation
```bash
python3 RecSimulator.py --start 1 --end 50 \
    --recombProb 0.01 \
    --dualInfectProb 0.1 \
    --genCount 20000 \
    --threads 16
```

### Example 3: Resume BEAST Only on Existing Replicates
```bash
python3 RecSimulator.py --start 1 --end 50 \
    --resume-beast-only \
    --bacter_chain 500000 \
    --beast_mem 16 \
    --threads 32
```

### Example 4: Integrate RDP Results (after manual analysis)
```bash
# 1. Run Santa-Sim + IQ-TREE + ClonalFrameML + Bacter
python3 RecSimulator.py --start 1 --end 100 --threads 8

# 2. FASTA files are in results/RDP/*.fasta
# 3. Analyze them with RDP GUI
# 4. Save results as *_Recombination events.csv in results/RDP/

# 5. Integrate RDP results
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 100
```

### Example 5: Extract Results from a Specific Replicate
```bash
python3 RecSimulator.py --extract-results results/15_50_NC_045512.2
```

---

## Pipeline Outputs

### 1. Per-Replicate Files

Each replicate creates a folder `results/<idx>_<num_seqs>_<ref_name>/` containing:

#### **Santa-Sim (`santa/`)**
- `final_sample.fasta` - Sampled sequences
- `final_sample_withroot.fasta` - With reference genome added
- `final_tree.newick` - True genealogical tree
- `rec_events.txt` - Recombination events with fragments and tags
- `sequence_events_map_final_sample.txt` - Sequence-to-event mapping
- `recombination_events_final_sample.txt` - Events with ancestral sequences
- `parental_tags_by_event.txt` -  Unused
- `cfg.xml` - Simulation configuration

#### **IQ-TREE (`iqtree/`)**
- `<tag>.treefile` - Maximum-likelihood tree (rooted)
- `<tag>_unrooted.tree` - Unrooted tree (for ClonalFrameML)
- `<tag>.iqtree` - Detailed report
- Other unused remnant files

#### **ClonalFrameML (`cfml/`)**
- `<tag>.importation_status.txt` - Detected events with breakpoints
- `<tag>.labelled_tree.newick` - Tree annotated with recombination nodes
- Other unused remnant files

#### **BEAST-Bacter (`bacter/`)**
- `bacter.log` - MCMC chain log
- `bacter.trees` - Sampled trees (extended NEXUS format)
- `summary.tree` - Consensus tree with recombination annotations (if converged)
- `<tag>.xml` - BEAST configuration
- Other unused remnant files

#### **RDP (`results/RDP/`)**
- `<tag>.fasta` - Sequences for manual RDP analysis

---

### 2. Global Results Files

#### **`recombination_summary.xlsx`** (Sheet: `events_by_tool`)
Contains **all detected events** from each tool, with columns:
- `tag` - Replicate identifier
- Simulation parameters (`popSize`, `mutationRate`, etc.)
- `tool` - Tool name (Santa, CFML, Bacter, RDP)
- `event_id` - Event ID in the external tool
- Tool-specific data

#### **`santa_vs_tools.xlsx`** (Sheet: `santa_vs_tools`)
Detailed **event-by-event comparisons** between Santa-Sim (ground truth) and each tool:

| Column | Description |
|--------|-------------|
| `tag` | Replicate identifier (`<idx>_<nSeq>_<ref>`) |
| `event` | Santa-Sim event number (ground truth) |
| Simulation parameters | `popSize`, `mutationRate`, etc.|
| `tool` | Comparison tool (CFML, Bacter, RDP) |
| `ext_id` | Event ID in the external tool |
| `par1_match_n` | Number of shared Parent 1 sequences (∣P1 ∩ T∣) |
| `par1_missing_santa` | Sequences in Santa P1 but not in tool (∣P1 \ T∣) |
| `par1_missing_tool` | Sequences in tool but not in Santa P1 (∣T \ P1∣) |
| `par1_match_pct` | Parent 1 recovery percentage: 100·∣P1∩T∣ / ∣P1∣ |
| `par2_match_n` | Similar for Parent 2 |
| `par2_missing_santa` | Similar for Parent 2 |
| `par2_missing_tool` | Similar for Parent 2 |
| `par2_match_pct` | Similar for Parent 2 |
| `rec_match_n` | Number of detected recombinants |
| `rec_missing_santa` | Recombinants in Santa but not in tool |
| `rec_missing_tool` | Tips flagged as recombinants by tool but not by Santa |
| `rec_match_pct` | Percentage: 100·∣REC∩T∣ / ∣REC∣ |
| `nonmatch_sim_prod_p1` | Similarity product for EXTRA sequences in P1 (false positive penalty) |
| `nonmatch_sim_prod_p2` | Similar for P2 |
| `break_overlap_nt` | Breakpoint intersection size (nucleotides) |
| `break_overlap_pct` | Dice overlap: 100 × (2·∣∩∣ / (∣S∣+∣T∣)) |
| `break_bp_sizes` | Sizes of (Santa_breakpoint, tool_breakpoint) |

**Interpretation:**
- **High `par1_match_pct` and `par2_match_pct`** → Tool correctly identified parental lineages
- **High `break_overlap_pct`** → Tool correctly localized recombination breakpoint
- **Low `nonmatch_sim_prod_p1/p2`** → Tool included incorrect sequences (penalty)

---

## RDP Integration

**RDP (Recombination Detection Program)** is a GUI-based software package that cannot be directly integrated into an automated pipeline. The workflow is:

### Step 1: Generate FASTA Files for RDP
Run the pipeline normally:
```bash
python3 RecSimulator.py --start 1 --end 100
```

FASTA files are automatically copied to `results/RDP/<tag>.fasta`.

### Step 2: Analyze with RDP (Manual)
1. Open RDP GUI
2. Load each `<tag>.fasta` file – You can analyze them all at once in series by selecting them together when opening the file.
3. Run detection methods – Automatic when you select more than one file / if not press "Run" and save results files
4. Export results as `<tag>_Recombination events.csv` in the same `results/RDP/` directory – Automatic name when you select more than one file

### Step 3: Integrate Results
```bash
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 100
```

This command:
- Reads each `*_Recombination events.csv` file
- Extracts recombination events with breakpoints and parentals (minor/major)
- Compares with Santa-Sim events
- Updates Excel files with RDP metrics

### Step 4 (Optional): Exclude Uncertain Breakpoint Flags

RDP marks some breakpoints with special flags [(see RDP5 documentation)](https://web.cbio.uct.ac.za/~darren/RDP5Manual.pdf):
- `~`
- `*` 
- `$` 
- `^` 

To ignore events with these marks:
```bash
python3 RecSimulator.py --rdp-dir results/RDP --exclude-RDP "~,*,$" --start 1 --end 100
```

---
## New Features (Under Review: [RecSimulator_v2.py](RecSimulator_v2.py))

**⚠️ Status: This version includes new experimental features that are currently under review and testing.**

### New Command-Line Options

#### `--only-RDP`
When processing RDP results or running simulations, extract and compare only Santa and RDP events, ignoring ClonalFrameML and Bacter directories.

**Behavior:**
- During simulation: Executes Santa-Sim and skips IQ-TREE, ClonalFrameML, and Bacter execution
- When processing RDP results: Extracts events from Santa and RDP only, ignoring CFML and Bacter directories

**Use case:** Useful when you only want to compare Santa-Sim ground truth with RDP detections, without requiring the computationally intensive ClonalFrameML and Bacter analyses.

#### `--write-also-no-detection`
Writes `santa_vs_tools.xlsx` even when some tools have no detections, leaving corresponding columns empty for missing tools.

**Behavior:**
- If Santa has events but no external tools: Writes Santa events with `tool='ND'` and empty comparison columns
- If RDP/CFML/Bacter have events but no Santa: Writes those tool events without comparisons
- Works with `--only-RDP`: Only considers Santa and RDP when active

**Use case:** Ensures comprehensive data collection even when detection tools fail or produce no results, maintaining complete records for analysis.

### Usage Examples

#### Example 1: Run simulations with Santa-Sim only (skip CFML and Bacter)
```bash
python3 RecSimulator.py --start 1 --end 5 --only-RDP \
    --num_seqs 50 \
    --ref NC_045512.2.fna
```

This command:
- Executes 5 simulations with Santa-Sim
- Copies FASTA files to `results/RDP/` for RDP analysis
- Skips IQ-TREE, ClonalFrameML, and Bacter execution
- Extracts only Santa events from results

#### Example 2: Process RDP results with Santa-RDP comparison only
```bash
python3 RecSimulator.py --rdp-dir results/RDP \
    --start 1 --end 5 \
    --only-RDP
```

This command:
- Reads RDP CSV files from `results/RDP/`
- Extracts Santa events from simulation directories
- Compares only Santa vs RDP (ignores CFML and Bacter)
- Updates Excel files with Santa-RDP comparisons

#### Example 3: Write Excel even when tools have no detections
```bash
python3 RecSimulator.py --rdp-dir results/RDP \
    --start 1 --end 10 \
    --only-RDP \
    --write-also-no-detection
```

This command:
- Processes RDP results with Santa-RDP comparison
- Writes `santa_vs_tools.xlsx` even if Santa or RDP have no detections
- Leaves columns empty for missing data instead of skipping rows

#### Example 4: Full pipeline with no-detection flag
```bash
python3 RecSimulator.py --start 1 --end 20 \
    --write-also-no-detection \
    --threads 8
```

This command:
- Runs full simulation pipeline (all tools)
- Ensures Excel files are written even when some tools detect no recombination
- Useful for comprehensive benchmarking studies

**Note:** These features are actively being tested. If you encounter issues, please report them through the project's issue tracker or luisdaniel.gonzalez@uvigo.es.

---

## SLURM Cluster Execution

To run 500 replicates in parallel on an HPC cluster:

### 1. Generate Submission Script
```bash
python3 RecSimulator.py --script_cluster \
    --cluster_array 1-500 \
    --cluster_cpus 16 \
    --cluster_mem 32G \
    --cluster_time 8:00:00 \
    --cluster_jobname RecSim_run
```

This creates `submit_cluster.sh` and a `logs/` directory.

### 2. Submit Job
```bash
sbatch submit_cluster.sh
```

Each array task executes:
```bash
python3 RecSimulator.py --start $SLURM_ARRAY_TASK_ID --end $SLURM_ARRAY_TASK_ID --threads 16
```

### 3. Monitor Progress
```bash
squeue -u $USER
cat logs/RecSim_run_<jobid>_<taskid>.out
```

### 4. Process RDP Results
```bash
# After manual RDP analysis
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 500
```

---

## Authors and Contact

**Developed by:**
- **Luis Daniel González-Vázquez**  
  Email: luisdaniel.gonzalez@uvigo.es / luisde365@gmail.com
- In collaboration with **Darren P. Martin**

**Development Period:**  
April - June 2025

**Institution:**  
University of Vigo (Faculty of Sciences and CINBIO) and University of Cape Town

**Funding:**  
Xunta de Galicia ED481A-2023/089, programa de axudas á etapa predoutoral da Xunta de Galicia (Consellería de Cultura, Educación, Formación Profesional e Universidades) cofinanciado pola Unión Europea no marco do Programa FSE+ Galicia 2021-2027.

**License:**  
MIT License (see [`LICENSE`](LICENSE) file)

---

## References

[These are just the main references, remember to also cite the methodologies included in these programs that you have used during your analysis and that are available in the documentation of each program]

- **Santa-Sim:** Jariani A, Warth C, Deforche K, Libin P, Drummond AJ, Rambaut A, Matsen Iv FA, Theys K. SANTA-SIM: simulating viral sequence evolution dynamics under selection and recombination. Virus Evol. 2019 Mar 8;5(1):vez003. doi: 10.1093/ve/vez003. PMID: 30863552; PMCID: PMC6407609.

- **IQ-TREE:** Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol Biol Evol. 2020 May 1;37(5):1530-1534. doi: 10.1093/molbev/msaa015. Erratum in: Mol Biol Evol. 2020 Aug 1;37(8):2461. doi: 10.1093/molbev/msaa131. PMID: 32011700; PMCID: PMC7182206.

- **ClonalFrameML:** Didelot X, Wilson DJ. ClonalFrameML: efficient inference of recombination in whole bacterial genomes. PLoS Comput Biol. 2015 Feb 12;11(2):e1004041. doi: 10.1371/journal.pcbi.1004041. PMID: 25675341; PMCID: PMC4326465.

- **BEAST-Bacter:** Vaughan TG, Welch D, Drummond AJ, Biggs PJ, George T, French NP. Inferring Ancestral Recombination Graphs from Bacterial Genomic Data. Genetics. 2017 Feb;205(2):857-870. doi: 10.1534/genetics.116.193425. Epub 2016 Dec 22. PMID: 28007885; PMCID: PMC5289856.

- **RDP:** Martin DP, Varsani A, Roumagnac P, Botha G, Maslamoney S, Schwab T, Kelz Z, Kumar V, Murrell B. RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. Virus Evol. 2020 Apr 12;7(1):veaa087. doi: 10.1093/ve/veaa087. PMID: 33936774; PMCID: PMC8062008.
