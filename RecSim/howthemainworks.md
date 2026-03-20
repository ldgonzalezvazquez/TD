# How the Main Function Works in RecSimulator.py

This document explains the execution flow of the [`RecSimulator.py`](RecSimulator.py) pipeline in a simplified manner.

---

## General Structure of Main

The entry point is the `main()` function (line (2657)[RecSimulator.py#L2657]-3024), which:

1. **Parses command-line arguments** using `argparse`
2. **Executes special modes** (if enabled):
   - `--download-software` → Downloads tools and exits
   - `--script_cluster` → Generates SLURM script and exits
   - `--extract-results DIR` → Extracts events from a folder and exits
   - `--rdp-dir DIR` → Integrates RDP results and exits
3. **Executes the full pipeline** for each replicate in the range `--start` to `--end`:
   - If `--resume-beast-only` → Jumps directly to BEAST
   - Otherwise → Runs Santa-Sim → IQ-TREE → ClonalFrameML → BEAST → Results extraction

---

## Full Pipeline Flow

### **Main Loop** (line [2804](RecSimulator.py#L2804))

```python
for idx in range(args.start, args.end + 1):
    tag = f"{idx}_{args.num_seqs}_{ref_name}"
    run_dir = root / tag
    # ... subdirectories santa/, iqtree/, cfml/, bacter/ are created
```

Each replicate `idx` generates its own folder `results/<tag>/`.

---

### **Step 1: Simulation with Santa-Sim** (lines [2898](RecSimulator.py#L2898)–2936)

```python
# 1. Create XML configuration file
create_santa_config_xml(
    args.ref, cfg_xml,
    args.popSize, args.mutationRate,
    args.recombProb, args.dualInfectProb,
    args.genCount, args.num_seqs
)

# 2. Execute Santa-Sim
shell([
    "java", "-jar", str(santa_jar),
    f"-seed={idx}",
    str(xml_abs)
], cwd=santa_dir)
```

**Outputs:** (see [README](README.md))
- `final_sample.fasta`
- `rec_events.txt`
- `sequence_events_map_final_sample.txt`
- `recombination_events_final_sample.txt`

**Note:** The reference genome is added to the FASTA as "root" to root the tree.

---

### **Step 2: Tree Reconstruction with IQ-TREE** (lines [2940](RecSimulator.py#L2940)–2948)

```python
# Execute IQ-TREE with model selection
shell([str(iqtree_bin), "-T", str(args.threads), 
       "-s", str(withroot_fa), "-m", "TEST", 
       "-nt", "AUTO", "-o", ref_id, 
       "--seed", str(idx), "-pre", str(mlpref)])

# Remove root and save unrooted tree
t = Phylo.read(rooted_tree, "newick")
t.rooted = False
t.prune("root")
Phylo.write(t, unroot_tree, "newick")
```

**Outputs:** (see [README](README.md))
- `<tag>.treefile`
- `<tag>_unrooted.tree`

---

### **Step 3: Detection with ClonalFrameML** (lines [2950](RecSimulator.py#L2950)–2953)

```python
cfml_prefix = cfml_dir / tag
shell([str(cfml_bin), str(unroot_tree), str(fasta), str(cfml_prefix)])
```

**Outputs:** (see [README](README.md))
- `<tag>.importation_status.txt`
- `<tag>.labelled_tree.newick`

---

### **Step 4: Detection with BEAST-Bacter** (lines [2954](RecSimulator.py#L2954)–3002)

#### **4.1 Create XML Configuration**
```python
build_bacter_xml(fasta, bacter_xml, tag,
                 args.bacter_chain, args.bacter_delta, args.bacter_rho)
```

Generates an XML with:
- GTR model
- Selected priors for δ and ρ
- Bacter-specific MCMC operators

#### **4.2 Execute BEAST with ESS Control**
```python
# First execution
run_beast([str(bacter_bin), "-threads", str(args.threads),
           "-seed", str(idx), str(bacter_xml.resolve())])

# Resume loop until ESS ≥ 200
while True:
    df = read_beast_log(log_file)
    ess = compute_ess(df, args.burnin)
    r, d = ess.get("rho"), ess.get("delta")
    if r >= 200 and d >= 200:
        break
    # Resume
    run_beast([..., "-resume", ...])
```

#### **4.3 Condense Tree**
```python
shell([str(condensetree_bin), "ACGAnnotator", "bacter.trees"], cwd=bacter_dir)
```

**Outputs:** (see [README](README.md))
- `bacter.log`
- `bacter.trees`
- `summary.tree`

---

### **Step 5: Results Extraction and Comparison** (lines [2998](RecSimulator.py#L2998)–3024)

```python
# 1. Extract events from all three tools
santa_dict, cfml_dict, bacter_dict = extract_recombination_results(
    run_dir, par_threshold=args.par_threshold)

# 2. Parse parameters from cfg.xml (Santa config)
run_params = parse_cfg_params(cfg_xml)

# 3. Generate Excel spreadsheets
write_event_summary(tag, run_params, santa_dict, cfml_dict, bacter_dict, rdp=None)
write_tool_metrics(df_summary, Path("results/performance_summary.csv"))
```

**Key Functions:**
- [`extract_recombination_results()`](RecSimulator.py#L1416) - Compute reading of Santa-Sim, ClonalFrameML and Bacter result files
   - [`parse_santasim()`](RecSimulator.py#L815) - Reads `rec_events.txt`, `sequence_events_map_filtered.txt`, computes detectability, ancestral similarities
   - [`parse_clonalframeml_pair()`](RecSimulator.py#L1209) - Reads `.importation_status.txt` and `.labelled_tree.newick`
   - [`parse_bacter_nexus()`](RecSimulator.py#L1295) - Parses `summary.tree` (extended NEXUS format)
- [`write_event_summary()`](RecSimulator.py#L2641) - Write comparison excel of the tools against Santa-Sim
   - [`build_event_rows_v3()`](RecSimulator.py#L1953) - Makes comparisons of the tools against Santa-Sim
- [`write_events_excel()`](RecSimulator.py#L2125) - Write raw output from methods for all tests -> **run it via [--extract-results](RecSimulator.py#L2754)**
- [`write_comparisons_excel()`](RecSimulator.py#L2231) - Write comparisons output of the tools against Santa-Sim for all tests -> **run it via [--extract-results](RecSimulator.py#L2754)**

**Outputs:** (see [README](README.md))
- `recombination_summary.xlsx` - Raw events by tool
- `santa_vs_tools.xlsx` - Detailed Santa vs. tools comparisons

---

## Modifications to Run Only Santa-Sim and RDP

If you want to **skip ClonalFrameML and BEAST** to run only Santa-Sim (simulation) and then integrate RDP, follow these steps:

### **Step 1: Comment Code Blocks** (Simple Method)

Edit [`RecSimulator.py` main](RecSimulator.py#L2657) and comment out the following sections:

#### **1. Comment IQ-TREE** (lines [2940](RecSimulator.py#L2940)–2949)
```python
# # ───── 2 · IQ‑TREE ─────
# shell([str(iqtree_bin), "-T", str(args.threads), "-s", str(withroot_fa), 
#        "-m", "TEST", "-nt", "AUTO", "-o", ref_id, 
#        "--seed", str(idx),"-pre", str(mlpref)])
# 
# rooted_tree   = mlpref.with_name(mlpref.name + ".treefile")
# unroot_tree = mlpref.with_name(mlpref.name + "_unrooted.tree")
# t = Phylo.read(rooted_tree, "newick")
# t.rooted = False
# t.prune("root")
# Phylo.write(t, unroot_tree, "newick")
```

#### **2. Comment ClonalFrameML** (lines [2950](RecSimulator.py#L2950)–2953)
```python
# # ───── 3 · ClonalFrameML ─────
# cfml_prefix = cfml_dir / tag
# shell([str(cfml_bin), str(unroot_tree), str(fasta), str(cfml_prefix)])
```

#### **3. Comment BEAST-Bacter** (lines [2954](RecSimulator.py#L2950)–2997)
```python
# # ───── 4 · Bacter ────────────────────────────────────────────────────
# bacter_xml = bacter_dir / f"{tag}.xml"
# build_bacter_xml(fasta, bacter_xml, tag,
#                  args.bacter_chain, args.bacter_delta, args.bacter_rho)
# 
# env_beast = os.environ.copy()
# env_beast["BEAST_OPTS"] = f"-Xmx{args.beast_mem}g -Xss256m"
# 
# def run_beast(cmd):
#     rc = shell(cmd, cwd=bacter_dir, env=env_beast, ignore_err=True)
#     if rc not in (0, 1, 137):
#         sys.exit(f"BEAST terminated with unexpected code {rc}")
#     return rc
# 
# run_beast([str(bacter_bin), "-threads", str(args.threads),"-seed", str(idx), str(bacter_xml.resolve())])
# 
# log_file = bacter_dir / "bacter.log"
# 
# r = d = float("nan")
# while True:
#     if not log_file.exists():
#         print("   (bacter.log not showing up yet, waiting…)")
#     else:
#         df  = read_beast_log(log_file)
#         ess = compute_ess(df, args.burnin)
#         r, d = ess.get("rho"), ess.get("delta")
#         if r is not None and d is not None:
#             print(f"   ➤ ESS rho.t={r:.1f}, delta.t={d:.1f}")
#             if r >= 200 and d >= 200:
#                 print(f"Convergence has been achieved with Rho = {r} and Delta = {d}", flush=True)
#                 break
#     print(f"   ↻ Insufficient ESS (rho = {r}; delta = {d}) or BEAST interrupted; resuming…", flush=True)
#     run_beast([str(bacter_bin), "-threads", str(args.threads),"-seed", str(idx), "-resume", str(bacter_xml.resolve())])
# 
# bater_trees_file = bacter_dir / "bacter.trees"
# condense_tree = bacter_dir / "summary.tree"
# 
# shell([ str(condensetree_bin), "ACGAnnotator", "bacter.trees" ], cwd=bacter_dir)
```
               
### **Step 2: Complete Santa-Sim → RDP Workflow**

If you want to run Santa-Sim and then integrate RDP:

#### **Step 1: Run Only Santa-Sim**
```bash
# After commenting the IQ-TREE, ClonalFrameML, and BEAST sections
python3 RecSimulator.py --start 1 --end 100 --num_seqs 50
```

This generates:
- `results/<tag>/santa/` with all Santa-Sim files
- `results/RDP/<tag>.fasta` (copied on line [2938](RecSimulator.py#L2938))

#### **Step 2: Analyze with RDP (Manual)**
1. Open RDP (GUI)
2. Load and analyze all `results/RDP/<tag>.fasta`
3. Move `<tag>_Recombination events.csv` to `results/RDP/`

#### **Step 3: Integrate Results**
```bash
python3 RecSimulator.py --rdp-dir results/RDP --start 1 --end 100
```

This:
- Reads RDP CSVs
- Compares with Santa-Sim (ground truth)
- Generates `recombination_summary.xlsx` and `santa_vs_tools.xlsx` with RDP metrics

---

## Pipeline Summary

```
┌────────────────────────────────────────────────────────────────┐
│  main()                                                        │
│    ├─ Parse arguments                                          │
│    ├─ Special modes (--download-software, --rdp-dir, etc.)     │
│    └─ Loop for idx in range(start, end+1):                     │
│         ├─ Create directories (santa, iqtree, cfml, bacter)    │
│         ├─ [1] Santa-Sim                                       │
│         │     ├─ create_santa_config_xml()                     │
│         │     └─ java -jar santa.jar cfg.xml                   │
│         │           → final_sample.fasta                       │
│         │           → rec_events.txt                           │
│         │           → sequence_events_map_final_sample.txt     │
│         │           → recombination_events_final_sample.txt    │
│         ├─ [2] IQ-TREE                                         │
│         │     └─ iqtree2 -s final_sample_withroot.fasta        │
│         │           → <tag>.treefile                           │
│         │           → <tag>_unrooted.tree                      │
│         ├─ [3] ClonalFrameML                                   │
│         │     └─ ClonalFrameML <tree> <fasta> <prefix>         │
│         │           → <tag>.importation_status.txt             │
│         │           → <tag>.labelled_tree.newick               │
│         ├─ [4] BEAST-Bacter                                    │
│         │     ├─ build_bacter_xml()                            │
│         │     ├─ beast -threads N <tag>.xml                    │
│         │     ├─ Loop while ESS < 200:                         │
│         │     │     └─ beast -resume <tag>.xml                 │
│         │     └─ ACGAnnotator bacter.trees                     │
│         │           → summary.tree                             │
│         ├─ [5] Results Extraction                              │
│         │     ├─ extract_recombination_results()               │
│         │     │     ├─ parse_santasim()                        │
│         │     │     ├─ parse_clonalframeml_pair()              │
│         │     │     └─ parse_bacter_nexus()                    │
│         │     ├─ write_events_excel()                          │
│         │     ├─ write_comparisons_excel()                     │
│         │     └─ write_tool_metrics()                          │
│         └─ Next replicate (idx+1)                              │
└────────────────────────────────────────────────────────────────┘
```

---

**Final Note:** This document is a simplified guide. For complete details on each function, consult the comments in the [`RecSimulator.py`](RecSimulator.py) source code and [README](README.md) or mail luisdaniel.gonzalez@uvigo.es.
