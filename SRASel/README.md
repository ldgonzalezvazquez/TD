# SRASel – Methods

**1) RG tagging and BAM preparation**
Addition of read groups (RG: ID, LB, PL, PU, SM), coordinate sorting, and indexing using Picard (AddOrReplaceReadGroups). 

[run_add_rg.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_add_rg.sbatch)

**2) Per-position variant profiling**
Generation, for each genomic position, of counts for A/C/G/T/N and indels with bam-readcount using:
- -q 20 (minimum mapping quality),
- -b 20 (minimum base quality),
- -w 0 (indel window; reports only at the start nucleotide).

[run_bam_readcount.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_bam_readcount.sbatch)

**3) Integrity control and case selection**
Reorganization of outputs and safety checks by number of positions (how many passed bam-readcount; >29,000 positions). Identification of time points that have different tokens but share a date (1,2; not 1a, 1b). All samples/patients proceed to downstream analyses.

[run_integrity_control.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_integrity_control.sbatch)

**4) Position×time matrices (per patient)**
Execution of [built_time_mutation_tsv.py](https://github.com/ldgonzalezvazquez/SRASel/blob/main/built_time_mutation_tsv.py), which, from bam-readcount, builds tables where each row is a position and each column (in increasing day order) encodes: days since time 0, coverage, and counts of A/C/T/G/N/ins/del.

[run_built_tables.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_built_tables.sbatch)

<img width="842" height="595" alt="PosTables" src="https://github.com/user-attachments/assets/73f39f67-57d3-4530-b1b8-06e2441351b4" />

**5) Gene/codon/AA annotation**
Execution of [annotate_patient_tables.py](https://github.com/ldgonzalezvazquez/SRASel/blob/main/annotate_patient_tables.py). Reading NC_045512.2.gb (GenBank), each table row (position) is annotated with gene, codon position (1,2,3; allowing simultaneous positions in template jumps), and amino-acid index. This is specifically adapted for SARS-CoV-2 and that reference due to the template jump and the +ssRNA genome; _it must be modified to run on other viruses_.

[run_annotate_tables.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_annotate_tables.sbatch)

**6) Selection by variant×position trajectories**
Execution of [SRASel.py](https://github.com/ldgonzalezvazquez/SRASel/blob/main/SRASel.py) // [run_srasel.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_srasel.sbatch). For each position and variant (A, C, T, G, N, INS, DEL; different indels are treated as distinct variants) a trajectory is plotted. A depth-weighted linear fit across days is performed (accounting for days since day 0) to obtain a slope β (change in frequency per day). A permutation p-value is computed (reshuffling dates 2000 times), and positions are labeled as selected if the frequency across the entire time range changed by ≥3% and p-value < 0.05; they are labeled positive or negative if β>0 or β<0, respectively. Additionally, for each such position, the codon and amino acid are inferred by consensus from the two neighboring positions—only that are clear are accepted: frequency ≥95% for the variant in the neighbor; it must not be an indel, and there must be only one selected position in the codon. If the codon/amino acid is not clear it is annotated as ambiguous.

<img width="842" height="595" alt="aaChangeConst" src="https://github.com/user-attachments/assets/00b4ee6f-104e-44b5-9695-6b05cb2a2912" />

<img width="842" height="595" alt="p-values" src="https://github.com/user-attachments/assets/5fb5fcb9-0404-41e9-86fc-d8ef29ed1393" />

Result plots are created (openable from the patients’ plot pop-ups in HTML) for all positions with a selected variant—showing the fitted slope and time points of all variants (point size encodes total and variant coverage at each time). Selected variants are colored red (β>0), blue (β<0), or gray (p-value > 0.05). Two additional guide lines are drawn: the mean β over all selected slopes for that patient (black dashed) and the mean β over all slopes, including non-selected (gray dashed).


Genomic plots are also created (openable from the pop-ups in HTML). These show a scatter of selected sites along the genome (red = β>0; blue = β<0) against their slope. In addition, sign-specific density curves are drawn and scaled to visualization (axis values are not meaningful): the curve reflects only the existence of selected sites, not slope magnitude, and is obtained by summing Gaussian kernels with σ = 350 bp per point (≈68% of each kernel’s mass within ±350 bp). The same two guide lines (black dashed for selected slopes’ mean; gray dashed for all slopes’ mean) are shown. A gene track based on NC_045512.2.gb is included. This plot visually identifies regions under stronger selection in specific patients without averaging.

Summary files required by the next script are also produced.

**7) Interactive panel and complementary metrics – [Selection.py](https://github.com/ldgonzalezvazquez/SRASel/blob/main/Selection_HTML_construction.py)**

[run_html_construction.sbatch](https://github.com/ldgonzalezvazquez/SRASel/blob/main/run_html_construction.sbatch)

An HTML report with two additional parameters is produced:
- Temporal variability per position [classifying non-variable positions – TVD]: Weighted average of the absolute frequency subtraction of all possible pairs of each iSNV at each position – the weighting is by the coverage of the time with the lowest coverage in each pair and TVDs are normalized by time difference. After that, positions below the 5th percentile of variability are marked non-variable; the rest are “neutral” or selected (if they have selected slopes). If a position is both non-variable and selected, selected prevails.

<img width="842" height="595" alt="TVD" src="https://github.com/user-attachments/assets/e18b9e30-e107-4c63-a0b7-ae68d72deb63" />

- UPGMA showing relationships between patients taking into account the intersections of selected positions [differentiating the variant and the direction of selection].


**References:**

- Bam-readcount: Khanna A, Larson DE, Srivatsan SN, Mosior M, Abbott TE, Kiwala S, Ley TJ, Duncavage EJ, Walter MJ, Walker JR, Griffith OL, Griffith M, Miller CA. Bam-readcount - rapid generation of basepair-resolution sequence metrics. ArXiv [Preprint]. 2021 Jul 27:arXiv:2107.12817v1. PMID: 34341766; PMCID: PMC8328062.

- Picard: no citation paper: https://broadinstitute.github.io/picard/

- Biopython: Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.

- NumPy: Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Del Río JF, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE. Array programming with NumPy. Nature. 2020 Sep;585(7825):357-362. doi: 10.1038/s41586-020-2649-2. Epub 2020 Sep 16. PMID: 32939066; PMCID: PMC7759461.

- McKinney, Wes. (2010). Data Structures for Statistical Computing in Python. 56-61. 10.25080/Majora-92bf1922-00a.

- Plotly (htmls): Plotly Technologies Inc. Title: Collaborative data science. Publisher: Plotly Technologies Inc. Place of publication: Montréal, QC. Date of publication: 2015. URL: https://plot.ly

---
