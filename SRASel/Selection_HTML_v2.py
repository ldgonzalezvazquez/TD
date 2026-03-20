#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# interactive_heatmap_upset.py · v2025-07-31  (UPGMA edition)
# link plots + genome PDF  ✚  *clean* protein view (uses /product= only)
# ---------------------------------------------------------------------------
#  ▸ Heat-map β mean per patient × gene / protein
#  ▸ Stacked-bar UpSet per category (click -> table with links)
#  ▸ Two views toggled with the button: gene | protein (mat_peptide products)
# ---------------------------------------------------------------------------

### ───────── PREAMBLE: FIX LOCALE ─────────
import os, locale
# Evita nl_langinfo errors en forks / spawn
os.environ.setdefault('LC_ALL', 'en_US.UTF-8')
os.environ.setdefault('LANG',   'en_US.UTF-8')
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    # Entornos HPC (p. ej. CESGA) pueden no tener en_US.UTF-8.
    # Fallback robusto sin abortar ejecución.
    for _loc in ('C.UTF-8', 'UTF-8', 'C'):
        try:
            locale.setlocale(locale.LC_ALL, _loc)
            break
        except locale.Error:
            continue
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

import argparse, glob, json, re, sys
from collections import defaultdict, OrderedDict
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.utils import PlotlyJSONEncoder
import itertools
from tqdm import tqdm
from multiprocessing import Pool, get_start_method
from datetime import datetime
import numpy as _np

# Debug flag for worker-level diagnostics
DEBUG_WORKER = False


def _json_dumps_compact(obj, **kwargs):
    """Compact JSON serializer for HTML payloads."""
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=False, **kwargs)


def _minify_html_inline(html_text):
    """Conservative inline minification preserving script/style content."""
    try:
        import minify_html  # optional dependency
        return minify_html.minify(
            html_text,
            minify_css=True,
            minify_js=False,
            remove_processing_instructions=True,
        )
    except Exception:
        pass

    parts = re.split(r"(<script.*?</script>|<style.*?</style>)", html_text, flags=re.S | re.I)
    out = []
    for chunk in parts:
        if not chunk:
            continue
        if re.match(r"^<script|^<style", chunk, flags=re.I):
            out.append(chunk)
            continue
        compact = re.sub(r">\s+<", "><", chunk)
        compact = re.sub(r"\s{2,}", " ", compact)
        out.append(compact)
    return "".join(out)

# ───────── WORKER FUNCTION ─────────

def worker(args):
    pat, pth, idx = args
    df_long = pd.read_csv(pth, sep="\t", dtype={"pos": int})
    cov_df = df_long[['pos','day','n']].copy()
    cov_df['pos'] = cov_df['pos'].astype(int)
    cov_df['day'] = cov_df['day'].astype(int)
    cov_df = cov_df.groupby(['pos','day'], sort=False, as_index=False)['n'].max()
    cov = cov_df.set_index(['pos','day'])['n'].sort_index()
    records = []
    for pos, grp in tqdm(df_long.groupby("pos"),
                         desc=pat, position=idx+1, leave=False):
        days = sorted(grp['day'].unique())
        num = den = 0.0
        for d0, d1 in itertools.combinations(days, 2):
            def build_freq(day):
                sub = grp[grp['day'] == day][['var','k','n']]
                if DEBUG_WORKER:
                    dup_vars = sub['var'].duplicated(keep=False)
                    if dup_vars.any():
                        dup_count = dup_vars.sum()
                        examples = sub['var'][dup_vars].unique()[:20]
                        print(f"[DEBUG] Duplicate var values in f0/f1 for {pat} pos={pos} day={day}: {dup_count} rows, examples={examples}, sub_size={len(sub)}")
                g = sub.groupby('var', sort=False, observed=True)[['k','n']].sum()
                return g['k'] / g['n']

            f0 = build_freq(d0)
            f1 = build_freq(d1)

            if DEBUG_WORKER:
                dup_f0 = f0.index.duplicated(keep=False)
                dup_f1 = f1.index.duplicated(keep=False)
                if dup_f0.any():
                    print(f"[DEBUG] f0 has duplicate index values for {pat} pos={pos} day={d0}: {dup_f0.sum()} duplicates, examples={f0.index[dup_f0][:20]}, size={len(f0)}")
                if dup_f1.any():
                    print(f"[DEBUG] f1 has duplicate index values for {pat} pos={pos} day={d1}: {dup_f1.sum()} duplicates, examples={f1.index[dup_f1][:20]}, size={len(f1)}")

            allv = list(f0.index) + list(f1.index)
            if DEBUG_WORKER:
                dup_allv = pd.Index(allv).duplicated()
                if dup_allv.any():
                    print(f"[DEBUG] allv has {dup_allv.sum()} duplicates for {pat} pos={pos} days=({d0},{d1})")
            allv = pd.Index(allv).unique()
            f0 = f0.reindex(allv, fill_value=0)
            f1 = f1.reindex(allv, fill_value=0)

            if DEBUG_WORKER:
                try:
                    assert not f0.index.duplicated().any()
                    assert not f1.index.duplicated().any()
                except AssertionError:
                    print(f"[DEBUG] Duplicate indices detected after reindex for {pat} pos={pos} days=({d0},{d1})")

            if DEBUG_WORKER:
                if f0.index.duplicated().any() or f1.index.duplicated().any():
                    print(f"[DEBUG] Duplicate indices remain after collapse for {pat} pos={pos} days=({d0},{d1})")
            tvd = 0.5 * (f0.subtract(f1).abs().sum())
            
            # Normalizar por tiempo transcurrido entre días
            time_diff = abs(d1 - d0)
            tvd_per_day = tvd / time_diff if time_diff > 0 else tvd
            
            key0 = (int(pos), int(d0))
            key1 = (int(pos), int(d1))
            if DEBUG_WORKER:
                a = cov.get(key0, 0)
                b = cov.get(key1, 0)
                if isinstance(a, pd.Series) or isinstance(b, pd.Series):
                    print(f"[DEBUG] cov lookup returned Series for {pat} file={pth} pos={pos} d0={d0} d1={d1}")
                    print(f"[DEBUG] type(a)={type(a)} type(b)={type(b)}")
                    if isinstance(a, pd.Series):
                        print(f"[DEBUG] a len={len(a)} idx={list(a.index[:5])} vals={list(a.values[:5])}")
                    if isinstance(b, pd.Series):
                        print(f"[DEBUG] b len={len(b)} idx={list(b.index[:5])} vals={list(b.values[:5])}")
                if cov.index.has_duplicates:
                    dup = cov.index.duplicated(keep=False)
                    print(f"[DEBUG] cov index has duplicates: {dup.sum()}")
            a = cov.get(key0, 0)
            b = cov.get(key1, 0)
            if isinstance(a, pd.Series):
                a = a.max()
            if isinstance(b, pd.Series):
                b = b.max()
            w = min(int(a), int(b))
            num += tvd_per_day * w
            den += w
        variability = num/den if den>0 else 0
        records.append({"__patient__": pat, "pos": pos, "variability": variability})
    return pd.DataFrame(records)


def main():
    # ───────── CLI ─────────
    ap = argparse.ArgumentParser()
    ap.add_argument("-d","--dir",required=True,help="folder *_selected_ranked.tsv")
    ap.add_argument("-g","--gb",default="sequence.gb",help="GenBank file")
    ap.add_argument("--min",type=int,default=1,help="min patients per site")
    ap.add_argument("-o","--out",default="selection_UPGMA.html")
    ap.add_argument("--cache", default=None, help="Ruta opcional a caché de variabilidad (.pkl, .txt o .pyc)")
    ap.add_argument("--debug",action="store_true")
    ap.add_argument("-w","--workers",type=int,default=os.cpu_count(),
                   help="número de procesos para calcular variabilidad")
    ap.add_argument(
        "--variability-decimals",
        type=int,
        default=4,
        help=(
            "Decimales para serializar variabilidad por posición (var-data y caché). "
            "Menos decimales reduce tamaño HTML/cache; más decimales conserva precisión."
        ),
    )
    if hasattr(argparse, "BooleanOptionalAction"):
        ap.add_argument(
            "--minify-html-inline",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Minifica el HTML inline final para reducir tamaño (sin externalizar datos).",
        )
    else:
        # Compatibilidad con Python < 3.9
        ap.add_argument(
            "--minify-html-inline",
            action="store_true",
            default=True,
            help="Minifica el HTML inline final para reducir tamaño (sin externalizar datos).",
        )
        ap.add_argument(
            "--no-minify-html-inline",
            dest="minify_html_inline",
            action="store_false",
            help="Desactiva minificación inline del HTML final.",
        )
    args = ap.parse_args()
    if args.variability_decimals < 0:
        sys.exit("❌ --variability-decimals debe ser >= 0")
    global DEBUG_WORKER
    DEBUG_WORKER = args.debug

    methodology_content = """SELECTION ANALYSIS

1) RG tagging and BAM preparation
Addition of read groups (RG: ID, LB, PL, PU, SM), coordinate sorting, and indexing using Picard (AddOrReplaceReadGroups).

2) Per-position variant profiling
Generation, for each genomic position, of counts for A/C/G/T/N and indels with bam-readcount using:
-q 20 (minimum mapping quality),
-b 20 (minimum base quality),
-w 0 (indel window; reports only at the start nucleotide).

3) Integrity control and case selection
Reorganization of outputs and safety checks by number of positions (how many passed bam-readcount; >29,000 positions). Identification of time points that have different tokens but share a date (1,2; not 1a, 1b). All samples/patients proceed to downstream analyses.

4) Position×time matrices (per patient)
Execution of built_time_mutation_tsv.py, which, from bam-readcount, builds tables where each row is a position and each column (in increasing day order) encodes: days since time 0, coverage, and counts of A/C/T/G/N/ins/del.

5) Gene/codon/AA annotation
Execution of annotate_patient_tables.py. Reading NC_045512.2.gb (GenBank), each table row (position) is annotated with gene, codon position (1,2,3; allowing simultaneous positions in template jumps), and amino-acid index. This is specifically adapted for SARS-CoV-2 and that reference due to the template jump and the +ssRNA genome; it must be modified to run on other viruses.

6) Selection by variant×position trajectories
Execution of SRASel.py. For each position and variant (A, C, T, G, N, INS, DEL; different indels are treated as distinct variants) a trajectory is plotted. A depth-weighted linear fit across days is performed (accounting for days since day 0) to obtain a slope β (change in frequency per day). A permutation p-value is computed (reshuffling dates 2000 times), and positions are labeled as selected if the frequency across the entire time range changed by ≥3% and p-value < 0.05; they are labeled positive or negative if β>0 or β<0, respectively. Additionally, for each such position, the codon and amino acid are inferred by consensus from the two neighboring positions—only that are clear are accepted: frequency ≥95% for the variant in the neighbor; it must not be an indel, and there must be only one selected position in the codon. If the codon/amino-acid is not clear it is annotated as ambiguous.

Result plots are created (openable from the patients’ plot pop-ups in HTML) for all positions with a selected variant—showing the fitted slope and time points of all variants (point size encodes total and variant coverage at each time). Selected variants are colored red (β>0), blue (β<0), or gray (p-value > 0.05). Two additional guide lines are drawn: the mean β over all selected slopes for that patient (black dashed) and the mean β over all slopes, including non-selected (gray dashed).

Genomic plots are also created (openable from the pop-ups in HTML). These show a scatter of selected sites along the genome (red = β>0; blue = β<0) against their slope. In addition, sign-specific density curves are drawn and scaled to visualization (axis values are not meaningful): the curve reflects only the existence of selected sites, not slope magnitude, and is obtained by summing Gaussian kernels with σ = 350 bp per point (≈68% of each kernel’s mass within ±350 bp). The same two guide lines (black dashed for selected slopes’ mean; gray dashed for all slopes’ mean) are shown. A gene track based on NC_045512.2.gb is included. This plot visually identifies regions under stronger selection in specific patients without averaging.

Summary files required by the next script are also produced.

7) Interactive panel and complementary metrics (Selection.py)
An HTML report with two additional parameters is produced:

– Temporal variability per position [classifying non-variable positions – TVD]: Weighted average of the absolute frequency subtraction of all possible pairs of each iSNV at each position – the weighting is by the coverage of the time with the lowest coverage in each pair and TVDs are normalized by time difference. After that, positions below the 5th percentile of variability are marked non-variable; the rest are “neutral” or selected (if they have selected slopes). If a position is both non-variable and selected, selected prevails.
– UPGMA showing relationships between patients taking into account the intersections of selected positions [differentiating the variant and the direction of selection].

Figures in this HTML:

Figure 1 (heatmap): the number of selected positions (pull the ratio toward 1) versus non-variable positions (pull the ratio toward 0) is shown. The ratio is calculated as 0.5 + 0.5×(#selected/total) and 0.5 − 0.5×(#non_variable/total). This is computed over nucleotide positions and values (both in protein and gene views). There are caveats here: when selected sites exist, non-selected no longer count, but they should. An average per gene/protein of all selected ratios is also shown in the HTML. A hover is provided with gene/protein information per patient; clicking tiles opens position-level (nucleotide and amino-acid) information and the plots generated in the previous step for each selected position.

Figure 2 (heatmap): the mean of slope values for the selected sites is displayed. The hover shows how many positions contribute to that mean and the number of slopes (a position can have >1 selected variant). An average per gene/protein of all ratios is also shown, scoring as 0 those genes/proteins for patients without selected positions.

Figure 3: the mean variability (TVD) per site among the selected patients is shown. Selected positions are highlighted with points, with color encoding the number of patients in which each position was selected. A genes/proteins track (per the active view) is shown at the bottom from NC_045512.2.gb. The hover of a selected position provides detailed information; clicking opens a comprehensive table about selection at that site and gives access to the tables/plots from the previous step.

Figure 4: UPGMA tree constructed to visualize the relationships among patients according to their profiles of selected genomic sites. Patient subsets can be defined either manually or through metadata-based filters, and the corresponding patients are highlighted (green) in the tree. This visualization facilitates the exploration of similarity patterns and the identification of potential clusters of patients sharing common evolutionary or genomic features.

Figure 5: stacked bars summarize changes found in selected positions for the chosen patients and genes/proteins. Each bar condenses distinguishable amino-acid changes by color. The hover expands details; clicking opens a summary table of selection per color segment per bar, with links to the plots generated in the previous step.

Additionally, the HTML provides filters for patients (manuals or by immunological characteristics or sample location), genes/proteins (as discussed), and a minimum number of patients per site (for convenience, points and bars are hidden according to how many patients share them).

Bibliography:

Bam-readcount: Khanna A, Larson DE, Srivatsan SN, Mosior M, Abbott TE, Kiwala S, Ley TJ, Duncavage EJ, Walter MJ, Walker JR, Griffith OL, Griffith M, Miller CA. Bam-readcount - rapid generation of basepair-resolution sequence metrics. ArXiv [Preprint]. 2021 Jul 27:arXiv:2107.12817v1. PMID: 34341766; PMCID: PMC8328062.

Picard: no citation paper: https://broadinstitute.github.io/picard/

Biopython: Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.

NumPy: Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, Wieser E, Taylor J, Berg S, Smith NJ, Kern R, Picus M, Hoyer S, van Kerkwijk MH, Brett M, Haldane A, Del Río JF, Wiebe M, Peterson P, Gérard-Marchant P, Sheppard K, Reddy T, Weckesser W, Abbasi H, Gohlke C, Oliphant TE. Array programming with NumPy. Nature. 2020 Sep;585(7825):357-362. doi: 10.1038/s41586-020-2649-2. Epub 2020 Sep 16. PMID: 32939066; PMCID: PMC7759461.

McKinney, Wes. (2010). Data Structures for Statistical Computing in Python. 56-61. 10.25080/Majora-92bf1922-00a.

Plotly (htmls): Plotly Technologies Inc. Title: Collaborative data science. Publisher: Plotly Technologies Inc. Place of publication: Montréal, QC. Date of publication: 2015. URL: https://plot.ly"""

    # ───────── TSVs (selected_ranked) ─────────
    tsv_paths = sorted(glob.glob(os.path.join(args.dir,"*_selected_ranked.tsv")))
    if not tsv_paths:
        sys.exit("❌  No TSVs found")

    dfs, patients = [], []
    for pth in tsv_paths:
        pat = os.path.basename(pth).split("_selected_ranked.tsv")[0]
        df  = pd.read_csv(
            pth, sep="\t",
            dtype={"pos": int},
            low_memory=False
        )
        df["__file__"] = pth
        if "pos" not in df.columns:
            continue
        if "aa_change" not in df.columns:
            df["aa_change"] = "X"
        df["__patient__"] = pat
        dfs.append(df.dropna(subset=["pos"]))
        patients.append(pat)
    if not dfs:
        sys.exit("❌  TSVs empty / bad columns")

    df_all = pd.concat(dfs, ignore_index=True)
    df_all["aa_idx"] = pd.to_numeric(df_all["aa_idx"], errors="coerce")
    df_all['beta'] = pd.to_numeric(df_all['beta'], errors='coerce')

    # ───────── Pacientes reales presentes en los datos cargados ─────────
    # Derivar SIEMPRE de df_all para que el heatmap muestre a todos los pacientes con datos.
    patients_all = sorted(df_all["__patient__"].unique().tolist())
    # A partir de aquí, usar siempre patients_all
    patients = patients_all #Se conservan solo los pacientes con selección, pero después se representan todos porque se suman con los de variabilidad que salen del long.tsv

    # ───────── Variabilidad con cache (acepta .pkl o .txt) ─────────
    output_html = args.out or "selection_UPGMA.html"
    cache_id = os.path.splitext(os.path.basename(output_html))[0]
    cache_id = cache_id.replace(" ", "_")
    cache_id = re.sub(r"[^A-Za-z0-9_-]", "", cache_id)
    cache_base_dir = os.path.join(args.dir, "cache_runs", cache_id)
    os.makedirs(cache_base_dir, exist_ok=True)
    cache_txt = os.path.join(cache_base_dir, "variability_cache.txt")
    cache_pkl = os.path.join(cache_base_dir, "variability_cache.pkl")
    # Permitir ruta arbitraria via --cache
    explicit_cache = args.cache
    cache_file = None
    if explicit_cache and os.path.exists(explicit_cache) and not args.debug:
        cache_file = explicit_cache
    elif os.path.exists(cache_txt) and not args.debug:
        # Preferir TXT (JSON) para reproducir exactamente la lógica de _direct_table
        cache_file = cache_txt
    elif os.path.exists(cache_pkl) and not args.debug:
        cache_file = cache_pkl

    if args.debug:
        print(f"[CACHE] cache_id = {cache_id}")
        print(f"[CACHE] cache_base_dir = {cache_base_dir}")
        print(f"[CACHE] cache_txt = {cache_txt}")
        print(f"[CACHE] cache_pkl = {cache_pkl}")
        print(f"[CACHE] cache_txt exists = {os.path.exists(cache_txt)}")

    def _round_var_data_values(data, ndigits):
        if not isinstance(data, dict):
            return data
        out = {}
        for patient, pos_map in data.items():
            if not isinstance(pos_map, dict):
                out[patient] = pos_map
                continue
            new_pos_map = {}
            for pos, val in pos_map.items():
                key = str(pos)
                try:
                    new_pos_map[key] = round(float(val), ndigits)
                except (TypeError, ValueError):
                    new_pos_map[key] = val
            out[patient] = new_pos_map
        return out

    if cache_file:
        print(f"[CACHE] Loading variability from {cache_file}")
        var_data = None
        # Intento 1: pickle binario
        try:
            import pickle
            with open(cache_file, "rb") as cf:
                var_data = pickle.load(cf)
        except Exception:
            pass
        # Intento 2: JSON texto
        if var_data is None:
            try:
                with open(cache_file, "r") as cf:
                    var_data = json.load(cf)
            except Exception:
                pass
        if var_data is None:
            sys.exit(f"❌ No se pudo leer el caché: {cache_file}. Usa .pkl (pickle) o .txt (JSON).")
        var_data = _round_var_data_values(var_data, args.variability_decimals)
    else:
        # calcular y luego guardar cache
        long_paths = sorted(glob.glob(os.path.join(args.dir, "*_long.tsv")))
        tasks = [
          (os.path.basename(p).split("_long.tsv")[0], p, idx)
          for idx, p in enumerate(long_paths)
        ]
        if args.debug and tasks:
            for pat, pth, _ in tasks[:5]:
                print(f"[DEBUG] cache_file_for_run={cache_txt} exists={os.path.exists(cache_txt)} pid={os.getpid()} pat={pat} file={pth}")
        if args.debug and tasks:
            print("[DEBUG] Running worker self-check on first task")
            _ = worker(tasks[0])
        with Pool(processes=args.workers) as pool:
            dfs_var = []
            for df in tqdm(pool.imap(worker, tasks),
                           total=len(tasks), desc="General", position=0):
                dfs_var.append(df)
        variability_df = pd.concat(dfs_var, ignore_index=True)
        var_data = (
            variability_df
            .groupby("__patient__")
            .apply(lambda g: g.set_index("pos")["variability"].to_dict())
            .to_dict()
        )
        var_data = _round_var_data_values(var_data, args.variability_decimals)
        # guardar cache por defecto en TXT (JSON)
        print(f"[CACHE] Saving variability to {cache_txt}")
        tmp_path = f"{cache_txt}.tmp.{os.getpid()}"
        with open(tmp_path, "w") as cf:
            cf.write(_json_dumps_compact(var_data))
        os.replace(tmp_path, cache_txt)

    # ───────── Combinar pacientes de TSV y caché de variabilidad ─────────
    # Ahora que var_data está disponible, combinar ambas fuentes
    patients_from_tsv = sorted(df_all["__patient__"].unique().tolist())
    patients_from_cache = sorted(var_data.keys()) if var_data else []
    
    # Unir ambas listas para tener todos los pacientes disponibles
    patients_all = sorted(set(patients_from_tsv + patients_from_cache))
    
    # Avisos útiles si hay discrepancias
    missing_in_df = sorted(set(patients) - set(patients_from_tsv))
    missing_in_cache = sorted(set(patients) - set(patients_from_cache))
    
    if missing_in_df:
        sys.stderr.write(
            "⚠️  Advertencia: estos pacientes se detectaron por nombre de archivo "
            "pero no tienen filas válidas en df_all (pos/gene NaN o TSV vacío tras filtros): "
            f"{', '.join(missing_in_df)}\n"
        )
    
    if missing_in_cache:
        sys.stderr.write(
            "⚠️  Advertencia: estos pacientes están en los TSV pero no en el caché de variabilidad: "
            f"{', '.join(missing_in_cache)}\n"
        )
    
    # A partir de aquí, usar la lista combinada
    patients = patients_all
    
    print(f"[PATIENTS] Total patients available: {len(patients)}")
    print(f"[PATIENTS] From TSV: {len(patients_from_tsv)}")
    print(f"[PATIENTS] From cache: {len(patients_from_cache)}")
    print(f"[PATIENTS] Combined: {patients}")

    # ───────── Leer datos de Inmuno_location_data.tsv para filtros ─────────
    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    immuno_data_path = os.path.join(script_dir, "Inmuno_location_data.tsv")
    if not os.path.exists(immuno_data_path):
        # Buscar en el directorio actual
        immuno_data_path = os.path.join(os.getcwd(), "Inmuno_location_data.tsv")
    
    patient_filter_data = {}
    patient_filter_rows = {}
    filter_columns = []
    filter_unique_values = {}
    patient_name_mapping = {}
    
    if os.path.exists(immuno_data_path):
        print(f"[FILTERS] Loading immuno data from: {immuno_data_path}")
        try:
            immuno_df = pd.read_csv(immuno_data_path, sep="\t", low_memory=False)
            
            # Normalizar nombres de pacientes (008-A -> 8-A, etc.)
            def normalize_patient(patient_str):
                if pd.isna(patient_str):
                    return None
                parts = str(patient_str).split("-")
                if len(parts) >= 2:
                    try:
                        num_part = parts[0]
                        if num_part.isdigit() or (num_part.startswith("0") and num_part[1:].isdigit()):
                            if not any(c.isalpha() for c in num_part[1:]):
                                normalized_num = str(int(num_part))
                            else:
                                normalized_num = num_part
                            return f"{normalized_num}-{parts[1]}"
                    except ValueError:
                        pass
                return str(patient_str)
            
            immuno_df["Pacient_normalized"] = immuno_df["Pacient"].apply(normalize_patient)
            
            # Obtener columnas de filtro (todas excepto Time_point y Pacient)
            filter_columns = [col for col in immuno_df.columns if col not in ["Time_point", "Pacient", "Pacient_normalized"]]
            
            # Crear estructura de datos: paciente -> {columna -> valor}
            for _, row in immuno_df.iterrows():
                patient = row["Pacient_normalized"]
                if pd.isna(patient):
                    continue
                
                if patient not in patient_filter_data:
                    patient_filter_data[patient] = {}
                if patient not in patient_filter_rows:
                    patient_filter_rows[patient] = []

                row_entry = {}
                if "Time_point" in row:
                    row_entry["_time_point"] = row["Time_point"]
                
                for col in filter_columns:
                    value = row[col]
                    if pd.isna(value) or str(value).strip() == "":
                        if col not in patient_filter_data[patient]:
                            patient_filter_data[patient][col] = []
                        if "__NO_DATA__" not in patient_filter_data[patient][col]:
                            patient_filter_data[patient][col].append("__NO_DATA__")
                        row_entry[col] = ["__NO_DATA__"]
                    else:
                        if col not in patient_filter_data[patient]:
                            patient_filter_data[patient][col] = []
                        values = [v.strip() for v in str(value).split(";") if v.strip()]
                        for v in values:
                            if v not in patient_filter_data[patient][col]:
                                patient_filter_data[patient][col].append(v)
                        row_entry[col] = values
                
                patient_filter_rows[patient].append(row_entry)
            
            # Generar valores únicos por columna (sin duplicados)
            for col in filter_columns:
                unique_vals = set()
                for patient_data in patient_filter_data.values():
                    if col in patient_data:
                        for val in patient_data[col]:
                            if val != "__NO_DATA__":
                                unique_vals.add(val)
                filter_unique_values[col] = sorted(list(unique_vals))
            
            # Crear mapeo bidireccional de nombres de pacientes (para manejar variaciones como C03-C vs 030A-C)
            # Mapear tanto el nombre normalizado como el original
            for _, row in immuno_df.iterrows():
                original = str(row["Pacient"])
                normalized = row["Pacient_normalized"]
                if pd.notna(normalized):
                    # Mapear en ambas direcciones (usar listas para que sea JSON serializable)
                    if normalized not in patient_name_mapping:
                        patient_name_mapping[normalized] = []
                    if original not in patient_name_mapping[normalized]:
                        patient_name_mapping[normalized].append(original)
                    if normalized not in patient_name_mapping[normalized]:
                        patient_name_mapping[normalized].append(normalized)
                    # También mapear el original a sí mismo y al normalizado
                    if original not in patient_name_mapping:
                        patient_name_mapping[original] = []
                    if original not in patient_name_mapping[original]:
                        patient_name_mapping[original].append(original)
                    if normalized not in patient_name_mapping[original]:
                        patient_name_mapping[original].append(normalized)
            
            print(f"[FILTERS] Loaded data for {len(patient_filter_data)} patients")
            print(f"[FILTERS] Filter columns: {filter_columns}")
            print(f"[FILTERS] Created patient name mapping with {len(patient_name_mapping)} entries")
            
        except Exception as e:
            print(f"[FILTERS] Error loading immuno data: {e}")
            patient_filter_rows = {}
    else:
        print(f"[FILTERS] Inmuno_location_data.tsv not found at {immuno_data_path}")

    # ───────── GenBank → mat_peptide products ─────────
    def parse_gb_products(gb_path):
        prod = OrderedDict()
        cur_s = cur_e = None
        for ln in open(gb_path):
            if ln.startswith("     ") and re.search(r"\d+\.\.\d+", ln):
                nums = re.findall(r"\d+", ln)
                cur_s, cur_e = int(nums[0]), int(nums[-1])
                continue
            if "/product=" in ln:
                txt = ln.split("=",1)[1].strip().strip('"').lower()
                if "polyprotein" in txt:
                    continue
                alias = re.search(r"nsp\d+", txt)
                name  = alias.group(0) if alias else re.sub(
                            r"\b(protein|peptide)\b", "", txt).strip()
                if name and name not in prod and cur_s:
                    prod[name] = (cur_s, cur_e)
        return sorted([(s, e, n) for n,(s,e) in prod.items()], key=lambda x: x[0])

    # ───────── GenBank → genes ─────────
    def parse_gb_genes(gb_path):
        genes = OrderedDict()
        cur_s = cur_e = None
        for ln in open(gb_path):
            if ln.startswith("     ") and re.search(r"\d+\.\.\d+", ln):
                nums = re.findall(r"\d+", ln)
                cur_s, cur_e = int(nums[0]), int(nums[-1])
                continue
            if "/gene=" in ln:
                gene_name = ln.split("=",1)[1].strip().strip('"')
                if gene_name and gene_name not in genes and cur_s:
                    genes[gene_name] = (cur_s, cur_e)
        return sorted([(s, e, n) for n,(s,e) in genes.items()], key=lambda x: x[0]) #S,N,E -> Start, end, name

    CDS_LIST = parse_gb_products(args.gb)
    GENE_LIST = parse_gb_genes(args.gb)
    print(f"CDS_LIST= {CDS_LIST}")
    print(f"GENE_LIST= {GENE_LIST}")
    
    def nt2prod(nt):
        for s,e,n in CDS_LIST: #Start, end, name
            if s <= nt <= e:
                return n
        return "unknown"
    
    def nt2gene(nt):
        for s,e,n in GENE_LIST: #Start, end, name
            if s <= nt <= e:
                return n
        return "unknown"
    
    def color_for_site(site: str, n_pos: int = 0, n_neg: int = 0) -> str:
        """
        site = 'K484 / 23012', '484E / 23012', 'amb484 / 23012', 'K484E / 23012'
        Colorea según la lógica de p_pos y p_neg:
        - Mixtos (positivos y negativos): naranja
        - Solo positivos: rojo
        - Solo negativos: azul
        - Ambiguos: gris
        """   
        # Nueva lógica basada en n_pos y n_neg
        if n_pos > 0 and n_neg > 0:
            return '#ff7f0e'  # naranja para mixtos
        elif n_pos > 0 and n_neg == 0:
            return '#E74C3C'  # rojo para solo positivos
        elif n_pos == 0 and n_neg > 0:
            return '#1f77b4'  # azul para solo negativos
        else:
            return '#8f8f8f'  # gris para ambiguos
    
    def combine_aa_changes(aa_changes):
        """Combina múltiples cambios aminoacídicos en una sola etiqueta
        Siempre genera formato: letra + número + sufijo
        """
        if len(aa_changes) <= 1:
            return aa_changes[0] if aa_changes else ""
        
        unique_changes = list(set(aa_changes))
        
        # Definir función de parsing (usada tanto para 2 como para 3+ cambios)
        def parse_change(change):
            """Extrae letra, número y sufijo de un cambio"""
            # Patrones posibles:
            # [A-Z]\d+ -> letra + número
            # \d+[A-Z] -> número + letra
            # [A-Z*]\d+ -> letra/stop + número
            # \d+[A-Z*] -> número + letra/stop
            # (indel)amb\d+ -> indel antes, amb después
            # \d+amb(indel) -> número, amb, indel después
            # amb\d+amb -> amb + número + amb
            
            # Caso: (indel)amb\d+ (formato: (indel)amb614)
            if re.match(r'^\(indel\)amb\d+$', change):
                number_match = re.search(r'\d+', change)
                if number_match:
                    return "(indel)amb", number_match.group(), ""
                else:
                    return change, "", ""
            # Caso: \d+amb(indel) (formato: 614amb(indel))
            elif re.match(r'^\d+amb\(indel\)$', change):
                number_match = re.search(r'\d+', change)
                if number_match:
                    return "", number_match.group(), "amb(indel)"
                else:
                    return change, "", ""
            # Caso: (indel)[A-Z]\d+ (formato: (indel)K614)
            elif re.match(r'^\(indel\)[A-Z]\d+$', change):
                letter_match = re.search(r'[A-Z]', change)
                number_match = re.search(r'\d+', change)
                if letter_match and number_match:
                    return f"(indel){letter_match.group()}", number_match.group(), ""
                else:
                    return change, "", ""
            # Caso: (indel)[A-Z*]\d+ (formato: (indel)*614 para stop codons)
            elif re.match(r'^\(indel\)[A-Z\*]\d+$', change):
                letter_match = re.search(r'[A-Z\*]', change)
                number_match = re.search(r'\d+', change)
                if letter_match and number_match:
                    return f"(indel){letter_match.group()}", number_match.group(), ""
                else:
                    return change, "", ""
            # Caso: \d+[A-Z](indel) (formato: 614K(indel))
            elif re.match(r'^\d+[A-Z]\(indel\)$', change):
                number_match = re.search(r'\d+', change)
                letter_match = re.search(r'[A-Z]', change)
                if number_match and letter_match:
                    return "", number_match.group(), f"{letter_match.group()}(indel)"
                else:
                    return "", "", change
            # Caso: \d+[A-Z*](indel) (formato: 614*(indel) para stop codons)
            elif re.match(r'^\d+[A-Z\*]\(indel\)$', change):
                number_match = re.search(r'\d+', change)
                letter_match = re.search(r'[A-Z\*]', change)
                if number_match and letter_match:
                    return "", number_match.group(), f"{letter_match.group()}(indel)"
                else:
                    return "", "", change
            # Caso: [A-Z*]\d+ (formato: K614 o *614)
            elif re.match(r'^[A-Z\*]\d+$', change):
                letter_match = re.search(r'[A-Z\*]', change)
                number_match = re.search(r'\d+', change)
                if letter_match and number_match:
                    return letter_match.group(), number_match.group(), ""
                else:
                    return "", "", change
            # Caso: \d+[A-Z*] (formato: 614K o 614*)
            elif re.match(r'^\d+[A-Z\*]$', change):
                number_match = re.search(r'\d+', change)
                letter_match = re.search(r'[A-Z\*]', change)
                if number_match and letter_match:
                    return "", number_match.group(), letter_match.group()
                else:
                    return "", "", change
            # Caso: amb\d+amb (formato: amb614amb) - MÁS ESPECÍFICO PRIMERO
            elif re.match(r'^amb\d+amb$', change):
                number_match = re.search(r'\d+', change)
                if number_match:
                    return "amb", number_match.group(), "amb"
                else:
                    return "", "", change
            # Caso: amb\d+ (formato: amb614)
            elif re.match(r'^amb\d+$', change):
                number_match = re.search(r'\d+', change)
                if number_match:
                    return "amb", number_match.group(), ""
                else:
                    return "", "", change
            # Caso: \d+amb (formato: 614amb)
            elif re.match(r'^\d+amb$', change):
                number_match = re.search(r'\d+', change)
                if number_match:
                    return "", number_match.group(), "amb"
                else:
                    return "", "", change
            elif change == 'ambND':
                # Casos especiales con ND
                return "amb", "ND", ""
            elif change == 'NDamb':
                # Casos especiales con ND
                return "", "ND", "amb"
            elif change == '(indel)ambND':
                # Caso: (indel)ambND
                return "(indel)amb", "ND", ""
            elif change == 'NDamb(indel)':
                # Caso: NDamb(indel)
                return "", "ND", "amb(indel)"
            elif change == '(indel)NDamb':
                # Caso: (indel)NDamb
                return "(indel)", "ND", "amb"
            elif change == 'ambND(indel)':
                # Caso: ambND(indel)
                return "amb", "ND", "(indel)"
        
        # Lógica unificada para 2 o más cambios
        if len(unique_changes) >= 2:
            # Parsear TODOS los cambios
            parsed = [parse_change(ch) for ch in unique_changes]
            
            # Extraer componentes
            letters = [p[0] for p in parsed]
            numbers = [p[1] for p in parsed]
            suffixes = [p[2] for p in parsed]
            
            # Verificar si todos comparten el mismo número
            unique_numbers = set(n for n in numbers if n)
            if len(unique_numbers) != 1:
                # Números diferentes → concatenar ordenado
                return "".join(sorted(unique_changes))
            
            num = unique_numbers.pop()
            
            # Limpiar letras y sufijos de (indel) para análisis
            letters_clean = [l.replace("(indel)", "") if l else "" for l in letters]
            suffixes_clean = [s.replace("(indel)", "") if s else "" for s in suffixes]
            
            # Detectar presencia de indels
            has_indel_start = any("(indel)" in str(l) for l in letters)
            has_indel_end = any("(indel)" in str(s) for s in suffixes)
            
            # Determinar letra inicial (antes del número)
            unique_letters = set(letters_clean) - {"", "amb"}
            if len(unique_letters) > 1:
                # Múltiples letras diferentes → amb
                letter_start = "amb"
            elif len(unique_letters) == 1:
                # Una sola letra específica
                letter_start = unique_letters.pop()
            elif "amb" in letters_clean:
                # Solo amb
                letter_start = "amb"
            else:
                # Sin letras → nada por defecto
                letter_start = ""
            
            # Determinar letra final (después del número)
            unique_suffixes = set(suffixes_clean) - {"", "amb"}
            if len(unique_suffixes) > 1:
                # Múltiples sufijos diferentes → amb
                letter_end = "amb"
            elif len(unique_suffixes) == 1:
                # Un solo sufijo específico
                letter_end = unique_suffixes.pop()
            elif "amb" in suffixes_clean:
                # Solo amb
                letter_end = "amb"
            else:
                # Sin sufijos → nada por defecto
                letter_end = ""
            
            # Construir resultado con indels si aplica
            if has_indel_start and has_indel_end:
                return f"(indel){letter_start}{num}{letter_end}(indel)"
            elif has_indel_start:
                return f"(indel){letter_start}{num}{letter_end}"
            elif has_indel_end:
                return f"{letter_start}{num}{letter_end}(indel)"
            else:
                return f"{letter_start}{num}{letter_end}"
        else:
            return "".join(unique_changes)
    
    df_all["product"] = df_all["pos"].map(nt2prod)
    df_all["gene"] = df_all["pos"].map(nt2gene)

    if args.debug:
        print("[protein view] products =", ", ".join(n for _,_,n in CDS_LIST))

    palette = (go.Figure().layout.template.layout.colorway or
               ["#636EFA","#EF553B","#00CC96","#AB63FA","#FFA15A",
                "#19D3F3","#FF6692","#B6E880","#FF97FF","#FECB52"])

    def wrap(t, w=15):
        t = str(t).replace("_", " ")
        return "<br>".join(re.findall(f".{{1,{w}}}", t)) if len(t) > w else t

    def create_ratio_heatmap(mat_ratio, title, colorscale, zmin=0.0, zmax=None, zmid=None, detailed_data=None):
        """Crea un heatmap de proporciones con la misma estructura que los heatmaps originales"""
        # Calcular promedio por columna (como en los heatmaps originales)
        avg_row = mat_ratio.mean(axis=0)
        z_with_avg = [avg_row.values.tolist()] + mat_ratio.values.tolist()
        y_with_avg = ["Average"] + mat_ratio.index.tolist()
        
        # Calcular el rango real de los datos
        all_values = []
        for row in z_with_avg:
            all_values.extend([v for v in row if not pd.isna(v)])
        
        if all_values:
            real_max = max(all_values)
            real_min = min(all_values)
        else:
            real_max = 1.0
            real_min = 0.0
        
        # Usar el valor máximo real si no se especifica zmax
        if zmax is None:
            zmax = real_max
        
        # Calcular zmid basado en los datos reales
        if zmid is None:
            zmid = (real_min + real_max) / 2
        
        # Crear hover text similar al original
        hov_with_avg = []
        
        # Hover para la fila Average
        avg_hover_row = []
        for cat in mat_ratio.columns:
            avg_val = avg_row[cat]
            if not pd.isna(avg_val):
                # Calcular totales para Average
                total_selected = sum(detailed_data[patient][cat]['selected_count'] for patient in mat_ratio.index if cat in detailed_data[patient])
                total_non_variable = sum(detailed_data[patient][cat]['non_variable_count'] for patient in mat_ratio.index if cat in detailed_data[patient])
                total_positive = sum(detailed_data[patient][cat]['positive_count'] for patient in mat_ratio.index if cat in detailed_data[patient])
                total_negative = sum(detailed_data[patient][cat]['negative_count'] for patient in mat_ratio.index if cat in detailed_data[patient])
                total_positions = sum(detailed_data[patient][cat]['total_count'] for patient in mat_ratio.index if cat in detailed_data[patient])
                
                if title == "Selected/Total":
                    txt = f"<b>Average</b> × <b>{cat}</b><br>{title} = {avg_val:.4f}<br><b>Nucleotides:</b><br>• Selected: {total_selected}<br>• Total: {total_positions}<br><i>Calculated from {len(mat_ratio.index)} patients</i>"
                elif title == "Non-Variable/Total":
                    txt = f"<b>Average</b> × <b>{cat}</b><br>{title} = {avg_val:.4f}<br><b>Nucleotides:</b><br>• Non-Variable: {total_non_variable}<br>• Total: {total_positions}<br><i>Calculated from {len(mat_ratio.index)} patients</i>"
                elif title == "Positive/Total":
                    txt = f"<b>Average</b> × <b>{cat}</b><br>{title} = {avg_val:.4f}<br><b>Nucleotides:</b><br>• Positive: {total_positive}<br>• Total: {total_positions}<br><i>Calculated from {len(mat_ratio.index)} patients</i>"
                elif title == "Negative/Total":
                    txt = f"<b>Average</b> × <b>{cat}</b><br>{title} = {avg_val:.4f}<br><b>Nucleotides:</b><br>• Negative: {total_negative}<br>• Total: {total_positions}<br><i>Calculated from {len(mat_ratio.index)} patients</i>"
                else:
                    txt = f"<b>Average</b> × <b>{cat}</b><br>{title} = {avg_val:.4f}<br><i>Calculated from {len(mat_ratio.index)} patients</i>"
            else:
                txt = f"<b>Average</b> × <b>{cat}</b><br>No data"
            avg_hover_row.append(txt)
        hov_with_avg.append(avg_hover_row)
        
        # Hover para cada paciente
        for patient in mat_ratio.index:
            row_hov = []
            for cat in mat_ratio.columns:
                val = mat_ratio.loc[patient, cat]
                if not pd.isna(val) and detailed_data and patient in detailed_data and cat in detailed_data[patient]:
                    data = detailed_data[patient][cat]
                    if title == "Selected/Total":
                        txt = f"<b>{patient}</b> × <b>{cat}</b><br>{title} = {val:.4f}<br><b>Nucleotides:</b><br>• Selected: {data['selected_count']}<br>• Total: {data['total_count']}"
                    elif title == "Non-Variable/Total":
                        txt = f"<b>{patient}</b> × <b>{cat}</b><br>{title} = {val:.4f}<br><b>Nucleotides:</b><br>• Non-Variable: {data['non_variable_count']}<br>• Total: {data['total_count']}"
                    elif title == "Positive/Total":
                        txt = f"<b>{patient}</b> × <b>{cat}</b><br>{title} = {val:.4f}<br><b>Nucleotides:</b><br>• Positive: {data['positive_count']}<br>• Total: {data['total_count']}"
                    elif title == "Negative/Total":
                        txt = f"<b>{patient}</b> × <b>{cat}</b><br>{title} = {val:.4f}<br><b>Nucleotides:</b><br>• Negative: {data['negative_count']}<br>• Total: {data['total_count']}"
                    else:
                        txt = f"<b>{patient}</b> × <b>{cat}</b><br>{title} = {val:.4f}"
                else:
                    txt = f"<b>{patient}</b> × <b>{cat}</b><br>No data"
                row_hov.append(txt)
            hov_with_avg.append(row_hov)
        
        # Crear el heatmap con la misma estructura que los originales
        heatmap = go.Figure(go.Heatmap(
            z=z_with_avg,
            x=list(mat_ratio.columns),
            y=y_with_avg,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            zmid=zmid,
            hovertext=hov_with_avg,
            hovertemplate="%{hovertext}<extra></extra>",
            hoverongaps=False,
            showscale=True,
            colorbar=dict(title=title, len=0.87)
        ))
        
        # Aplicar el mismo layout que los heatmaps originales (copiado exactamente del var_heat)
        heatmap.update_layout(
            height=max(440, len(mat_ratio.index)*22+210),
            xaxis=dict(side="top", tickangle=-45, tickfont_size=9,
                       tickmode="array", tickvals=mat_ratio.columns,
                       ticktext=[wrap(c) for c in mat_ratio.columns], constrain="domain"),
            yaxis=dict(autorange="reversed"), plot_bgcolor="rgba(255,255,255,0.8)",
            paper_bgcolor="rgba(255,255,255,0)",
            margin=dict(l=190, r=50, t=20, b=80),
            legend=dict(
                font=dict(size=10),
                x=1.02,
                y=1,
                xanchor="left",
                yanchor="top",
                orientation="v",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                borderwidth=1
            )
        )
        
        return heatmap.to_plotly_json()

    def build(col):
        print(f"[DEBUG] === BUILD START for {col} ===")
        # Usar el mapeo pos_to_aa ya creado (ahora con posiciones compartidas)
        nt_to_aa_mapping = {}
        aa_to_nt_mapping = {}
        for nt_pos, aa_data in pos_to_aa.items():
            # Para compatibilidad, usar la primera posición aminoacídica
            if aa_data:
                nt_to_aa_mapping[nt_pos] = aa_data[0][1]  # (gene, aa_pos, frame) -> aa_pos
                for gene, aa_pos, frame in aa_data:
                    if aa_pos not in aa_to_nt_mapping:
                        aa_to_nt_mapping[aa_pos] = []
                    aa_to_nt_mapping[aa_pos].append(nt_pos)
        
        # ───────── EJE X (orden) robusto por vista ─────────
        if col == "gene":
            # 1) GenBank
            all_genes = [n for _,_,n in GENE_LIST]
            # 2) Fallback: lo que realmente haya en df_all['gene'] (excluye 'unknown'/NaN)
            if not all_genes:
                if "gene" in df_all.columns:
                    all_genes = [g for g in df_all["gene"].dropna().unique().tolist()
                                 if g and str(g).lower() != "unknown"]
            # 3) Último recurso: no dejarlo vacío
            if not all_genes:
                all_genes = ["ORF1ab","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10"]
            order = list(dict.fromkeys(all_genes))
        else:
            # product / proteínas
            all_prot = [n for _,_,n in CDS_LIST]
            if not all_prot and "product" in df_all.columns:
                all_prot = [p for p in df_all["product"].dropna().unique().tolist()
                            if p and str(p).lower() != "unknown"]
            if not all_prot:
                all_prot = ["nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10",
                            "nsp12","nsp13","nsp14","nsp15","nsp16","S","E","M","N"]
            order = list(dict.fromkeys(all_prot))

        # ───────── Nuevas visualizaciones para heatmaps ─────────
        # Función para calcular proporciones (se definirá después de var_detailed_data)
        def calculate_ratios_for_region(patient, region, col):
            """Calcula las proporciones para una región específica usando los datos ya calculados en var_detailed_data"""
            # Obtener posiciones nucleotídicas para esta región
            if col == "gene":
                region_nt_positions = []
                for start, end, gene_name in GENE_LIST:
                    if gene_name == region:
                        region_nt_positions = list(range(start, end + 1))
                        break
            else:  # col == "product"
                region_nt_positions = []
                for start, end, protein_name in CDS_LIST:
                    if protein_name == region:
                        region_nt_positions = list(range(start, end + 1))
                        break
            
            if not region_nt_positions:
                return 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0  # selected, non_variable, positive, negative, selected_count, non_variable_count, positive_count, negative_count, total_count
            
            total_nt_positions = len(region_nt_positions)
            
            # Usar los datos ya calculados en var_detailed_data
            if patient in var_detailed_data and region in var_detailed_data[patient]:
                detailed_data = var_detailed_data[patient][region]
                # var_detailed_data almacena listas de posiciones, necesitamos el conteo
                selected_nt_positions = len(detailed_data['selected_nt'])
                neutral_nt_positions = len(detailed_data['neutral_nt'])
                non_variable_nt_positions = len(detailed_data['non_variable_nt'])
            else:
                # Fallback si no hay datos (no debería pasar)
                selected_nt_positions = 0
                neutral_nt_positions = 0
                non_variable_nt_positions = 0
            
            # Contar posiciones positivas y negativas (en df_all) para heatmap 2
            # NUEVA LÓGICA: Una posición puede contar para ambos si tiene tanto positivos como negativos
            positive_positions = 0
            negative_positions = 0
            
            for nt_pos in region_nt_positions:
                pos_beta_data_all = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt_pos)]
                if not pos_beta_data_all.empty and pos_beta_data_all['beta'].notna().any():
                    # Verificar si hay valores positivos Y negativos en esta posición
                    has_positive = False
                    has_negative = False
                    
                    for _, row in pos_beta_data_all.iterrows():
                        beta_val = row['beta']
                        if pd.notna(beta_val):
                            if beta_val > 0:
                                has_positive = True
                            elif beta_val < 0:
                                has_negative = True
                    
                    # Contar independientemente: una posición puede ser tanto positiva como negativa
                    if has_positive:
                        positive_positions += 1
                    if has_negative:
                        negative_positions += 1
            
            # Calcular proporciones
            selected_ratio = selected_nt_positions / total_nt_positions if total_nt_positions > 0 else 0.0
            non_variable_ratio = non_variable_nt_positions / total_nt_positions if total_nt_positions > 0 else 0.0
            positive_ratio = positive_positions / total_nt_positions if total_nt_positions > 0 else 0.0
            negative_ratio = negative_positions / total_nt_positions if total_nt_positions > 0 else 0.0
            
            return selected_ratio, non_variable_ratio, positive_ratio, negative_ratio, selected_nt_positions, non_variable_nt_positions, positive_positions, negative_positions, total_nt_positions
        
        # Las matrices de proporciones se crearán después de var_detailed_data


        print(f"[DEBUG] === BUILD: Antes de construir heatmaps ===")
        
        # ───────── Calcular umbral dinámico de variabilidad (percentil 5) ─────────
        all_variability_values = []
        for patient in patients:
            patient_var_data = var_data.get(patient, {})
            for pos_str, var_val in patient_var_data.items():
                if var_val is not None:
                    all_variability_values.append(float(var_val))
        
        if all_variability_values:
            import numpy as np
            variability_threshold = np.percentile(all_variability_values, 5)
        else:
            variability_threshold = 0.001
        
        print(f"[DEBUG] Umbral dinámico de variabilidad (percentil 5): {variability_threshold:.6f}")
        
        # ───────── Heatmap 1 (β media) usando SOLO sitios variables ─────────
        # Conjuntos de posiciones variables por paciente (nt y aa)
        variable_nt_positions_by_patient = {}
        variable_aa_positions_by_patient = {}
        for _pat in patients:
            nt_vars = {int(p) for p, v in var_data.get(_pat, {}).items()
                       if v is not None and float(v) >= variability_threshold}
            variable_nt_positions_by_patient[_pat] = nt_vars
            # mapear a aa usando nt_to_aa_mapping (si existe)
            aa_vars = set()
            if nt_vars:
                for ntp in nt_vars:
                    aa_idx = nt_to_aa_mapping.get(ntp)
                    if aa_idx is not None and not pd.isna(aa_idx):
                        aa_vars.add(int(aa_idx))
            variable_aa_positions_by_patient[_pat] = aa_vars

        # 3) Filtrar df_all a filas que estén en posiciones variables (SIEMPRE por nucleótidos)
        is_variable_mask = df_all.apply(
            lambda r: int(r.pos) in variable_nt_positions_by_patient.get(r.__patient__, set()), axis=1
        )

        df_var = df_all[is_variable_mask]
        # Fallback si el filtro dejó todo vacío (evita heatmap totalmente blanco)
        if df_var.empty:
            print("[WARN] df_var vacío con umbral", variability_threshold,
                  "→ usando todos los sitios para el heatmap β.")
            df_var = df_all.copy()

        # 4) Construir matriz β̄ solo con sitios variables (para conteos) y matriz con TODOS los datos para pintar betas
        mat = (df_var.groupby(["__patient__", col])["beta"]
               .mean()
               .unstack(col)
               .reindex(index=patients, columns=order))
        # Matriz con TODOS los datos (sin filtrar variabilidad) para el heatmap 2 (betas ±)
        mat_all = (df_all.groupby(["__patient__", col])["beta"]
                   .mean()
                   .unstack(col)
                   .reindex(index=patients, columns=order))
        # Asegurar numérico (evita dtype object que Plotly no pinta)
        mat = mat.apply(pd.to_numeric, errors="coerce")
        mat_all = mat_all.apply(pd.to_numeric, errors="coerce")
        # Si faltan columnas (None), forzar presencia
        if mat.columns.isnull().any():
            mat = mat.rename(columns={None: "unknown"}).reindex(columns=order, fill_value=float('nan'))
        
        print(f"[MATRICES] mat shape: {mat.shape}, patients: {list(mat.index)}")
        print(f"[MATRICES] mat_all shape: {mat_all.shape}, patients: {list(mat_all.index)}")
        


        # Calcular media por posición (para la gráfica de variabilidad agregada)
        # Para Heatmap 1, no se usa df_avg; mantenemos para otras vistas
        df_avg = df_all.groupby("pos")["beta"].mean().dropna().round(4)
        var_data["__avg__"] = df_avg.to_dict()

        # β⁺ / β⁻ lists per cell
        def posneg(g):
            """
            Construye dos sets de posiciones aminoacídicas (p = beta>0, n = beta<0)
            usando siempre r.aa_idx. Descarta (y avisa) si aa_idx es NaN.
            """
            p, n = set(), set()
            for aa_pos, b in g[["aa_idx", "beta"]].values:
                # saltar betas no numéricas
                if pd.isna(b):
                    continue
                # si por alguna razón faltara aa_idx, avisamos y descartamos
                if pd.isna(aa_pos):
                    sys.stderr.write(
                        f"⚠️  ¡Atención! en posneg faltó aa_idx para fila:\n{g.head(1)}\n"
                    )
                    continue
                num = int(aa_pos)
                (p if b > 0 else n).add(num)
            return sorted(p), sorted(n)

        lst = (df_all.groupby(["__patient__", col]).apply(posneg)
               .unstack(col).reindex_like(mat_all))

        hov = []
        for pat in mat_all.index:
            row = []
            for cat in mat_all.columns:
                if pd.notna(mat_all.at[pat, cat]):
                    # Calcular sitios seleccionados y número de slopes
                    selected_sites = 0
                    number_of_slopes = 0
                    
                    if col == "gene":
                        # gene: contar nt seleccionados en ese gen
                        if cat in [n for _,_,n in GENE_LIST]:
                            for s,e,n in GENE_LIST:
                                if n == cat:
                                    nt_range = range(s, e+1)
                                    break
                            # Contar posiciones seleccionadas (con beta) en este gen
                            for nt_pos in nt_range:
                                pos_beta_data = df_all[(df_all['__patient__'] == pat) & (df_all['pos'] == nt_pos)]
                                if not pos_beta_data.empty and pos_beta_data['beta'].notna().any():
                                    selected_sites += 1
                                    number_of_slopes += len(pos_beta_data[pos_beta_data['beta'].notna()])
                    else:
                        # protein: contar aa seleccionados en esa proteína
                        if cat in [n for _,_,n in CDS_LIST]:
                            for s,e,n in CDS_LIST:
                                if n == cat:
                                    nt_range = range(s, e+1)
                                    break
                            # Contar posiciones seleccionadas (con beta) en esta proteína
                            for nt_pos in nt_range:
                                pos_beta_data = df_all[(df_all['__patient__'] == pat) & (df_all['pos'] == nt_pos)]
                                if not pos_beta_data.empty and pos_beta_data['beta'].notna().any():
                                    selected_sites += 1
                                    number_of_slopes += len(pos_beta_data[pos_beta_data['beta'].notna()])
                    
                    # Debug print para casos específicos
                    #if pat == "016-B" and cat in ["ORF8", "N"]:
                    #    print(f"[DEBUG HOVER] {pat} × {cat}: selected_sites={selected_sites}, number_of_slopes={number_of_slopes}")
                    
                    txt = (
                        f"<b>{pat}</b> × <b>{cat}</b><br>β̄ (all sites) = {mat_all.at[pat,cat]:.4f}"
                        f"<br><i>Selected sites used</i>: {selected_sites}<br>"
                        f"<i>Number of slopes</i>: {number_of_slopes}"
                    )
                    row.append(txt)
                else:
                    row.append("")
            hov.append(row)

        # Escala de colores robusta
        import numpy as np
        # mat_all ya tiene el índice correcto, solo llenar NaN con 0.0
        safe_mat = mat_all.fillna(0.0)
        vmax = np.nanmax(np.abs(safe_mat.values))
        vmax = max(vmax, 0.001)
        
        # Debug
        #print("β matrix min/max:", mat_all.min().min(), mat_all.max().max())
        
        # Construir matriz de customdata para el popup
        def build_custom_matrix():
            cmat = []
            for pat in mat_all.index:
                row = []
                for cat in mat_all.columns:
                    if pd.isna(mat_all.at[pat,cat]) or pd.isna(lst.at[pat,cat]):
                        row.append("{}")
                    else:
                        p_pos, p_neg = lst.at[pat,cat]
                        row.append(json.dumps({
                            "p_pos": p_pos,
                            "p_neg": p_neg
                        }))
                cmat.append(row)
            return cmat

        custom_mat = build_custom_matrix()

        # ── Prepend average row (across patients) that is always visible
        # Calcular la media de todos los pacientes para cada columna
        avg_row = []
        for col_name in safe_mat.columns:
            col_values = safe_mat[col_name].dropna()  # Ignorar NaN
            if len(col_values) > 0:
                avg_val = col_values.mean()
            else:
                avg_val = 0.0
            avg_row.append(avg_val)
        
        # Agregar fila Average al principio
        z_with_avg = [avg_row] + safe_mat.values.tolist()
        y_with_avg = ["Average"] + list(mat_all.index)
        
        # Agregar fila Average al hover
        avg_hover_row = []
        for j, cat in enumerate(mat_all.columns):
            if j < len(avg_row):
                avg_val = avg_row[j]
                txt = f"<b>Average</b> × <b>{cat}</b><br>β̄ = {avg_val:.6f}<br><i>Calculated from {len(patients)} patients</i>"
            else:
                txt = f"<b>Average</b> × <b>{cat}</b><br>No data"
            avg_hover_row.append(txt)
        
        hov_with_avg = [avg_hover_row] + hov
        
        print(f"[HEATMAP] z_with_avg rows: {len(z_with_avg)}")
        print(f"[HEATMAP] y_with_avg labels: {y_with_avg}")
        
        # Ajustar colores para β: azul negativo, blanco 0, rojo positivo
        heat = go.Figure(go.Heatmap(
            z=z_with_avg,
            x=list(mat_all.columns),
            y=y_with_avg,
            colorscale=[[0.0, "#2166ac"], [0.5, "#f7f7f7"], [1.0, "#b2182b"]],
            zmid=0.0,
            zmin=-vmax,
            zmax=vmax,
            hovertext=hov_with_avg,
            hovertemplate="%{hovertext}<extra></extra>",
            customdata=custom_mat,
            showscale=True,
            colorbar=dict(title="β", len=0.87)))
        heat.update_layout(
            height=max(440, len(patients)*22+210),
            xaxis=dict(side="top", tickangle=-45, tickfont_size=9,
                       tickmode="array", tickvals=mat_all.columns,
                       ticktext=[wrap(c) for c in mat_all.columns], constrain="domain"),
            yaxis=dict(autorange="reversed"), plot_bgcolor="rgba(255,255,255,0.8)",
            paper_bgcolor="rgba(255,255,255,0)",
            margin=dict(l=190, r=50, t=20, b=80),
            legend=dict(
                font=dict(size=10),
                x=1.02,
                y=1,
                xanchor="left",
                yanchor="top",
                orientation="v",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                borderwidth=1
            ))
        


        # Crear heatmap de variabilidad con datos reales
        var_mat = pd.DataFrame(index=patients, columns=order, dtype=float)
        print(f"[DEBUG INIT] var_mat initialized with shape: {var_mat.shape}")
        print(f"[DEBUG INIT] var_mat initial values (first 3x3): {var_mat.iloc[:3, :3].values.tolist()}")
        
        # Función auxiliar para calcular ratio de selección según especificaciones
        def calculate_selection_ratio(selected_nt_positions, neutral_nt_positions, variable_nt_positions, total_nt_positions):
            # Balancear entre posiciones seleccionadas vs no variables, normalizado por total de posiciones
            # Escala: 0.0-0.5 = verde (no variables), 0.5 = blanco (neutral), 0.5-1.0 = naranja (seleccionadas)
            
            non_variable_positions = total_nt_positions - variable_nt_positions
            
            if total_nt_positions == 0:
                # Sin datos
                return 0.5  # Blanco (neutral)
            
            # Normalizar por el total de posiciones de la región
            # Esto permite comparar regiones de diferentes tamaños
            selected_ratio = selected_nt_positions / total_nt_positions
            non_variable_ratio = non_variable_positions / total_nt_positions
            
            # Calcular balance: positivo si más seleccionadas, negativo si más no variables
            balance = selected_ratio - non_variable_ratio
            
            # Convertir a escala 0.0-1.0 donde 0.5 es neutral
            # balance va de -1.0 (todas no variables) a +1.0 (todas seleccionadas)
            # Transformar a 0.0-1.0: (balance + 1.0) / 2.0
            selection_score = (balance + 1.0) / 2.0
            
            return selection_score
        
        # Usar el umbral dinámico ya calculado arriba (no recalcular)
        # El variability_threshold ya se calculó en las líneas 720-731
        
        # Calcular datos de variabilidad reales (basados en nucleótidos)
        for patient in patients:
            for region in order:
                # Obtener posiciones nucleotídicas para esta región
                if col == "gene":
                    # Para genes, usar EXACTAMENTE la misma lógica que proteínas
                    # Usar GENE_LIST para obtener el rango completo del gen
                    region_nt_positions = []
                    for start, end, name in GENE_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                else:
                    # Para proteínas, usar el mismo enfoque que genes pero con CDS_LIST
                    region_nt_positions = []
                    for start, end, name in CDS_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                
                if not region_nt_positions:
                    print(f"[DEBUG NO POSITIONS] {patient} × {region}: No positions found in {'GENE_LIST' if col == 'gene' else 'CDS_LIST'}")
                    var_mat.at[patient, region] = 0.5
                    print(f"[DEBUG ASSIGN 0.5] {patient} × {region}: assigned 0.5 (no positions)")
                    continue
                
                # Obtener datos de variabilidad para este paciente
                patient_var_data = var_data.get(patient, {})
                
                # Denominador: TODAS las posiciones nt de la región (coincidir con hover)
                total_nt_positions = len(region_nt_positions)
                variable_nt_positions = 0
                selected_nt_positions = 0
                neutral_nt_positions = 0

                for nt_pos in region_nt_positions:
                    nt_key = str(nt_pos)
                    
                    # PRIORIDAD 1: Verificar si está seleccionada (tiene beta) - SIEMPRE
                    pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt_pos)]
                    beta_has_value = (not pos_beta_data.empty) and (pos_beta_data['beta'].notna().any())
                    
                    if beta_has_value:
                        # Seleccionada → clasificar como SELECTED (independiente de variabilidad)
                        selected_nt_positions += 1
                        # También contar como variable si supera el threshold
                        if nt_key in patient_var_data:
                            variability_value = patient_var_data[nt_key]
                            if variability_value >= variability_threshold:
                                variable_nt_positions += 1
                    elif nt_key in patient_var_data:
                        # No seleccionada → clasificar por variabilidad
                        variability_value = patient_var_data[nt_key]
                        if variability_value >= variability_threshold:
                            # Variable pero no seleccionada → NEUTRAL
                            variable_nt_positions += 1
                            neutral_nt_positions += 1
                        # Si < threshold → implícitamente NON-VARIABLE (no se cuenta)


                # Debug: verificar qué valores se están pasando a calculate_selection_ratio
                #if region in ['ORF6', 'N', 'ORF8']:  # Solo para genes con datos
                #    print(f"[DEBUG RATIO] {patient} × {region}: selected={selected_nt_positions}, neutral={neutral_nt_positions}, variable={variable_nt_positions}, total={total_nt_positions}")
                
                # Usar exactamente el mismo ratio que mostramos en el hover
                selection_ratio = calculate_selection_ratio(
                    selected_nt_positions, neutral_nt_positions,
                    variable_nt_positions, total_nt_positions
                )
                
                #if region in ['ORF6', 'N', 'ORF8']:  # Solo para genes con datos
                #    print(f"[DEBUG RATIO] {patient} × {region}: ratio = {selection_ratio}")
                
                var_mat.at[patient, region] = float(selection_ratio)
                #print(f"[DEBUG ASSIGN RATIO] {patient} × {region}: assigned {selection_ratio}")
        
        # Solo reemplazar NaN con 0.5 cuando sea apropiado, forzar float y acotar a [0,1]
        # Solo acceder a columnas que existen
        #if 'ORF6' in var_mat.columns:
        #    print(f"[DEBUG BEFORE PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist()}")
        #if 'N' in var_mat.columns:
        #    print(f"[DEBUG BEFORE PROCESSING] var_mat values for N: {var_mat['N'].tolist()}")
        #if 'ORF8' in var_mat.columns:
        #    print(f"[DEBUG BEFORE PROCESSING] var_mat values for ORF8: {var_mat['ORF8'].tolist()}")
        
        var_mat = (
            var_mat.apply(pd.to_numeric, errors="coerce")
                   .astype(float)
                   .clip(lower=0.0, upper=1.0)
        )
        
        # Solo acceder a columnas que existen
        #if 'ORF6' in var_mat.columns:
        #    print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist()}")
        #if 'N' in var_mat.columns:
        #    print(f"[DEBUG AFTER PROCESSING] var_mat values for N: {var_mat['N'].tolist()}")
        #if 'ORF8' in var_mat.columns:
        #    print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF8: {var_mat['ORF8'].tolist()}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat min/max: {var_mat.min().min()}, {var_mat.max().max()}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG AFTER PROCESSING] var_mat values for ORF6: {var_mat['ORF6'].tolist() if 'ORF6' in var_mat.columns else 'Column not found'}")
        
        # Crear hover text detallado para el heatmap de variabilidad
        # Construir hover en el MISMO ORDEN que los datos finales: [Average] + [pacientes]
        var_hov_with_avg = []
        
        # 1. PRIMERO agregar hover de pacientes (en el mismo orden que var_mat.index)
        for patient in var_mat.index:  # Usar var_mat.index para mantener el orden correcto
            row_hov = []
            for region in order:
                # Obtener posiciones nucleotídicas para esta región
                if col == "gene":
                    # Para genes, usar EXACTAMENTE la misma lógica que proteínas
                    # Usar GENE_LIST para obtener el rango completo del gen
                    region_nt_positions = []
                    for start, end, name in GENE_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                else:
                    # Para proteínas, usar el mismo enfoque que genes pero con CDS_LIST
                    region_nt_positions = []
                    for start, end, name in CDS_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                
                if not region_nt_positions:
                    row_hov.append(f"<b>{patient}</b> × <b>{region}</b><br>No data")
                    continue
                
                # Obtener datos de variabilidad para este paciente
                patient_var_data = var_data.get(patient, {})
                
                # Contar posiciones nucleotídicas usando el umbral
                total_nt_positions = len(region_nt_positions)
                variable_nt_positions = 0
                selected_nt_positions = 0
                neutral_nt_positions = 0
                
                for nt_pos in region_nt_positions:
                    nt_pos_str = str(nt_pos)
                    
                    # PRIORIDAD 1: Verificar si está seleccionada (tiene beta) - SIEMPRE
                    pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt_pos)]
                    beta_has_value = (not pos_beta_data.empty) and pos_beta_data['beta'].notna().any()
                    
                    if beta_has_value:
                        # Seleccionada → clasificar como SELECTED (independiente de variabilidad)
                        selected_nt_positions += 1
                        # También contar como variable si supera el threshold
                        if nt_pos_str in patient_var_data:
                            variability_value = patient_var_data[nt_pos_str]
                            if variability_value >= variability_threshold:
                                variable_nt_positions += 1
                    elif nt_pos_str in patient_var_data:
                        # No seleccionada → clasificar por variabilidad
                        variability_value = patient_var_data[nt_pos_str]
                        if variability_value >= variability_threshold:
                            # Variable pero no seleccionada → NEUTRAL
                            variable_nt_positions += 1
                            neutral_nt_positions += 1
                        # Si < threshold → implícitamente NON-VARIABLE (no se cuenta)
                
                # Amino acid site classification by ANY nucleotide in its codon (independent categories)
                # Usar el nuevo mapeo pos_to_aa en lugar de nt_to_aa_mapping
                aa_positions_in_region = sorted({get_first_aa_position(nt_pos) for nt_pos in region_nt_positions if get_first_aa_position(nt_pos) is not None})
                aa_positions_in_region = [aa for aa in aa_positions_in_region if aa is not None]

                aa_selected_set = set()
                aa_neutral_set = set()
                aa_nonvar_set = set()

                for aa_pos in aa_positions_in_region:
                    # Construir nt_list usando el mapeo inverso de pos_to_aa
                    nt_list = [nt for nt in region_nt_positions if get_first_aa_position(nt) == aa_pos]
                    
                    # Contar nucleótidos por categoría
                    selected_count = 0
                    neutral_count = 0
                    nonvar_count = 0
                    
                    # Inicializar variables de clasificación
                    any_selected = False
                    any_neutral = False
                    any_nonvar = False
                    
                    for nt in nt_list:
                        key = str(nt)
                        
                        # PRIORIDAD 1: Verificar si está seleccionada (tiene beta) - SIEMPRE
                        pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt)]
                        beta_has_value = (not pos_beta_data.empty) and (pos_beta_data['beta'].notna().any())
                        
                        if beta_has_value:
                            # Seleccionada → clasificar como SELECTED (independiente de variabilidad)
                            selected_count += 1
                        elif key in patient_var_data:
                            # No seleccionada → clasificar por variabilidad
                            if patient_var_data[key] >= variability_threshold:
                                # Variable pero no seleccionada → NEUTRAL
                                neutral_count += 1
                            else:
                                # No variable y no seleccionada → NON-VARIABLE
                                nonvar_count += 1
                        else:
                            # Sin datos de variabilidad y no seleccionada → NON-VARIABLE
                            nonvar_count += 1
                    
                    # Clasificar aminoácido: SELECTED > NEUTRAL > NON_VARIABLE
                    total_nt = len(nt_list)
                    if selected_count > 0:
                        any_selected = True
                    elif neutral_count > 0:
                        any_neutral = True
                    else:
                        any_nonvar = True
                    
                    if any_selected:
                        aa_selected_set.add(aa_pos)
                    elif any_neutral:
                        aa_neutral_set.add(aa_pos)
                    else:
                        aa_nonvar_set.add(aa_pos)

                total_aa_positions = len(aa_positions_in_region)
                variable_aa_positions = len({aa for aa in aa_positions_in_region if (aa in aa_selected_set or aa in aa_neutral_set)})
                selected_aa_positions = len(aa_selected_set)
                neutral_aa_positions = len(aa_neutral_set)
                non_variable_aa_positions = len(aa_nonvar_set)
                
                # Calcular ratio de selección usando la función auxiliar
                selection_ratio = calculate_selection_ratio(
                    selected_nt_positions, neutral_nt_positions, 
                    variable_nt_positions, total_nt_positions
                )
                
                # IMPRIMIR VALORES EXACTOS DEL HEATMAP
                print(f"[HEATMAP HOVER] {patient} × {region}: SELECTED={selected_nt_positions}, NEUTRAL={neutral_nt_positions}, NON_VARIABLE={total_nt_positions - variable_nt_positions}, TOTAL={total_nt_positions}")
                
                hover_text = (
                    f"<b>{patient}</b> × <b>{region}</b><br>"
                    f"<b>Nucleotides:</b><br>"
                    f"• Total: {total_nt_positions}<br>"
                    f"• Variable: {variable_nt_positions} ({variable_nt_positions/max(total_nt_positions,1)*100:.1f}%)<br>"
                    f"• Selected: {selected_nt_positions} ({selected_nt_positions/max(total_nt_positions,1)*100:.1f}%)<br>"
                    f"• Neutral: {neutral_nt_positions} ({neutral_nt_positions/max(total_nt_positions,1)*100:.1f}%)<br>"
                    f"• Non-variable: {total_nt_positions - variable_nt_positions} ({(total_nt_positions - variable_nt_positions)/max(total_nt_positions,1)*100:.1f}%)<br>"
                    f"<b>Amino acids:</b><br>"
                    f"• Total: {total_aa_positions}<br>"
                    f"• Variable: {variable_aa_positions} ({variable_aa_positions/max(total_aa_positions,1)*100:.1f}%)<br>"
                    f"• Selected: {selected_aa_positions} ({selected_aa_positions/max(total_aa_positions,1)*100:.1f}%)<br>"
                    f"• Neutral: {neutral_aa_positions} ({neutral_aa_positions/max(total_aa_positions,1)*100:.1f}%)<br>"
                    f"• Non-variable: {non_variable_aa_positions} ({non_variable_aa_positions/max(total_aa_positions,1)*100:.1f}%)<br>"
                    f"Selection ratio: {selection_ratio:.3f}"
                )
                row_hov.append(hover_text)
            var_hov_with_avg.append(row_hov)
        
        # Crear datos detallados para cada celda del heatmap de variabilidad
        var_detailed_data = {}
        for i, patient in enumerate(patients):
            var_detailed_data[patient] = {}
            for j, region in enumerate(order):
                # Obtener posiciones nucleotídicas para esta región
                if col == "gene":
                    # Para genes, usar posiciones nucleotídicas del gen
                    region_nt_positions = []
                    for start, end, name in GENE_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                else:
                    # Para proteínas, usar el mismo enfoque que genes pero con CDS_LIST
                    region_nt_positions = []
                    for start, end, name in CDS_LIST:
                        if name == region:
                            region_nt_positions = list(range(start, end + 1))
                            break
                
                patient_var_data = var_data.get(patient, {})
                
                # Clasificar posiciones nucleotídicas
                selected_nt_positions = []
                neutral_nt_positions = []
                non_variable_nt_positions = []
                
                for nt_pos in region_nt_positions:
                    nt_pos_str = str(nt_pos)
                    
                    # PRIORIDAD 1: Verificar si está seleccionada (tiene beta) - SIEMPRE
                    pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt_pos)]
                    beta_has_value = (not pos_beta_data.empty) and (pos_beta_data['beta'].notna().any())
                    
                    if beta_has_value:
                        # Seleccionada → clasificar como SELECTED (independiente de variabilidad)
                        selected_nt_positions.append(nt_pos)
                    elif nt_pos_str in patient_var_data:
                        # No seleccionada → clasificar por variabilidad
                        variability_value = patient_var_data[nt_pos_str]
                        if variability_value >= variability_threshold:
                            # Variable pero no seleccionada → NEUTRAL
                            neutral_nt_positions.append(nt_pos)
                        else:
                            # No variable y no seleccionada → NON-VARIABLE
                            non_variable_nt_positions.append(nt_pos)
                    else:
                        # Sin datos de variabilidad y no seleccionada → NON-VARIABLE
                        non_variable_nt_positions.append(nt_pos)
                
                # Amino acid lists by ANY nucleotide status (independent categories)
                aa_selected_list = []
                aa_neutral_list = []
                aa_nonvar_list = []
                # Build sets per AA via codon nts - usar el nuevo mapeo pos_to_aa
                aa_positions_in_region = sorted({get_first_aa_position(nt_pos) for nt_pos in region_nt_positions if get_first_aa_position(nt_pos) is not None})
                aa_positions_in_region = [aa for aa in aa_positions_in_region if aa is not None]
                for aa_pos in aa_positions_in_region:
                    nt_list = [nt for nt in region_nt_positions if get_first_aa_position(nt) == aa_pos]
                    # Contar nucleótidos por categoría
                    selected_count = 0
                    neutral_count = 0
                    nonvar_count = 0
                    
                    for nt in nt_list:
                        key = str(nt)
                        
                        # PRIORIDAD 1: Verificar si está seleccionada (tiene beta) - SIEMPRE
                        pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt)]
                        beta_has_value = (not pos_beta_data.empty) and (pos_beta_data['beta'].notna().any())
                        
                        if beta_has_value:
                            # Seleccionada → clasificar como SELECTED (independiente de variabilidad)
                            selected_count += 1
                        elif key in patient_var_data:
                            # No seleccionada → clasificar por variabilidad
                            if patient_var_data[key] >= variability_threshold:
                                # Variable pero no seleccionada → NEUTRAL
                                neutral_count += 1
                            else:
                                # No variable y no seleccionada → NON-VARIABLE
                                nonvar_count += 1
                        else:
                            # Sin datos de variabilidad y no seleccionada → NON-VARIABLE
                            nonvar_count += 1
                    
                    # Clasificar aminoácido: SELECTED > NEUTRAL > NON_VARIABLE
                    if selected_count > 0:
                        aa_selected_list.append(aa_pos)
                    elif neutral_count > 0:
                        aa_neutral_list.append(aa_pos)
                    else:
                        aa_nonvar_list.append(aa_pos)


                
                var_detailed_data[patient][region] = {
                    'selected_nt': selected_nt_positions,
                    'neutral_nt': neutral_nt_positions,
                    'non_variable_nt': non_variable_nt_positions,
                    'selected_aa': sorted(set(aa_selected_list)),
                    'neutral_aa': sorted(set(aa_neutral_list)),
                    'non_variable_aa': sorted(set(aa_nonvar_list))
                }
        
        # Crear matrices de proporciones usando los datos ya calculados en var_detailed_data
        mat_selected_ratio = pd.DataFrame(index=patients, columns=order, dtype=float)
        mat_non_variable_ratio = pd.DataFrame(index=patients, columns=order, dtype=float)
        mat_positive_ratio = pd.DataFrame(index=patients, columns=order, dtype=float)
        mat_negative_ratio = pd.DataFrame(index=patients, columns=order, dtype=float)
        
        # Almacenar datos detallados para hover
        detailed_data = {}
        
        for patient in patients:
            detailed_data[patient] = {}
            for region in order:
                selected_ratio, non_variable_ratio, positive_ratio, negative_ratio, selected_count, non_variable_count, positive_count, negative_count, total_count = calculate_ratios_for_region(patient, region, col)
                mat_selected_ratio.loc[patient, region] = selected_ratio
                mat_non_variable_ratio.loc[patient, region] = non_variable_ratio
                mat_positive_ratio.loc[patient, region] = positive_ratio
                mat_negative_ratio.loc[patient, region] = negative_ratio
                
                # Guardar datos detallados para hover
                detailed_data[patient][region] = {
                    'selected_count': selected_count,
                    'non_variable_count': non_variable_count,
                    'positive_count': positive_count,
                    'negative_count': negative_count,
                    'total_count': total_count
                }
        
        # Agregar fila Average al heatmap de variabilidad
        # Calcular la media de todos los pacientes para cada columna
        var_avg_row = []
        for col_name in var_mat.columns:
            col_values = var_mat[col_name].dropna()  # Ignorar NaN
            if len(col_values) > 0:
                avg_val = col_values.mean()
            else:
                avg_val = 0.5  # Neutral (blanco)
            var_avg_row.append(avg_val)
        
        # ⬇️ Colocar "Average" AL PRINCIPIO (mismo criterio que en el heatmap de β)
        # Nota: con yaxis.autorange="reversed", la PRIMERA categoría de `y`
        # se dibuja ARRIBA. Por eso "Average" debe ir primero en z, y y hover.
        z_var_with_avg = [var_avg_row] + var_mat.values.tolist()
        y_var_with_avg = ["Average"] + var_mat.index.tolist()
        
        # 2. CREAR HOVER DE AVERAGE CON INFORMACIÓN DETALLADA
        var_avg_hover_row = []
        for j, cat in enumerate(var_mat.columns):
            if j < len(var_avg_row):
                avg_val = var_avg_row[j]
                
                # Calcular estadísticas agregadas para esta categoría
                total_nt_all_patients = 0
                variable_nt_all_patients = 0
                selected_nt_all_patients = 0
                neutral_nt_all_patients = 0
                non_variable_nt_all_patients = 0
                
                total_aa_all_patients = 0
                variable_aa_all_patients = 0
                selected_aa_all_patients = 0
                neutral_aa_all_patients = 0
                non_variable_aa_all_patients = 0
                
                # Obtener posiciones para esta categoría
                if col == "gene":
                    region_nt_positions = []
                    for start, end, name in GENE_LIST:
                        if name == cat:
                            region_nt_positions = list(range(start, end + 1))
                            break
                else:
                    region_nt_positions = []
                    for start, end, name in CDS_LIST:
                        if name == cat:
                            region_nt_positions = list(range(start, end + 1))
                            break
                
                if region_nt_positions:
                    # Calcular estadísticas agregadas
                    for patient in patients:
                        patient_var_data = var_data.get(patient, {})
                        
                        # Nucleótidos
                        total_nt_all_patients += len(region_nt_positions)
                        for nt_pos in region_nt_positions:
                            nt_pos_str = str(nt_pos)
                            if nt_pos_str in patient_var_data:
                                if patient_var_data[nt_pos_str] >= variability_threshold:
                                    variable_nt_all_patients += 1
                                    # Verificar selección
                                    pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt_pos)]
                                    if not pos_beta_data.empty and pos_beta_data['beta'].notna().any():
                                        selected_nt_all_patients += 1
                                    else:
                                        neutral_nt_all_patients += 1
                                else:
                                    non_variable_nt_all_patients += 1
                            else:
                                non_variable_nt_all_patients += 1
                        
                        # Aminoácidos - usar el nuevo mapeo pos_to_aa
                        aa_positions_in_region = sorted({get_first_aa_position(nt_pos) for nt_pos in region_nt_positions if get_first_aa_position(nt_pos) is not None})
                        aa_positions_in_region = [aa for aa in aa_positions_in_region if aa is not None]
                        total_aa_all_patients += len(aa_positions_in_region)
                        
                        for aa_pos in aa_positions_in_region:
                            nt_list = [nt for nt in region_nt_positions if get_first_aa_position(nt) == aa_pos]
                            any_selected = False
                            any_neutral = False
                            any_nonvar = False
                            
                            for nt in nt_list:
                                key = str(nt)
                                if key in patient_var_data:
                                    if patient_var_data[key] >= variability_threshold:
                                        pos_beta_data = df_all[(df_all['__patient__'] == patient) & (df_all['pos'] == nt)]
                                        if not pos_beta_data.empty and pos_beta_data['beta'].notna().any():
                                            any_selected = True
                                        else:
                                            any_neutral = True
                                    else:
                                        any_nonvar = True
                                else:
                                    any_nonvar = True
                            
                            if any_selected:
                                variable_aa_all_patients += 1
                                selected_aa_all_patients += 1
                            elif any_neutral:
                                neutral_aa_all_patients += 1
                                variable_aa_all_patients += 1
                            elif any_nonvar:
                                non_variable_aa_all_patients += 1
                
                # Crear hover detallado para Average
                txt = (
                    f"<b>Average</b> × <b>{cat}</b><br>"
                    f"<b>Selection ratio = {avg_val:.3f}</b><br>"
                    f"<i>Calculated from {len(patients)} patients</i><br><br>"
                    f"<b>Aggregated Nucleotides:</b><br>"
                    f"• Total: {total_nt_all_patients}<br>"
                    f"• Variable: {variable_nt_all_patients} ({variable_nt_all_patients/max(total_nt_all_patients,1)*100:.1f}%)<br>"
                    f"• Selected: {selected_nt_all_patients} ({selected_nt_all_patients/max(total_nt_all_patients,1)*100:.1f}%)<br>"
                    f"• Neutral: {neutral_nt_all_patients} ({neutral_nt_all_patients/max(total_nt_all_patients,1)*100:.1f}%)<br>"
                    f"• Non-variable: {non_variable_nt_all_patients} ({non_variable_nt_all_patients/max(total_nt_all_patients,1)*100:.1f}%)<br><br>"
                    f"<b>Aggregated Amino acids:</b><br>"
                    f"• Total: {total_aa_all_patients}<br>"
                    f"• Variable: {variable_aa_all_patients} ({variable_aa_all_patients/max(total_aa_all_patients,1)*100:.1f}%)<br>"
                    f"• Selected: {selected_aa_all_patients} ({selected_aa_all_patients/max(total_aa_all_patients,1)*100:.1f}%)<br>"
                    f"• Neutral: {neutral_aa_all_patients} ({neutral_aa_all_patients/max(total_aa_all_patients,1)*100:.1f}%)<br>"
                    f"• Non-variable: {non_variable_aa_all_patients} ({non_variable_aa_all_patients/max(total_aa_all_patients,1)*100:.1f}%)"
                )
            else:
                txt = f"<b>Average</b> × <b>{cat}</b><br>No data"
            
            var_avg_hover_row.append(txt)
        
        # 3. CONSTRUIR ARRAY FINAL DE HOVER EN ORDEN CORRECTO
        # El orden debe ser: [Average, 016-B, 019-A] para que coincida con z/y
        var_hov_with_avg = [var_avg_hover_row] + var_hov_with_avg
        
        # 4. REORDENAR TODO PARA QUE COINCIDA CON EL ORDEN VISUAL (autorange="reversed")
        # Con autorange="reversed", el orden visual es: [019-A, 016-B, Average] (de abajo hacia arriba)
        # Pero queremos que Average esté arriba, así que reordenamos TODO
        z_for_visual_order = list(reversed(z_var_with_avg))
        y_for_visual_order = list(reversed(y_var_with_avg))
        hover_for_visual_order = list(reversed(var_hov_with_avg))




        # Debug: verificar qué valores se están pasando al heatmap de variabilidad
        #print(f"[DEBUG VARIABILIDAD] var_mat shape: {var_mat.shape}")
        #print(f"[DEBUG VARIABILIDAD] var_mat min/max: {var_mat.min().min()}, {var_mat.max().max()}")
        #print(f"[DEBUG VARIABILIDAD] z_var_with_avg[0] (Average): {z_var_with_avg[0][:5]}...")
        #print(f"[DEBUG VARIABILIDAD] var_avg_row values: {var_avg_row[:5]}...")
        #print(f"[DEBUG VARIABILIDAD] z_var_with_avg[1] (First patient): {z_var_with_avg[1][:5]}...")
        
        # Debug: verificar si hay valores diferentes a 0.5 en var_mat
        # Solo acceder a columnas que existen
        #if 'ORF6' in var_mat.columns:
            #print(f"[DEBUG VARIABILIDAD] var_mat values for ORF6: {var_mat['ORF6'].tolist()}")
        #if 'N' in var_mat.columns:
            #print(f"[DEBUG VARIABILIDAD] var_mat values for N: {var_mat['N'].tolist()}")
        #if 'ORF8' in var_mat.columns:
            #print(f"[DEBUG VARIABILIDAD] var_mat values for ORF8: {var_mat['ORF8'].tolist()}")
        
        # Debug: verificar las primeras columnas disponibles
        #print(f"[DEBUG VARIABILIDAD] var_mat columns (first 5): {list(var_mat.columns[:5])}")
        #print(f"[DEBUG VARIABILIDAD] var_mat values for first column: {var_mat.iloc[:, 0].tolist()}")
        
        # Debug: verificar si hay algún valor diferente a 0.5 en toda la matriz
        #unique_values = var_mat.values.flatten()
        #unique_values = [v for v in unique_values if not pd.isna(v)]
        #print(f"[DEBUG VARIABILIDAD] Unique values in var_mat: {sorted(set(unique_values))}")
        #print(f"[DEBUG VARIABILIDAD] Count of values = 0.5: {sum(1 for v in unique_values if v == 0.5)}")
        #print(f"[DEBUG VARIABILIDAD] Count of values != 0.5: {sum(1 for v in unique_values if v != 0.5)}")
        
        # Debug: verificar si se está llamando build múltiples veces
        #print(f"[DEBUG BUILD CALL] === BUILD END for {col} ===")
        
        # Debug: verificar si se está sobrescribiendo var_mat en algún lugar
        #print(f"[DEBUG BEFORE RETURN] var_mat min/max: {var_mat.min().min()}, {var_mat.max().max()}")

        # ───────── Heatmap de variabilidad (hover alineado + HTML + Average con valor) ─────────
        # Alineación directa: ya construiste z/y/hover como:
        #   z_var_with_avg[0]   = Average
        #   y_var_with_avg[0]   = "Average"
        #   var_hov_with_avg[0] = hover de Average (con el valor medio)
        _cols = var_mat.columns.tolist()
        # Verificaciones de forma para prevenir desalineos
        assert len(z_var_with_avg) == len(y_var_with_avg) == len(var_hov_with_avg), \
            "Desalineo: filas de z/y/hover no coinciden"
        assert all(len(row) == len(_cols) for row in z_var_with_avg), \
            "Desalineo: columnas de z no coinciden con columnas del heatmap"
        assert all(len(row) == len(_cols) for row in var_hov_with_avg), \
            "Desalineo: columnas de hover no coinciden con columnas del heatmap"

        # 3) Construir figura con hover HTML completo y puntería perfecta
        # Usar la misma estrategia que funciona en el heatmap de beta slope:
        # NO reordenar nada, mantener el orden original
        z_for_visual_order = z_var_with_avg
        y_for_visual_order = y_var_with_avg
        hover_for_visual_order = var_hov_with_avg
        

        
        var_heat = go.Figure(go.Heatmap(
            z=z_for_visual_order,                         # Datos reordenados para orden visual
            x=_cols,
            y=y_for_visual_order,                         # Etiquetas reordenadas para orden visual
            colorscale=[[0.0, "rgb(34,139,34)"], [0.5, "white"], [1.0, "rgb(255,140,0)"]],
            zmin=0.0, zmax=1.0, zmid=0.5,
            hovertext=hover_for_visual_order,              # Hover reordenado para orden visual
            hovertemplate="%{hovertext}<extra></extra>",  # renderiza <b>, <br>, etc.
            hoverongaps=False,
            showscale=True,
            colorbar=dict(title="Selection Ratio (0=neg, 1=pos)", len=0.87)))
        
        var_heat.update_layout(
            height=max(440, len(patients)*22+210),
            xaxis=dict(side="top", tickangle=-45, tickfont_size=9,
                       tickmode="array", tickvals=var_mat.columns,
                       ticktext=[wrap(c) for c in var_mat.columns], constrain="domain"),
            yaxis=dict(autorange="reversed"), plot_bgcolor="rgba(255,255,255,0.8)",
            paper_bgcolor="rgba(255,255,255,0)",
            margin=dict(l=190, r=50, t=20, b=80),
            legend=dict(
                font=dict(size=10),
                x=1.02,
                y=1,
                xanchor="left",
                yanchor="top",
                orientation="v",
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                borderwidth=1
            ))
        




        # ─── stacked bars ───
        df_all["conf"] = df_all["aa_change"].str.replace(
            ".*ambiguous.*", "ambiguous", regex=True, case=False)
        bars = {}
        
        # Filtrar order para incluir solo categorías que tienen datos
        if col in df_all.columns:
            #print(f"[DEBUG] df_all['{col}'].unique(): {df_all[col].unique()}")
            #print(f"[DEBUG] order: {order}")
            available_categories = [cat for cat in order if cat in df_all[col].unique()]
            #print(f"[DEBUG] available_categories: {available_categories}")
            #print(f"[BARS] Processing {len(available_categories)} categories with data out of {len(order)} total")
        else:
            # Si la columna no existe, usar solo las categorías que sabemos que tienen datos
            available_categories = []
            #print(f"[BARS] Column '{col}' not found in df_all. Available columns: {list(df_all.columns)}")
            #print(f"[BARS] Skipping stacked bars for {col}")
        
        for cat in available_categories:
            sub = df_all[df_all[col] == cat]
            
            # ───────── NUEVA LÓGICA: FUSIONAR POSICIONES CON BETA POS/NEG ─────────
            print(f"[FUSION] Processing category {cat} with {len(sub)} records")
            
            # Agrupar por posición (gene/protein, aa_idx, nt) - SIN aa_change para combinar
            position_groups = defaultdict(list)
            for _, r in sub.iterrows():
                if pd.isna(r.pos) or pd.isna(r[col]) or pd.isna(r.aa_idx):
                    continue
                nt = int(r.pos)
                position_groups[(r[col], int(r.aa_idx), nt)].append(r)
            
            print(f"[FUSION] Found {len(position_groups)} unique positions")
            
            site2 = defaultdict(lambda: defaultdict(
                lambda: {"detail": {}, "p_pos": set(), "p_neg": set()}))
            
            for (gene_prot, aa_idx, nt), records in position_groups.items():
                print(f"[FUSION] Processing position {gene_prot} aa={aa_idx} nt={nt} with {len(records)} records")

                # Agrupar por paciente en esta posición
                by_patient = defaultdict(list)
                for r in records:
                    by_patient[r.__patient__].append(r)

                # Para cada paciente, combinar los cambios aminoacídicos
                for pat, recs in by_patient.items():
                    # Obtener todos los cambios aminoacídicos para este paciente en esta posición
                    aa_changes = []
                    for r in recs:
                        aa_change = str(r.get("aa_change", "")) if pd.notna(r.get("aa_change")) else "ambiguous"
                        if aa_change != "ambiguous":
                            aa_changes.append(aa_change)
                    
                    # Aplicar lógica de combinación
                    aa_label = combine_aa_changes(aa_changes) if aa_changes else f"{aa_idx}"

                    site = f"{aa_label} / {nt}"
                    st = site2[site]["by_variant"]  # agrupamos por variante en el mismo x

                    # Un único β y p por paciente si sólo hay uno; si hay varios, concatenamos (sin duplicar altura)
                    betas = [r.beta for r in recs if pd.notna(r.beta)]
                    pvals = [r.get("pval") for r in recs if pd.notna(r.get("pval"))]

                    beta_str = "<br/>".join(f"{b:.4f}" for b in betas) if betas else ""
                    pval_str = "<br/>".join(f"{p:.4f}" for p in pvals) if pvals else ""

                    plot_pdf   = f"{args.dir}/{pat}_plots/pos{nt:05d}_ALL.pdf"
                    genome_pdf = f"{args.dir}/{pat}_genome_selection.pdf"
                    
                    # Obtener datos de covariación específicos por variante
                    covarying_data_list = []
                    covarying_count_list = []
                    covarying_with_list = []
                    date_list = []
                    
                    for r in recs:
                        var_value = r['var'] if 'var' in r else (r.var if hasattr(r, 'var') else None)
                        if var_value is not None and pd.notna(var_value):
                            variant_key = f"{nt}_{var_value}"
                            
                            if variant_key in pos_to_covarying:
                                variant_data = pos_to_covarying[variant_key]
                                covarying_data_list.append(variant_data.get('covarying', False))
                                covarying_count_list.append(variant_data.get('covarying_count', 0))
                                covarying_with_list.append(variant_data.get('covarying_with', ''))
                                date_list.append(variant_data.get('date', ''))
                                #print(f"[COVARYING DEBUG] Position {nt}, variant {var_value}: Found specific covarying data: {variant_data}")
                            else:
                                # Fallback a datos por posición
                                if nt in pos_to_covarying:
                                    position_data = pos_to_covarying[nt]
                                    covarying_data_list.append(position_data.get('covarying', False))
                                    covarying_count_list.append(position_data.get('covarying_count', 0))
                                    covarying_with_list.append(position_data.get('covarying_with', ''))
                                    date_list.append(position_data.get('date', ''))
                                    #print(f"[COVARYING DEBUG] Position {nt}, variant {var_value}: Using fallback covarying data: {position_data}")
                                else:
                                    covarying_data_list.append(False)
                                    covarying_count_list.append(0)
                                    covarying_with_list.append('')
                                    date_list.append('')
                                    #print(f"[COVARYING DEBUG] Position {nt}, variant {var_value}: NO covarying data found")
                    
                    # Si no hay datos específicos por variante, usar datos por posición
                    if not covarying_data_list and nt in pos_to_covarying:
                        position_data = pos_to_covarying[nt]
                        covarying_data_list.append(position_data.get('covarying', False))
                        covarying_count_list.append(position_data.get('covarying_count', 0))
                        covarying_with_list.append(position_data.get('covarying_with', ''))
                        date_list.append(position_data.get('date', ''))
                        #print(f"[COVARYING DEBUG] Position {nt}: Using fallback covarying data: {position_data}")
                    elif not covarying_data_list:
                        covarying_data_list.append(False)
                        covarying_count_list.append(0)
                        covarying_with_list.append('')
                        date_list.append('')
                        #print(f"[COVARYING DEBUG] Position {nt}: NO covarying data found")
                    
                    # Crear strings separados por | para múltiples variantes
                    covarying_str = ' | '.join(str(x) for x in covarying_data_list) if len(covarying_data_list) > 1 else str(covarying_data_list[0]) if covarying_data_list else 'False'
                    # Asegurar que covarying_count_str no esté vacío si hay datos
                    if covarying_count_list:
                        covarying_count_str = ' | '.join(str(x) for x in covarying_count_list) if len(covarying_count_list) > 1 else str(covarying_count_list[0])
                    else:
                        covarying_count_str = '0'
                    # Para covarying_with, mantener las covariaciones separadas por variante usando un separador especial
                    covarying_with_str = ' ||| '.join(str(x) for x in covarying_with_list) if len(covarying_with_list) > 1 else str(covarying_with_list[0]) if covarying_with_list else ''
                    # Asegurar que date_str no esté vacío si hay fechas
                    date_str = '<br>'.join(str(x).strip() for x in date_list if x and str(x).strip() != '') if date_list else ''
                    if not date_str and date_list:
                        date_str = str(date_list[0]).strip() if date_list[0] and str(date_list[0]).strip() != '' else ''
                    
                    st["detail"][pat] = {
                        "beta": beta_str or None,
                        "pval": pval_str or None,
                        "plot": plot_pdf,
                        "genome": genome_pdf,
                        "covarying": covarying_str,
                        "covarying_count": covarying_count_str,
                        "covarying_with": covarying_with_str,
                        "date": date_str,
                        "aa_change": aa_label
                    }
                    
                    # Llenar conjuntos p_pos y p_neg basándose en los valores de beta
                    for beta in betas:
                        if beta > 0:
                            st["p_pos"].add(pat)
                        elif beta < 0:
                            st["p_neg"].add(pat)
                    
                #print(f"[FUSION] Done position {gene_prot} aa={aa_idx} nt={nt}")

            #print(f"[DEBUG] site2 keys: {list(site2.keys())}")
            #print(f"[DEBUG] site2 values sample: {list(site2.values())[:2] if site2 else 'empty'}")
            keep = [s for s in site2 if len({p for d in site2[s].values() for p in d["detail"]}) >= args.min]
            # Ordenar por número de pacientes (descendente) y luego por posición nucleotídica (ascendente)
            def sort_key(s):
                num_patients = len({p for d in site2[s].values() for p in d["detail"]})
                # Extraer posición nucleotídica para ordenar como segundo criterio
                nt_pos = 0
                if ' / ' in s:
                    try:
                        nt_pos = int(s.split(' / ')[1])
                    except (ValueError, IndexError):
                        nt_pos = 0
                # Devolver tupla: (-num_patients, nt_pos) para orden descendente por pacientes, ascendente por posición
                return (-num_patients, nt_pos)
            
            # Debug: mostrar información de ordenamiento
            keep.sort(key=sort_key)
            
            # DEBUG: Verificar ordenamiento después del sort
            #print(f"[DEBUG FINAL ORDER] keep after sort: {keep}")
            #for i, s in enumerate(keep):
            #    num_patients = len({p for d in site2[s].values() for p in d["detail"]})
            #    nt_pos = 0
            #    if ' / ' in s:
            #        try:
            #            nt_pos = int(s.split(' / ')[1])
            #        except (ValueError, IndexError):
            #            nt_pos = 0
            #        print(f"  [{i}] {s} -> patients: {num_patients}, nt_pos: {nt_pos}")
            
            # Crear etiquetas para el eje X (posición nucleotídica y aminoacídica sin letras)
            x_labels = []
            for s in keep:
                if ' / ' in s:
                    # Extraer posición nucleotídica y aminoacídica (ej: "L9Q / 28299" -> "28299 (9)")
                    parts = s.split(' / ')
                    if len(parts) == 2:
                        nt_pos = parts[1]  # "28299"
                        aa_part = parts[0]  # "L9Q"
                        # Extraer solo el número de la posición aminoacídica
                        import re
                        aa_match = re.search(r'(\d+)', aa_part)
                        if aa_match:
                            aa_pos = aa_match.group(1)  # "9"
                            x_labels.append(f"{nt_pos} ({aa_pos})")
                        else:
                            x_labels.append(nt_pos)
                    else:
                        x_labels.append(s)
                else:
                    x_labels.append(s)
            if not keep:
                print(f"[DEBUG] No keep for {cat}, creating empty chart")
                # Mostrar un gráfico vacío si no hay sitios suficientes
                traces = []
                x_labels = []
                bars[cat] = dict(
                    data=traces,
                    layout=dict(
                        barmode="stack",
                        showlegend=False,
                        plot_bgcolor="rgba(255,255,255,0.8)",
                        paper_bgcolor="rgba(255,255,255,0)",
                        xaxis=dict(type="category", categoryorder="array",
                                   categoryarray=x_labels,
                                   tickangle=-45, tickfont_size=9),
                        yaxis=dict(title="Number of patients", tickmode="linear", tick0=0, dtick=1),
                        xaxis_title="aa / nt",
                        margin=dict(l=70, r=20, t=60, b=80),
                        height=300
                    )
                )
                continue
            confs = sorted({c for d in site2.values() for c in d})
            #colour = {c: palette[i % len(palette)] for i, c in enumerate(confs)} #No se usa
            traces = []
            for conf in confs:
                y, hov2, cdat, x_positions = [], [], [], []
                for i, site in enumerate(keep):
                    det = site2[site].get(conf, {})
                    if det and det.get("detail") and len(det["detail"]) > 0:
                        n_pos = len(det["p_pos"]); n_neg = len(det["p_neg"])
                        pos_list = " · ".join(sorted(det["p_pos"])) or "–"
                        neg_list = " · ".join(sorted(det["p_neg"])) or "–"
                        table = "<br>".join(
                            f"{p} | β: {d['beta'].replace('<br/>', ' / ')} | p: {d['pval'].replace('<br/>', ' / ')}"
                            for p, d in sorted(det["detail"].items()))
                        hov2.append(
                            f"<b>{site}</b><br>{conf}" +
                            f"<br><span style='color:red'>β⁺ ({n_pos})</span>: {pos_list}" +
                            f"<br><span style='color:blue'>β⁻ ({n_neg})</span>: {neg_list}" +
                            f"<br>--<br>{table}")
                        cdat.append(json.dumps(det["detail"], cls=PlotlyJSONEncoder))
                        y.append(len(det["detail"])) #Number of patients
                        x_positions.append(x_labels[i])  # Etiquetas simplificadas para eje X
                if max(y) == 0:
                    continue
                colors = []
                for i, site in enumerate(keep):
                    det = site2[site].get(conf, {})
                    if det and det.get("detail") and len(det["detail"]) > 0:
                        n_pos = len(det["p_pos"])
                        n_neg = len(det["p_neg"])
                        colors.append(color_for_site(site, n_pos, n_neg))
                    else:
                        colors.append(color_for_site(site, 0, 0))
                traces.append(go.Bar(
                    x=x_positions, y=y, name=conf,
                    marker=dict(color=colors, line=dict(color='black', width=1)),
                    hovertext=hov2, customdata=cdat,
                    hovertemplate="%{hovertext}<extra></extra>").to_plotly_json())

            # Asegurar que todas las y sean numéricas (evita object)
            for tr in traces:
                if "y" in tr and tr["y"] is not None:
                    tr["y"] = [float(y) if y is not None else 0.0 for y in tr["y"]]
            # DEBUG: Verificar categoryarray antes de crear el layout
            #print(f"[DEBUG CATEGORYARRAY] x_labels for {cat}: {x_labels}")
            #print(f"[DEBUG CATEGORYARRAY] x_labels length: {len(x_labels)}")
            
            bars[cat] = dict(
                data=traces,
                layout=dict(barmode="stack", showlegend=False, 
                            plot_bgcolor="rgba(255,255,255,0.8)",
                            paper_bgcolor="rgba(255,255,255,0)",
                            xaxis=dict(type="category", categoryorder="array",
                                       categoryarray=x_labels,
                                       tickangle=-45, tickfont_size=9),
                            yaxis=dict(title="Number of patients", tickmode="linear", tick0=0, dtick=1), xaxis_title="aa / nt",
                            margin=dict(l=70, r=20, t=60, b=80),
                            height=max(480, 290 + len(keep)*6)),
                patients=patients,  # Incluir lista de pacientes para filtrado
                keep=keep  # Incluir sitios para referencia
            )
        btns = "".join(f'<button class="cdsbtn" data-cds="{c}">{c}</button>' for c in order)
        print(f"[DEBUG] === BUILD END for {col} ===")
        print(f"[DEBUG] bars keys: {list(bars.keys())}")
        print(f"[DEBUG] order: {order}")
        # Devolvemos también el orden para poblar selects de forma consistente
        return {
            "heat": heat.to_plotly_json(),
            "var_heat": var_heat.to_plotly_json(),
            "bars": bars,
            "btns": btns,
            "var_detailed": var_detailed_data,
            "order": order,
            # Nuevas visualizaciones para heatmaps (usando la misma estructura que los originales)
            "heat_selected_ratio": create_ratio_heatmap(mat_selected_ratio, "Selected/Total", 
                                                       [[0.0, "white"], [1.0, "rgb(255,140,0)"]], 
                                                       zmin=0.0, detailed_data=detailed_data),
            "heat_non_variable_ratio": create_ratio_heatmap(mat_non_variable_ratio, "Non-Variable/Total", 
                                                           [[0.0, "white"], [1.0, "rgb(34,139,34)"]], 
                                                           zmin=0.0, detailed_data=detailed_data),
            "heat_positive_ratio": create_ratio_heatmap(mat_positive_ratio, "Positive/Total", 
                                                       [[0.0, "white"], [1.0, "rgb(220,20,60)"]], 
                                                       zmin=0.0, detailed_data=detailed_data),
            "heat_negative_ratio": create_ratio_heatmap(mat_negative_ratio, "Negative/Total", 
                                                       [[0.0, "white"], [1.0, "rgb(30,144,255)"]], 
                                                       zmin=0.0, detailed_data=detailed_data)
        }
    # ───────── Crear mapeo de posiciones a genes y proteínas ─────────
    def create_position_mapping():
        """Crea un mapeo de posiciones nucleotídicas a genes, proteínas y posiciones aminoacídicas leyendo desde archivos TSV."""
        pos_to_gene = {}
        pos_to_protein = {}
        pos_to_aa = {} # Almacenará una lista de (gene, aa_pos, frame) para posiciones compartidas
        pos_to_covarying = {} # Almacenará datos de covariación por posición

        # Buscar archivos *_genes.tsv en el directorio de datos
        genes_files = glob.glob(os.path.join(args.dir, "*_genes.tsv"))
        
        # Buscar archivos *_selected_ranked.tsv para datos de covariación
        tsv_files = glob.glob(os.path.join(args.dir, "*_selected_ranked.tsv"))
        
        if not tsv_files:
            print(f"Error: No se encontraron archivos *_selected_ranked.tsv en '{args.dir}'")
            return {}, {}, {}, {}
        
        if not genes_files:
            print(f"Error: No se encontraron archivos *_genes.tsv en '{args.dir}'")
            return {}, {}, {}, {}
        
        # ───────── PROCESAR ARCHIVOS *_genes.tsv PARA MAPEO NT->AA ─────────
        # Usar el primer archivo encontrado (todos deberían tener la misma estructura)
        tsv_file = genes_files[0]
        print(f"[DEBUG] Cargando correspondencias desde: {tsv_file}")
        
        try:
            with open(tsv_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line_num == 1:  # Saltar header
                        continue
                        
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                        
                    position = int(parts[0])
                    gene = parts[1]
                    codon_pos = parts[2]
                    aa_idx = parts[3]
                    
                    # Saltar posiciones sin gen
                    if gene == '' or gene == 'NA':
                        continue
                    
                    # Inicializar posición
                    if position not in pos_to_aa:
                        pos_to_aa[position] = []
                    
                    # Debug específico para la posición 13468
                    #if position == 13468:
                    #    print(f"[DEBUG] Procesando pos 13468: gene='{gene}', aa_idx='{aa_idx}'")
                    
                    # Procesar genes múltiples (separados por coma)
                    genes = [g.strip() for g in gene.split(',')]
                    aa_indices = [int(a.strip()) for a in aa_idx.split(',')]
                    
                    # Debug específico para la posición 13468
                    #if position == 13468:
                    #    print(f"[DEBUG] pos 13468: genes={genes}, aa_indices={aa_indices}")
                    
                    # Mapear cada gen/posición aminoacídica
                    # Caso especial: si hay un solo gen pero múltiples posiciones AA (frameshift)
                    if len(genes) == 1 and len(aa_indices) > 1:
                        g = genes[0]
                        if g and g != 'NA':
                            for aa_pos in aa_indices:
                                pos_to_aa[position].append((g, aa_pos, "single"))
                                # Debug específico para la posición 13468
                                #if position == 13468:
                                #    print(f"[DEBUG] pos 13468: agregado ({g}, {aa_pos}, 'single')")
                            pos_to_gene[position] = g
                            # Obtener proteína específica basada en la posición nucleotídica
                            protein_name = get_specific_protein_for_position(position)
                            pos_to_protein[position] = protein_name
                    else:
                        # Caso normal: mapear genes con sus posiciones AA correspondientes
                        for i, (g, aa_pos) in enumerate(zip(genes, aa_indices)):
                            if g and g != 'NA':
                                pos_to_aa[position].append((g, aa_pos, "single"))
                                pos_to_gene[position] = g
                                
                                # Obtener proteína específica basada en la posición nucleotídica
                                protein_name = get_specific_protein_for_position(position)
                                pos_to_protein[position] = protein_name
                                
                                # Debug específico para la posición 13468
                                #if position == 13468:
                                #    print(f"[DEBUG] pos 13468: agregado ({g}, {aa_pos}, 'single')")
                            
        except Exception as e:
            print(f"Error leyendo archivo TSV: {e}")
            return {}, {}, {}, {}

        # Debugging para la posición 13468
        #if 13468 in pos_to_aa:
        #    print(f"[DEBUG] Pos 13468: {pos_to_aa[13468]}")
        #else:
        #    print("[DEBUG] Pos 13468 no encontrada en pos_to_aa")

        #print(f"[DEBUG] Mapeo nt→aa creado: {len(pos_to_aa)} posiciones mapeadas")
        
        # ───────── PROCESAR ARCHIVOS *_selected_ranked.tsv PARA DATOS DE COVARIACIÓN ─────────
        # Funciones auxiliares para procesamiento de datos (misma lógica que versión 1)
        def _to_bool(x):
            if x is None: return False
            xs = str(x).strip().lower()
            return xs in ("true","1","yes","y","✓","si","sí")
        def _norm_int(x, default=0):
            try:
                if pd.isna(x): return default
            except Exception:
                pass
            try: return int(x)
            except Exception:
                try: return int(float(x))
                except Exception: return default
        def _norm_str(x):
            return '' if x is None or (isinstance(x,float) and pd.isna(x)) else str(x)
        
        try:
            for tsv_file in tsv_files:
                print(f"[DEBUG] Cargando datos de covariación desde: {tsv_file}")
                lines_processed = 0
                with open(tsv_file, 'r') as f:
                    for line_num, line in enumerate(f, 1):
                        if line_num == 1:  # Saltar header
                            continue
                            
                        parts = line.strip().split('\t')
                        if len(parts) < 20:  # Necesitamos al menos las columnas hasta date
                            print(f"[COVARYING SKIP] Skipping line {line_num} due to insufficient columns: {len(parts)}")
                            continue
                        
                        lines_processed += 1
                            
                        gene = parts[0]
                        aa_idx = parts[1]
                        try:
                            position = int(parts[2])
                        except ValueError:
                            print(f"[COVARYING SKIP] Skipping line {line_num} - invalid position: '{parts[2]}'")
                            continue
                        var = parts[3]
                        beta = parts[4]
                        pval = parts[5]
                        score = parts[6]
                        tag = parts[7]
                        delta_f = parts[8]
                        codon = parts[9]
                        aa_change = parts[10]
                        ambiguous = parts[11]
                        annot_detail = parts[12]
                        depth = parts[13]
                        n_days = parts[14]
                        variant_suffix = parts[15]
                        date = parts[16]
                        covarying = parts[17]
                        covarying_count = parts[18]  # covarying_co
                        covarying_with = parts[19]   # covarying_wi
                        days_from_tc = parts[20]
                        
                        # índice de covariación por variante y por posición (misma lógica que versión 1)
                        cov = _to_bool(covarying) if covarying else False
                        cnt = _norm_int(covarying_count, 0)
                        with_ = covarying_with if covarying_with else ""
                        # Asegurar que date_val no esté vacío si hay una fecha válida
                        date_val = date.strip() if date and str(date).strip() else ""
                        
                        # Crear entradas por variante
                        variant_key = f"{position}_{var}"
                        pos_to_covarying[variant_key] = {
                            "covarying": cov, 
                            "covarying_count": cnt, 
                            "covarying_with": with_, 
                            "date": date_val,
                            "position": position,
                            "variant": var
                        }
                        
                        # También crear entrada por posición (para fallback cuando no se encuentra por variante)
                        # Solo si no existe ya una entrada para esta posición, o si esta tiene más información
                        if position not in pos_to_covarying:
                            pos_to_covarying[position] = {
                                "covarying": cov, 
                                "covarying_count": cnt, 
                                "covarying_with": with_, 
                                "date": date_val,
                                "position": position,
                                "variant": var
                            }
                        else:
                            # Si ya existe, actualizar solo si esta variante tiene más información (mayor count o fecha)
                            existing = pos_to_covarying[position]
                            if cnt > existing.get("covarying_count", 0) or (date_val and not existing.get("date")):
                                pos_to_covarying[position] = {
                                    "covarying": cov, 
                                    "covarying_count": cnt, 
                                    "covarying_with": with_, 
                                    "date": date_val,
                                    "position": position,
                                    "variant": var
                                }
                        
                        # Debug para entender qué se está cargando
                        #print(f"[COVARYING LOAD] Loading position {position}, variant {var}: covarying={cov}, count={cnt}")
        except Exception as e:
            print(f"Error leyendo archivo TSV: {e}")
            return {}, {}, {}, {}
        
        #print(f"[COVARYING DEBUG] Final pos_to_covarying keys: {list(pos_to_covarying.keys())}")
        #print(f"[DEBUG] Datos de covariación cargados: {len(pos_to_covarying)} posiciones")
        return pos_to_gene, pos_to_protein, pos_to_aa, pos_to_covarying

    def format_covarying_with(covarying_with_str):
        """Formatea el string de covariación para mejor legibilidad."""
        if not covarying_with_str or covarying_with_str.strip() == '':
            return 'No covariation data'
        
        # Dividir por " | " para separar cada covariación
        covarying_items = covarying_with_str.split(' | ')
        formatted_items = []
        
        for item in covarying_items:
            item = item.strip()
            if not item:
                continue
                
            # Parsear el formato: nt=25511(C); aa=ORF3a:40S [Δβ=0.xx, Δday=X]
            if 'nt=' in item and 'aa=' in item:
                # Extraer posición nucleotídica
                nt_match = item.split('nt=')[1].split(';')[0] if 'nt=' in item else ''
                # Extraer información aminoacídica
                aa_part = item.split('aa=')[1].split(' [')[0] if 'aa=' in item else ''
                # Extraer diferencia de beta y día
                delta_beta_match = ''
                delta_day_match = ''
                if '[Δβ=' in item:
                    bracket_content = item.split('[')[1].split(']')[0] if '[' in item and ']' in item else ''
                    # Puede ser "Δβ=0.xx, Δday=X" o solo "Δβ=0.xx"
                    if ', Δday=' in bracket_content:
                        parts = bracket_content.split(', ')
                        delta_beta_match = parts[0].replace('Δβ=', '') if len(parts) > 0 else ''
                        delta_day_match = parts[1].replace('Δday=', '') if len(parts) > 1 else ''
                    else:
                        delta_beta_match = bracket_content.replace('Δβ=', '')
                
                if nt_match and aa_part:
                    # Formatear como: "Nucleotide: 25511(C) | Amino acid: ORF3a:40S | Δβ: 0.0001 | Δday: 3"
                    formatted_item = f"<strong>Nucleotide:</strong> {nt_match}<br/>"
                    formatted_item += f"<strong>Amino acid:</strong> {aa_part}<br/>"
                    if delta_beta_match:
                        formatted_item += f"<strong>Δβ:</strong> {delta_beta_match}"
                    if delta_day_match:
                        formatted_item += f" | <strong>Δday:</strong> {delta_day_match}"
                    formatted_items.append(formatted_item)
            else:
                # Si no sigue el formato esperado, mostrar tal como está
                formatted_items.append(item)
        
        return '<br/><br/>'.join(formatted_items) if formatted_items else 'No covariation data'

    def get_protein_name_from_gene(gene_name):
        """Obtiene el nombre de la proteína a partir del nombre del gen."""
        gene_to_protein = {
            'ORF1ab': 'pp1ab',
            'ORF1a': 'pp1a', 
            'ORF1b': 'pp1b',
            'S': 'S',
            'ORF3a': 'ORF3a',
            'E': 'E',
            'M': 'M',
            'ORF6': 'ORF6',
            'ORF7a': 'ORF7a',
            'ORF7b': 'ORF7b',
            'ORF8': 'ORF8',
            'N': 'N',
            'ORF10': 'ORF10'
        }
        return gene_to_protein.get(gene_name, gene_name)
    
    def get_specific_protein_for_position(nt_position):
        """Obtiene la proteína específica para una posición nucleotídica usando CDS_LIST."""
        for start, end, protein in CDS_LIST:
            if start <= nt_position <= end:
                return protein
        return 'Unknown'
    
    # Funciones auxiliares para manejar posiciones compartidas
    def get_aa_positions(nt_pos):
        """Obtiene todas las posiciones aminoacídicas para una posición nucleotídica"""
        if nt_pos in pos_to_aa:
            return [aa_pos for gene, aa_pos, frame in pos_to_aa[nt_pos]]
        return []
    
    def get_aa_positions_with_genes(nt_pos):
        """Obtiene posiciones aminoacídicas con información de genes"""
        if nt_pos in pos_to_aa:
            return [(gene, aa_pos) for gene, aa_pos, frame in pos_to_aa[nt_pos]]
        return []
    
    def get_unique_aa_positions(nt_positions):
        """Obtiene posiciones aminoacídicas únicas para una lista de posiciones nucleotídicas"""
        unique_aa = set()
        for nt_pos in nt_positions:
            unique_aa.update(get_aa_positions(nt_pos))
        return sorted(unique_aa)
    
    def format_aa_position_for_hover(nt_pos):
        """Formatea la posición aminoacídica para el hover, mostrando gen solo si hay múltiples frames/genes."""
        if nt_pos not in pos_to_aa:
            return ""
        
        aa_data_list = pos_to_aa[nt_pos]
        
        # Si solo hay una entrada, mostrar solo la posición AA
        if len(aa_data_list) == 1:
            return str(aa_data_list[0][1])
        
        # Si hay múltiples, formatear como "pos-gene; pos-gene"
        formatted_parts = []
        for gene, aa_pos, frame in aa_data_list:
            formatted_parts.append(f"{aa_pos}-{gene}")
        
        return "; ".join(sorted(list(set(formatted_parts))))
    
    def get_first_aa_position(nt_pos):
        """Obtiene la primera posición aminoacídica para una posición nucleotídica (para compatibilidad)"""
        if nt_pos in pos_to_aa and pos_to_aa[nt_pos]:
            return pos_to_aa[nt_pos][0][1]  # (gene, aa_pos, frame) -> aa_pos
        return None
    
    # Crear mapeo de posiciones antes de build()
    pos_to_gene, pos_to_protein, pos_to_aa, pos_to_covarying = create_position_mapping()
    
    # Crear una copia de pos_to_covarying para preservar los datos específicos por variante
    pos_to_covarying_backup = {k: v.copy() if isinstance(v, dict) else v for k, v in pos_to_covarying.items()}
    
    VIEWS = {"gene": build("gene"), "protein": build("product")}
    print("[GENES]", VIEWS["gene"]["order"])
    print("[PROTEINS]", VIEWS["protein"]["order"])
    #print(f"[DEBUG] VIEWS['gene']['order'] length: {len(VIEWS['gene']['order'])}")
    #print(f"[DEBUG] VIEWS['protein']['order'] length: {len(VIEWS['protein']['order'])}")
    
    # GENERAR OPCIONES DE GEN Y NSP PARA LOS SELECTORES (coherentes con los ejes)
    gene_options = "".join(f'<option value="{g}" selected>{g}</option>'
                           for g in VIEWS["gene"]["order"])
    nsp_options  = "".join(f'<option value="{p}" selected>{p}</option>'
                           for p in VIEWS["protein"]["order"])
    
    #print(f"[DEBUG] gene_options length: {len(gene_options)}")
    #print(f"[DEBUG] nsp_options length: {len(nsp_options)}")
    #print(f"[DEBUG] gene_options sample: {gene_options[:200]}...")
    #print(f"[DEBUG] nsp_options sample: {nsp_options[:200]}...")
    
    # ───────── Función para leer posiciones seleccionadas ─────────
    def read_selected_positions():
        """Lee automáticamente los archivos *_selected_ranked.tsv para obtener posiciones seleccionadas"""
        selected_positions = {}
        
        # Buscar todos los archivos *_selected_ranked.tsv
        selected_files = glob.glob(os.path.join(args.dir, "*_selected_ranked.tsv"))
        
        for file_path in selected_files:
            try:
                # Extraer nombre del paciente del nombre del archivo
                patient_name = os.path.basename(file_path).replace("_selected_ranked.tsv", "")
                
                # Leer el archivo
                df_selected = pd.read_csv(file_path, sep='\t')
                
                # Procesar cada posición seleccionada
                for _, row in df_selected.iterrows():
                    if pd.notna(row['pos']):
                        pos = int(row['pos'])
                        
                        # Inicializar la posición si no existe
                        if pos not in selected_positions:
                            selected_positions[pos] = {'count': 0, 'patients': set()}
                        
                        # Agregar el paciente a esta posición
                        selected_positions[pos]['patients'].add(patient_name)
                
            except Exception as e:
                print(f"[ERROR] Error reading {file_path}: {e}")
        
        # Convertir sets a listas para serialización JSON
        for pos in selected_positions:
            selected_positions[pos]['count'] = len(selected_positions[pos]['patients'])
            selected_positions[pos]['patients'] = list(selected_positions[pos]['patients'])
        
        return selected_positions
    
    def read_selected_positions_with_flags():
        """Lee automáticamente los archivos *_selected_ranked.tsv para obtener posiciones seleccionadas con flags de existencia de PDFs"""
        selected_positions = {}
        
        # Buscar todos los archivos *_selected_ranked.tsv
        selected_files = glob.glob(os.path.join(args.dir, "*_selected_ranked.tsv"))
        
        for file_path in selected_files:
            try:
                # Extraer nombre del paciente del nombre del archivo
                patient_name = os.path.basename(file_path).replace("_selected_ranked.tsv", "")
                
                # Leer el archivo
                df_selected = pd.read_csv(file_path, sep='\t')
                
                # Inicializar el paciente si no existe
                if patient_name not in selected_positions:
                    selected_positions[patient_name] = {}
                
                # Procesar cada posición seleccionada
                for _, row in df_selected.iterrows():
                    if pd.notna(row['pos']):
                        pos = int(row['pos'])
                        
                        # Verificar existencia de PDFs
                        plot_path = os.path.join(args.dir, f"{patient_name}_plots", f"pos{pos:05d}_ALL.pdf")
                        genome_path = os.path.join(args.dir, f"{patient_name}_genome_selection.pdf")
                        
                        plot_exists = os.path.exists(plot_path)
                        genome_exists = os.path.exists(genome_path)
                        
                        # Buscar datos de covariación para esta variante específica
                        covarying_data = None
                        
                        # Primero intentar buscar por variante específica (posición + variante)
                        variant_value = row.get('var', '')
                        if variant_value:
                            variant_key = f"{pos}_{variant_value}"
                            if variant_key in pos_to_covarying_backup:
                                covarying_data = pos_to_covarying_backup[variant_key]
                        
                        # Si no se encontró por variante específica, buscar por posición general
                        if covarying_data is None and pos in pos_to_covarying_backup:
                            covarying_data = pos_to_covarying_backup[pos]
                        
                        # Crear entrada de datos para esta posición
                        position_entry = {
                            'gene': row.get('gene', ''),
                            'product': row.get('product', ''),
                            'aa_idx': row.get('aa_idx', ''),
                            'var': row.get('var', ''),
                            'aa_change': row.get('aa_change', ''),
                            'beta': row.get('beta', None),
                            'pval': row.get('pval', None),
                            'plot': plot_path,
                            'genome': genome_path,
                            'plot_exists': plot_exists,
                            'genome_exists': genome_exists
                        }
                        
                        # Agregar datos de covariación si están disponibles
                        if covarying_data:
                            position_entry['covarying'] = covarying_data.get('covarying', False)
                            position_entry['covarying_count'] = covarying_data.get('covarying_count', 0)
                            position_entry['covarying_with'] = covarying_data.get('covarying_with', '')
                            position_entry['date'] = covarying_data.get('date', '')
                        else:
                            position_entry['covarying'] = False
                            position_entry['covarying_count'] = 0
                            position_entry['covarying_with'] = ''
                            position_entry['date'] = ''
                        
                        
                        # Si ya existe una entrada para esta posición, convertir a lista
                        if pos in selected_positions[patient_name]:
                            if not isinstance(selected_positions[patient_name][pos], list):
                                # Convertir la entrada existente a lista
                                selected_positions[patient_name][pos] = [selected_positions[patient_name][pos]]
                            # Agregar la nueva entrada
                            selected_positions[patient_name][pos].append(position_entry)
                        else:
                            # Primera entrada para esta posición
                            selected_positions[patient_name][pos] = position_entry
                
            except Exception as e:
                print(f"[ERROR] Error reading {file_path}: {e}")
        
        # Aplicar combinación de aminoácidos para posiciones con múltiples entradas
        for patient_name in selected_positions:
            for pos in selected_positions[patient_name]:
                if isinstance(selected_positions[patient_name][pos], list) and len(selected_positions[patient_name][pos]) > 1:
                    # Obtener todos los cambios aminoacídicos
                    aa_changes = [entry['aa_change'] for entry in selected_positions[patient_name][pos] if entry['aa_change']]
                    
                    if aa_changes:
                        # Aplicar la función de combinación
                        combined_aa = combine_aa_changes(aa_changes)
                        
                        # Actualizar el primer entry con el aminoácido combinado
                        selected_positions[patient_name][pos][0]['aa_change'] = combined_aa
                        
                        # NO combinar los datos de covariación - mantener las entradas separadas
                        # El frontend debe mostrar cada variante con sus propios datos de covariación
                        # Las entradas se mantienen como una lista para que cada una conserve sus datos únicos
        
        return selected_positions

    # Leer posiciones seleccionadas ANTES de crear VIEWS para preservar pos_to_covarying
    selected_positions_data = read_selected_positions()
    selected_positions_with_flags = read_selected_positions_with_flags()

    def _public_view_payload(view_payload):
        if not isinstance(view_payload, dict):
            return view_payload
        payload = dict(view_payload)
        payload.pop("fig1_counts", None)
        payload.pop("fig5_counts", None)
        payload.pop("btns", None)
        return payload

    def blob(i, o):
        return f'<script id="{i}" type="application/json">{_json_dumps_compact(o, cls=PlotlyJSONEncoder)}</script>'

    BLOBS = "__VIEW_BLOBS__"
    var_blob = (
        f'<script id="var-data" type="application/json">'
        f'{_json_dumps_compact(var_data, cls=PlotlyJSONEncoder)}'
        f'</script>'
    )
    selected_positions_blob = (
        f'<script id="selected-positions" type="application/json">'
        f'{_json_dumps_compact(selected_positions_data, cls=PlotlyJSONEncoder)}'
        f'</script>'
    )
    selected_positions_with_flags_blob = (
        f'<script id="selected-positions-flags" type="application/json">'
        f'{_json_dumps_compact(selected_positions_with_flags, cls=PlotlyJSONEncoder)}'
        f'</script>'
    )
    # var_detailed ya va embebido dentro de cada vista (gene/protein)
    # Crear JSON para JavaScript con formateo correcto para hover
    aa_positions_formatted = {}
    for nt_pos, aa_data in pos_to_aa.items():
        aa_positions_formatted[nt_pos] = format_aa_position_for_hover(nt_pos)
    
    position_mapping = {
        'genes': pos_to_gene,
        'proteins': pos_to_protein,
        'aa_positions': aa_positions_formatted
    }
    
    mapping_blob = (
        f'<script id="position-mapping" type="application/json">'
        f'{_json_dumps_compact(position_mapping, cls=PlotlyJSONEncoder)}'
        f'</script>'
    )
    
    # Blob para datos de filtros
    filter_data_blob = (
        f'<script id="filter-data" type="application/json">'
        f'{_json_dumps_compact({"patient_data": patient_filter_data, "patient_rows": patient_filter_rows, "columns": filter_columns, "unique_values": filter_unique_values, "patient_name_mapping": patient_name_mapping}, cls=PlotlyJSONEncoder)}'
        f'</script>'
    )

    # ───────── UPGMA (per-gene and per-protein, Jaccard binario sin pesos) ─────────
    def _build_pos_to_regions(segments):
        pos2names = {}
        for s, e, name in segments:
            for nt in range(int(s), int(e)+1):
                pos2names.setdefault(nt, set()).add(name)
        return pos2names

    def _build_patient_tokens_by_region(selected_positions_with_flags, pos2names):
        # dict: region_name -> dict: patient -> set(tokens)
        data = { "__ALL__": {} }
        # inicializa por región
        for nameset in set().union(*pos2names.values()) if pos2names else set():
            data[nameset] = {}
        # recorrer pacientes
        for patient, posmap in selected_positions_with_flags.items():
            # asegúrate de inicializar en "__ALL__"
            data["__ALL__"].setdefault(patient, set())
            for pos, entry in posmap.items():
                entries = entry if isinstance(entry, list) else [entry]
                for ent in entries:
                    beta = ent.get('beta')
                    if beta is None or not isinstance(beta, (int,float)) or beta == 0:
                        continue
                    sign = 'pos' if beta > 0 else 'neg'
                    varv = ent.get('var') or ent.get('aa_change') or ''
                    token = f"{int(pos)}|{varv}|{sign}"
                    # ALL
                    data["__ALL__"].setdefault(patient, set()).add(token)
                    # por regiones (si la pos pertenece)
                    for name in pos2names.get(int(pos), []):
                        data.setdefault(name, {}).setdefault(patient, set()).add(token)
        return data

    def _upgma_from_token_sets(patients, token_sets_map):
        import numpy as _np
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, dendrogram
        
        # Filtrar pacientes que NO tienen tokens (sin posiciones seleccionadas)
        patients_with_tokens = [p for p in patients if token_sets_map.get(p, set())]
        
        # Si no hay suficientes pacientes con tokens, retornar None
        if len(patients_with_tokens) < 2:
            return None, None
        
        # construir vocabulario solo de pacientes con tokens
        vocab = sorted({t for p in patients_with_tokens for t in token_sets_map.get(p, set())})
        if len(vocab) == 0:
            return None, None
        
        idx = {t:i for i,t in enumerate(vocab)}
        X = _np.zeros((len(patients_with_tokens), len(vocab)), dtype=bool)
        for r, p in enumerate(patients_with_tokens):
            for t in token_sets_map.get(p, set()):
                j = idx.get(t)
                if j is not None:
                    X[r, j] = True
        if X.shape[1] == 0:
            return None, None
        d = pdist(X, metric="jaccard")
        if d.size == 0:
            return None, None
        Z = linkage(d, method="average")
        leaves = dendrogram(Z, labels=list(patients_with_tokens), no_plot=True)["leaves"]
        # redondear alturas para aligerar
        Z_list = Z.tolist()
        for row in Z_list:
            row[2] = float(f"{row[2]:.6f}")
        return Z_list, leaves

    # Construir mapeos pos->gen y pos->proteína
    POS2GENE = _build_pos_to_regions(GENE_LIST)
    POS2PROT = _build_pos_to_regions(CDS_LIST)

    # Construir tokens por región
    tokens_by_gene = _build_patient_tokens_by_region(selected_positions_with_flags, POS2GENE)
    tokens_by_prot = _build_patient_tokens_by_region(selected_positions_with_flags, POS2PROT)

    # Ensamblar UPGMA por cada clave
    def _pack_upgma(tokens_by_region, patients):
        out = {}
        for region, by_patient in tokens_by_region.items():
            # Filtrar pacientes que tienen tokens en esta región
            patients_with_tokens = [p for p in patients if by_patient.get(p, set())]
            
            Z_list, leaves = _upgma_from_token_sets(patients, by_patient)
            if Z_list is None:
                # Solo incluir pacientes con tokens en trivial
                out[region] = {"trivial": True, "patients": patients_with_tokens}
            else:
                # Solo incluir pacientes con tokens en resultado
                out[region] = {"patients": patients_with_tokens, "leaves_order": list(leaves), "Z": Z_list}
        return out

    upgma_genes = _pack_upgma(tokens_by_gene, patients)
    upgma_proteins = _pack_upgma(tokens_by_prot, patients)

    upgma_blob = "__UPGMA_BLOB__"

    # Segments for genome ribbon (genes & proteins)
    segments = {
        "genes": [{"name": n, "start": s, "end": e} for s,e,n in GENE_LIST],
        "proteins": [{"name": n, "start": s, "end": e} for s,e,n in CDS_LIST]
    }
    segments_blob = "__SEGMENTS_BLOB__"
    # ───────── JavaScript (links included) ─────────
    JS = """
         console.log('=== JAVASCRIPT LOADED ===');
         // Comprobar Plotly y cargar fallback si fuese necesario
          // === CARGA PLOTLY DE FORMA SEGURA ===
          (function ensurePlotly(){
            if (typeof Plotly === 'undefined'){
              const s = document.createElement('script');
              s.src = 'https://cdn.plot.ly/plotly-2.32.0.min.js';
              s.defer = true;
              s.onload = () => { console.log('[Plotly] loaded'); whenPlotlyReady(() => window.__do_initial_render__ && window.__do_initial_render__()); };
              s.onerror = () => { console.error('[Plotly] failed to load'); };
              document.head.appendChild(s);
            }
          })();

          // Botón de código
          const codeBtn = document.getElementById('code-btn');
          if (codeBtn) {
            codeBtn.onclick = function() {
              window.open('https://github.com/ldgonzalezvazquez/SRASel', '_blank');
            };
          }

          // Modal de metodología
          const methodologyBtn = document.getElementById('methodology-btn');
          const methodologyModal = document.getElementById('methodology-modal');
          const closeMethodologyBtn = document.getElementById('close-methodology-modal');

          if (methodologyBtn && methodologyModal && closeMethodologyBtn) {
            methodologyBtn.onclick = function() {
              methodologyModal.style.display = 'block';
            };

            closeMethodologyBtn.onclick = function() {
              methodologyModal.style.display = 'none';
            };

            // Close modal when clicking outside of it
            window.onclick = function(event) {
              if (event.target === methodologyModal) {
                methodologyModal.style.display = 'none';
              }
            };
          }

          function whenPlotlyReady(cb){
            const start = () => { try { cb(); } catch(e){ console.error(e); } };
            if (document.readyState === 'loading'){
              document.addEventListener('DOMContentLoaded', () => {
                if (window.Plotly) return start();
                const t = setInterval(() => { if (window.Plotly){ clearInterval(t); start(); } }, 50);
              });
            } else {
              if (window.Plotly) return start();
              const t = setInterval(() => { if (window.Plotly){ clearInterval(t); start(); } }, 50);
            }
          }
         function parseJsonScript(id, fallback){
           const el = document.getElementById(id);
           if (!el) return fallback;
           try {
             return JSON.parse(el.textContent);
           } catch (e) {
             console.error(`Error parsing JSON script #${id}:`, e);
             return fallback;
           }
         }
         const DATA = {
           gene    : parseJsonScript('gene', null),
           protein : parseJsonScript('protein', null)
         };
        const requestedInitialMode = parseJsonScript('initial-mode', 'gene');
        const toggleTarget = parseJsonScript('toggle-target', '');
        let mode = (requestedInitialMode && DATA[requestedInitialMode])
          ? requestedInitialMode
          : (DATA.gene ? 'gene' : 'protein');
        const fig        = document.getElementById('fig');
        const varFig     = document.getElementById('var-fig');
        const back       = document.getElementById('back');
        const toggle     = document.getElementById('toggle');
        const patientSel = document.getElementById('patient-selector');
        const varPlotDiv = document.getElementById('var-plot');
        const varData    = parseJsonScript('var-data', {});
         const positionMapping = parseJsonScript('position-mapping', { genes:{}, proteins:{}, aa_positions:{} });
         const UPGMA_GENES = parseJsonScript('upgma-genes', null);
         const UPGMA_PROT  = parseJsonScript('upgma-proteins', null);
         const SEGMENTS    = parseJsonScript('segments', { genes:[], proteins:[] });
        const selectedPositionsData = parseJsonScript('selected-positions', {});
        const SELECTED_POSITIONS = selectedPositionsData; // Alias para compatibilidad
        function getVarDetailed() {
          return (DATA[mode] && DATA[mode].var_detailed) ? DATA[mode].var_detailed : {};
        }

        // nuevos selectores de listas y contenedores de barras
        const geneList     = document.getElementById('gene-list');
        const nspList      = document.getElementById('nsp-list');
        const geneChartDiv = document.getElementById('gene-chart');
        const nspChartDiv  = document.getElementById('nsp-chart');
        const geneControls = document.getElementById('gene-controls');
        const proteinControls = document.getElementById('protein-controls');
        const applyFiltersBtn = document.getElementById('apply-filters');
        const resetFiltersBtn = document.getElementById('reset-filters');
        const minPatientsInput = document.getElementById('min-patients');

        function syncToggleLabel() {
          if (!toggle) return;
          toggle.textContent = (mode === 'gene') ? '🔄 Protein View' : '🔄 Gene View';
        }

        // Estado de filtros de visibilidad SOLO para UPGMA
        let upgmaFilterState = {
          categories: {}, // {column: true/false}
          includePatientsNotInFilterData: false
        };

        function getMinPatientsThreshold(){
          console.log('getMinPatientsThreshold called');
          console.log('minPatientsInput:', minPatientsInput);
          console.log('minPatientsInput value:', minPatientsInput ? minPatientsInput.value : 'undefined');
          const v = parseInt(minPatientsInput && minPatientsInput.value || '1', 10);
          const result = isNaN(v) ? 1 : Math.max(1, v);
          console.log('getMinPatientsThreshold returning:', result);
          return result;
        }


        // ===== UPGMA dendrogram (prune-only) =====
        function buildTreeFromLinkage(Z, n) {
          const left = new Int32Array(n-1), right = new Int32Array(n-1), height = new Float64Array(n-1);
          for (let r=0; r<n-1; r++){
            const row = Z[r];
            left[r]  = row[0]|0;
            right[r] = row[1]|0;
            height[r]= +row[2];
          }
          return {left, right, height, n, root: n+(n-1)-1};
        }
        function findPatientKeyForUpgma(patientName) {
          if (!filterData || !filterData.patient_data) {
            return null;
          }
          if (patientName in filterData.patient_data) {
            return patientName;
          }
          const patientNameMapping = filterData.patient_name_mapping || {};
          if (patientNameMapping && typeof patientNameMapping === 'object') {
            for (const [key, variations] of Object.entries(patientNameMapping)) {
              if (Array.isArray(variations) && variations.includes(patientName)) {
                for (const variation of variations) {
                  if (variation in filterData.patient_data) {
                    return variation;
                  }
                }
                return key;
              }
            }
          }
          return null;
        }
        function patientHasCategoryData(patientKey, category) {
          if (!filterData || !filterData.patient_data || !patientKey) {
            return false;
          }
          const patientData = filterData.patient_data[patientKey] || {};
          const values = patientData[category] || [];
          return Array.isArray(values) && values.some(v => v !== "__NO_DATA__" && String(v).trim() !== "");
        }
        function getUpgmaVisiblePatientSet(patients) {
          if (!patients || patients.length === 0) {
            return new Set();
          }
          if (!filterData || !filterData.columns || filterData.columns.length === 0) {
            return new Set(patients);
          }
          const selectedCategories = Object.keys(upgmaFilterState.categories || {}).filter(
            c => upgmaFilterState.categories[c]
          );
          if (selectedCategories.length === 0) {
            return new Set(patients);
          }
          const visible = new Set();
          for (const patient of patients) {
            const patientKey = findPatientKeyForUpgma(patient);
            const inFilterData = patientKey && filterData.patient_data && patientKey in filterData.patient_data;
            if (!inFilterData) {
              if (upgmaFilterState.includePatientsNotInFilterData) {
                visible.add(patient);
              }
              continue;
            }
            let hasAll = true;
            for (const category of selectedCategories) {
              if (!patientHasCategoryData(patientKey, category)) {
                hasAll = false;
                break;
              }
            }
            if (hasAll) {
              visible.add(patient);
            }
          }
          return visible;
        }
        function getUpgmaHighlightSet(selectedPatients) {
          const highlight = new Set();
          if (!selectedPatients || selectedPatients.length === 0) {
            return highlight;
          }
          if (!filterData || !filterData.patient_data) {
            return highlight;
          }
          for (const patient of selectedPatients) {
            const patientKey = findPatientKeyForUpgma(patient);
            if (patientKey && patientKey in filterData.patient_data) {
              highlight.add(patient);
            }
          }
          return highlight;
        }
        function renderDendrogramPruned(containerId, upgmaDataMap, selectedItemKey, selectedPatients){
          const data = upgmaDataMap[selectedItemKey] || upgmaDataMap["__ALL__"];
          if (!data || data.trivial){ 
            const el = document.getElementById(containerId);
            if (el) el.innerHTML = '<div class="warn">Not enough signal to build this tree.</div>';
            return;
          }
          const patients = data.patients;
          const leavesOrder = data.leaves_order;
          const n = patients.length;
          const tree = buildTreeFromLinkage(data.Z, n);
          // active leaves in order (mostrar siempre todos los pacientes)
          const highlightSet = getUpgmaHighlightSet(selectedPatients);
          const visibleSet = getUpgmaVisiblePatientSet(patients);
          const activeLeaves = leavesOrder.filter(i => visibleSet.has(patients[i]));
          if (n < 2){
            const el = document.getElementById(containerId);
            if (el) el.innerHTML = '<div class="warn">Not enough patients to build this tree.</div>';
            return;
          }
          if (activeLeaves.length < 2){
            const el = document.getElementById(containerId);
            if (el) el.innerHTML = '<div class="warn">Not enough visible patients to build this tree.</div>';
            return;
          }
          const xOfLeaf = new Map();
          activeLeaves.forEach((leaf, i) => xOfLeaf.set(leaf, i));
          const xs=[], ys=[];
          function layout(idx){
            if (idx < n){
              const active = xOfLeaf.has(idx);
              const x = active ? xOfLeaf.get(idx) : null;
              return {active, x, xmin:x, xmax:x, h:0};
            }
            const r = idx - n;
            const L = tree.left[r], R = tree.right[r], H = tree.height[r];
            const A = layout(L), B = layout(R);
            if (A.active && B.active){
              const yL = (A.xmin + A.xmax)/2, yR = (B.xmin + B.xmax)/2;
              xs.push(A.h, H, null); ys.push(yL, yL, null);
              xs.push(B.h, H, null); ys.push(yR, yR, null);
              xs.push(H,  H, null); ys.push(yL, yR, null);
              return {active:true, xmin:Math.min(A.xmin,B.xmin), xmax:Math.max(A.xmax,B.xmax), h:H};
            } else if (A.active){ return A; }
              else if (B.active){ return B; }
              else { return {active:false, xmin:null, xmax:null, h:H}; }
          }
          layout(tree.root);
          const unselectedLabels = [];
          const unselectedY = [];
          const selectedLabels = [];
          const selectedY = [];
          activeLeaves.forEach(i => {
            const name = patients[i];
            const y = xOfLeaf.get(i);
            if (highlightSet.has(name)) {
              selectedLabels.push(name);
              selectedY.push(y);
            } else {
              unselectedLabels.push(name);
              unselectedY.push(y);
            }
          });
          const trace = { 
            x: xs, y: ys, mode: "lines", type: "scatter", hoverinfo:"skip",
            line: { color: '#34495e', width: 2 }
          };
          const traces = [trace];
          if (unselectedLabels.length > 0) {
            traces.push({
              x: unselectedLabels.map(_ => -0.05),
              y: unselectedY,
              text: unselectedLabels, textposition: "middle left",
              mode: "text", type: "scatter", hoverinfo:"skip",
              textfont: { size: 10, color: '#2c3e50' }
            });
          }
          const highlightAnnotations = selectedLabels.map((label, idx) => ({
            x: -0.05,
            y: selectedY[idx],
            xref: 'x',
            yref: 'y',
            text: label,
            showarrow: false,
            xanchor: 'left',
            align: 'left',
            bgcolor: 'rgba(46, 204, 113, 0.28)',
            bordercolor: 'rgba(46, 204, 113, 0.45)',
            borderpad: 2,
            font: { size: 10, color: '#2c3e50' }
          }));
          Plotly.newPlot(containerId, traces, {
            xaxis: { showline:false, showticklabels:false, showgrid:false, zeroline:false },
            yaxis: { showline:false, showticklabels:false, showgrid:false, zeroline:false },
            margin: {l:20, r:20, b:20, t:20},
            showlegend: false,
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            annotations: highlightAnnotations
          }, {responsive:true, displayModeBar:false});
        }

        function drawUPGMA(){
          const selectedPatients = getSelectedPatients();
          const selectedItems = getSelectedItems();
          
          // Determinar si todas las proteínas/genes están seleccionados
          const allItems = mode === 'gene' ? 
            Array.from(geneList.options).map(o => o.value) : 
            Array.from(nspList.options).map(o => o.value);
          
          const isAllSelected = selectedItems.length === allItems.length && 
            allItems.every(item => selectedItems.includes(item));
          
          // Usar "__ALL__" si todas están seleccionadas, sino usar la primera seleccionada
          const key = isAllSelected ? "__ALL__" : (selectedItems && selectedItems.length ? selectedItems[0] : "__ALL__");
          
          console.log('UPGMA Debug:', {
            mode: mode,
            selectedItems: selectedItems,
            allItems: allItems,
            isAllSelected: isAllSelected,
            key: key
          });
          
          if (mode === 'gene') {
            if (document.getElementById('upgma-genes-div')) renderDendrogramPruned('upgma-genes-div', UPGMA_GENES, key, selectedPatients);
            const other = document.getElementById('upgma-proteins-div'); if (other) other.innerHTML = '';
          } else {
            if (document.getElementById('upgma-proteins-div')) renderDendrogramPruned('upgma-proteins-div', UPGMA_PROT, key, selectedPatients);
            const other = document.getElementById('upgma-genes-div'); if (other) other.innerHTML = '';
          }
        }
        // Variables para almacenar las selecciones actuales
        let currentPatients = [];
        let currentGenes = [];
        let currentProteins = [];

        // Función para obtener genes/proteínas seleccionados
        function getSelectedItems() {
          if (mode === 'gene') {
            return Array.from(geneList.selectedOptions).map(o => o.value);
          } else {
            return Array.from(nspList.selectedOptions).map(o => o.value);
          }
        }

        // Función para obtener pacientes seleccionados
        function getSelectedPatients() {
          const selected = Array.from(patientSel.selectedOptions).map(o => o.value);
          const allOptions = Array.from(patientSel.options).map(o => o.value);
          
          console.log('getSelectedPatients - selected options:', selected);
          console.log('getSelectedPatients - all options:', allOptions);
          
          // Verificar si hay filtros avanzados activos (panel de resumen visible)
          const filtersSummary = document.getElementById('filters-summary');
          const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
          
          // Si hay filtros avanzados activos y no hay pacientes seleccionados, devolver una lista vacía
          if (advancedFiltersActive && selected.length === 0) {
            console.log('getSelectedPatients - advanced filters active and no patients selected, returning empty array.');
            return [];
          }
          
          // Si no hay selección específica (y no hay filtros avanzados activos), devolver todos los pacientes disponibles
          if (selected.length === 0) {
            console.log('getSelectedPatients - no selection and no advanced filters, returning all:', allOptions);
            return allOptions;
          }
          
          console.log('getSelectedPatients - returning selected:', selected);
          return selected;
        }

        // Función para aplicar filtros
        function applyFilters() {
          console.log('=== APPLY FILTERS CALLED ===');
          
          const selectedPatients = getSelectedPatients();
          const selectedItems = getSelectedItems();
          
          console.log('Current selections:', {
            patients: selectedPatients,
            items: selectedItems,
            mode: mode
          });
          
          // Redibujar todas las gráficas
          drawVarHeat();
          drawHeat();
          drawVarPlot();
          drawUPGMA();
          if (mode === 'gene') {
            drawGeneBars();
          } else {
            drawNspBars();
          }
        }

        // Función para resetear filtros
        function resetFilters() {
          // Seleccionar todos los pacientes
          for (let i = 0; i < patientSel.options.length; i++) {
            patientSel.options[i].selected = true;
          }
          
          // Seleccionar todos los genes
          for (let i = 0; i < geneList.options.length; i++) {
            geneList.options[i].selected = true;
          }
          
          // Seleccionar todas las proteínas
          for (let i = 0; i < nspList.options.length; i++) {
            nspList.options[i].selected = true;
          }
          
          // Aplicar filtros reseteados
          applyFilters();
        }

        // Función para filtrar el heatmap según selección
        function filterHeatmap() {
          console.log('=== FILTERHEATMAP CALLED ===');
          // Heatmap principal = β (slope) mean analysis (Figure 2)
          // Filtrado real (pacientes y items) preservando fila Average (index 0)
          const selectedPatients = getSelectedPatients();
          const selectedItems = getSelectedItems();
          
          console.log('Selected patients:', selectedPatients);
          console.log('Selected items:', selectedItems);
          console.log('Mode:', mode);
          
          // Si no hay datos disponibles, mostrar error
          if (!DATA[mode] || !DATA[mode].heat) {
            console.error('No heatmap data available for mode:', mode);
            return;
          }
          
          const heatData = DATA[mode].heat.data[0];
          const layout = DATA[mode].heat.layout;
          
          console.log('Original heatmap data:', {
            x: heatData.x ? heatData.x.length : 'no x',
            y: heatData.y ? heatData.y.length : 'no y',
            z: heatData.z ? heatData.z.length : 'no z'
          });
          
          // --- Selección de filas (pacientes) EXCLUYENDO SIEMPRE la fila 'Average' original ---
          const yToIdx = new Map(heatData.y.map((lab,i)=>[lab,i]));
          // Verificar si hay filtros avanzados activos
          const filtersSummary = document.getElementById('filters-summary');
          const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
          
          const patientIndices = (selectedPatients.length > 0)
            ? selectedPatients
                .map(p => yToIdx.get(p))
                .filter(i => typeof i === 'number' && i !== 0)   // nunca incluir 'Average' original
            : (advancedFiltersActive 
                ? [] // Si hay filtros avanzados activos y no hay pacientes, no mostrar ninguno (solo Average)
                : Array.from({length: heatData.y.length}, (_,i)=>i)
                    .filter(i => heatData.y[i] !== 'Average'));        // sin selección → todos menos 'Average'
          
          // Etiquetas iniciales (sin Average) y datos Z de pacientes
          let filteredY = patientIndices.map(i => heatData.y[i]);
          let filteredZ = patientIndices.map((pi, k) => {
            const row = heatData.z[pi];
            if (!row || !Array.isArray(row)) {
              console.error('Invalid row at index', k, ':', row);
              return new Array(heatData.x.length).fill(0);
            }
            return row.slice();
          });
          // --- Filtrado por columnas (genes/proteínas) si hay selección ---
          // Usamos colIdx para evitar conflictos con variables existentes
          const colIdx = (selectedItems.length > 0)
            ? heatData.x.map((lab, i) => selectedItems.includes(lab) ? i : -1).filter(i => i !== -1)
            : Array.from({length: heatData.x.length}, (_, i) => i);

          // Etiquetas de columnas filtradas
          let filteredX = colIdx.map(i => heatData.x[i]);

          // Filtrar Z por columnas según colIdx
          filteredZ = filteredZ.map(row => {
            if (!row || !Array.isArray(row)) {
              console.error('Invalid row data:', row);
              return new Array(colIdx.length).fill(0);
            }
            return colIdx.map(i => {
              if (i >= row.length) {
                console.error(`Index ${i} out of bounds for row length ${row.length}`);
                return 0;
              }
              return (row[i] !== undefined && row[i] !== null) ? row[i] : 0;
            });
          });
          // ================== Average DINÁMICO (encabezado) ==================
          // Calcular la media por columna con las filas actualmente filtradas
          const nSel = filteredZ.length; // nº de pacientes usados
          const avgRow = (filteredX.length > 0)
            ? filteredX.map((_, j) => {
                let sum = 0, n = 0;
                for (let r = 0; r < filteredZ.length; r++) {
                  const v = filteredZ[r][j];
                  if (typeof v === 'number' && !Number.isNaN(v)) { sum += v; n++; }
                }
                return n ? +(sum / n).toFixed(4) : null;
              })
            : [];
          // Insertar Average al principio de Y y Z
          if (avgRow.length > 0) {
            filteredY.unshift("Average");
            filteredZ.unshift(avgRow);
          }
          
          console.log('Filtered data sizes:', {
            x: filteredX.length,
            y: filteredY.length,
            z: filteredZ.length
          });
          console.log('Filtered Z data sample:', filteredZ.slice(0, 2));
          console.log('Filtered X:', filteredX);
          console.log('Filtered Y:', filteredY);
          
            // Verificar si hay valores no-cero en los datos
          const hasNonZeroValues = filteredZ.some(row => Array.isArray(row) && row.some(val => val !== 0 && val !== null && val !== undefined));
          console.log('Has non-zero values:', hasNonZeroValues);
          console.log('Z value range:', {
            min: Math.min(...filteredZ.flat().filter(v => v !== null && v !== undefined)),
            max: Math.max(...filteredZ.flat().filter(v => v !== null && v !== undefined))
          });
          
          // --- hovertext y customdata ALINEADOS incluyendo Average dinámico ---
          const baseHT = Array.isArray(heatData.hovertext) ? heatData.hovertext : [];
          const baseCD = Array.isArray(heatData.customdata) ? heatData.customdata : [];
          // Pacientes (sin Average) ya filtrados por columnas
          let filteredHover = patientIndices.map(r => {
            const hrow = Array.isArray(baseHT[r]) ? baseHT[r] : [];
            return colIdx.map(c => (hrow[c] !== undefined ? hrow[c] : ''));
          });
          let filteredCustom2 = patientIndices.map(r => {
            const crow = Array.isArray(baseCD[r]) ? baseCD[r] : [];
            return colIdx.map(c => (crow[c] !== undefined ? crow[c] : '{}'));
          });
          // Fila de hover para el Average dinámico (muestra valor y n pacientes)
          if (avgRow.length > 0) {
            const avgHoverRow = filteredX.map((lab, j) => {
              const v = avgRow[j];
              return (v === null)
                ? `Average | ${lab}: NA (n=${nSel})`
                : `Average | ${lab}: ${v.toFixed(4)} (n=${nSel})`;
            });
            filteredHover.unshift(avgHoverRow);
            filteredCustom2.unshift(filteredX.map(_ => '{}'));
          }

          const filteredData = [{
            type: 'heatmap',
            x: filteredX,
            y: filteredY,
            z: filteredZ,
            colorscale: heatData.colorscale,
            zmin: heatData.zmin,
            zmax: heatData.zmax,
            zmid: heatData.zmid,
            hovertext: filteredHover,
            customdata: filteredCustom2,
            hovertemplate: heatData.hovertemplate || '%{x}, %{y}: %{z}<extra></extra>',
            showscale: (typeof heatData.showscale !== 'undefined') ? heatData.showscale : true,
            colorbar: heatData.colorbar || {title: 'β', len: 0.87}
          }];

          const filteredLayout = {
            ...layout,
            xaxis: {
              ...layout.xaxis,
              tickvals: filteredX,
              ticktext: filteredX.map(c => c.length > 15 ? c.substring(0, 15) + '...' : c)
            },
            yaxis: {
              ...layout.yaxis,
              tickvals: filteredY,
              ticktext: filteredY
            },
            height: getHeatmapHeight(filteredY.length)
          };
          
          console.log('Applying filtered heatmap with data:', filteredData);
          Plotly.newPlot(fig, filteredData, filteredLayout)
            .then(() => {
              fig.on('plotly_click', ev => {
                const region = ev.points[0].x;
                const patient = ev.points[0].y;
                if (patient === 'Average') {
                  console.log('[HEATMAP1] Click en Average: se ignora.');
                  return;
                }
                console.log('[HEATMAP1] Click:', { patient, region, mode });

                // Obtener el valor Z del punto clickeado
                const zValue = ev.points[0].z || 0;
                
                // Obtener filtros actuales
                const minEvt = +(document.getElementById('min-events')?.value || 0);
                const minMr = +(document.getElementById('min-reads')?.value || 0);
                const minMf = +(document.getElementById('min-freq')?.value || 0) / 100;

                // Buscar eventos que contribuyen a esta posición
                const events = [];
                if (DATA.events) {
                  for (const evt of DATA.events) {
                    if (evt.gene === region && evt.sample === patient) {
                      const meetsFilters = eventPasses(evt, minMr, minMf, minEvt);
                      events.push({
                        event: {
                          id: `${evt.gene}-${evt.sample}`,
                          gene: evt.gene,
                          sample: evt.sample,
                          region: `${Math.min(evt.pA, evt.pB)}-${Math.max(evt.pA, evt.pB)}`
                        },
                        minEvents: getEventSupport(evt),
                        mireads: evt.mireads || 0,
                        minfreq: evt.minfreq || 0,
                        meetsFilters: meetsFilters
                      });
                    }
                  }
                }

                // Mostrar popup con eventos detallados - DESHABILITADO
                // showEventsPopup(region, patient, zValue, events, minEvt, minMr, minMf);
              });
            });
        }

        function drawVarHeat(){
          try {
            const selectedPatients = getSelectedPatients();
            const selectedItems = getSelectedItems();
            
            console.log('Drawing variability heatmap with filters:', {
              patients: selectedPatients,
              items: selectedItems,
              mode: mode
            });
            console.log('selectedPatients length:', selectedPatients.length);
            console.log('selectedPatients values:', selectedPatients);
            
            // Obtener datos originales
            const varHeatData = DATA[mode].var_heat.data[0];
            const varHeatLayout = DATA[mode].var_heat.layout;
            
            // Índices de filas/columnas según selección (respetando orden de selección)
            const yToIdx = new Map(varHeatData.y.map((label, i) => [label, i]));
            // Verificar si hay filtros avanzados activos
            const filtersSummary = document.getElementById('filters-summary');
            const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
            
            const patientIndices = (selectedPatients.length > 0)
              ? selectedPatients
                  .map(p => yToIdx.get(p))
                  .filter(i => typeof i === 'number' && i !== 0) // nunca incluir "Average" original
              : (advancedFiltersActive 
                  ? [] // Si hay filtros avanzados activos y no hay pacientes, no mostrar ninguno (solo Average)
                  : Array.from({length: varHeatData.y.length}, (_, i) => i)
                      .filter(i => varHeatData.y[i] !== 'Average')); // excluir "Average" original

            const itemIndices = (selectedItems.length > 0)
              ? varHeatData.x.map((lab, i) => selectedItems.includes(lab) ? i : -1).filter(i => i !== -1)
              : Array.from({length: varHeatData.x.length}, (_, i) => i);
            
                        // Definir colIdx primero para poder usarlo en el cálculo del average
            const colIdx = (selectedItems.length > 0)
              ? itemIndices
              : Array.from({length: varHeatData.x.length}, (_,i)=>i);

            // Filas a usar en Heatmap 1 (NO hay fila 'Average' original nunca)
            const rowIdx = patientIndices;
            
            console.log('Row indices:', rowIdx);
            console.log('Original Y labels:', varHeatData.y);
            console.log('Selected patients:', selectedPatients);
            // console.log('allIdx:', allIdx);  // removed: allIdx undefined here
            console.log('Filtered rowIdx:', rowIdx);
            console.log('varHeatData.y values:', varHeatData.y);
            
            // Calcular media dinámica basada en pacientes seleccionados
            let avgRow = [];
            if (colIdx.length > 0) {
              console.log('Calculating average for columns:', colIdx);
              avgRow = colIdx.map((colIdx, colIndex) => {
                let sum = 0;
                let count = 0;
                // Solo considerar pacientes seleccionados
                const patientData = rowIdx.map(rowIdx => {
                  if (varHeatData.z && varHeatData.z[rowIdx] && varHeatData.z[rowIdx][colIdx] !== undefined) {
                    return varHeatData.z[rowIdx][colIdx];
                  }
                  return 0.5; // valor neutro para datos faltantes
                });
                
                if (patientData.length > 0) {
                  sum = patientData.reduce((a, b) => a + b, 0);
                  count = patientData.length;
                  const avg = count > 0 ? sum / count : 0.5;
                  if (colIndex < 3) console.log(`Column ${colIdx}: sum=${sum}, count=${count}, avg=${avg}`);
                  return avg;
                }
                return 0.5;
              });
              console.log('Average row calculated:', avgRow.slice(0, 5));
            }

            // Crear datos filtrados coherentes (robustos a filas/columnas indefinidas)
            let filteredY = rowIdx.map(i => (Array.isArray(varHeatData.y) && varHeatData.y[i] !== undefined) ? varHeatData.y[i] : i);
            let filteredX = colIdx.map(i => (Array.isArray(varHeatData.x) && varHeatData.x[i] !== undefined) ? varHeatData.x[i] : i);

            // Construir HOVER y CUSTOMDATA filtrados (sin Average por ahora)
            const baseHT = Array.isArray(varHeatData.hovertext) ? varHeatData.hovertext : [];
            const baseCD = Array.isArray(varHeatData.customdata) ? varHeatData.customdata : [];
            let filteredHover = rowIdx.map(r => {
              const hrow = Array.isArray(baseHT[r]) ? baseHT[r] : [];
              return colIdx.map(c => (hrow[c] !== undefined ? hrow[c] : ''));
            });
            let filteredCustom2 = rowIdx.map(r => {
              const crow = Array.isArray(baseCD[r]) ? baseCD[r] : [];
              return colIdx.map(c => (crow[c] !== undefined ? crow[c] : '{}'));
            });
            
            // Agregar etiqueta "Average" para la primera fila
            if (filteredY.length > 0 && avgRow.length > 0) {
              filteredY.unshift("Average");
              // Construir hover dinámico para "Average" usando la selección actual
              const detailed = DATA[mode].var_detailed || {};
              const selectedNames = rowIdx.map(i => varHeatData.y[i]).filter(n => n && n !== "Average");
              const pct = (num, den) => den > 0 ? (num/den*100).toFixed(1) : "0.0";
              const avgHoverRow = colIdx.map((c, j) => {
                const region = filteredX[j];
                let nt_tot=0, nt_var=0, nt_sel=0, nt_neu=0, nt_non=0;
                let aa_tot=0, aa_var=0, aa_sel=0, aa_neu=0, aa_non=0;
                let n_pat = 0;
                for (const pname of selectedNames) {
                  const det = (detailed[pname] && detailed[pname][region]) ? detailed[pname][region] : null;
                  if (!det) continue;
                  const s_nt = Array.isArray(det.selected_nt) ? det.selected_nt.length : (det.selected_nt || 0);
                  const n_nt = Array.isArray(det.neutral_nt) ? det.neutral_nt.length : (det.neutral_nt || 0);
                  const nv_nt = Array.isArray(det.non_variable_nt) ? det.non_variable_nt.length : (det.non_variable_nt || 0);
                  const tot_nt = s_nt + n_nt + nv_nt;
                  const s_aa = Array.isArray(det.selected_aa) ? det.selected_aa.length : (det.selected_aa || 0);
                  const n_aa = Array.isArray(det.neutral_aa) ? det.neutral_aa.length : (det.neutral_aa || 0);
                  const nv_aa = Array.isArray(det.non_variable_aa) ? det.non_variable_aa.length : (det.non_variable_aa || 0);
                  const tot_aa = s_aa + n_aa + nv_aa;
                  nt_tot += tot_nt; nt_sel += s_nt; nt_neu += n_nt; nt_non += nv_nt; nt_var += (s_nt + n_nt);
                  aa_tot += tot_aa; aa_sel += s_aa; aa_neu += n_aa; aa_non += nv_aa; aa_var += (s_aa + n_aa);
                  n_pat++;
                }
                const m_nt_tot = n_pat ? nt_tot / n_pat : 0;
                const m_nt_var = n_pat ? nt_var / n_pat : 0;
                const m_nt_sel = n_pat ? nt_sel / n_pat : 0;
                const m_nt_neu = n_pat ? nt_neu / n_pat : 0;
                const m_nt_non = n_pat ? nt_non / n_pat : 0;
                const m_aa_tot = n_pat ? aa_tot / n_pat : 0;
                const m_aa_var = n_pat ? aa_var / n_pat : 0;
                const m_aa_sel = n_pat ? aa_sel / n_pat : 0;
                const m_aa_neu = n_pat ? aa_neu / n_pat : 0;
                const m_aa_non = n_pat ? aa_non / n_pat : 0;
                const selRatio = (avgRow && typeof avgRow[j] === 'number') ? avgRow[j] : 0.5;
                return (
                  `<b>Average</b> × <b>${region}</b><br>` +
                  `<b>Selection ratio = ${selRatio.toFixed(3)}</b><br>` +
                  `<i>Calculated from ${n_pat} patients</i><br><br>` +
                  `<b>Average Nucleotides:</b><br>` +
                  `• Total: ${m_nt_tot.toFixed(1)}<br>` +
                  `• Variable: ${m_nt_var.toFixed(1)} (${pct(m_nt_var, m_nt_tot)}%)<br>` +
                  `• Selected: ${m_nt_sel.toFixed(1)} (${pct(m_nt_sel, m_nt_tot)}%)<br>` +
                  `• Neutral: ${m_nt_neu.toFixed(1)} (${pct(m_nt_neu, m_nt_tot)}%)<br>` +
                  `• Non-variable: ${m_nt_non.toFixed(1)} (${pct(m_nt_non, m_nt_tot)}%)<br><br>` +
                  `<b>Average Amino acids:</b><br>` +
                  `• Total: ${m_aa_tot.toFixed(1)}<br>` +
                  `• Variable: ${m_aa_var.toFixed(1)} (${pct(m_aa_var, m_aa_tot)}%)<br>` +
                  `• Selected: ${m_aa_sel.toFixed(1)} (${pct(m_aa_sel, m_aa_tot)}%)<br>` +
                  `• Neutral: ${m_aa_neu.toFixed(1)} (${pct(m_aa_neu, m_aa_tot)}%)<br>` +
                  `• Non-variable: ${m_aa_non.toFixed(1)} (${pct(m_aa_non, m_aa_tot)}%)`
                );
              });
              filteredHover.unshift(avgHoverRow);
              filteredCustom2.unshift(colIdx.map(_ => '{}'));
            }
            
            console.log('Filtered Y before processing:', filteredY);
            console.log('Filtered X before processing:', filteredX);
            
            // Construir matriz Z con la fila de media dinámica
            let filteredZ = [];
            
            // Primera fila: media dinámica
            if (avgRow.length > 0) {
              filteredZ.push(avgRow);
              console.log('Average row added:', avgRow);
            }
            
            // Resto de filas: datos de pacientes
            for (let i = 0; i < rowIdx.length; i++) {
              let row = [];
              if (Array.isArray(varHeatData.z)) {
                if (Array.isArray(varHeatData.z[rowIdx[i]])) {
                  row = varHeatData.z[rowIdx[i]];
                } else if (Array.isArray(varHeatData.z[0])) {
                  // si el índice no existe, usa fila 0 como fallback de tamaño
                  row = varHeatData.z[0].map(_ => 0);
                } else {
                  row = colIdx.map(_ => 0);
                }
              } else {
                row = colIdx.map(_ => 0);
              }
              filteredZ.push(colIdx.map(j => (row[j] !== undefined && row[j] !== null) ? row[j] : 0));
            }
            
            console.log('Filtered Z matrix dimensions:', filteredZ.length, 'x', filteredZ[0] ? filteredZ[0].length : 0);
            console.log('Filtered Y (labels):', filteredY);
            console.log('Average row values:', avgRow);
            
            // ---- colores y rango centrados en 0.5 (blanco) ----
            const base = DATA[mode].var_heat.data[0];
            const cs = base && base.colorscale
              ? base.colorscale
              : [[0.0, "rgb(34,139,34)"], [0.5, "white"], [1.0, "rgb(255,140,0)"]];
            const zmin = (base && typeof base.zmin !== "undefined") ? base.zmin : 0.0;
            const zmax = (base && typeof base.zmax !== "undefined") ? base.zmax : 1.0;
            const zmid = (base && typeof base.zmid !== "undefined") ? base.zmid : 0.5;

            const filteredData = [{
              type: 'heatmap',
              x: filteredX,
              y: filteredY,
              z: filteredZ,
              colorscale: cs,
              zmin: zmin,
              zmax: zmax,
              zmid: zmid,
              hovertemplate: (varHeatData.hovertemplate !== undefined
                             ? varHeatData.hovertemplate
                             : '%{x}, %{y}: %{z}<extra></extra>'),
              hovertext: filteredHover,
              customdata: filteredCustom2,
              showscale: (varHeatData.showscale !== undefined) ? varHeatData.showscale : true,
              colorbar: varHeatData.colorbar || {title: 'Selection Ratio (0=neg, 1=pos)', len: 0.87}
            }];
            
            const filteredLayout = {
              ...varHeatLayout,
              xaxis: {
                ...(varHeatLayout && varHeatLayout.xaxis ? varHeatLayout.xaxis : {}),
                tickvals: filteredX,
                ticktext: filteredX.map(c => (c+"").length > 15 ? (c+"").substring(0, 15) + '...' : c)
              },
              yaxis: {
                ...(varHeatLayout && varHeatLayout.yaxis ? varHeatLayout.yaxis : {}),
                tickvals: filteredY,
                ticktext: filteredY
              },
              height: getHeatmapHeight(filteredY.length)
            };
            
            Plotly.newPlot(varFig, filteredData, filteredLayout)
              .then(() => {
                // Agregar evento de click después de dibujar
                varFig.on('plotly_click', function(data) {
              const point = data.points[0];
              const patient = point.y;
              const region = point.x;
              
              console.log('Clicked on:', patient, region);
              
                            // Obtener datos detallados
              console.log('Click detected on:', patient, region);
              console.log('Available data structure:', DATA[mode].var_detailed);
              
              const detailedData = DATA[mode].var_detailed[patient] && DATA[mode].var_detailed[patient][region];
              
              console.log('Detailed data found:', detailedData);
              
              if (!detailedData) {
                console.error('No detailed data found for:', patient, region);
                console.log('Available patients:', Object.keys(DATA[mode].var_detailed));
                if (DATA[mode].var_detailed[patient]) {
                  console.log('Available regions for patient:', Object.keys(DATA[mode].var_detailed[patient]));
                }
                return;
              }
              
              // Crear tabla con información fusionada de las barras
              let tableHTML = `
                <div class="popup-overlay" onclick="closePopup()">
                  <div class="popup-content" onclick="event.stopPropagation()">
                    <div class="popup-header">
                      <h3>Positions in ${region} - ${patient}</h3>
                      <button onclick="closePopup()" class="close-btn">×</button>
                    </div>
                    <div class="popup-body">
                      <div class="position-section">
                        <h4>🟠 Selected Positions (${detailedData.selected_nt.length}) - Nucleotide</h4>
                        <div class="position-list">
                          ${detailedData.selected_nt.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                        </div>
                      </div>
                      <div class="position-section">
                        <h4>🟠 Selected Positions (${detailedData.selected_aa.length}) - Amino Acid</h4>
                        <div class="position-list">
                          ${detailedData.selected_aa.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                        </div>
                      </div>
                      <div class="position-section">
                        <h4>⚪ Neutral Positions (${detailedData.neutral_nt.length}) - Nucleotide</h4>
                        <div class="position-list">
                          ${detailedData.neutral_nt.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                        </div>
                      </div>
                      <div class="position-section">
                        <h4>⚪ Neutral Positions (${detailedData.neutral_aa.length}) - Amino Acid</h4>
                        <div class="position-list">
                          ${detailedData.neutral_aa.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                        </div>
                      </div>
                      <div class="position-section">
                        <h4>🟢 Non-Variable Positions (${detailedData.non_variable_nt.length}) - Nucleotide</h4>
                        <div class="position-list">
                          ${detailedData.non_variable_nt.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                        </div>
                      </div>
                      <div class="position-section">
                        <h4>🟢 Non-Variable Positions (${detailedData.non_variable_aa.length}) - Amino Acid</h4>
                        <div class="position-list">
                          ${detailedData.non_variable_aa.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                        </div>
                      </div>
                      
                      <!-- Nueva sección: Información fusionada de las barras -->
                      <div class="position-section">
                        <h4>📊 Detailed Selection Data</h4>
                        <table class="detailed-table" border="1" cellpadding="8" style="width: 100%; min-width: 600px;">
                          <thead>
                            <tr><th>Position</th><th>β</th><th>p-value</th><th>Amino Acid</th><th>Covarying</th><th>Date</th><th>Plot</th><th>Genome</th></tr>
                          </thead>
                          <tbody>
                            ${(() => {
                              // Buscar información fusionada en los datos de las barras
                              const barsData = DATA[mode].bars;
                              let rows = '';
                              
                              console.log('DEBUG: mode =', mode);
                              console.log('DEBUG: barsData =', barsData);
                              console.log('DEBUG: region =', region);
                              console.log('DEBUG: patient =', patient);
                              
                              if (barsData && barsData[region]) {
                                const barData = barsData[region];
                                if (barData.data && barData.data.length > 0) {
                                  for (const trace of barData.data) {
                                    if (trace.customdata) {
                                      for (let i = 0; i < trace.customdata.length; i++) {
                                        const customData = JSON.parse(trace.customdata[i] || '{}');
                                        if (customData[patient]) {
                                          // Usar nombres completos fusionados de barData.keep en lugar de trace.x
                                          const fullNames = barData.keep || [];
                                          const pos = fullNames[i] || trace.x[i]; // Fallback a trace.x si no hay keep
                                          const data = customData[patient];
                                          const aaDisplay = data.aa_change && data.aa_change.trim() !== '' ? data.aa_change : (pos.includes(' / ') ? pos.split(' / ')[0] : 'N/A');
                                          rows += `
                                            <tr>
                                              <td>${pos}</td>
                                              <td>${data.beta || 'N/A'}</td>
                                              <td>${data.pval || 'N/A'}</td>
                                              <td><strong>${aaDisplay}</strong></td>
                                              <td>${renderCovaryingCell(data, pos)}</td>
                                              <td>${data.date && data.date.trim() !== '' ? data.date : 'N/A'}</td>
                                              <td><a href="${data.plot}" target="_blank">📊 View</a></td>
                                              <td><a href="${data.genome}" target="_blank">🧬 View</a></td>
                                            </tr>`;
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                              
                              return rows || '<tr><td colspan="8">No selected positions</td></tr>';
                            })()}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  </div>
                </div>
              `;
              
              // Limpiar popups anteriores antes de mostrar el nuevo
              const existingPopups = document.querySelectorAll('.popup-overlay');
              existingPopups.forEach(popup => popup.remove());
              
              // Mostrar popup
              document.body.insertAdjacentHTML('beforeend', tableHTML);
            });
              });
            
          } catch (error) {
            console.error('Error drawing variability heatmap:', error);
            Plotly.newPlot(varFig, DATA[mode].var_heat.data, DATA[mode].var_heat.layout)
              .then(() => {
                // Agregar evento de click después de dibujar
                varFig.on('plotly_click', function(data) {
                  const point = data.points[0];
                  const patient = point.y;
                  const region = point.x;
                  
                  console.log('Clicked on:', patient, region);
                  
                  // Obtener datos detallados
                  console.log('Click detected on:', patient, region);
                  console.log('Available data structure:', DATA[mode].var_detailed);
                  
                  const detailedData = DATA[mode].var_detailed[patient] && DATA[mode].var_detailed[patient][region];
                  
                  console.log('Detailed data found:', detailedData);
                  
                  if (!detailedData) {
                    console.error('No detailed data found for:', patient, region);
                    console.log('Available patients:', Object.keys(DATA[mode].var_detailed));
                    if (DATA[mode].var_detailed[patient]) {
                      console.log('Available regions for patient:', Object.keys(DATA[mode].var_detailed[patient]));
                    }
                    return;
                  }
                  
                  // Crear tabla con información fusionada de las barras
                  let tableHTML = `
                    <div class="popup-overlay" onclick="closePopup()">
                      <div class="popup-content" onclick="event.stopPropagation()">
                        <div class="popup-header">
                          <h3>Positions in ${region} - ${patient}</h3>
                          <button onclick="closePopup()" class="close-btn">×</button>
                        </div>
                        <div class="popup-body">
                          <div class="position-section">
                            <h4>🟠 Selected Positions (${detailedData.selected_nt.length}) - Nucleotide</h4>
                            <div class="position-list">
                              ${detailedData.selected_nt.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                            </div>
                          </div>
                          <div class="position-section">
                            <h4>🟠 Selected Positions (${detailedData.selected_aa.length}) - Amino Acid</h4>
                            <div class="position-list">
                              ${detailedData.selected_aa.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                            </div>
                          </div>
                          <div class="position-section">
                            <h4>⚪ Neutral Positions (${detailedData.neutral_nt.length}) - Nucleotide</h4>
                            <div class="position-list">
                              ${detailedData.neutral_nt.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                            </div>
                          </div>
                          <div class="position-section">
                            <h4>⚪ Neutral Positions (${detailedData.neutral_aa.length}) - Amino Acid</h4>
                            <div class="position-list">
                              ${detailedData.neutral_aa.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                            </div>
                          </div>
                          <div class="position-section">
                            <h4>🟢 Non-Variable Positions (${detailedData.non_variable_nt.length}) - Nucleotide</h4>
                            <div class="position-list">
                              ${detailedData.non_variable_nt.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                            </div>
                          </div>
                          <div class="position-section">
                            <h4>🟢 Non-Variable Positions (${detailedData.non_variable_aa.length}) - Amino Acid</h4>
                            <div class="position-list">
                              ${detailedData.non_variable_aa.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                            </div>
                          </div>
                          
                          <!-- Nueva sección: Información fusionada de las barras -->
                          <div class="position-section">
                            <h4>📊 Detailed Selection Data</h4>
                            <table class="detailed-table" border="1" cellpadding="8" style="width: 100%; min-width: 600px;">
                              <thead>
                                <tr><th>Position</th><th>β</th><th>p-value</th><th>Covarying</th><th>Plot</th><th>Genome</th></tr>
                              </thead>
                              <tbody>
                                ${(() => {
                                  // Buscar información fusionada en los datos de las barras
                                  const barsData = DATA[mode].bars;
                                  let rows = '';
                                  
                                  console.log('DEBUG: mode =', mode);
                                  console.log('DEBUG: barsData =', barsData);
                                  console.log('DEBUG: region =', region);
                                  console.log('DEBUG: patient =', patient);
                                  
                                  if (barsData && barsData[region]) {
                                    const barData = barsData[region];
                                    if (barData.data && barData.data.length > 0) {
                                      for (const trace of barData.data) {
                                        if (trace.customdata) {
                                          for (let i = 0; i < trace.customdata.length; i++) {
                                            const customData = JSON.parse(trace.customdata[i] || '{}');
                                            if (customData[patient]) {
                                              // Usar nombres completos fusionados de barData.keep en lugar de trace.x
                                              const fullNames = barData.keep || [];
                                              const pos = fullNames[i] || trace.x[i]; // Fallback a trace.x si no hay keep
                                              const data = customData[patient];
                                              rows += `
                                                <tr>
                                                  <td>${pos}</td>
                                                  <td>${data.beta || 'N/A'}</td>
                                                  <td>${data.pval || 'N/A'}</td>
                                                  <td>${renderCovaryingCell(data, pos)}</td>
                                                  <td><a href="${data.plot}" target="_blank">📊 View</a></td>
                                                  <td><a href="${data.genome}" target="_blank">🧬 View</a></td>
                                                </tr>`;
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                  
                                  return rows || '<tr><td colspan="6">No selected positions</td></tr>';
                                })()}
                              </tbody>
                            </table>
                          </div>
                        </div>
                      </div>
                    </div>
                  `;
                  
                  // Limpiar popups anteriores antes de mostrar el nuevo
                  const existingPopups = document.querySelectorAll('.popup-overlay');
                  existingPopups.forEach(popup => popup.remove());
                  
                  // Mostrar popup
                  document.body.insertAdjacentHTML('beforeend', tableHTML);
                });
              });
          }
        }

        function drawHeat(){
          try {
            filterHeatmap();
          } catch (error) {
            console.error('Error drawing heatmap:', error);
            // Fallback: mostrar datos originales
          Plotly.newPlot(fig, DATA[mode].heat.data, DATA[mode].heat.layout);
          }
          back.style.display = 'none';
        }

        function getHeatmapHeight(rowCount) {
          return Math.max(450, rowCount * 20);
        }

        function drawBar(cat){
          const B = DATA[mode].bars[cat];
          if (!B) return;
          Plotly.newPlot(fig, B.data, B.layout);
          back.style.display = 'inline';
          fig.on('plotly_click', ev => {
            const det = JSON.parse(ev.points[0].customdata || '{}');
            if (Object.keys(det).length) {
              // Pasar información adicional sobre la traza clickeada
              const traceName = ev.points[0].data.name || '';
              const pointIndex = ev.points[0].pointIndex || 0;
              popupGenePlot(cat, ev.points[0].x, det, traceName, pointIndex);
            }
          });
        }

        function renderCellPopup(patient, region, container) {
          // COPIAR EXACTAMENTE EL CÓDIGO DEL HEATMAP 1 ORIGINAL
          console.log('Clicked on:', patient, region);
          
          // Obtener datos detallados
          console.log('Click detected on:', patient, region);
          console.log('Available data structure:', DATA[mode].var_detailed);
          
          const detailedData = DATA[mode].var_detailed[patient] && DATA[mode].var_detailed[patient][region];
          
          console.log('Detailed data found:', detailedData);
          
          if (!detailedData) {
            console.error('No detailed data found for:', patient, region);
            console.log('Available patients:', Object.keys(DATA[mode].var_detailed));
            if (DATA[mode].var_detailed[patient]) {
              console.log('Available regions for patient:', Object.keys(DATA[mode].var_detailed[patient]));
            }
            return;
          }
          
          // Crear tabla con información fusionada de las barras
          let tableHTML = `
            <div class="popup-overlay" onclick="closePopup()">
              <div class="popup-content" onclick="event.stopPropagation()">
                <div class="popup-header">
                  <h3>Positions in ${region} - ${patient}</h3>
                  <button onclick="closePopup()" class="close-btn">×</button>
                </div>
                <div class="popup-body">
                  <div class="position-section">
                    <h4>🟠 Selected Positions (${detailedData.selected_nt.length}) - Nucleotide</h4>
                    <div class="position-list">
                      ${detailedData.selected_nt.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                    </div>
                  </div>
                  <div class="position-section">
                    <h4>🟠 Selected Positions (${detailedData.selected_aa.length}) - Amino Acid</h4>
                    <div class="position-list">
                      ${detailedData.selected_aa.map(pos => `<span class="position-tag selected">${pos}</span>`).join('')}
                    </div>
                  </div>
                  <div class="position-section">
                    <h4>⚪ Neutral Positions (${detailedData.neutral_nt.length}) - Nucleotide</h4>
                    <div class="position-list">
                      ${detailedData.neutral_nt.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                    </div>
                  </div>
                  <div class="position-section">
                    <h4>⚪ Neutral Positions (${detailedData.neutral_aa.length}) - Amino Acid</h4>
                    <div class="position-list">
                      ${detailedData.neutral_aa.map(pos => `<span class="position-tag neutral">${pos}</span>`).join('')}
                    </div>
                  </div>
                  <div class="position-section">
                    <h4>🟢 Non-Variable Positions (${detailedData.non_variable_nt.length}) - Nucleotide</h4>
                    <div class="position-list">
                      ${detailedData.non_variable_nt.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                    </div>
                  </div>
                  <div class="position-section">
                    <h4>🟢 Non-Variable Positions (${detailedData.non_variable_aa.length}) - Amino Acid</h4>
                    <div class="position-list">
                      ${detailedData.non_variable_aa.map(pos => `<span class="position-tag non-variable">${pos}</span>`).join('')}
                    </div>
                  </div>
                  
                  <!-- Nueva sección: Información fusionada de las barras -->
                  <div class="position-section">
                    <h4>📊 Detailed Selection Data</h4>
                    <table class="detailed-table" border="1" cellpadding="8" style="width: 100%; min-width: 600px;">
                      <thead>
                        <tr><th>Position</th><th>β</th><th>p-value</th><th>Amino Acid</th><th>Covarying</th><th>Date</th><th>Plot</th><th>Genome</th></tr>
                      </thead>
                      <tbody>
                        ${(() => {
                          // Buscar información fusionada en los datos de las barras
                          const barsData = DATA[mode].bars;
                          let rows = '';
                          
                          console.log('DEBUG: mode =', mode);
                          console.log('DEBUG: barsData =', barsData);
                          console.log('DEBUG: region =', region);
                          console.log('DEBUG: patient =', patient);
                          
                          if (barsData && barsData[region]) {
                            const barData = barsData[region];
                            if (barData.data && barData.data.length > 0) {
                              for (const trace of barData.data) {
                                if (trace.customdata) {
                                  for (let i = 0; i < trace.customdata.length; i++) {
                                    const customData = JSON.parse(trace.customdata[i] || '{}');
                                    if (customData[patient]) {
                                      // Usar nombres completos fusionados de barData.keep en lugar de trace.x
                                      const fullNames = barData.keep || [];
                                      const pos = fullNames[i] || trace.x[i]; // Fallback a trace.x si no hay keep
                                      const data = customData[patient];
                                      const aaDisplay = data.aa_change && data.aa_change.trim() !== '' ? data.aa_change : (pos.includes(' / ') ? pos.split(' / ')[0] : 'N/A');
                                      rows += `
                                        <tr>
                                          <td>${pos}</td>
                                          <td>${data.beta || 'N/A'}</td>
                                          <td>${data.pval || 'N/A'}</td>
                                          <td><strong>${aaDisplay}</strong></td>
                                          <td>${renderCovaryingCell(data, pos)}</td>
                                          <td>${data.date && data.date.trim() !== '' ? data.date : 'N/A'}</td>
                                          <td><a href="${data.plot}" target="_blank">📊 View</a></td>
                                          <td><a href="${data.genome}" target="_blank">🧬 View</a></td>
                                        </tr>`;
                                    }
                                  }
                                }
                              }
                            }
                          }
                          
                          return rows || '<tr><td colspan="8">No selected positions</td></tr>';
                        })()}
                      </tbody>
                    </table>
                  </div>
                </div>
              </div>
            </div>
          `;
          
          // Limpiar popups anteriores antes de mostrar el nuevo
          const existingPopups = document.querySelectorAll('.popup-overlay');
          existingPopups.forEach(popup => popup.remove());
          
          // Mostrar popup
          document.body.insertAdjacentHTML('beforeend', tableHTML);
        }

        // Función mejorada para mostrar popup con botón de cerrar
        function popupGenePlot(cat, x, det, traceName, pointIndex, container){
          // DEBUG: Log de entrada
          console.log('DEBUG popupGenePlot:', {cat, x, traceName, pointIndex, container});
          
          // Buscar el nombre completo fusionado en barData.keep
          let fullName = x; // Fallback a x si no se encuentra
          if (DATA[mode] && DATA[mode].bars && DATA[mode].bars[cat]) {
            const barData = DATA[mode].bars[cat];
            const fullNames = barData.keep || [];
            console.log('DEBUG barData available:', {fullNames, cat, mode, traceName, pointIndex});
            
            // Si solo hay un elemento en fullNames, usarlo directamente
            if (fullNames.length === 1) {
              fullName = fullNames[0];
              console.log('DEBUG using single fullName:', fullName);
            } else {
              // Usar el índice del punto para obtener la conformación específica
              if (pointIndex !== undefined && pointIndex < fullNames.length) {
                fullName = fullNames[pointIndex];
                console.log('DEBUG using fullName by pointIndex:', fullName);
              } else {
                // Fallback: buscar por coincidencia parcial en la posición nucleotídica
              const ntPos = x.split(' - ')[1]; // ej: "28226 (111)"
              const matchingFullName = fullNames.find(fn => fn.includes(ntPos.split(' ')[0]));
              if (matchingFullName) {
                fullName = matchingFullName;
                console.log('DEBUG found matching fullName:', fullName);
              } else {
                console.log('DEBUG no matching fullName found, using x as fallback');
                }
              }
            }
          }
          
          // Usar exactamente la etiqueta clicada para evitar desalineación por índices/filtros.
          let title = (typeof x === 'string' && x.trim()) ? x.trim() : String(cat || '').trim();
          if (title && cat && !title.startsWith(`${cat} - `) && !title.startsWith(`${cat}:`)) {
            title = `${cat} - ${title}`;
          }
          console.log('DEBUG final title (from clicked label):', title);
          
          // Obtener datos de posiciones seleccionadas con flags
          const selectedPositionsWithFlags = parseJsonScript('selected-positions-flags', {});
          
          // Crear tabla detallada de posiciones seleccionadas
          let detailedTable = '';
          if (Object.keys(det).length > 0) {
            const rows = [];
            Object.keys(det).sort().forEach(p => {
              const data = det[p];
              
              // Manejar múltiples betas y p-values como en la Figura 3
              let beta, pval;
              if (data.beta != null && data.pval != null) {
                // Mantener <br/> para mostrar múltiples valores en líneas separadas
                beta = data.beta;
                pval = data.pval;
              } else {
                beta = (data.beta != null) ? parseFloat(data.beta).toFixed(4) : 'N/A';
                pval = (data.pval != null) ? parseFloat(data.pval).toFixed(4) : 'N/A';
              }
              
              // Construir rutas dinámicamente como en la Figura 3
              const dataDir = parseJsonScript('data-directory', '');
              const extractNtPos = () => {
                if (data && data.pos !== undefined && data.pos !== null && String(data.pos).trim() !== '') {
                  return String(data.pos).trim();
                }
                if (typeof x === 'string') {
                  const mX = x.match(/-\s*(\d+)\s*(?:\(|$)/);
                  if (mX) return mX[1];
                }
                if (typeof fullName === 'string') {
                  const mFull = fullName.match(/\/\s*(\d+)\s*$/);
                  if (mFull) return mFull[1];
                }
                return '';
              };
              const pos = extractNtPos();
              const plotPath = pos
                ? `${dataDir}/${p}_plots/pos${pos.padStart(5, '0')}_ALL.pdf`
                : '#';
              const genomePath = `${dataDir}/${p}_genome_selection.pdf`;
              
              // Determinar color según el beta (lógica como en Figura 3)
              let rowColor = '#ffffff'; // blanco por defecto (para β mixtos)
              if (data.beta != null && typeof data.beta === 'string') {
                // Manejar diferentes formatos de beta
                let hasPositive = false;
                let hasNegative = false;
                
                // Verificar si hay múltiples betas separados por <br/>
                if (data.beta.includes('<br/>')) {
                  const betaValues = data.beta.split('<br/>');
                  betaValues.forEach(betaStr => {
                    const betaNum = parseFloat(betaStr.trim());
                    if (!isNaN(betaNum)) {
                      if (betaNum > 0) hasPositive = true;
                      if (betaNum < 0) hasNegative = true;
                    }
                  });
                } else {
                  // Un solo valor
                const betaNum = parseFloat(data.beta);
                  if (!isNaN(betaNum)) {
                    if (betaNum > 0) hasPositive = true;
                    if (betaNum < 0) hasNegative = true;
                  }
                }
                
                // Asignar color según el tipo de β
                if (hasPositive && hasNegative) {
                  rowColor = '#fff3e0'; // naranja claro para β mixtos
                } else if (hasPositive) {
                  rowColor = '#ffebee'; // rojo claro para β positivo
                } else if (hasNegative) {
                  rowColor = '#e3f2fd'; // azul claro para β negativo
                }
              }
              
              const posLabel = data.aa_change ? `${data.aa_change}/${pos}` : `${pos}`;
              // Mostrar aminoácidos correctamente
              const aaDisplay = data.aa_change && data.aa_change.trim() !== '' ? data.aa_change : (posLabel.includes('/') ? posLabel.split('/')[0] : 'N/A');
              // Mostrar enlaces directamente como en la Figura 3
              const plotLink = `<a href="${plotPath}" target="_blank">📊 View</a>`;
              const genomeLink = `<a href="${genomePath}" target="_blank">🧬 View</a>`;
              
              rows.push(`
                <tr style="background-color: ${rowColor};">
                  <td><strong>${p}</strong></td>
                  <td><strong>${beta}</strong></td>
                  <td><strong>${pval}</strong></td>
                  <td><strong>${aaDisplay}</strong></td>
                  <td>${renderCovaryingCell(data, pos)}</td>
                  <td>${data.date && data.date.trim() !== '' ? data.date : 'N/A'}</td>
                  <td>${plotLink}</td>
                  <td>${genomeLink}</td>
                </tr>
              `);
            });
            
            if (rows.length > 0) {
              detailedTable = `
                <div class="position-section">
                  <h4>📊 Detailed Selection Data</h4>
                  <table class="detailed-table" border="1" cellpadding="8" style="width:100%;min-width:600px">
                    <thead>
                      <tr><th>Patient</th><th>β</th><th>p-value</th><th>Amino Acid</th><th>Covarying</th><th>Date</th><th>Position PDF</th><th>Genome PDF</th></tr>
                    </thead>
                    <tbody>${rows.join('')}</tbody>
                  </table>
                </div>`;
            }
          }
          
          let html = `<div class="popup-table">
            <button class="close-btn" onclick="closePopup('${container}')">×</button>
            <h3>${title}</h3>
            ${detailedTable}
          </div>`;

          // Mostrar en el contenedor correcto
          if (container === 'gene') {
            geneChartDiv.innerHTML = html;
          } else {
            nspChartDiv.innerHTML = html;
          }
        }

        function renderCovaryingCell(data, pos) {
          // Función unificada para renderizar la celda de covariación
          if (!data || !data.covarying) {
            return '<span style="color: #95a5a6;">✗ No</span>';
          }
          
          // Verificar si hay covariación (nuevo formato con | o formato simple)
          // Aceptar True, 'True', o strings que contengan 'True'
          const hasCovarying = data.covarying === true || 
                             data.covarying === 'True' || 
                             (typeof data.covarying === 'string' && 
                              (data.covarying.includes('True') || data.covarying.includes('true')));
          
          if (!hasCovarying) {
            return '<span style="color: #95a5a6;">✗ No</span>';
          }
          
          // Si hay múltiples betas (separados por <br/>), mostrar múltiples enlaces
          if (data.beta && data.beta.includes('<br/>')) {
            const betaCount = (data.beta.match(/<br\/>/g) || []).length + 1;
            const covaryingLinks = [];
            
            // Dividir los datos de covariación por variante usando el separador especial |||
            const covaryingWithValues = data.covarying_with ? data.covarying_with.split(' ||| ') : [''];
            const covaryingCountValues = data.covarying_count ? String(data.covarying_count).split(' | ') : ['0'];
            
            // Dividir los aa_change por variante (separados por <br/>)
            const aaChangeValues = data.aa_change ? data.aa_change.split('<br/>') : [''];
            
            console.log(`[RENDER COVARYING CELL] Position ${pos}:`, {
              covaryingWith: data.covarying_with,
              covaryingWithValues: covaryingWithValues,
              covaryingCountValues: covaryingCountValues,
              aaChangeValues: aaChangeValues,
              betaCount: betaCount
            });
            
            for (let i = 0; i < betaCount; i++) {
              const currentCovaryingWith = covaryingWithValues[i] || '';
              const currentCovaryingCount = covaryingCountValues[i] ? parseInt(covaryingCountValues[i]) : 0;
              const currentAaChange = aaChangeValues[i] || '';
              
              // DEBUG: Verificar qué se está pasando al modal
              console.log(`[RENDER COVARYING CELL] Variant ${i}:`, {
                position: pos,
                covaryingWith: currentCovaryingWith,
                covaryingCount: currentCovaryingCount,
                aaChange: currentAaChange
              });
              
              // Si no hay covariaciones (count = 0), mostrar "No"
              if (currentCovaryingCount === 0) {
                covaryingLinks.push(`<span style="color: #95a5a6;">✗ No</span>`);
              } else {
                covaryingLinks.push(`<span style="color: #e74c3c; font-weight: bold; cursor: pointer;" onclick="showCovaryingDetails('${pos}', '${currentCovaryingWith}', ${i}, '${currentAaChange}')">✓ Yes (${currentCovaryingCount})</span>`);
              }
            }
            
            return covaryingLinks.join('<br/>');
          } else {
            // Caso simple: una sola variante
            const covaryingCount = data.covarying_count ? parseInt(data.covarying_count) : 0;
            const aaChange = data.aa_change || '';
            if (covaryingCount === 0) {
              return `<span style="color: #95a5a6;">✗ No</span>`;
            } else {
              return `<span style="color: #e74c3c; font-weight: bold; cursor: pointer;" onclick="showCovaryingDetails('${pos}', '${data.covarying_with || ''}', 0, '${aaChange}')">✓ Yes (${covaryingCount})</span>`;
            }
          }
        }

        function formatCovaryingWith(covaryingWithStr) {
          if (!covaryingWithStr || covaryingWithStr.trim() === '') {
            return 'To keep this HTML lightweight for large patient cohorts, detailed covariation partners are shown only in Figure 5. Please check Figure 5 for the complete list.';
          }
          
          console.log('=== FORMAT COVARYING WITH DEBUG ===');
          console.log('Input string:', covaryingWithStr);
          console.log('Input length:', covaryingWithStr.length);
          
          // Dividir por " | " para separar cada covariación
          const covaryingItems = covaryingWithStr.split(' | ');
          console.log('Split result:', covaryingItems);
          console.log('Number of items:', covaryingItems.length);
          
          const formattedItems = [];
          
          for (let i = 0; i < covaryingItems.length; i++) {
            const item = covaryingItems[i];
            const trimmedItem = item.trim();
            console.log(`Processing item ${i}:`, trimmedItem);
            
            if (!trimmedItem) continue;
            
            // Parsear el formato: nt=25511(C); aa=ORF3a:40S [Δβ=0.xx, Δday=X]
            if (trimmedItem.includes('nt=') && trimmedItem.includes('aa=')) {
              const ntMatch = trimmedItem.split('nt=')[1].split(';')[0];
              const aaPart = trimmedItem.split('aa=')[1].split(' [')[0];
              
              // Extraer diferencia de beta y día del contenido entre corchetes
              let deltaBetaMatch = '';
              let deltaDayMatch = '';
              if (trimmedItem.includes('[') && trimmedItem.includes(']')) {
                const bracketContent = trimmedItem.split('[')[1].split(']')[0];
                // Puede ser "Δβ=0.xx, Δday=X" o solo "Δβ=0.xx"
                if (bracketContent.includes(', Δday=')) {
                  const parts = bracketContent.split(', ');
                  deltaBetaMatch = parts[0].replace('Δβ=', '');
                  deltaDayMatch = parts[1].replace('Δday=', '');
                } else if (bracketContent.includes('Δβ=')) {
                  deltaBetaMatch = bracketContent.replace('Δβ=', '');
                }
              }
              
              console.log(`Item ${i} parsed:`, {ntMatch, aaPart, deltaBetaMatch, deltaDayMatch});
              
              if (ntMatch && aaPart) {
                let formattedItem = '<div style="margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 5px; border-left: 3px solid #3498db;">';
                formattedItem += '<strong>Nucleotide:</strong> ' + ntMatch + '<br/>';
                formattedItem += '<strong>Amino acid:</strong> ' + aaPart + '<br/>';
                if (deltaBetaMatch) {
                  formattedItem += '<strong>Δβ:</strong> ' + deltaBetaMatch;
                }
                if (deltaDayMatch) {
                  formattedItem += ' | <strong>Δday:</strong> ' + deltaDayMatch + ' days';
                }
                formattedItem += '</div>';
                formattedItems.push(formattedItem);
                console.log(`Item ${i} formatted successfully`);
              }
            } else {
              // Si no sigue el formato esperado, mostrar tal como está
              formattedItems.push('<div style="margin-bottom: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;">' + trimmedItem + '</div>');
              console.log(`Item ${i} added as raw text`);
            }
          }
          
          console.log('Final formatted items count:', formattedItems.length);
          console.log('Final result:', formattedItems.join(''));
          
          return formattedItems.length > 0 ? formattedItems.join('') : 'To keep this HTML lightweight for large patient cohorts, detailed covariation partners are shown only in Figure 5. Please check Figure 5 for the complete list.';
        }

        function getFallbackCovaryingWith(position) {
          try {
            const posStr = String(position);
            const posNum = Number(position);
            const byPatient = selectedPositionsWithFlags || {};
            for (const patientData of Object.values(byPatient)) {
              if (!patientData || typeof patientData !== 'object') continue;
              const rawEntry = patientData[posStr] ?? patientData[posNum];
              if (!rawEntry) continue;
              const entries = Array.isArray(rawEntry) ? rawEntry : [rawEntry];
              for (const e of entries) {
                const cw = (e && e.covarying_with) ? String(e.covarying_with).trim() : '';
                if (cw) return cw;
              }
            }
          } catch (err) {
            console.warn('Fallback covarying lookup error:', err);
          }
          return '';
        }

        function showCovaryingDetails(position, covaryingWith, variantIndex = 0, aaChange = '') {
          if (!covaryingWith || String(covaryingWith).trim() === '') {
            covaryingWith = getFallbackCovaryingWith(position);
          }
          // DEBUG: Imprimir datos recibidos
          console.log('=== SHOW COVARYING DETAILS DEBUG ===');
          console.log('Position:', position);
          console.log('Variant Index:', variantIndex);
          console.log('Covarying With:', covaryingWith);
          console.log('AA Change:', aaChange);
          console.log('Covarying With length:', covaryingWith ? covaryingWith.length : 0);
          console.log('Covarying With split by | :', covaryingWith ? covaryingWith.split(' | ').length : 0);
          console.log('Covarying With split result:', covaryingWith ? covaryingWith.split(' | ') : []);
          
          // DEBUG: Verificar si hay covariaciones
          if (covaryingWith && covaryingWith.trim() !== '') {
            const covaryingItems = covaryingWith.split(' | ');
            console.log('Number of covarying items:', covaryingItems.length);
            console.log('First few items:', covaryingItems.slice(0, 3));
          }
          
          // Crear título específico para la variante usando el aa_change
          let title = `🔗 Covarying Details for Position ${position}`;
          if (aaChange && aaChange.trim() !== '') {
            title = `🔗 Covarying Details for Position ${position}`;
          }
          
          console.log('Final title:', title);
          
          // Crear modal para mostrar detalles de covariación
          const modal = document.createElement('div');
          modal.id = 'covarying-modal';
          modal.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.5);
            z-index: 10000;
            display: flex;
            justify-content: center;
            align-items: center;
          `;
          
          modal.innerHTML = `
            <div style="
              background-color: white;
              padding: 20px;
              border-radius: 10px;
              max-width: 80%;
              max-height: 80%;
              overflow-y: auto;
              box-shadow: 0 4px 20px rgba(0,0,0,0.3);
            ">
              <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
                <h3 style="margin: 0; color: #2c3e50;">🔗 Covarying Details for Position ${position}</h3>
                <button onclick="closeCovaryingModal()" style="
                  background: #e74c3c;
                  color: white;
                  border: none;
                  border-radius: 50%;
                  width: 30px;
                  height: 30px;
                  font-size: 16px;
                  cursor: pointer;
                  display: flex;
                  align-items: center;
                  justify-content: center;
                  font-weight: bold;
                ">×</button>
              </div>
              <div style="line-height: 1.6; color: #34495e;">
                <p><strong>Covarying with:</strong></p>
                <div style="
                  max-height: 500px;
                  overflow-y: auto;
                ">${formatCovaryingWith(covaryingWith)}</div>
              </div>
            </div>
          `;
          
          document.body.appendChild(modal);
        }

        function closeCovaryingModal() {
          const modal = document.getElementById('covarying-modal');
          if (modal) {
            modal.remove();
          }
        }

        // Funciones auxiliares para el manejo de eventos del heatmap
        function getEventSupport(ev) {
          // Detecta diferentes claves de soporte
          return ev.support ?? ev.events ?? ev.n ?? ev.n_events ?? ev.support_count ?? ev.bigger ?? ev.big ?? 0;
        }

        function eventPasses(ev, minMr, minMf, minEvt) {
          if (minMr && ev.mireads < minMr) return false;
          if (minMf && ev.minfreq < minMf) return false;
          if (minEvt && getEventSupport(ev) < minEvt) return false;
          return true;
        }

        function showEventsPopup(gene, sample, zValue, events, minEvt, minMr, minMf) {
          const modal = document.getElementById('events-modal');
          const content = document.getElementById('events-content');
          // Construir contribuyentes (paciente -> lista de tiempos) como en _past_
          const byPatient = new Map();
          for (const e of events) {
            if (!e.meetsFilters) continue;
            const smp = (e.event && e.event.sample) ? e.event.sample : '';
            if (!smp) continue;
            const [pat, t] = smp.split(' | ');
            if (!byPatient.has(pat)) byPatient.set(pat, new Set());
            if (t) byPatient.get(pat).add(String(t));
          }
          const contribRows = [...byPatient.entries()]
              .map(([pat, times]) => `<div><b>${pat}</b>: ${[...times].sort().join(', ')}</div>`)
              .join('') || '<i>No contributors passing filters.</i>';

          let html = `
            <div style="margin-bottom: 15px;">
              <strong>Gene:</strong> ${gene}<br>
              <strong>Sample:</strong> ${sample}<br>
              <strong>Z-score:</strong> ${zValue.toFixed(4)}<br>
              <strong>Current Filters:</strong> Min Events: ${minEvt}, Min Reads: ${minMr}, Min Freq: ${minMf}
            </div>
            <div style="margin: 10px 0; padding: 8px; background:#f7f9fc; border-radius:8px;">
              <div style="font-weight:600; margin-bottom:6px;">Contributors (patient → times)</div>
              ${contribRows}
            </div>
          `;
          
          if (events.length === 0) {
            html += '<p>No events found for this combination.</p>';
          } else {
            // Separar eventos que cumplen y no cumplen filtros
            const meetsFilters = events.filter(e => e.meetsFilters);
            const filteredOut = events.filter(e => !e.meetsFilters);
            
            html += `<h4>Events (${events.length} total)</h4>`;
            
            if (meetsFilters.length > 0) {
              html += `<h5 style="color: #28a745;">✓ Meets Current Filters (${meetsFilters.length})</h5>`;
              html += '<table style="width: 100%; border-collapse: collapse; margin-bottom: 20px;">';
              html += '<tr style="background: #f8f9fa;"><th style="border: 1px solid #ddd; padding: 8px;">Event ID</th><th style="border: 1px solid #ddd; padding: 8px;">Min Events</th><th style="border: 1px solid #ddd; padding: 8px;">Min Reads</th><th style="border: 1px solid #ddd; padding: 8px;">Min Freq</th></tr>';
              
              meetsFilters.forEach(eventData => {
                const event = eventData.event;
                html += `<tr>
                  <td style="border: 1px solid #ddd; padding: 8px;">${event.id || 'N/A'}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${eventData.minEvents || 1}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${eventData.mireads ?? 0}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${(eventData.minfreq ?? 0).toFixed(4)}</td>
                </tr>`;
              });
              html += '</table>';
            }
            
            if (filteredOut.length > 0) {
              html += `<h5 style="color: #dc3545;">✗ Filtered Out (${filteredOut.length})</h5>`;
              html += '<table style="width: 100%; border-collapse: collapse;">';
              html += '<tr style="background: #f8f9fa;"><th style="border: 1px solid #ddd; padding: 8px;">Event ID</th><th style="border: 1px solid #ddd; padding: 8px;">Min Events</th><th style="border: 1px solid #ddd; padding: 8px;">Min Reads</th><th style="border: 1px solid #ddd; padding: 8px;">Min Freq</th><th style="border: 1px solid #ddd; padding: 8px;">Reason</th></tr>';
              
              filteredOut.forEach(eventData => {
                const event = eventData.event;
                let reason = [];
                if ((eventData.minEvents || 1) < minEvt) reason.push(`Min Events < ${minEvt}`);
                if ((eventData.mireads ?? 0) < minMr) reason.push(`Reads < ${minMr}`);
                if ((eventData.minfreq ?? 0) < minMf) reason.push(`Freq < ${minMf}`);
                
                html += `<tr>
                  <td style="border: 1px solid #ddd; padding: 8px;">${event.id || 'N/A'}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${eventData.minEvents || 1}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${eventData.mireads ?? 0}</td>
                  <td style="border: 1px solid #ddd; padding: 8px;">${(eventData.minfreq ?? 0).toFixed(4)}</td>
                  <td style="border: 1px solid #ddd; padding: 8px; color: #dc3545;">${reason.join(', ')}</td>
                </tr>`;
              });
              html += '</table>';
            }
          }
          
          content.innerHTML = html;
          modal.style.display = 'flex';
        }

        function closePopup(container) {
          if (container === 'gene') {
            drawGeneBars();
          } else if (container === 'protein') {
            drawNspBars();
          } else {
            // Cerrar popup general (para el heatmap de variabilidad)
            const popup = document.querySelector('.popup-overlay');
            if (popup) {
              popup.remove();
            }
          }
        }

        const BAR_COLOR_POSITIVE = '#e74c3c'; // rojo: solo beta positiva
        const BAR_COLOR_NEGATIVE = '#3498db'; // azul: solo beta negativa
        const BAR_COLOR_MIXED = '#f39c12';    // naranja: mixto
        const BAR_COLOR_NEUTRAL = '#95a5a6';  // gris: sin beta interpretable

        function inferSelectionDirectionColorFromCustomData(customDataObj) {
          try {
            const rawBeta = customDataObj ? customDataObj.beta : null;
            if (rawBeta === null || rawBeta === undefined || rawBeta === '') return BAR_COLOR_NEUTRAL;
            const parts = String(rawBeta).split('<br/>').map(v => parseFloat(String(v).trim())).filter(v => !Number.isNaN(v));
            if (!parts.length) return BAR_COLOR_NEUTRAL;
            const hasPos = parts.some(v => v > 0);
            const hasNeg = parts.some(v => v < 0);
            if (hasPos && hasNeg) return BAR_COLOR_MIXED;
            if (hasPos) return BAR_COLOR_POSITIVE;
            if (hasNeg) return BAR_COLOR_NEGATIVE;
            return BAR_COLOR_NEUTRAL;
          } catch (_e) {
            return BAR_COLOR_NEUTRAL;
          }
        }

        function drawGeneBars(){
          if (!DATA.gene || !DATA.gene.bars) {
            geneChartDiv.innerHTML = '<div class="loading">Gene data not available in this view.</div>';
            return;
          }
          let sels = Array.from(geneList.selectedOptions).map(o => o.value);
          
          // Si no hay genes seleccionados, mostrar todos los genes
          if (sels.length === 0){
            sels = Array.from(geneList.options).map(o => o.value);
          }
          
          // Combinar datos de todos los genes seleccionados
          let allTraces = [];
          let allXLabels = [];
          let allKeep = [];
          let maxHeight = 300;
          
          for (const cat of sels) {
            const B = DATA['gene'].bars[cat];
            if (!B || !B.data || B.data.length === 0) {
              continue;
            }
            
            // Agregar prefijo de gen a las etiquetas X
            const geneTraces = B.data.map(trace => {
              const newTrace = {...trace};
              if (newTrace.x && newTrace.x.length > 0) {
                newTrace.x = newTrace.x.map(x => `${cat} - ${x}`);
              }
              return newTrace;
            });
            
            allTraces = allTraces.concat(geneTraces);
            
            // Agregar etiquetas X con prefijo de gen
            if (B.keep && B.keep.length > 0) {
              const geneXLabels = B.keep.map(site => {
                if (site && site.includes(' / ')) {
                  const parts = site.split(' / ');
                  if (parts.length === 2) {
                    const nt_pos = parts[1];
                    const aa_part = parts[0];
                    const aa_match = aa_part.match(/(\d+)/);
                    if (aa_match) {
                      const aa_pos = aa_match[1];
                      return `${cat} - ${nt_pos} (${aa_pos})`;
                    }
                  }
                }
                return `${cat} - ${site}`;
              });
              allXLabels = allXLabels.concat(geneXLabels);
              allKeep = allKeep.concat(B.keep);
            }
            
            maxHeight = Math.max(maxHeight, 290 + B.keep.length * 6);
          }
          
          if (allTraces.length === 0) {
            geneChartDiv.innerHTML = '<div class="loading">No valid sites for selected genes.</div>';
            return;
          }
          
          // Filtrar por pacientes seleccionados
          const selectedPatients = getSelectedPatients();
          const minPatients = getMinPatientsThreshold();
          
          // Primero, calcular pacientes únicos por posición (sin prefijo de gen)
          const positionPatientCount = {};
          const positionCovaryingData = {};
          const positionDateData = {};
          
          for (const cat of sels) {
            const B = DATA['gene'].bars[cat];
            if (!B || !B.data || B.data.length === 0) continue;
            
            for (const trace of B.data) {
              if (trace.customdata) {
                for (let i = 0; i < trace.customdata.length; i++) {
                  const customData = JSON.parse(trace.customdata[i] || '{}');
                  // Verificar si hay filtros avanzados activos
                  const filtersSummary = document.getElementById('filters-summary');
                  const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
                  
                  const relevantPatients = Object.keys(customData).filter(p => {
                    // Si hay filtros avanzados activos y no hay pacientes seleccionados, no incluir ningún paciente
                    if (advancedFiltersActive && selectedPatients.length === 0) {
                      return false;
                    }
                    // Si no hay selección específica (y no hay filtros avanzados), incluir todos
                    if (selectedPatients.length === 0) {
                      return true;
                    }
                    // Incluir solo los pacientes seleccionados
                    return selectedPatients.includes(p);
                  });
                  
                  // Extraer posición original (sin prefijo de gen)
                  const originalX = trace.x ? trace.x[i] : '';
                  
                  if (!positionPatientCount[originalX]) {
                    positionPatientCount[originalX] = new Set();
                    positionCovaryingData[originalX] = {
                      covarying: false,
                      covarying_count: 0,
                      covarying_with: ''
                    };
                    positionDateData[originalX] = [];
                  }
                  
                  // Agregar pacientes únicos a esta posición
                  relevantPatients.forEach(p => positionPatientCount[originalX].add(p));
                  
                  // Recopilar datos de covariación del primer paciente relevante
                  if (relevantPatients.length > 0 && !positionCovaryingData[originalX].covarying) {
                    const firstPatient = relevantPatients[0];
                    const patientData = customData[firstPatient];
                    console.log('=== COVARYING DEBUG ===');
                    console.log('Position:', originalX);
                    console.log('Patient:', firstPatient);
                    console.log('PatientData:', patientData);
                    console.log('covarying value:', patientData?.covarying);
                    console.log('covarying type:', typeof patientData?.covarying);
                    console.log('covarying_count:', patientData?.covarying_count);
                    console.log('covarying_with:', patientData?.covarying_with);
                    
                    if (patientData && patientData.covarying !== undefined) {
                      const isCovarying = patientData.covarying === true || patientData.covarying === 'True' || patientData.covarying === 'true';
                      positionCovaryingData[originalX].covarying = isCovarying;
                      positionCovaryingData[originalX].covarying_count = parseInt(patientData.covarying_count || 0);
                      positionCovaryingData[originalX].covarying_with = patientData.covarying_with || '';
                      console.log('Final covarying result:', isCovarying);
                    }
                  }
                  
                  // Recopilar datos de fechas de todos los pacientes relevantes
                  relevantPatients.forEach(p => {
                    const patientData = customData[p];
                    if (patientData && patientData.date) {
                      const dates = patientData.date.split("<br>");
                      dates.forEach(date => {
                        if (date && date.trim() !== "") {
                          positionDateData[originalX].push(date.trim());
                        }
                      });
                    }
                  });
                }
              }
            }
          }
          
          // Filtrar posiciones que cumplan el criterio mínimo de pacientes
          console.log('=== FILTER DEBUG ===');
          console.log('minPatients:', minPatients);
          console.log('positionPatientCount keys:', Object.keys(positionPatientCount));
          console.log('positionCovaryingData:', positionCovaryingData);
          
          const validPositions = Object.keys(positionPatientCount).filter(pos => {
            const patientCount = positionPatientCount[pos].size;
            const covaryingData = positionCovaryingData[pos];
            
            console.log(`Position ${pos}: patients=${patientCount}, covarying=${covaryingData.covarying}, count=${covaryingData.covarying_count}`);
            
            // Filtro de pacientes mínimos
            if (patientCount < minPatients) {
              console.log(`Position ${pos}: REJECTED - insufficient patients (${patientCount} < ${minPatients})`);
              return false;
            }
            
            console.log(`Position ${pos}: filter - ACCEPTED`);
            return true;
          });
          
          console.log('Valid positions after filtering:', validPositions);
          
          // Aplicar filtros a todas las trazas combinadas
          let filteredTraces = allTraces.map(trace => {
            const filteredTrace = { ...trace };
            if (trace.customdata) {
              const filteredCustomData = [];
              const filteredY = [];
              const filteredHoverText = [];
              const filteredX = [];
              
              for (let i = 0; i < trace.customdata.length; i++) {
                const xLabel = trace.x ? trace.x[i] : '';
                
                // Extraer posición original (sin prefijo de gen)
                const originalX = xLabel.replace(/^[^-]+ - /, '');
                
                // Solo incluir si la posición cumple el criterio de pacientes mínimos
                if (validPositions.includes(originalX)) {
                const customData = JSON.parse(trace.customdata[i] || '{}');
                // Verificar si hay filtros avanzados activos
                const filtersSummary = document.getElementById('filters-summary');
                const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
                
                const relevantPatients = Object.keys(customData).filter(p => {
                  // Si hay filtros avanzados activos y no hay pacientes seleccionados, no incluir ningún paciente
                  if (advancedFiltersActive && selectedPatients.length === 0) {
                    return false;
                  }
                  // Si no hay selección específica (y no hay filtros avanzados), incluir todos
                  if (selectedPatients.length === 0) {
                    return true;
                  }
                  // Incluir solo los pacientes seleccionados
                  return selectedPatients.includes(p);
                });
                
                  // Crear nuevo customdata solo con pacientes seleccionados
                  const newCustomData = {};
                  relevantPatients.forEach(p => {
                    newCustomData[p] = customData[p];
                  });
                  
                  filteredCustomData.push(JSON.stringify(newCustomData));
                  filteredY.push(relevantPatients.length);
                  filteredHoverText.push(trace.hovertext ? trace.hovertext[i] : '');
                  filteredX.push(trace.x ? trace.x[i] : '');
                }
              }
              
              filteredTrace.customdata = filteredCustomData;
              filteredTrace.y = filteredY;
              filteredTrace.x = filteredX;
              const directionColors = filteredCustomData.map(cd => {
                try {
                  const parsed = JSON.parse(cd || '{}');
                  return inferSelectionDirectionColorFromCustomData(parsed);
                } catch (_e) {
                  return BAR_COLOR_NEUTRAL;
                }
              });
              filteredTrace.marker = {
                ...(trace.marker || {}),
                color: directionColors,
                line: { color: '#2c3e50', width: 0.4 }
              };
              if (trace.hovertext) {
                filteredTrace.hovertext = filteredHoverText;
              }
            }
            return filteredTrace;
          });
          
          // Filtrar trazas vacías
          filteredTraces = filteredTraces.filter(trace => trace.y && trace.y.length > 0);
          
          if (filteredTraces.length === 0) {
            geneChartDiv.innerHTML = '<div class="loading">No valid sites for selected genes after filtering.</div>';
            return;
          }
          
          // Build filtered category list ordered by nucleotide position
          const categorySet = new Set();
          filteredTraces.forEach(trace => (trace.x || []).forEach(xVal => categorySet.add(xVal)));
          const getNtPosFromLabel = (label) => {
            const s = String(label || '');
            const m = s.match(/-\s*(\d+)\s*(?:\(|$)/);
            return m ? parseInt(m[1], 10) : Number.MAX_SAFE_INTEGER;
          };
          const filteredCategoryOrder = Array.from(categorySet).sort((a, b) => {
            const na = getNtPosFromLabel(a);
            const nb = getNtPosFromLabel(b);
            if (na !== nb) return na - nb;
            return String(a).localeCompare(String(b));
          });

          // Crear layout combinado
          const combinedLayout = {
            barmode: "stack",
            showlegend: false,
            plot_bgcolor: "rgba(255,255,255,0.8)",
            paper_bgcolor: "rgba(255,255,255,0)",
            xaxis: {
              type: "category",
              categoryorder: "array",
              categoryarray: filteredCategoryOrder,
              tickangle: -45,
              tickfont_size: 9
            },
            yaxis_title: "# pacientes",
            yaxis: { tickmode: 'linear', dtick: 1, tickformat: 'd' },
            xaxis_title: "Gene - aa / nt",
            margin: {l: 90, r: 20, t: 60, b: 180},
            height: Math.min(850, maxHeight + 50),
            title: { 
              text: `Mutation frequency – Genes: ${sels.join(', ')}`, 
              x: 0.01, 
              xanchor: "left" 
            }
          };
          
          geneChartDiv.innerHTML = "";
          Plotly.newPlot(geneChartDiv, filteredTraces, combinedLayout);
          geneChartDiv.on('plotly_click', ev => {
            const det = JSON.parse(ev.points[0].customdata || '{}');
            if (Object.keys(det).length) {
              // Extraer el gen del título de la barra (ej: "S - 22424 (8)" -> "S")
              const barTitle = ev.points[0].x;
              const geneMatch = barTitle.match(/^([^-]+) - /);
              const gene = geneMatch ? geneMatch[1] : 'Unknown';
              const traceName = ev.points[0].data.name || '';
              const pointIndex = ev.points[0].pointIndex || 0;
              popupGenePlot(gene, barTitle, det, traceName, pointIndex, 'gene');
            }
          });
        }

        function drawNspBars(){
          if (!DATA.protein || !DATA.protein.bars) {
            nspChartDiv.innerHTML = '<div class="loading">Protein data not available in this view.</div>';
            return;
          }
          let sels = Array.from(nspList.selectedOptions).map(o => o.value);
          
          // Si no hay proteínas seleccionadas, mostrar todas las proteínas
          if (sels.length === 0){
            sels = Array.from(nspList.options).map(o => o.value);
          }
          
          // Combinar datos de todas las proteínas seleccionadas
          let allTraces = [];
          let allXLabels = [];
          let allKeep = [];
          let maxHeight = 300;
          
          for (const cat of sels) {
            const B = DATA['protein'].bars[cat];
            if (!B || !B.data || B.data.length === 0) {
              continue;
            }
            
            // Agregar prefijo de proteína a las etiquetas X
            const proteinTraces = B.data.map(trace => {
              const newTrace = {...trace};
              if (newTrace.x && newTrace.x.length > 0) {
                newTrace.x = newTrace.x.map(x => `${cat} - ${x}`);
              }
              return newTrace;
            });
            
            allTraces = allTraces.concat(proteinTraces);
            
            // Agregar etiquetas X con prefijo de proteína
            if (B.keep && B.keep.length > 0) {
              const proteinXLabels = B.keep.map(site => {
                if (site && site.includes(' / ')) {
                  const parts = site.split(' / ');
                  if (parts.length === 2) {
                    const nt_pos = parts[1];
                    const aa_part = parts[0];
                    const aa_match = aa_part.match(/(\d+)/);
                    if (aa_match) {
                      const aa_pos = aa_match[1];
                      return `${cat} - ${nt_pos} (${aa_pos})`;
                    }
                  }
                }
                return `${cat} - ${site}`;
              });
              allXLabels = allXLabels.concat(proteinXLabels);
              allKeep = allKeep.concat(B.keep);
            }
            
            maxHeight = Math.max(maxHeight, 290 + B.keep.length * 6);
          }
          
          if (allTraces.length === 0) {
            nspChartDiv.innerHTML = '<div class="loading">No valid sites for selected proteins.</div>';
            return;
          }
          
          // Filtrar por pacientes seleccionados
          const selectedPatients = getSelectedPatients();
          const minPatients = getMinPatientsThreshold();
          
          // Primero, calcular pacientes únicos por posición (sin prefijo de proteína)
          const positionPatientCount = {};
          const positionDateData = {};
          const positionCovaryingData = {};
          
          for (const cat of sels) {
            const B = DATA['protein'].bars[cat];
            if (!B || !B.data || B.data.length === 0) continue;
            
            for (const trace of B.data) {
              if (trace.customdata) {
                for (let i = 0; i < trace.customdata.length; i++) {
                  const customData = JSON.parse(trace.customdata[i] || '{}');
                  // Verificar si hay filtros avanzados activos
                  const filtersSummary = document.getElementById('filters-summary');
                  const advancedFiltersActive = filtersSummary && filtersSummary.style.display !== 'none';
                  
                  const relevantPatients = Object.keys(customData).filter(p => {
                    // Si hay filtros avanzados activos y no hay pacientes seleccionados, no incluir ningún paciente
                    if (advancedFiltersActive && selectedPatients.length === 0) {
                      return false;
                    }
                    // Si no hay selección específica (y no hay filtros avanzados), incluir todos
                    if (selectedPatients.length === 0) {
                      return true;
                    }
                    // Incluir solo los pacientes seleccionados
                    return selectedPatients.includes(p);
                  });
                  
                  // Extraer posición original (sin prefijo de proteína)
                  const originalX = trace.x ? trace.x[i] : '';
                  
                  if (!positionPatientCount[originalX]) {
                    positionPatientCount[originalX] = new Set();
                    positionCovaryingData[originalX] = {
                      covarying: false,
                      covarying_count: 0,
                      covarying_with: ''
                    };
                    positionDateData[originalX] = [];
                  }
                  
                  // Agregar pacientes únicos a esta posición
                  relevantPatients.forEach(p => positionPatientCount[originalX].add(p));
                  
                  // Recopilar datos de covariación del primer paciente relevante
                  if (relevantPatients.length > 0 && !positionCovaryingData[originalX].covarying) {
                    const firstPatient = relevantPatients[0];
                    const patientData = customData[firstPatient];
                    console.log('=== COVARYING DEBUG ===');
                    console.log('Position:', originalX);
                    console.log('Patient:', firstPatient);
                    console.log('PatientData:', patientData);
                    console.log('covarying value:', patientData?.covarying);
                    console.log('covarying type:', typeof patientData?.covarying);
                    console.log('covarying_count:', patientData?.covarying_count);
                    console.log('covarying_with:', patientData?.covarying_with);
                    
                    if (patientData && patientData.covarying !== undefined) {
                      const isCovarying = patientData.covarying === true || patientData.covarying === 'True' || patientData.covarying === 'true';
                      positionCovaryingData[originalX].covarying = isCovarying;
                      positionCovaryingData[originalX].covarying_count = parseInt(patientData.covarying_count || 0);
                      positionCovaryingData[originalX].covarying_with = patientData.covarying_with || '';
                  
                  // Recopilar datos de fechas de todos los pacientes relevantes
                  relevantPatients.forEach(p => {
                    const patientData = customData[p];
                    if (patientData && patientData.date) {
                      const dates = patientData.date.split("<br>");
                      dates.forEach(date => {
                        if (date && date.trim() !== "") {
                          positionDateData[originalX].push(date.trim());
                        }
                      });
                    }
                  });
                      console.log('Final covarying result:', isCovarying);
                    }
                  }
                }
              }
            }
          }
          
          // Filtrar posiciones que cumplan el criterio mínimo de pacientes
          console.log('=== FILTER DEBUG ===');
          console.log('minPatients:', minPatients);
          console.log('positionPatientCount keys:', Object.keys(positionPatientCount));
          console.log('positionCovaryingData:', positionCovaryingData);
          
          const validPositions = Object.keys(positionPatientCount).filter(pos => {
            const patientCount = positionPatientCount[pos].size;
            const covaryingData = positionCovaryingData[pos];
            
            console.log(`Position ${pos}: patients=${patientCount}, covarying=${covaryingData.covarying}, count=${covaryingData.covarying_count}`);
            
            // Filtro de pacientes mínimos
            
            if (patientCount < minPatients) {
              console.log(`Position ${pos}: REJECTED - insufficient patients (${patientCount} < ${minPatients})`);
              return false;
            }
            
            console.log(`Position ${pos}: filter - ACCEPTED`);
            return true;
          });
          
          console.log('Valid positions after filtering:', validPositions);
          console.log("=== NSP BARS DEBUG ===");
          console.log("selectedPatients:", selectedPatients);
          console.log("minPatients:", minPatients);
          console.log("positionPatientCount:", positionPatientCount);
          console.log("validPositions:", validPositions);
          
          // Aplicar filtros a todas las trazas combinadas
          let filteredTraces = allTraces.map(trace => {
            const filteredTrace = { ...trace };
            if (trace.customdata) {
              const filteredCustomData = [];
              const filteredY = [];
              const filteredHoverText = [];
              const filteredX = [];
              
              for (let i = 0; i < trace.customdata.length; i++) {
                const xLabel = trace.x ? trace.x[i] : '';
                
                // Extraer posición original (sin prefijo de proteína)
                // Usar split para manejar nombres de proteínas con guiones
                const originalX = xLabel.includes(' - ') ? xLabel.split(' - ').slice(1).join(' - ') : xLabel;
                
                // Solo incluir si la posición cumple el criterio de pacientes mínimos
                if (validPositions.includes(originalX)) {
                const customData = JSON.parse(trace.customdata[i] || '{}');
                const relevantPatients = Object.keys(customData).filter(p => 
                  selectedPatients.length === 0 || selectedPatients.includes(p)
                );
                
                  // Solo agregar si hay pacientes relevantes
                  if (relevantPatients.length > 0) {
                  // Crear nuevo customdata solo con pacientes seleccionados
                  const newCustomData = {};
                  relevantPatients.forEach(p => {
                    newCustomData[p] = customData[p];
                  });
                  
                  filteredCustomData.push(JSON.stringify(newCustomData));
                  filteredY.push(relevantPatients.length);
                  filteredHoverText.push(trace.hovertext ? trace.hovertext[i] : '');
                  filteredX.push(trace.x ? trace.x[i] : '');
                }
                }
              }
              
              filteredTrace.customdata = filteredCustomData;
              filteredTrace.y = filteredY;
              filteredTrace.x = filteredX;
              const directionColors = filteredCustomData.map(cd => {
                try {
                  const parsed = JSON.parse(cd || '{}');
                  return inferSelectionDirectionColorFromCustomData(parsed);
                } catch (_e) {
                  return BAR_COLOR_NEUTRAL;
                }
              });
              filteredTrace.marker = {
                ...(trace.marker || {}),
                color: directionColors,
                line: { color: '#2c3e50', width: 0.4 }
              };
              if (trace.hovertext) {
                filteredTrace.hovertext = filteredHoverText;
              }
            }
            return filteredTrace;
          });
          
          // Filtrar trazas vacías
          filteredTraces = filteredTraces.filter(trace => trace.y && trace.y.length > 0);
          
          if (filteredTraces.length === 0) {
            nspChartDiv.innerHTML = '<div class="loading">No valid sites for selected proteins after filtering.</div>';
            return;
          }
          
          // Build filtered category list ordered by nucleotide position
          const categorySet = new Set();
          filteredTraces.forEach(trace => (trace.x || []).forEach(xVal => categorySet.add(xVal)));
          const getNtPosFromLabel = (label) => {
            const s = String(label || '');
            const m = s.match(/-\s*(\d+)\s*(?:\(|$)/);
            return m ? parseInt(m[1], 10) : Number.MAX_SAFE_INTEGER;
          };
          const filteredCategoryOrder = Array.from(categorySet).sort((a, b) => {
            const na = getNtPosFromLabel(a);
            const nb = getNtPosFromLabel(b);
            if (na !== nb) return na - nb;
            return String(a).localeCompare(String(b));
          });

          // Crear layout combinado
          const combinedLayout = {
            barmode: "stack",
            showlegend: false,
            plot_bgcolor: "rgba(255,255,255,0.8)",
            paper_bgcolor: "rgba(255,255,255,0)",
            xaxis: {
              type: "category",
              categoryorder: "array",
              categoryarray: filteredCategoryOrder,
              tickangle: -45,
              tickfont_size: 9
            },
            yaxis_title: "# pacientes",
            yaxis: { tickmode: 'linear', dtick: 1, tickformat: 'd' },
            xaxis_title: "Protein - aa / nt",
            margin: {l: 90, r: 20, t: 60, b: 280},
            height: Math.min(950, maxHeight + 150),
            title: { 
              text: `Mutation frequency – Proteins: ${sels.join(', ')}`, 
              x: 0.01, 
              xanchor: "left" 
            }
          };
          
          nspChartDiv.innerHTML = "";
          Plotly.newPlot(nspChartDiv, filteredTraces, combinedLayout);
          nspChartDiv.on('plotly_click', ev => {
            const det = JSON.parse(ev.points[0].customdata || '{}');
            if (Object.keys(det).length) {
              // Extraer la proteína del título de la barra (ej: "S - 22424 (8)" -> "S")
              const barTitle = ev.points[0].x;
              const proteinMatch = barTitle.match(/^([^-]+) - /);
              const protein = proteinMatch ? proteinMatch[1] : 'Unknown';
              const traceName = ev.points[0].data.name || '';
              const pointIndex = ev.points[0].pointIndex || 0;
              popupGenePlot(protein, barTitle, det, traceName, pointIndex, 'protein');
            }
          });
        }

        function updateSelectors() {
          syncToggleLabel();
          // Seleccionar todos los pacientes por defecto
          console.log('updateSelectors - selecting all patients by default...');
          for (let i = 0; i < patientSel.options.length; i++) {
            patientSel.options[i].selected = true;
          }
          console.log('updateSelectors - all patients selected. Options count:', patientSel.options.length);
          
          if (mode === 'gene') {
            geneControls.classList.remove('hidden');
            proteinControls.classList.add('hidden');
            // Mostrar solo la sección de genes
            document.getElementById('gene-section').style.display = 'block';
            document.getElementById('protein-section').style.display = 'none';
            // UPGMA visible solo para genes
            const pcaG = document.getElementById('upgma-genes-div');
            const pcaP = document.getElementById('upgma-proteins-div');
            if (pcaG) pcaG.parentElement.style.display = 'block';
            if (pcaP) pcaP.parentElement.style.display = 'none';
          } else {
            geneControls.classList.add('hidden');
            proteinControls.classList.remove('hidden');
            // Mostrar solo la sección de proteínas
            document.getElementById('gene-section').style.display = 'none';
            document.getElementById('protein-section').style.display = 'block';
            // UPGMA visible solo para proteínas
            const pcaG = document.getElementById('upgma-genes-div');
            const pcaP = document.getElementById('upgma-proteins-div');
            if (pcaG) pcaG.parentElement.style.display = 'none';
            if (pcaP) pcaP.parentElement.style.display = 'block';
          }
          
          // Actualizar las visualizaciones de heatmaps al cambiar de modo
          updateHeatmap1View();
          updateHeatmap2View();
        }

        // Funciones para cambiar visualizaciones de heatmaps
        function updateHeatmap1View() {
          const view = document.getElementById('heatmap1-view').value;
          let heatmapData;
          
          switch(view) {
            case 'mean':
              heatmapData = DATA[mode].var_heat;
              break;
            case 'selected':
              heatmapData = DATA[mode].heat_selected_ratio;
              break;
            case 'non-variable':
              heatmapData = DATA[mode].heat_non_variable_ratio;
              break;
            default:
              heatmapData = DATA[mode].var_heat;
          }
          
          try {
            // Aplicar filtros antes de mostrar los datos
            const selectedPatients = getSelectedPatients();
            const selectedItems = getSelectedItems();
            
            // Filtrar datos según selección actual
            if (selectedPatients.length > 0 || selectedItems.length > 0) {
              // Aplicar filtrado similar a drawVarHeat()
              const varHeatData = heatmapData.data[0];
              const varHeatLayout = heatmapData.layout;
              
              // Índices de filas/columnas según selección
              const yToIdx = new Map(varHeatData.y.map((label, i) => [label, i]));
              const patientIndices = (selectedPatients.length > 0)
                ? selectedPatients
                    .map(p => yToIdx.get(p))
                    .filter(i => typeof i === 'number' && i !== 0) // nunca incluir "Average" original
                : Array.from({length: varHeatData.y.length}, (_, i) => i)
                    .filter(i => varHeatData.y[i] !== 'Average'); // excluir "Average" original

              const itemIndices = (selectedItems.length > 0)
                ? varHeatData.x.map((lab, i) => selectedItems.includes(lab) ? i : -1).filter(i => i !== -1)
                : Array.from({length: varHeatData.x.length}, (_, i) => i);
              
              // Filtrar datos
              let filteredY = patientIndices.map(i => varHeatData.y[i]);
              let filteredZ = patientIndices.map((pi, k) => {
                const row = varHeatData.z[pi];
                if (!row || !Array.isArray(row)) {
                  return new Array(varHeatData.x.length).fill(0);
                }
                return row.slice();
              });
              
              // Filtrar por columnas
              const colIdx = (selectedItems.length > 0)
                ? varHeatData.x.map((lab, i) => selectedItems.includes(lab) ? i : -1).filter(i => i !== -1)
                : Array.from({length: varHeatData.x.length}, (_, i) => i);

              let filteredX = colIdx.map(i => varHeatData.x[i]);
              filteredZ = filteredZ.map(row => {
                if (!row || !Array.isArray(row)) {
                  return new Array(colIdx.length).fill(0);
                }
                return colIdx.map(i => {
                  if (i >= row.length) return 0;
                  return (row[i] !== undefined && row[i] !== null) ? row[i] : 0;
                });
              });
              
              // Recalcular Average
              const nSel = filteredZ.length;
              const avgRow = (filteredX.length > 0)
                ? filteredX.map((_, j) => {
                    let sum = 0, n = 0;
                    for (let r = 0; r < filteredZ.length; r++) {
                      const v = filteredZ[r][j];
                      if (typeof v === 'number' && !Number.isNaN(v)) { sum += v; n++; }
                    }
                    return n ? +(sum / n).toFixed(4) : null;
                  })
                : [];
              
              // Insertar Average al principio
              if (avgRow.length > 0) {
                filteredY.unshift("Average");
                filteredZ.unshift(avgRow);
              }
              
              // Crear datos filtrados
              const filteredData = [{
                ...varHeatData,
                x: filteredX,
                y: filteredY,
                z: filteredZ
              }];
              
              const heatHeight = getHeatmapHeight(filteredY.length);
              Plotly.newPlot(varFig, filteredData, { ...varHeatLayout, height: heatHeight })
                .then(() => {
                  // Agregar el mismo manejador de clicks que el heatmap 1 original
                  varFig.on('plotly_click', ev => {
                    const region = ev.points[0].x;
                    const patient = ev.points[0].y;
                    if (patient === 'Average') {
                      console.log('[HEATMAP1] Click en Average: se ignora.');
                      return;
                    }
                    console.log('[HEATMAP1] Click:', { patient, region, mode });

                  // SIEMPRE usar los datos del heatmap original (var_heat) para la ventana emergente
                  // independientemente de la visualización actual
                  renderCellPopup(patient, region, /*container*/ 'gene');
                });
              });
            } else {
              // Sin filtros, mostrar datos originales
              Plotly.newPlot(varFig, heatmapData.data, heatmapData.layout)
                .then(() => {
                  // Agregar el mismo manejador de clicks que el heatmap 1 original
                  varFig.on('plotly_click', ev => {
                    const region = ev.points[0].x;
                    const patient = ev.points[0].y;
                    if (patient === 'Average') {
                      console.log('[HEATMAP1] Click en Average: se ignora.');
                      return;
                    }
                    console.log('[HEATMAP1] Click:', { patient, region, mode });
                    // SIEMPRE usar los datos del heatmap original (var_heat) para la ventana emergente
                    // independientemente de la visualización actual
                    renderCellPopup(patient, region, /*container*/ 'gene');
                  });
                });
            }
            
            console.log('Heatmap 1 view updated to:', view, 'with filters applied');
          } catch(e) {
            console.error('Error updating heatmap 1:', e);
          }
        }

        function updateHeatmap2View() {
          const view = document.getElementById('heatmap2-view').value;
          let heatmapData;
          
          switch(view) {
            case 'mean':
              heatmapData = DATA[mode].heat;
              break;
            case 'positive':
              heatmapData = DATA[mode].heat_positive_ratio;
              break;
            case 'negative':
              heatmapData = DATA[mode].heat_negative_ratio;
              break;
            default:
              heatmapData = DATA[mode].heat;
          }
          
          try {
            // Aplicar filtros antes de mostrar los datos
            const selectedPatients = getSelectedPatients();
            const selectedItems = getSelectedItems();
            
            // Filtrar datos según selección actual
            if (selectedPatients.length > 0 || selectedItems.length > 0) {
              // Aplicar filtrado similar a filterHeatmap()
              const heatData = heatmapData.data[0];
              const layout = heatmapData.layout;
              
              // Índices de filas/columnas según selección
              const yToIdx = new Map(heatData.y.map((lab,i)=>[lab,i]));
              const patientIndices = (selectedPatients.length > 0)
                ? selectedPatients
                    .map(p => yToIdx.get(p))
                    .filter(i => typeof i === 'number' && i !== 0)   // nunca incluir 'Average' original
                : Array.from({length: heatData.y.length}, (_,i)=>i)
                    .filter(i => heatData.y[i] !== 'Average');        // sin selección → todos menos 'Average'
              
              // Etiquetas iniciales (sin Average) y datos Z de pacientes
              let filteredY = patientIndices.map(i => heatData.y[i]);
              let filteredZ = patientIndices.map((pi, k) => {
                const row = heatData.z[pi];
                if (!row || !Array.isArray(row)) {
                  return new Array(heatData.x.length).fill(0);
                }
                return row.slice();
              });
              
              // Filtrar por columnas (genes/proteínas) si hay selección
              const colIdx = (selectedItems.length > 0)
                ? heatData.x.map((lab, i) => selectedItems.includes(lab) ? i : -1).filter(i => i !== -1)
                : Array.from({length: heatData.x.length}, (_, i) => i);

              // Etiquetas de columnas filtradas
              let filteredX = colIdx.map(i => heatData.x[i]);

              // Filtrar Z por columnas según colIdx
              filteredZ = filteredZ.map(row => {
                if (!row || !Array.isArray(row)) {
                  return new Array(colIdx.length).fill(0);
                }
                return colIdx.map(i => {
                  if (i >= row.length) {
                    return 0;
                  }
                  return (row[i] !== undefined && row[i] !== null) ? row[i] : 0;
                });
              });
              
              // Recalcular Average
              const avgRow = (filteredX.length > 0)
                ? filteredX.map((_, j) => {
                    let sum = 0, n = 0;
                    for (let r = 0; r < filteredZ.length; r++) {
                      const v = filteredZ[r][j];
                      if (typeof v === 'number' && !Number.isNaN(v)) { sum += v; n++; }
                    }
                    return n ? +(sum / n).toFixed(4) : null;
                  })
                : [];
              
              // Insertar Average al principio de Y y Z
              if (avgRow.length > 0) {
                filteredY.unshift("Average");
                filteredZ.unshift(avgRow);
              }
              
              // Crear datos filtrados
              const filteredData = [{
                ...heatData,
                x: filteredX,
                y: filteredY,
                z: filteredZ
              }];
              
              const heatHeight = getHeatmapHeight(filteredY.length);
              Plotly.newPlot(fig, filteredData, { ...layout, height: heatHeight })
                .then(() => {
                  // Agregar el mismo manejador de clicks que el heatmap 2 original
                  fig.on('plotly_click', ev => {
                    const region = ev.points[0].x;
                    const patient = ev.points[0].y;
                    if (patient === 'Average') {
                      console.log('[HEATMAP2] Click en Average: se ignora.');
                      return;
                    }
                    console.log('[HEATMAP2] Click:', { patient, region, mode });
                    // SIEMPRE usar los datos del heatmap original (var_heat) para la ventana emergente
                    // independientemente de la visualización actual
                    renderCellPopup(patient, region, /*container*/ 'gene');
                  });
                });
            } else {
              // Sin filtros, mostrar datos originales
              Plotly.newPlot(fig, heatmapData.data, heatmapData.layout)
                .then(() => {
                  // Agregar el mismo manejador de clicks que el heatmap 2 original
                  fig.on('plotly_click', ev => {
                    const region = ev.points[0].x;
                    const patient = ev.points[0].y;
                    if (patient === 'Average') {
                      console.log('[HEATMAP2] Click en Average: se ignora.');
                      return;
                    }
                    console.log('[HEATMAP2] Click:', { patient, region, mode });
                    // SIEMPRE usar los datos del heatmap original (var_heat) para la ventana emergente
                    // independientemente de la visualización actual
                    renderCellPopup(patient, region, /*container*/ 'gene');
                  });
                });
            }
            
            console.log('Heatmap 2 view updated to:', view, 'with filters applied');
          } catch(e) {
            console.error('Error updating heatmap 2:', e);
          }
        }

        // Event listeners para los selectores de visualización
        document.getElementById('heatmap1-view').onchange = updateHeatmap1View;
        document.getElementById('heatmap2-view').onchange = updateHeatmap2View;

        back.onclick = drawHeat;
        toggle.onclick = () => {
          const directTarget = parseJsonScript('toggle-target', '') || toggleTarget;
          if (directTarget) {
            window.location.assign(directTarget);
          } else {
            console.warn('toggle-target no definido');
          }
        };
        patientSel.onchange = () => { /* No actualizar automáticamente */ };
        geneList.onchange = () => { /* No actualizar automáticamente */ };
        nspList.onchange  = () => { /* No actualizar automáticamente */ };
        if (minPatientsInput) {
          minPatientsInput.onchange = () => { /* No actualizar automáticamente */ };
        }
        
        
        applyFiltersBtn.onclick = function() {
          console.log('=== APPLY FILTERS BUTTON CLICKED ===');
          applyFilters();
        };
        resetFiltersBtn.onclick = function() {
          filtersApplied = false;
          resetFilters();
          resetAdvancedFiltersState(true);
          updateFiltersSummary();
        };

        // Modal functionality for Min Patients Per Site info
        const minPatientsInfoBtn = document.getElementById('min-patients-info');
        const minPatientsModal = document.getElementById('min-patients-modal');
        const closeModalBtn = document.getElementById('close-modal');

        if (minPatientsInfoBtn && minPatientsModal && closeModalBtn) {
          minPatientsInfoBtn.onclick = function() {
            minPatientsModal.style.display = 'block';
          };

          closeModalBtn.onclick = function() {
            minPatientsModal.style.display = 'none';
          };

          // Close modal when clicking outside of it
          window.onclick = function(event) {
            if (event.target === minPatientsModal) {
              minPatientsModal.style.display = 'none';
            }
          };
        }



        // al cargar la página…
        window.__do_initial_render__ = () => {
          updateSelectors();
          // Dibujo base sin filtros para asegurar render inicial
          try {
            Plotly.newPlot(varFig, DATA[mode].var_heat.data, DATA[mode].var_heat.layout);
          } catch(e) { console.error('base var_heat error', e); }
          try {
            Plotly.newPlot(fig, DATA[mode].heat.data, DATA[mode].heat.layout);
          } catch(e) { console.error('base heat error', e); }
          try {
            if (mode === 'gene' && DATA.gene) drawGeneBars();
            if (mode === 'protein' && DATA.protein) drawNspBars();
          } catch(e) { console.error('base bars error', e); }
          // Ahora aplicar filtros iniciales
          applyFilters();
          // UPGMA inicial
          try {
            if (mode === 'gene') {
              if (document.getElementById('upgma-genes-div')) drawUPGMA();
              const pcadiv = document.getElementById('upgma-proteins-div');
              if (pcadiv) pcadiv.innerHTML = '';
            } else {
              if (document.getElementById('upgma-proteins-div')) drawUPGMA();
              const pcadiv = document.getElementById('upgma-genes-div');
              if (pcadiv) pcadiv.innerHTML = '';
            }
          } catch(e) { console.error('upgma render error', e); }
        };
        whenPlotlyReady(() => window.__do_initial_render__ && window.__do_initial_render__());
        
        function drawVarPlot(){
          console.log('=== DRAWVARPLOT CALLED ===');
          try {
            // Verificar que el div existe
            if (!varPlotDiv) {
              console.error('varPlotDiv not found');
              return;
            }
            
            // Mostrar indicador de carga
            varPlotDiv.innerHTML = '<div class="loading">Loading variability data...</div>';
            
            const selected = getSelectedPatients();
            console.log('Selected patients in drawVarPlot:', selected);
            
            if (selected.length === 0){
              varPlotDiv.innerHTML = '<div class="loading">No patients selected.</div>';
              return;
            }
            
            // Verificar que todas las variables globales estén disponibles
            if (typeof varData === 'undefined') {
              console.error('varData is undefined');
              varPlotDiv.innerHTML = '<div class="loading">Error: varData not available</div>';
              return;
            }
            
            if (typeof positionMapping === 'undefined') {
              console.error('positionMapping is undefined');
              varPlotDiv.innerHTML = '<div class="loading">Error: positionMapping not available</div>';
              return;
            }
            
            if (typeof selectedPositionsData === 'undefined') {
              console.error('selectedPositionsData is undefined');
              varPlotDiv.innerHTML = '<div class="loading">Error: selectedPositionsData not available</div>';
              return;
            }

          const positionMap = new Map();
            const geneMap = new Map();
            const proteinMap = new Map();

            // Obtener datos de variabilidad y mapeo de genes/proteínas
            console.log('varData available:', Object.keys(varData));
            console.log('varData type:', typeof varData);
            console.log('varData sample:', Object.keys(varData).slice(0, 3));
            console.log('varData structure sample:', varData[Object.keys(varData)[0]]);
            console.log('varData first patient data type:', typeof varData[Object.keys(varData)[0]]);
            console.log('varData first patient data keys:', Object.keys(varData[Object.keys(varData)[0]] || {}));
            
            console.log('positionMapping available:', positionMapping);
            console.log('positionMapping type:', typeof positionMapping);
            console.log('positionMapping genes sample:', Object.keys(positionMapping.genes || {}).slice(0, 5));
            console.log('positionMapping proteins sample:', Object.keys(positionMapping.proteins || {}).slice(0, 5));
            
            console.log('selectedPositionsData available:', selectedPositionsData);
            console.log('selectedPositionsData type:', typeof selectedPositionsData);
            console.log('selectedPositionsData sample:', Object.keys(selectedPositionsData || {}).slice(0, 5));
            
            for (const p of selected) {
              console.log('Processing patient:', p);
              const data = varData[p];
              if (!data) {
                console.log('No data for patient:', p);
                continue;
              }
              console.log('Data for patient', p, ':', Object.keys(data).length, 'positions');
              console.log('Sample data for patient', p, ':', Object.entries(data).slice(0, 3));

              for (const [k, v] of Object.entries(data)) {
                const key = +k;
                if (!positionMap.has(key)) {
                  positionMap.set(key, []);
                  // Mapear posición a gen/proteína
                  try {
                    const gene = getGeneAtPosition(key);
                    const protein = getProteinAtPosition(key);
                    geneMap.set(key, gene);
                    proteinMap.set(key, protein);
                  } catch (e) {
                    console.error('Error getting gene/protein for position', key, ':', e);
                  }
                }
                positionMap.get(key).push(v);
              }
            }

          const pos = Array.from(positionMap.keys()).sort((a, b) => a - b);
          console.log('Positions found:', pos.length, 'positions');
          
          const val = pos.map(k => {
            const vals = positionMap.get(k);
            return vals.reduce((a, b) => a + b, 0) / vals.length;
          });

          if (pos.length === 0){
            console.log('No positions found, showing error message');
            varPlotDiv.innerHTML = '<div class="loading">No variability data available.</div>';
            return;
          }

            // Filtrar por gen/proteína seleccionado SOLO si hay una selección específica
            // Si todos los elementos están seleccionados (reset), mostrar todo el genoma
            const selectedItems = getSelectedItems();
            let filteredPos = pos;
            let filteredVal = val;
            
            // Solo filtrar si hay una selección específica (no todos los elementos)
            const allItems = mode === 'gene' ? 
              Array.from(geneList.options).map(o => o.value) : 
              Array.from(nspList.options).map(o => o.value);
            
            const isAllSelected = selectedItems.length === allItems.length && 
              allItems.every(item => selectedItems.includes(item));
            
            if (selectedItems.length > 0 && !isAllSelected) {
              const filteredIndices = [];
              
              console.log('=== FILTERING DEBUG ===');
              console.log('mode:', mode);
              console.log('selectedItems:', selectedItems);
              console.log('pos.length:', pos.length);
              
              for (let i = 0; i < pos.length; i++) {
                const position = pos[i];
                if (mode === 'gene') {
                  const gene = geneMap.get(position);
                  if (selectedItems.includes(gene)) {
                    filteredIndices.push(i);
                    console.log(`Gene match: pos ${position} -> gene ${gene}`);
                  }
                } else {
                  const protein = proteinMap.get(position);
                  if (selectedItems.includes(protein)) {
                    filteredIndices.push(i);
                    console.log(`Protein match: pos ${position} -> protein ${protein}`);
                  }
                }
              }
              
              console.log('filteredIndices.length:', filteredIndices.length);
              filteredPos = filteredIndices.map(i => pos[i]);
              filteredVal = filteredIndices.map(i => val[i]);
              console.log('filteredPos.length:', filteredPos.length);
            }

          // Limpiar primero el contenido
          varPlotDiv.innerHTML = "";

            // Crear hover text con información de gen/proteína y posición aminoacídica
            const hoverText = filteredPos.map(pos => {
              const gene = geneMap.get(pos);
              const protein = proteinMap.get(pos);
              const aaPos = getAminoAcidPosition(pos);
              
              let text = `Position: ${pos}`;
              if (gene) text += `<br>Gene: ${gene}`;
              if (protein) text += `<br>Protein: ${protein}`;
              if (aaPos !== null) text += `<br>Amino Acid Position: ${aaPos}`;
              
              console.log(`Hover text for position ${pos}:`, text);
              return text;
            });

            // Crear datos para la gráfica principal (línea)
            const mainTrace = {
              x: filteredPos,
              y: filteredVal,
            type: 'scatter',
              mode: 'lines',
              name: 'Mean Variability',
              line: {color: '#667eea', width: 2},
              text: hoverText,
              hoverinfo: 'text+y'
            };

            // Crear datos para puntos verdes en posiciones seleccionadas
            // Obtener posiciones seleccionadas agregadas (pos -> {count, patients:Set})
            console.log('Calling getSelectedPositionsFromBars...');
            const selectedPositions = getSelectedPositionsFromBars();
            console.log('getSelectedPositionsFromBars returned:', selectedPositions);
            console.log('selectedPositions type:', typeof selectedPositions);
            console.log('selectedPositions size:', selectedPositions ? selectedPositions.size : 'undefined');
            
            const minPatients = getMinPatientsThreshold();
            console.log('minPatients threshold:', minPatients);
            

            console.log('Selected positions found:', selectedPositions);
            const selectedPos = [];
            const selectedVal = [];
            const selectedHoverText = [];
            const selectedColors = [];
            
            let matchesFound = 0;
            for (let i = 0; i < filteredPos.length; i++) {
              const nt = filteredPos[i];
              if (!selectedPositions.has(nt)) {
                continue;
              }
              matchesFound++;

              // Obtener información de la posición
              const info = selectedPositions.get(nt);
              const patientsAll = Array.from(info.patients || []);
              // Filtrar por pacientes actualmente seleccionados
              const relevantPatients = (selected.length === 0)
                ? patientsAll
                : patientsAll.filter(p => selected.includes(p));
              const count = relevantPatients.length;
              // Aplicar umbral mínimo de pacientes
              if (count < minPatients) {
                continue;
              }

              selectedPos.push(nt);
              selectedVal.push(filteredVal[i]);
              selectedHoverText.push(
                hoverText[i] +
                `<br><b>Patients (${count})</b>: ${relevantPatients.join(', ')}`
              );
              selectedColors.push(getColorForPatientCount(count));
            }
            console.log('Selected points created:', selectedPos.length);

            const traces = [mainTrace];   // línea primero
            if (selectedPos.length) {
              traces.push({
                x: selectedPos,
                y: selectedVal,
                type: 'scatter',
                mode: 'markers',
                name: 'Selected',
                marker: {
                  color: selectedColors,
                  size: 12,
                  symbol: 'circle',
                  line: { color: 'black', width: 1 }
                },
                text: selectedHoverText,
                hoverinfo: 'text+y',
                cliponaxis: false        // evita recorte en bordes
              });
            }

            // Genome ribbon (genes/proteins) bajo la línea - cambiar según el modo
            const segs = (mode === 'gene')
              ? (Array.isArray(SEGMENTS.genes) ? SEGMENTS.genes : [])
              : (Array.isArray(SEGMENTS.proteins) ? SEGMENTS.proteins : []);
            const colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'];
            const name2color = {};
            segs.forEach((s, i)=>{ name2color[s.name] = colors[i % colors.length]; });

            const yMin = Math.min(...filteredVal);
            const yMax = Math.max(...filteredVal);
            const trackH = (yMax - yMin) * 0.08;
            const shapes = segs.map((s, i)=>({
              type: 'rect', xref: 'x', yref: 'y',
              x0: s.start, x1: s.end, y0: yMin - trackH, y1: yMin - trackH/2,
              fillcolor: name2color[s.name], line: {width: 0}, opacity: 0.8
            }));
            const legendTraces = segs.map((s,i)=>({
              x:[null], y:[null], mode:'lines', name: s.name, line:{color:name2color[s.name], width:8}, showlegend:true
            }));

            // Calcular rango de zoom si hay elementos seleccionados específicos
            let xRange = undefined;
            if (selectedItems.length > 0 && !isAllSelected && filteredPos.length > 0) {
              const minPos = Math.min(...filteredPos);
              const maxPos = Math.max(...filteredPos);
              const range = maxPos - minPos;
              // Añadir margen del 10% en cada lado
              const margin = range * 0.1;
              xRange = [minPos - margin, maxPos + margin];
            }
            
            Plotly.newPlot(varPlotDiv, [...traces, ...legendTraces], {
              title: "Per-site temporal variability (mean across patients)",
              xaxis: {
                title: "Genomic position",
                range: xRange
              },
              yaxis: {
                title: "Variability (TVD)", 
                range: [yMin - trackH*1.2, yMax],
                tickformat: '.3f'
              },
              margin: {l: 60, r: 20, t: 50, b: 80},
              height: 350,
              plot_bgcolor: 'rgba(255,255,255,0.8)',
              paper_bgcolor: 'rgba(255,255,255,0)',
              showlegend: true,
              shapes: shapes
            }).then(() => {
              // Agregar evento de click para los puntos seleccionados usando la forma correcta
              varPlotDiv.on('plotly_click', ev => {
                const point = ev.points[0];
                console.log('=== CLICK DETECTED ===');
                console.log('Point curveNumber:', point.curveNumber);
                console.log('Point data:', point);
                
                // Procesar solo puntos seleccionados (trace index 1)
                if (point.curveNumber === 1) {
                  const position = point.x;
                  console.log('=== CLICK ON VARIABILITY POINT ===');
                  console.log('Position clicked:', position);
                  console.log('Current mode:', mode);
                  console.log('Point data:', point);
                  

                  
                  // Buscar en los datos de barras del modo actual (como en las barras de la Figura 4)
                  const barsData = DATA[mode].bars;
                  let foundData = null;
                  let foundCategory = null;
                  let foundSiteLabel = null;
                  
                  // Buscar en todos los genes/proteínas
                  for (const [category, barData] of Object.entries(barsData)) {
                    if (barData.data && barData.data.length > 0) {
                      for (const trace of barData.data) {
                        if (trace.x && trace.customdata) {
                          for (let i = 0; i < trace.x.length; i++) {
                            const siteLabel = trace.x[i];
                            
                            // Extraer posición nucleotídica del formato "aa / nt"
                            const ntMatch = siteLabel.match(/\/(\d+)$/);
                            if (ntMatch) {
                              const ntPos = parseInt(ntMatch[1]);
                              if (ntPos === position) {
                                const customData = JSON.parse(trace.customdata[i] || '{}');
                                if (Object.keys(customData).length > 0) {
                                  foundData = customData;
                                  foundCategory = category;
                                  foundSiteLabel = siteLabel;
                                  break;
                                }
                              }
                            }
                          }
                          if (foundData) break;
                        }
                      }
                      if (foundData) break;
                    }
                  }
                  
                  // SOLUCIÓN DIRECTA: Crear tabla manualmente
                  const selectedPositionsWithFlags = JSON.parse(document.getElementById('selected-positions-flags').textContent);
                  console.log('selectedPositionsWithFlags:', selectedPositionsWithFlags);
                  
                  // Siempre mostrar tabla para todas las posiciones
                  console.log('Creating table for position:', position);
                  {
                    
                    // Obtener información del gen/proteína para el título (con sitio aminoacídico)
                    let titleInfo = `position ${position}`;
                    const barsData = DATA[mode].bars;
                    if (barsData) {
                      for (const [category, barData] of Object.entries(barsData)) {
                        if (barData.data && barData.data.length > 0) {
                          // Buscar en los nombres completos fusionados (keep) en lugar de trace.x
                          const fullNames = barData.keep || [];
                          for (const fullName of fullNames) {
                            const positionFromName = fullName.split(' / ')[1];
                            if (positionFromName === position.toString()) {
                              // Extraer el sitio aminoacídico (índice) del cambio aminoacídico
                              const aaPart = fullName.split(' / ')[0]; // ej: "L9Q", "I124T", "K484"
                              let aaSite = '';
                              
                              // Extraer el número del sitio aminoacídico
                              const aaMatch = aaPart.match(/(\d+)/);
                              if (aaMatch) {
                                aaSite = aaMatch[1];
                                titleInfo = `${category} - position ${position} (aa ${aaSite})`;
                              } else {
                                titleInfo = `${category} - position ${position}`;
                              }
                              break;
                            }
                          }
                        }
                      }
                    }
                    
                    // Crear tabla HTML con la nueva lógica mejorada
                    const rows = [];
                    
                    // Obtener todos los pacientes que tienen datos para esta posición
                    const allPatients = Object.keys(selectedPositionsWithFlags);
                    const patients = [];
                    allPatients.forEach(patient => {
                      // Verificar si este paciente tiene datos para esta posición
                      const patientData = selectedPositionsWithFlags[patient];
                      if (patientData && patientData[position]) {
                        patients.push(patient);
                      }
                    });
                    
                    // Si no hay pacientes específicos para esta posición, usar todos los pacientes disponibles
                    if (patients.length === 0) {
                      patients.push(...allPatients);
                    }
                    
                    patients.forEach(patient => {
                      // Obtener el directorio de datos
                      const dataDir = JSON.parse(document.getElementById('data-directory').textContent);
                      const plotPath = `${dataDir}/${patient}_plots/pos${position.toString().padStart(5, '0')}_ALL.pdf`;
                      const genomePath = `${dataDir}/${patient}_genome_selection.pdf`;
                      
                      // Usar directamente los datos de posiciones seleccionadas con flags
                      let positionInfo = position.toString();
                      let betaInfo = 'N/A';
                      let pvalInfo = 'N/A';
                      let aaChange = '';
                      
                      const selectedPositionsWithFlags = JSON.parse(document.getElementById('selected-positions-flags').textContent);
                      const patientSpecificData = selectedPositionsWithFlags[patient] && selectedPositionsWithFlags[patient][position];
                      
                      // Inicializar filteredPositionData en el scope correcto
                      let filteredPositionData = [];
                      
                      if (patientSpecificData) {
                        // Verificar si hay múltiples entradas para la misma posición
                        const allPatientData = selectedPositionsWithFlags[patient];
                        let positionData = allPatientData[position];
                        
                        // Si positionData no es un array, convertirlo a array para manejo uniforme
                        if (positionData && !Array.isArray(positionData)) {
                          positionData = [positionData];
                        }
                        
                        if (positionData && positionData.length > 0) {
                          filteredPositionData = positionData;
                          
                          // Separar betas positivos y negativos
                          const positiveEntries = filteredPositionData.filter(entry => entry.beta && parseFloat(entry.beta) > 0);
                          const negativeEntries = filteredPositionData.filter(entry => entry.beta && parseFloat(entry.beta) < 0);
                          
                          // Aplicar lógica de combinación de cambios aminoacídicos
                          let combinedAaChange = '';
                          let combinedBetaInfo = 'N/A';
                          let combinedPvalInfo = 'N/A';
                          
                          if (filteredPositionData.length === 1) {
                            // Caso simple: una sola entrada
                            const entry = filteredPositionData[0];
                            combinedAaChange = entry.aa_change || '';
                            combinedBetaInfo = entry.beta ? parseFloat(entry.beta).toFixed(4) : 'N/A';
                            combinedPvalInfo = entry.pval ? parseFloat(entry.pval).toFixed(4) : 'N/A';
                          } else if (positiveEntries.length > 0 && negativeEntries.length > 0) {
                            // Caso especial: hay tanto beta positivo como negativo
                            const posEntry = positiveEntries[0];
                            const negEntry = negativeEntries[0];
                            
                            // Combinar cambios aminoacídicos: aaneg+posaa+aapos
                            const posAa = posEntry.aa_change || '';
                            const negAa = negEntry.aa_change || '';
                            
                            if (posAa && negAa) {
                              // Extraer el aminoácido de la posición negativa (formato: I124 -> I)
                              const negAaMatch = negAa.match(/^([A-Z])\d+/);
                              const posAaMatch = posAa.match(/^([A-Z])\d+([A-Z])$/);
                              
                              if (negAaMatch && posAaMatch) {
                                const negAaChar = negAaMatch[1];
                                const posAaChar = posAaMatch[1];
                                const finalAaChar = posAaMatch[2];
                                const posNum = posAa.match(/\d+/)[0];
                                
                                // Formato: I124T (aaneg+posaa+aapos)
                                combinedAaChange = `${negAaChar}${posNum}${finalAaChar}`;
                              } else {
                                // Fallback: combinar directamente
                                combinedAaChange = `${negAa}+${posAa}`;
                              }
                            } else {
                              combinedAaChange = posAa || negAa || '';
                            }
                            
                            // Combinar beta y pval: mostrar ambos en líneas separadas
                            const posBeta = parseFloat(posEntry.beta).toFixed(4);
                            const negBeta = parseFloat(negEntry.beta).toFixed(4);
                            const posPval = parseFloat(posEntry.pval).toFixed(4);
                            const negPval = parseFloat(negEntry.pval).toFixed(4);
                            
                            combinedBetaInfo = `${negBeta}<br>${posBeta}`;
                            combinedPvalInfo = `${negPval}<br>${posPval}`;
                            
                            console.log(`[DEBUG VARIABILITY] Combined positive and negative betas for ${patient} at ${position}: ${combinedAaChange}, beta=${combinedBetaInfo}, pval=${combinedPvalInfo}`);
                          } else {
                            // Otros casos con múltiples entradas
                            const aaChanges = filteredPositionData.map(entry => entry.aa_change).filter(change => change);
                            const betas = filteredPositionData.map(entry => entry.beta).filter(beta => beta);
                            const pvals = filteredPositionData.map(entry => entry.pval).filter(pval => pval);
                            
                            if (aaChanges.length > 1) {
                              const uniqueChanges = [...new Set(aaChanges)];
                              if (uniqueChanges.length === 2) {
                                const [first, second] = uniqueChanges.sort();
                                
                                // Caso 1: Q3729 + 3729K = Q3729K
                                if (first.match(/^[A-Z]\d+$/) && second.match(/^\d+[A-Z]$/)) {
                                  const firstMatch = first.match(/^([A-Z])(\d+)$/);
                                  const secondMatch = second.match(/^(\d+)([A-Z])$/);
                                  
                                  if (firstMatch && secondMatch && firstMatch[2] === secondMatch[1]) {
                                    // Los números coinciden, combinar: Q + 3729 + K = Q3729K
                                    combinedAaChange = `${firstMatch[1]}${firstMatch[2]}${secondMatch[2]}`;
                                  } else {
                                    // Los números no coinciden, usar formato original
                                    combinedAaChange = `${first}${second}`;
                                  }
                                } else if (first.match(/^\d+[A-Z]$/) && second.match(/^\d+[A-Z]$/)) {
                                  // Formato: 484E + 484K = amb484E
                                  combinedAaChange = `amb${first}`;
                                } else {
                                  // Otros casos
                                  combinedAaChange = uniqueChanges.join('+');
                                }
                              } else {
                                combinedAaChange = uniqueChanges.join('+');
                              }
                            } else {
                              combinedAaChange = aaChanges[0] || '';
                            }
                            
                            // Para betas y pvals, mostrar el primero o combinar si hay múltiples
                            if (betas.length > 1) {
                              combinedBetaInfo = betas.map(b => parseFloat(b).toFixed(4)).join(' / ');
                            } else {
                              combinedBetaInfo = betas[0] || 'N/A';
                            }
                            
                            if (pvals.length > 1) {
                              combinedPvalInfo = pvals.map(p => parseFloat(p).toFixed(4)).join(' / ');
                            } else {
                              combinedPvalInfo = pvals[0] || 'N/A';
                            }
                          }
                          
                          positionInfo = `${combinedAaChange}/${position}`;
                          aaChange = combinedAaChange;
                          betaInfo = combinedBetaInfo || 'N/A';
                          pvalInfo = combinedPvalInfo || 'N/A';
                          
                          console.log(`[DEBUG VARIABILITY] Found patient-specific data: ${patient} -> ${positionInfo}, beta=${betaInfo}, pval=${pvalInfo}`);
                        } else {
                          console.log(`[DEBUG VARIABILITY] No patient-specific data found for ${patient} at position ${position}`);
                        }
                      } else {
                        console.log(`[DEBUG VARIABILITY] No patient-specific data found for ${patient} at position ${position}`);
                      }
                      
                      // Determinar color según el beta
                      let rowColor = '#ffffff'; // blanco por defecto (para β mixtos)
                      if (betaInfo && betaInfo !== 'N/A' && typeof betaInfo === 'string') {
                        // Manejar diferentes formatos de beta
                        let hasPositive = false;
                        let hasNegative = false;
                        
                        if (betaInfo.includes('<br>')) {
                          // Formato: "neg<br>pos" - β mixtos
                          const parts = betaInfo.split('<br>');
                          for (const part of parts) {
                            const num = parseFloat(part.trim());
                            if (!isNaN(num)) {
                              if (num > 0) hasPositive = true;
                              if (num < 0) hasNegative = true;
                            }
                          }
                        } else if (betaInfo.includes(' / ')) {
                          // Formato: "neg / pos" - β mixtos
                          const parts = betaInfo.split(' / ');
                          for (const part of parts) {
                            const num = parseFloat(part.trim());
                            if (!isNaN(num)) {
                              if (num > 0) hasPositive = true;
                              if (num < 0) hasNegative = true;
                            }
                          }
                        } else if (betaInfo.includes('(') && betaInfo.includes(')')) {
                          // Formato: "neg (pos)" - usar el negativo para el color
                          const negBetaMatch = betaInfo.match(/^([-\d.]+)/);
                          const num = negBetaMatch ? parseFloat(negBetaMatch[1]) : parseFloat(betaInfo);
                          if (!isNaN(num)) {
                            if (num > 0) hasPositive = true;
                            if (num < 0) hasNegative = true;
                          }
                        } else {
                          // Formato simple: un solo número
                          const num = parseFloat(betaInfo);
                          if (!isNaN(num)) {
                            if (num > 0) hasPositive = true;
                            if (num < 0) hasNegative = true;
                          }
                        }
                        
                        // Asignar color según el tipo de β
                        if (hasPositive && hasNegative) {
                          rowColor = '#fff3e0'; // naranja claro para β mixtos
                        } else if (hasPositive) {
                          rowColor = '#ffebee'; // rojo claro para β positivo
                        } else if (hasNegative) {
                          rowColor = '#e3f2fd'; // azul claro para β negativo
                        }
                      }
                      
                      // No redondear beta y pval ya que se procesaron anteriormente
                      const betaRounded = betaInfo;
                      const pvalRounded = pvalInfo;
                      
                      // Obtener datos de covariación - manejar múltiples variantes
                      let covaryingInfo = '✗ No';
                      if (patientSpecificData) {
                        const allPatientData = selectedPositionsWithFlags[patient];
                        let positionData = allPatientData[position];
                        
                        // Si positionData no es un array, convertirlo a array para manejo uniforme
                        if (positionData && !Array.isArray(positionData)) {
                          positionData = [positionData];
                        }
                        
                        if (filteredPositionData && filteredPositionData.length > 0) {
                          // Manejar múltiples variantes para covariación, igual que para betas y pvals
                          if (filteredPositionData.length === 1) {
                            // Caso simple: una sola variante
                            const entry = filteredPositionData[0];
                            if (entry.covarying === true || entry.covarying === 'True') {
                              const aaChange = entry.aa_change || '';
                              covaryingInfo = `<span style="color: #e74c3c; font-weight: bold; cursor: pointer;" onclick="showCovaryingDetails('${position}', '${entry.covarying_with || ''}', 0, '${aaChange}')">✓ Yes (${entry.covarying_count || 0})</span>`;
                            }
                          } else {
                            // Múltiples variantes: mostrar covariación para cada una
                            const covaryingEntries = [];
                            
                            for (let i = 0; i < filteredPositionData.length; i++) {
                              const entry = filteredPositionData[i];
                              if (entry.covarying === true || entry.covarying === 'True') {
                                const aaChange = entry.aa_change || '';
                                covaryingEntries.push(`<span style="color: #e74c3c; font-weight: bold; cursor: pointer;" onclick="showCovaryingDetails('${position}', '${entry.covarying_with || ''}', ${i}, '${aaChange}')">✓ Yes (${entry.covarying_count || 0})</span>`);
                              } else {
                                covaryingEntries.push('✗ No');
                              }
                            }
                            
                            // Combinar las covariaciones, similar a como se hace con betas y pvals
                            if (covaryingEntries.length > 1) {
                              covaryingInfo = covaryingEntries.join('<br>');
                            } else {
                              covaryingInfo = covaryingEntries[0] || '✗ No';
                            }
                          }
                        }
                      }
                      
                      // Obtener información de fecha
                      let dateInfo = 'N/A';
                      if (patientSpecificData) {
                        const allPatientData = selectedPositionsWithFlags[patient];
                        let positionData = allPatientData[position];
                        
                        // Si positionData no es un array, convertirlo a array para manejo uniforme
                        if (positionData && !Array.isArray(positionData)) {
                          positionData = [positionData];
                        }
                        
                        if (filteredPositionData && filteredPositionData.length > 0) {
                          // Manejar múltiples variantes para fecha, igual que para betas y pvals
                          if (filteredPositionData.length === 1) {
                            // Caso simple: una sola variante
                            const entry = filteredPositionData[0];
                            dateInfo = entry.date || 'N/A';
                          } else {
                            // Múltiples variantes: mostrar fechas para cada una
                            const dates = filteredPositionData.map(entry => entry.date || 'N/A');
                            if (dates.length > 1) {
                              dateInfo = dates.join('<br>');
                            } else {
                              dateInfo = dates[0] || 'N/A';
                            }
                          }
                        }
                      }
                      
                      rows.push(`
                        <tr style="background-color: ${rowColor};">
                          <td style="border: 1px solid black; padding: 8px;"><strong>${patient}</strong></td>
                          <td style="border: 1px solid black; padding: 8px;"><strong>${positionInfo}</strong></td>
                          <td style="border: 1px solid black; padding: 8px;"><strong>${betaRounded}</strong></td>
                          <td style="border: 1px solid black; padding: 8px;"><strong>${pvalRounded}</strong></td>
                          <td style="border: 1px solid black; padding: 8px;">${covaryingInfo}</td>
                          <td style="border: 1px solid black; padding: 8px;">${dateInfo}</td>
                          <td style="border: 1px solid black; padding: 8px;"><a href="${plotPath}" target="_blank">📊 View</a></td>
                          <td style="border: 1px solid black; padding: 8px;"><a href="${genomePath}" target="_blank">🧬 View</a></td>
                        </tr>
                      `);
                    });
                    
                    let tableHTML = `
                      <div class="popup-overlay" onclick="closePopup()">
                        <div class="popup-content" onclick="event.stopPropagation()">
                          <div class="popup-header">
                            <h3>PDFs for ${titleInfo}</h3>
                            <button onclick="closePopup()" class="close-btn">×</button>
                          </div>
                          <div class="popup-body">
                            <table class="pdf-table" style="width: 100%; min-width: 600px; border-collapse: collapse; border: 2px solid black;">
                              <thead>
                                <tr>
                                  <th style="border: 1px solid black; padding: 8px;">Patient</th>
                                  <th style="border: 1px solid black; padding: 8px;">Position</th>
                                  <th style="border: 1px solid black; padding: 8px;">β</th>
                                  <th style="border: 1px solid black; padding: 8px;">p-value</th>
                                  <th style="border: 1px solid black; padding: 8px;">Covarying</th>
                                  <th style="border: 1px solid black; padding: 8px;">Date</th>
                                  <th style="border: 1px solid black; padding: 8px;">Position PDF</th>
                                  <th style="border: 1px solid black; padding: 8px;">Genome PDF</th>
                                </tr>
                              </thead>
                              <tbody>${rows.join('')}</tbody>
                            </table>
                          </div>
                        </div>
                      </div>`;
                    
                    // Limpiar popups anteriores antes de mostrar el nuevo
                    const existingPopups = document.querySelectorAll('.popup-overlay');
                    existingPopups.forEach(popup => popup.remove());
                    
                    // Insertar tabla en el DOM
                    document.body.insertAdjacentHTML('beforeend', tableHTML);
                    console.log('Table inserted successfully');
                  }
                }
              });
            });
          } catch (error) {
            console.error('Error drawing variability plot:', error);
            varPlotDiv.innerHTML = '<div class="loading">Error loading variability data.</div>';
          }
        }

        // Función para obtener las posiciones seleccionadas desde los datos automáticos
        function getSelectedPositionsFromBars() {
          const selectedPositions = new Map(); // pos -> {count: number, patients: Set}
          
          // Usar los datos automáticos leídos desde los archivos *_selected_ranked.tsv
          for (const [pos, data] of Object.entries(selectedPositionsData)) {
            const posInt = parseInt(pos);
            selectedPositions.set(posInt, {
              count: data.count,
              patients: new Set(data.patients)
            });
          }
          
          return selectedPositions;
        }

        // Función para obtener el color según el número de pacientes
        function getColorForPatientCount(count) {
          // Escala de colores: verde (pocos pacientes) a rojo (muchos pacientes)
          if (count <= 1) return '#00FF00';      // Verde brillante
          if (count <= 2) return '#FFFF00';      // Amarillo
          if (count <= 3) return '#FFA500';      // Naranja
          if (count <= 5) return '#FF8C00';      // Naranja oscuro
          if (count <= 8) return '#FF4500';      // Rojo naranja
          return '#FF0000';                      // Rojo brillante
        }

        // Función simple para obtener datos de PDFs para una posición específica
        function getPDFDataForPosition(position) {
          console.log('=== SEARCHING FOR POSITION ===');
          console.log('Position to search:', position);
          console.log('Current mode:', mode);
          
          // Buscar en los datos de barras del modo actual
          const barsData = DATA[mode].bars;
          console.log('Bars data keys:', Object.keys(barsData));
          
          for (const [category, barData] of Object.entries(barsData)) {
            if (barData.data && barData.data.length > 0) {
              for (const trace of barData.data) {
                if (trace.x && trace.customdata) {
                  for (let i = 0; i < trace.x.length; i++) {
                    const siteLabel = trace.x[i];
                    
                    // Extraer posición nucleotídica del formato "aa / nt"
                    const ntMatch = siteLabel.match(/\/(\d+)$/);
                    if (ntMatch) {
                      const ntPos = parseInt(ntMatch[1]);
                      if (ntPos === position) {
                        const customData = JSON.parse(trace.customdata[i] || '{}');
                        if (Object.keys(customData).length > 0) {
                          return {
                            category: category,
                            siteLabel: siteLabel,
                            data: customData
                          };
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          return null;
        }

        // Función para obtener la posición aminoacídica
        function getAminoAcidPosition(pos) {
          const result = positionMapping.aa_positions[pos] || null;
          console.log(`getAminoAcidPosition(${pos}) = ${result}`);
          return result;
        }

        // Función para obtener gen en posición
        function getGeneAtPosition(pos) {
          const posKey = typeof pos === 'string' ? parseInt(pos, 10) : pos;
          return positionMapping.genes[posKey] || positionMapping.genes[pos] || positionMapping.genes[String(pos)] || 'Unknown';
        }

        // Función para obtener proteína en posición
        function getProteinAtPosition(pos) {
          const posKey = typeof pos === 'string' ? parseInt(pos, 10) : pos;
          return positionMapping.proteins[posKey] || positionMapping.proteins[pos] || positionMapping.proteins[String(pos)] || 'Unknown';
        }
        
        // ───────── FILTROS AVANZADOS (igual que en Recombinación) ─────────
        let filterState = {
          selections: {}, // {column: {value: true/false}}
          includeNoData: {}, // {column: true/false} - true por defecto
          includePatientsNotInFilterData: false, // Si true, incluye pacientes que no están en filterData
          patientSelectionMode: 'union' // strict | any | union | coverage
        };
        let filtersApplied = false;
        
        // Cargar datos de filtros
        let filterData = null;
        try {
          const filterDataEl = document.getElementById('filter-data');
          if (filterDataEl) {
            filterData = JSON.parse(filterDataEl.textContent);
          }
        } catch(e) {
          console.error('Error loading filter data:', e);
          filterData = {patient_data: {}, patient_rows: {}, columns: [], unique_values: {}, patient_name_mapping: {}};
        }
        
        // Inicializar todos los valores como seleccionados por defecto
        if (filterData && filterData.columns) {
          for (const col of filterData.columns) {
            filterState.selections[col] = {};
            filterState.includeNoData[col] = true;
            const uniqueVals = filterData.unique_values[col] || [];
            for (const val of uniqueVals) {
              filterState.selections[col][val] = true;
            }
          }
        }
        
        const moreFiltersBtn = document.getElementById('more-filters-btn');
        const moreFiltersModal = document.getElementById('more-filters-modal');
        const closeMoreFiltersBtn = document.getElementById('close-more-filters-modal');
        const applyFiltersModalBtn = document.getElementById('apply-filters-btn');
        const resetAdvancedFiltersBtn = document.getElementById('reset-advanced-filters-btn');
        const filtersTableContainer = document.getElementById('filters-table-container');
        
        function resetAdvancedFiltersState(updateUI) {
          if (!filterData || !filterData.columns) {
            return;
          }
          filterState.selections = {};
          filterState.includeNoData = {};
          for (const col of filterData.columns) {
            filterState.selections[col] = {};
            const uniqueVals = filterData.unique_values[col] || [];
            for (const val of uniqueVals) {
              filterState.selections[col][val] = true;
            }
            filterState.includeNoData[col] = true;
          }
          filterState.includePatientsNotInFilterData = false;
          filterState.patientSelectionMode = 'union';

          if (updateUI) {
            const strictRadio = document.getElementById('patient-mode-strict-sel');
            const anyRadio = document.getElementById('patient-mode-any-sel');
            const unionRadio = document.getElementById('patient-mode-union-sel');
            const coverageRadio = document.getElementById('patient-mode-coverage-sel');
            if (strictRadio) strictRadio.checked = false;
            if (anyRadio) anyRadio.checked = false;
            if (coverageRadio) coverageRadio.checked = false;
            if (unionRadio) unionRadio.checked = true;
            const includeNotInTsv = document.getElementById('patient-include-not-in-tsv-sel');
            if (includeNotInTsv) includeNotInTsv.checked = false;
            document.querySelectorAll('input[type="checkbox"][data-column]').forEach(cb => {
              cb.checked = true;
            });
          }
        }

        // Función para generar la tabla de filtros
        function generateFiltersTable() {
          if (!filterData || !filterData.columns || filterData.columns.length === 0) {
            filtersTableContainer.innerHTML = '<p>No filter data available.</p>';
            return;
          }
          
          // Inicializar includePatientsNotInFilterData si no existe
          if (filterState.includePatientsNotInFilterData === undefined) {
            filterState.includePatientsNotInFilterData = false;
          }
          if (!filterState.patientSelectionMode) {
            filterState.patientSelectionMode = 'union';
          }
          
          let html = '<div style="background-color: #f8f9fa; padding: 15px; margin-bottom: 20px; border: 1px solid #dee2e6; border-radius: 4px;">';
          html += '<label style="font-weight: bold; display: inline-block; margin-bottom: 10px;">Patient Selection Mode:</label>';
          html += '<button id="patient-selection-mode-info" style="margin-left: 6px; background: #3498db; color: white; border: none; border-radius: 50%; width: 16px; height: 16px; font-size: 10px; cursor: pointer; display: inline-flex; align-items: center; justify-content: center; font-weight: bold;" title="Click for information">ℹ</button>';
          html += '<div style="margin-left: 10px;">';
          html += '<label style="display: block; margin-bottom: 8px; cursor: pointer;">';
          html += '<input type="radio" name="patient-selection-mode-sel" id="patient-mode-strict-sel" value="strict" ' + (filterState.patientSelectionMode === 'strict' ? 'checked' : '') + ' style="margin-right: 8px;">';
          html += 'All time points must match (strict)';
          html += '</label>';
          html += '<label style="display: block; margin-bottom: 8px; cursor: pointer;">';
          html += '<input type="radio" name="patient-selection-mode-sel" id="patient-mode-any-sel" value="any" ' + (filterState.patientSelectionMode === 'any' ? 'checked' : '') + ' style="margin-right: 8px;">';
          html += 'Any time point matches (row-level match)';
          html += '</label>';
          html += '<label style="display: block; cursor: pointer;">';
          html += '<input type="radio" name="patient-selection-mode-sel" id="patient-mode-union-sel" value="union" ' + (filterState.patientSelectionMode === 'union' ? 'checked' : '') + ' style="margin-right: 8px;">';
          html += 'Patient-level aggregation (union across time points)';
          html += '</label>';
          html += '<label style="display: block; margin-top: 6px; cursor: pointer;">';
          html += '<input type="radio" name="patient-selection-mode-sel" id="patient-mode-coverage-sel" value="coverage" ' + (filterState.patientSelectionMode === 'coverage' ? 'checked' : '') + ' style="margin-right: 8px;">';
          html += 'All time points covered by filters (coverage mode)';
          html += '</label>';
          html += '</div>';
          html += '<div style="margin-top: 10px; margin-left: 10px;">';
          html += '<label style="display: block; cursor: pointer;">';
          html += '<input type="checkbox" id="patient-include-not-in-tsv-sel" ' + (filterState.includePatientsNotInFilterData ? 'checked' : '') + ' style="margin-right: 8px;">';
          html += 'Include patients not present in the filter TSV';
          html += '</label>';
          html += '</div>';
          html += '</div>';
          
          html += '<table style="width: 100%; border-collapse: collapse; margin-top: 20px;">';
          html += '<thead><tr style="background-color: #f8f9fa; border-bottom: 2px solid #dee2e6;">';
          html += '<th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Column</th>';
          html += '<th style="padding: 12px; text-align: left; border: 1px solid #dee2e6;">Values</th>';
          html += '<th style="padding: 12px; text-align: center; border: 1px solid #dee2e6;">Include No Data';
          html += '<button id="include-no-data-info" style="margin-left: 6px; background: #3498db; color: white; border: none; border-radius: 50%; width: 16px; height: 16px; font-size: 10px; cursor: pointer; display: inline-flex; align-items: center; justify-content: center; font-weight: bold;" title="Click for information">ℹ</button>';
          html += '</th>';
          html += '</tr></thead><tbody>';
          
          for (const col of filterData.columns) {
            html += '<tr style="border-bottom: 1px solid #dee2e6;">';
            html += '<td style="padding: 12px; font-weight: bold; vertical-align: top; border: 1px solid #dee2e6;">' + col + '</td>';
            html += '<td style="padding: 12px; border: 1px solid #dee2e6;">';
            
            // Botones de selección rápida
            html += '<div style="margin-bottom: 10px;">';
            html += '<button class="select-all-btn" data-column="' + col + '" style="background: #3498db; color: white; border: none; border-radius: 4px; padding: 5px 10px; cursor: pointer; margin-right: 5px; font-size: 12px;">Select All</button>';
            html += '<button class="deselect-all-btn" data-column="' + col + '" style="background: #95a5a6; color: white; border: none; border-radius: 4px; padding: 5px 10px; cursor: pointer; font-size: 12px;">Deselect All</button>';
            html += '</div>';
            
            // Valores únicos de esta columna
            const uniqueVals = filterData.unique_values[col] || [];
            html += '<div style="max-height: 200px; overflow-y: auto; border: 1px solid #dee2e6; padding: 10px; border-radius: 4px;">';
            
            // Inicializar estado si no existe - TODOS SELECCIONADOS POR DEFECTO
            if (!filterState.selections[col]) {
              filterState.selections[col] = {};
              for (const val of uniqueVals) {
                filterState.selections[col][val] = true;
              }
            }
            if (filterState.includeNoData[col] === undefined) {
              filterState.includeNoData[col] = true;
            }
            
            for (const val of uniqueVals) {
              const isChecked = filterState.selections[col][val] !== false;
              const checked = isChecked ? 'checked' : '';
              const checkboxId = ('filter-' + col + '-' + val).replace(/[^a-zA-Z0-9-]/g, '_');
              html += '<div style="margin: 5px 0;">';
              html += '<input type="checkbox" id="' + checkboxId + '" data-column="' + col + '" data-value="' + val.replace(/"/g, '&quot;') + '" ' + checked + ' style="margin-right: 8px;">';
              html += '<label for="' + checkboxId + '" style="cursor: pointer;">' + val + '</label>';
              html += '</div>';
            }
            
            html += '</div>';
            html += '</td>';
            
            // Checkbox para incluir pacientes sin datos
            const noDataChecked = filterState.includeNoData[col] !== false ? 'checked' : '';
            const noDataId = ('no-data-' + col).replace(/[^a-zA-Z0-9-]/g, '_');
            html += '<td style="padding: 12px; text-align: center; border: 1px solid #dee2e6;">';
            html += '<input type="checkbox" id="' + noDataId + '" data-column="' + col + '" ' + noDataChecked + ' style="cursor: pointer;">';
            html += '</td>';
            
            html += '</tr>';
          }
          
          html += '</tbody></table>';
          filtersTableContainer.innerHTML = html;
          
          // Añadir event listeners para los checkboxes
          document.querySelectorAll('input[type="checkbox"][data-column]').forEach(cb => {
            cb.addEventListener('change', function() {
              const col = this.getAttribute('data-column');
              const val = this.getAttribute('data-value');
              
              if (val) {
                if (!filterState.selections[col]) {
                  filterState.selections[col] = {};
                }
                filterState.selections[col][val] = this.checked;
              } else {
                filterState.includeNoData[col] = this.checked;
              }
            });
          });
          
          // Añadir event listeners para los botones de selección rápida
          document.querySelectorAll('.select-all-btn').forEach(btn => {
            btn.addEventListener('click', function() {
              const col = this.getAttribute('data-column');
              const checkboxes = document.querySelectorAll('input[type="checkbox"][data-column="' + col + '"][data-value]');
              checkboxes.forEach(cb => {
                cb.checked = true;
                const val = cb.getAttribute('data-value');
                if (!filterState.selections[col]) {
                  filterState.selections[col] = {};
                }
                filterState.selections[col][val] = true;
              });
            });
          });
          
          document.querySelectorAll('.deselect-all-btn').forEach(btn => {
            btn.addEventListener('click', function() {
              const col = this.getAttribute('data-column');
              const checkboxes = document.querySelectorAll('input[type="checkbox"][data-column="' + col + '"][data-value]');
              checkboxes.forEach(cb => {
                cb.checked = false;
                const val = cb.getAttribute('data-value');
                if (filterState.selections[col]) {
                  filterState.selections[col][val] = false;
                }
              });
            });
          });
          
          // Agregar event listeners para el modo de selección y la inclusión de pacientes fuera del TSV
          const strictModeRadio = document.getElementById('patient-mode-strict-sel');
          const anyModeRadio = document.getElementById('patient-mode-any-sel');
          const unionModeRadio = document.getElementById('patient-mode-union-sel');
          const coverageModeRadio = document.getElementById('patient-mode-coverage-sel');
          const includeNotInTsv = document.getElementById('patient-include-not-in-tsv-sel');
          if (strictModeRadio) {
            strictModeRadio.addEventListener('change', function() {
              if (this.checked) {
                filterState.patientSelectionMode = 'strict';
              }
            });
          }
          if (anyModeRadio) {
            anyModeRadio.addEventListener('change', function() {
              if (this.checked) {
                filterState.patientSelectionMode = 'any';
              }
            });
          }
          if (unionModeRadio) {
            unionModeRadio.addEventListener('change', function() {
              if (this.checked) {
                filterState.patientSelectionMode = 'union';
              }
            });
          }
          if (coverageModeRadio) {
            coverageModeRadio.addEventListener('change', function() {
              if (this.checked) {
                filterState.patientSelectionMode = 'coverage';
              }
            });
          }
          if (includeNotInTsv) {
            includeNotInTsv.addEventListener('change', function() {
              filterState.includePatientsNotInFilterData = this.checked;
            });
          }

          // Info buttons for Advanced Filters
          const patientModeInfoBtn = document.getElementById('patient-selection-mode-info');
          const includeNoDataInfoBtn = document.getElementById('include-no-data-info');
          const PATIENT_SELECTION_INFO_TEXT = `Patient Selection Mode help
- Mode A — All time points must match (strict): Selects a patient only if *every* time point satisfies the selected values for *each* active category. This is the most restrictive mode and is useful when you want a patient’s full longitudinal profile to belong to the same subgroup (e.g., always the same sampling compartment). Example: selecting Location = ETA keeps only patients whose samples are *always* ETA across all time points.
- Mode B — Any time point matches (row-level): Selects a patient if there is *at least one* time point where *all active categories* are satisfied simultaneously. Use this when you want to capture patients showing a pattern at a specific time (e.g., a treatment present at a particular visit together with a given sampling location). Example: Location = ETA selects patients who have ETA in at least one time point. With multiple categories, they must co-occur in the same time point.
- Mode C — Patient-level union across time points: Selects a patient if the patient’s history contains evidence for *each* active category at *some* time point, not necessarily the same one. This is useful when you want “ever observed” evidence across the follow-up. Example: Location = Blood AND Antiviral = X selects patients who had Blood at any time point and Antiviral X at any (possibly different) time point.
- Mode D — All time points covered by filters (coverage): Selects a patient only if *no time point contradicts* the selected values. Each time point must be compatible with the allowed values for every active category. This is useful when you want your selected filters to “explain” the full longitudinal trajectory. Example: if a patient has NP at one time point and ETA at another, selecting only ETA will exclude them unless NP is also allowed.`;

          const INCLUDE_NO_DATA_INFO_TEXT = `Include No Data help
- “Include No Data” is applied **per category**. When enabled for a category, missing values (empty/NA) in that category do not automatically exclude a patient/time point for that category. When disabled, missing values are treated as non-matching.
- This is different from “Include patients not present in the filter TSV”, which controls whether patients without any TSV metadata are kept.
- The impact depends on the selection mode: in strict/coverage modes, missing values can exclude time points unless “Include No Data” is enabled for that category.`;

          if (patientModeInfoBtn) {
            patientModeInfoBtn.addEventListener('click', function() {
              console.log('[INFO] Patient Selection Mode info clicked');
              openAdvancedFiltersInfo(PATIENT_SELECTION_INFO_TEXT);
            });
          }
          if (includeNoDataInfoBtn) {
            includeNoDataInfoBtn.addEventListener('click', function() {
              console.log('[INFO] Include No Data info clicked');
              openAdvancedFiltersInfo(INCLUDE_NO_DATA_INFO_TEXT);
            });
          }
        }

        // ===== UPGMA visibility filters (solo etiquetas del dendrograma) =====
        function buildUpgmaFilterModal() {
          const container = document.getElementById('upgma-filter-categories');
          if (!container) {
            return;
          }
          if (!filterData || !filterData.columns || filterData.columns.length === 0) {
            container.innerHTML = '<p>No filter data available.</p>';
            return;
          }
          let html = '';
          for (const col of filterData.columns) {
            const checked = upgmaFilterState.categories[col] === true ? 'checked' : '';
            const checkboxId = ('upgma-cat-' + col).replace(/[^a-zA-Z0-9-]/g, '_');
            html += '<div style="margin: 6px 0;">';
            html += '<input type="checkbox" id="' + checkboxId + '" data-upgma-category="' + col + '" ' + checked + ' style="margin-right: 8px;">';
            html += '<label for="' + checkboxId + '" style="cursor: pointer;">' + col + '</label>';
            html += '</div>';
          }
          container.innerHTML = html;
          const includeMissing = document.getElementById('upgma-include-missing');
          if (includeMissing) {
            includeMissing.checked = !!upgmaFilterState.includePatientsNotInFilterData;
          }
        }

        function applyUpgmaFilterModal() {
          upgmaFilterState.categories = {};
          document.querySelectorAll('input[type="checkbox"][data-upgma-category]').forEach(cb => {
            const col = cb.getAttribute('data-upgma-category');
            upgmaFilterState.categories[col] = cb.checked;
          });
          const includeMissing = document.getElementById('upgma-include-missing');
          upgmaFilterState.includePatientsNotInFilterData = includeMissing ? includeMissing.checked : false;
          drawUPGMA();
        }

        const upgmaFilterBtn = document.getElementById('upgma-filter-btn');
        const upgmaFilterModal = document.getElementById('upgma-filter-modal');
        const closeUpgmaFilterBtn = document.getElementById('close-upgma-filter-modal');
        const applyUpgmaFilterBtn = document.getElementById('apply-upgma-filter-btn');
        if (upgmaFilterBtn && upgmaFilterModal && closeUpgmaFilterBtn && applyUpgmaFilterBtn) {
          upgmaFilterBtn.onclick = function() {
            buildUpgmaFilterModal();
            upgmaFilterModal.style.display = 'block';
          };
          closeUpgmaFilterBtn.onclick = function() {
            upgmaFilterModal.style.display = 'none';
          };
          applyUpgmaFilterBtn.onclick = function() {
            applyUpgmaFilterModal();
            upgmaFilterModal.style.display = 'none';
          };
          window.addEventListener('click', function(event) {
            if (event.target === upgmaFilterModal) {
              upgmaFilterModal.style.display = 'none';
            }
          });
        }

        if (resetAdvancedFiltersBtn) {
          resetAdvancedFiltersBtn.onclick = function() {
            resetAdvancedFiltersState(true);
          };
        }

        function openAdvancedFiltersInfo(bodyText) {
          const overlay = document.getElementById('advanced-filters-info-overlay');
          const titleEl = document.getElementById('advanced-filters-info-title');
          const bodyEl = document.getElementById('advanced-filters-info-content');
          if (!overlay || !titleEl || !bodyEl) {
            return;
          }
          titleEl.textContent = '';
          titleEl.style.display = 'none';
          bodyEl.textContent = bodyText;
          overlay.style.display = 'block';
        }

        function closeAdvancedFiltersInfo() {
          const overlay = document.getElementById('advanced-filters-info-overlay');
          if (overlay) {
            overlay.style.display = 'none';
          }
        }

        const infoCloseBtn = document.getElementById('advanced-filters-info-close');
        const infoOverlay = document.getElementById('advanced-filters-info-overlay');
        if (infoCloseBtn) {
          infoCloseBtn.onclick = closeAdvancedFiltersInfo;
        }
        if (infoOverlay) {
          infoOverlay.addEventListener('click', function(e) {
            if (e.target === infoOverlay) {
              closeAdvancedFiltersInfo();
            }
          });
        }
        document.addEventListener('keydown', function(e) {
          if (e.key === 'Escape') {
            closeAdvancedFiltersInfo();
          }
        });
        
        // Función para actualizar el panel de resumen
        function updateFiltersSummary() {
          const summaryDiv = document.getElementById('filters-summary');
          const summaryContent = document.getElementById('filters-summary-content');
          const patientsListDiv = document.getElementById('selected-patients-list');
          const patientsCountSpan = document.getElementById('selected-patients-count');
          
          if (!summaryDiv || !summaryContent || !patientsListDiv || !patientsCountSpan) {
            return;
          }
          
          // Obtener filtros activos
          const activeFilters = [];
          if (filterData && filterData.columns) {
            for (const col of filterData.columns) {
              const colSelections = filterState.selections[col] || {};
              const selectedValues = Object.keys(colSelections).filter(k => colSelections[k] === true);
              const includeNoData = filterState.includeNoData[col] || false;
              
              if (selectedValues.length > 0 || includeNoData) {
                const filterParts = [];
                if (selectedValues.length > 0) {
                  filterParts.push(selectedValues.slice(0, 3).join(', ') + (selectedValues.length > 3 ? ' (+' + (selectedValues.length - 3) + ' more)' : ''));
                }
                if (includeNoData) {
                  filterParts.push('Include No Data');
                }
                activeFilters.push('<strong>' + col + ':</strong> ' + filterParts.join('; '));
              }
            }
          }
          
          // Obtener pacientes seleccionados (usar patientSel que está en el scope global)
          const selectedPatients = patientSel ? [...patientSel.selectedOptions].map(o => o.value) : [];
          
          // Mostrar el panel solo después de aplicar filtros
          if (filtersApplied && activeFilters.length > 0) {
            summaryDiv.style.display = 'block';
            
            // Actualizar lista de pacientes
            patientsCountSpan.textContent = selectedPatients.length;
            if (selectedPatients.length > 0) {
              if (selectedPatients.length <= 30) {
                patientsListDiv.innerHTML = selectedPatients.join(', ');
              } else {
                patientsListDiv.innerHTML = selectedPatients.slice(0, 30).join(', ') + 
                  ' <em style="opacity: 0.8;">(+' + (selectedPatients.length - 30) + ' more)</em>';
              }
            } else {
              patientsListDiv.innerHTML = '<em style="opacity: 0.7;">No patients match these filters</em>';
            }
            
            // Actualizar contenido de filtros
            summaryContent.innerHTML = activeFilters.join('<br>');
          } else {
            summaryDiv.style.display = 'none';
          }
        }
        
        // Función para aplicar los filtros (igual que en Recombinación, pero usa patientSel y applyFilters)
        function applyAdvancedFilters() {
          const DEBUG_PATIENT_SELECTION_MODE = false;
          if (!filterData || !filterData.patient_data) {
            console.warn('No filter data available');
            return;
          }
          
          const patientRowsMap = filterData.patient_rows || {};
          const matchingPatients = new Set();
          const allPatients = patientSel ? Array.from(patientSel.options).map(opt => opt.value) : [];
          const mode = filterState.patientSelectionMode || 'union';
          const includeNotInTsvEffective = filtersApplied ? filterState.includePatientsNotInFilterData : true;

          function findPatientKey(patientName) {
            if (patientName in filterData.patient_data) {
              return patientName;
            }
            const patientNameMapping = filterData.patient_name_mapping || {};
            if (patientNameMapping && typeof patientNameMapping === 'object') {
              for (const [key, variations] of Object.entries(patientNameMapping)) {
                if (Array.isArray(variations) && variations.includes(patientName)) {
                  for (const variation of variations) {
                    if (variation in filterData.patient_data) {
                      return variation;
                    }
                  }
                  return key;
                }
              }
            }
            return null;
          }

          function buildActiveFilters() {
            const active = [];
            for (const col of filterData.columns) {
              const colSelections = filterState.selections[col] || {};
              const selectedValues = Object.keys(colSelections).filter(k => colSelections[k] === true);
              const includeNoData = filterState.includeNoData[col] !== false;
              const uniqueVals = filterData.unique_values[col] || [];
              const allSelected = selectedValues.length === uniqueVals.length && includeNoData;
              active.push({ col, selectedValues, includeNoData, allSelected });
            }
            return active;
          }

          function hasAnyFilter(activeFilters) {
            for (const f of activeFilters) {
              if (f.selectedValues.length === 0 && !f.includeNoData) {
                return true;
              }
              if (!f.allSelected) {
                return true;
              }
            }
            return false;
          }

          function rowMatchesFilters(row, activeFilters) {
            for (const f of activeFilters) {
              const values = row[f.col] || [];
              const hasNoData = values.includes("__NO_DATA__") || values.length === 0;
              if (f.selectedValues.length === 0) {
                if (!f.includeNoData) {
                  return false;
                }
                continue;
              }
              if (hasNoData) {
                if (!f.includeNoData) {
                  return false;
                }
                continue;
              }
              if (!f.selectedValues.some(v => values.includes(v))) {
                return false;
              }
            }
            return true;
          }

          function evaluatePatientSelected(patientKey, activeFilters, mode) {
            const rows = patientRowsMap[patientKey] || [];
            if (rows.length === 0) {
              return false;
            }
            if (mode === 'strict') {
              return rows.every(row => rowMatchesFilters(row, activeFilters));
            }
            if (mode === 'any') {
              return rows.some(row => rowMatchesFilters(row, activeFilters));
            }
            if (mode === 'coverage') {
              return rows.every(row => rowMatchesFilters(row, activeFilters));
            }
            // union across time points
            const agg = {};
            for (const f of activeFilters) {
              agg[f.col] = { values: new Set(), hasNoData: false };
            }
            for (const row of rows) {
              for (const f of activeFilters) {
                const values = row[f.col] || [];
                if (values.length === 0 || values.includes("__NO_DATA__")) {
                  agg[f.col].hasNoData = true;
                }
                values.forEach(v => {
                  if (v !== "__NO_DATA__" && String(v).trim() !== "") {
                    agg[f.col].values.add(v);
                  }
                });
              }
            }
            for (const f of activeFilters) {
              if (f.selectedValues.length === 0) {
                if (!f.includeNoData) {
                  return false;
                }
                continue;
              }
              if (f.includeNoData && agg[f.col].hasNoData) {
                continue;
              }
              const intersects = f.selectedValues.some(v => agg[f.col].values.has(v));
              if (!intersects) {
                return false;
              }
            }
            return true;
          }

          const activeFilters = buildActiveFilters();
          const hasFilters = hasAnyFilter(activeFilters);
          
          if (!hasFilters) {
            // Todos los valores están seleccionados en todas las columnas = sin filtros efectivos
            if (includeNotInTsvEffective) {
              // Si está activado, incluir todos los pacientes
              allPatients.forEach(p => matchingPatients.add(p));
            } else {
              // Si no está activado, solo incluir pacientes que están en filterData
              for (const patient of allPatients) {
                const patientKey = findPatientKey(patient);
                if (patientKey !== null && patientKey in filterData.patient_data) {
                  matchingPatients.add(patient);
                }
              }
            }
          } else {
            for (const patient of allPatients) {
              const patientKey = findPatientKey(patient);
              const patientInFilterData = patientKey !== null && patientKey in filterData.patient_data;
              
              if (!patientInFilterData) {
                if (includeNotInTsvEffective) {
                  matchingPatients.add(patient);
                }
                continue;
              }
              const matches = evaluatePatientSelected(patientKey, activeFilters, mode);
              if (matches) {
                matchingPatients.add(patient);
              }
            }
          }

          if (DEBUG_PATIENT_SELECTION_MODE) {
            const tsvPatients = Object.keys(filterData.patient_data || {}).length;
            console.log('[FILTERS] Selection mode:', mode);
            console.log('[FILTERS] Patients in TSV:', tsvPatients);
            console.log('[FILTERS] Patients selected:', matchingPatients.size);
            const exampleCol = filterData.columns.includes('sample_location_condensed')
              ? 'sample_location_condensed'
              : (filterData.columns[0] || null);
            if (exampleCol) {
              const examplePatient = Object.keys(patientRowsMap).find(p => (patientRowsMap[p] || []).length > 1);
              if (examplePatient) {
                const rows = patientRowsMap[examplePatient] || [];
                const rowVals = rows.map(r => r[exampleCol] || []);
                console.log('[FILTERS] Example patient:', examplePatient, 'category:', exampleCol, 'rows:', rowVals);
                if (mode === 'coverage') {
                  const compatCount = rows.filter(r => rowMatchesFilters(r, activeFilters)).length;
                  console.log('[FILTERS] Coverage mode example:', examplePatient, 'compatible rows:', compatCount, '/', rows.length);
                }
              }
            }
          }
          
          // Actualizar el selector de pacientes
          if (patientSel) {
            for (let i = 0; i < patientSel.options.length; i++) {
              patientSel.options[i].selected = false;
            }
            
            const matchingArray = Array.from(matchingPatients);
            for (let i = 0; i < patientSel.options.length; i++) {
              const option = patientSel.options[i];
              if (matchingArray.includes(option.value)) {
                option.selected = true;
              }
            }
            
            console.log('Applied filters: ' + matchingArray.length + ' patients selected out of ' + patientSel.options.length);
          }
          
          // Cerrar la modal
          if (moreFiltersModal) {
            moreFiltersModal.style.display = 'none';
          }
          
          // Actualizar el panel de resumen ANTES de aplicar filtros
          updateFiltersSummary();
          
          // Aplicar los filtros (esto disparará el re-renderizado)
          // Usar setTimeout para asegurar que el DOM se actualice antes de applyFilters
          setTimeout(() => {
            applyFilters();
          }, 0);
        }
        
        // Event listeners para el modal de filtros
        const DEBUG_ADVANCED_FILTERS_APPLY = true;
        if (moreFiltersBtn && moreFiltersModal && closeMoreFiltersBtn && applyFiltersModalBtn) {
          moreFiltersBtn.onclick = function() {
            generateFiltersTable();
            moreFiltersModal.style.display = 'block';
          };
          
          closeMoreFiltersBtn.onclick = function() {
            moreFiltersModal.style.display = 'none';
          };
          
          applyFiltersModalBtn.onclick = function() {
            console.log('[ADV_FILTERS] Apply clicked');
            const beforeCount = patientSel ? patientSel.selectedOptions.length : 0;
            filtersApplied = true;
            applyAdvancedFilters();
            updateFiltersSummary();
            const afterCount = patientSel ? patientSel.selectedOptions.length : 0;
            if (DEBUG_ADVANCED_FILTERS_APPLY) {
              console.log('[ADV_FILTERS] selected patients before:', beforeCount);
              console.log('[ADV_FILTERS] selected patients after:', afterCount);
              console.log('[ADV_FILTERS] filtersApplied:', filtersApplied);
            }
          };
          
          // Close modal when clicking outside of it
          window.onclick = function(event) {
            if (event.target === moreFiltersModal) {
              moreFiltersModal.style.display = 'none';
            }
          };
        }
        
        // Actualizar resumen cuando cambian los selectores manualmente
        if (patientSel) {
          patientSel.addEventListener('change', function() {
            updateFiltersSummary();
          });
        }
        
        // Actualizar resumen inicial
        setTimeout(() => {
          updateFiltersSummary();
        }, 100);
    """
    # ───────── HTML ─────────
    # Usar order directamente en lugar de bars.keys() que está vacío
    gene_options = "".join(f'<option value="{g}" selected>{g}</option>'
                           for g in VIEWS["gene"]["order"])
    nsp_options  = "".join(f'<option value="{p}" selected>{p}</option>'
                           for p in VIEWS["protein"]["order"])

    html = f"""<!DOCTYPE html>
    <!-- Generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} -->
    <html><head><meta charset="utf-8">
    <title>SARS-CoV-2 Selection Analysis</title>
    <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
    <style>
     * {{
       margin: 0;
       padding: 0;
       box-sizing: border-box;
     }}
     
     body {{
       font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
       background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
       min-height: 100vh;
       color: #333;
       line-height: 1.6;
     }}
     
     .container {{
       max-width: 1400px;
       margin: 0 auto;
       padding: 20px;
     }}
     
     .header {{
       background: rgba(255, 255, 255, 0.95);
       backdrop-filter: blur(10px);
       border-radius: 15px;
       padding: 25px;
       margin-bottom: 30px;
       box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
       border: 1px solid rgba(255, 255, 255, 0.2);
       position: sticky;
       top: 20px;
       z-index: 1000;
     }}
     
     .header h1 {{
       color: #2c3e50;
       font-size: 2.2em;
       margin-bottom: 15px;
       text-align: center;
       font-weight: 300;
     }}
     
     .controls {{
       display: flex;
       flex-wrap: wrap;
       gap: 20px;
       align-items: center;
       justify-content: center;
     }}
     
     .control-group {{
       display: flex;
       flex-direction: column;
       align-items: center;
       gap: 8px;
     }}
     
     .control-group.hidden {{
       display: none;
     }}
     
     .control-group label {{
       font-weight: 600;
       color: #34495e;
       font-size: 0.9em;
       text-transform: uppercase;
       letter-spacing: 0.5px;
     }}
     
     select {{
       padding: 10px 15px;
       border: 2px solid #e0e6ed;
       border-radius: 8px;
       background: white;
       font-size: 14px;
       min-width: 150px;
       transition: all 0.3s ease;
       box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
     }}
     
     select:focus {{
       outline: none;
       border-color: #667eea;
       box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
     }}
     
     button {{
       padding: 12px 24px;
       border: none;
       border-radius: 8px;
       font-weight: 600;
       cursor: pointer;
       transition: all 0.3s ease;
       font-size: 14px;
       text-transform: uppercase;
       letter-spacing: 0.5px;
     }}
     
     #toggle {{
       background: linear-gradient(135deg, #667eea, #764ba2);
       color: white;
       box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
     }}
     
     #toggle:hover {{
       transform: translateY(-2px);
       box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
     }}
     
     #back {{
       background: linear-gradient(135deg, #e74c3c, #c0392b);
       color: white;
       box-shadow: 0 4px 15px rgba(231, 76, 60, 0.3);
     }}
     
     #back:hover {{
       transform: translateY(-2px);
       box-shadow: 0 6px 20px rgba(231, 76, 60, 0.4);
     }}
     
     .chart-section {{
       background: rgba(255, 255, 255, 0.95);
       backdrop-filter: blur(10px);
       border-radius: 15px;
       padding: 25px;
       margin-bottom: 30px;
       box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
       border: 1px solid rgba(255, 255, 255, 0.2);
     }}
     
     .chart-section h2 {{
       color: #2c3e50;
       font-size: 1.5em;
       margin-bottom: 20px;
       font-weight: 400;
       border-bottom: 2px solid #e0e6ed;
       padding-bottom: 10px;
     }}
     
     .chart-container {{
       margin-bottom: 30px;
       position: relative;
       overflow: hidden;
       min-height: 350px;
     }}
     
     .chart-container:last-child {{
       margin-bottom: 0;
     }}
     
     .chart-container > div {{
       width: 100% !important;
       height: auto !important;
       min-height: 350px;
     }}
     
     .footnote {{
       font-size: 12px;
       margin-top: 15px;
       color: #7f8c8d;
       font-style: italic;
       text-align: center;
       padding: 10px;
       background: rgba(236, 240, 241, 0.5);
       border-radius: 8px;
     }}
     
     .popup-table {{
       background: white;
       border-radius: 10px;
       padding: 20px;
       box-shadow: 0 10px 30px rgba(0, 0, 0, 0.2);
       max-height: 400px;
       overflow-y: auto;
     }}
     
     .popup-table h3 {{
       color: #2c3e50;
       margin-bottom: 15px;
       font-size: 1.3em;
     }}
     
     .popup-table table {{
       width: 100%;
       border-collapse: collapse;
       margin-top: 15px;
     }}
     
     .popup-table th, .popup-table td {{
       padding: 10px;
       text-align: left;
       border-bottom: 1px solid #ecf0f1;
     }}
     
     .popup-table th {{
       background: #f8f9fa;
       font-weight: 600;
       color: #2c3e50;
     }}
     
     .popup-table a {{
       color: #667eea;
       text-decoration: none;
       font-weight: 500;
       padding: 4px 8px;
       border-radius: 4px;
       background: rgba(102, 126, 234, 0.1);
       transition: all 0.3s ease;
     }}
     
     .popup-table a:hover {{
       background: rgba(102, 126, 234, 0.2);
       transform: translateY(-1px);
     }}
     
     .close-btn {{
       position: absolute;
       top: 10px;
       right: 15px;
       background: #e74c3c;
       color: white;
       border: none;
       border-radius: 50%;
       width: 30px;
       height: 30px;
       cursor: pointer;
       font-size: 16px;
       display: flex;
       align-items: center;
       justify-content: center;
       transition: all 0.3s ease;
     }}
     
     .close-btn:hover {{
       background: #c0392b;
       transform: scale(1.1);
     }}
     
     .loading {{
       text-align: center;
       padding: 40px;
       color: #7f8c8d;
       font-style: italic;
     }}
     
     .popup-overlay {{
       position: fixed;
       top: 0;
       left: 0;
       width: 100%;
       height: 100%;
       background: rgba(0, 0, 0, 0.5);
       display: flex;
       align-items: center;
       justify-content: center;
       z-index: 10000;
     }}
     
     .popup-content {{
       background: white;
       padding: 20px;
       border-radius: 10px;
       box-shadow: 0 10px 30px rgba(0,0,0,0.3);
       max-width: 90vw;
       max-height: 80vh;
       overflow-y: auto;
       position: relative;
     }}
     
     .popup-header {{
       display: flex;
       justify-content: space-between;
       align-items: center;
       margin-bottom: 20px;
       padding-bottom: 10px;
       border-bottom: 1px solid #ecf0f1;
     }}
     
     .popup-header h3 {{
       margin: 0;
       color: #2c3e50;
     }}
     
     .position-section {{
       margin-bottom: 20px;
     }}
     
     .position-section h4 {{
       margin: 0 0 10px 0;
       color: #333;
       font-size: 14px;
     }}
     
     .position-list {{
       display: flex;
       flex-wrap: wrap;
       gap: 5px;
       max-height: 200px;
       overflow-y: auto;
       padding: 10px;
       background: #f8f9fa;
       border-radius: 5px;
     }}
     
     .position-tag {{
       padding: 2px 8px;
       border-radius: 12px;
       font-size: 11px;
       font-weight: bold;
       display: inline-block;
     }}
     
             .position-tag.selected {{
          background: #fd7e14;
          color: white;
        }}
        
        .position-tag.neutral {{
          background: #6c757d;
          color: white;
        }}
        
        .position-tag.non-variable {{
          background: #198754;
          color: white;
        }}
     
     @media (max-width: 768px) {{
       .controls {{
         flex-direction: column;
         align-items: stretch;
       }}
       
       .control-group {{
         width: 100%;
       }}
       
       select {{
         width: 100%;
       }}
     }}
    </style></head><body>
    <div class="container">
      <!-- Panel de resumen de filtros activos (se mostrará cuando haya filtros) -->
      <div id="filters-summary" style="display: none; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px 20px; margin-bottom: 20px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.2);">
        <h3 style="margin: 0 0 15px 0; font-size: 18px; font-weight: bold; text-align: center; border-bottom: 2px solid rgba(255,255,255,0.3); padding-bottom: 10px;">
          👥 Patients with these filters (<span id="selected-patients-count">0</span>)
        </h3>
        <div id="selected-patients-list" style="font-size: 15px; line-height: 2; max-height: 120px; overflow-y: auto; background: rgba(255,255,255,0.15); padding: 15px; border-radius: 4px; margin-bottom: 15px; text-align: center; font-weight: 500;">
          <!-- Lista de pacientes generada dinámicamente -->
        </div>
        <div style="border-top: 1px solid rgba(255,255,255,0.2); padding-top: 15px;">
          <h4 style="margin: 0 0 10px 0; font-size: 16px; font-weight: bold;">📋 Active Filters:</h4>
          <div id="filters-summary-content" style="font-size: 14px; line-height: 1.8;">
            <!-- Contenido generado dinámicamente -->
          </div>
        </div>
      </div>
      
      <div class="header">
        <h1>🧬 SARS-CoV-2 Selection Analysis</h1>
        <div class="controls">
          <div class="control-group">
            <button id="code-btn" style="background: #2c3e50; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; margin-right: 10px;">💻 Code</button>
            <button id="methodology-btn" style="background: #3498db; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; margin-right: 10px;">📖 Methodology</button>
            <button id="more-filters-btn" style="background: #9b59b6; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; margin-right: 10px;">🔍 More Filters</button>
          </div>
          <div class="control-group">
            <label for="patient-selector">👥 Patients</label>
            <select id="patient-selector" multiple size="5">
        {"".join(f'<option value="{p}" selected>{p}</option>' for p in patients)}
      </select>
          </div>
          <div class="control-group" id="gene-controls">
            <label for="gene-list">🧬 Genes</label>
            <select id="gene-list" multiple size="5">
        {gene_options}
      </select>
          </div>
          <div class="control-group hidden" id="protein-controls">
            <label for="nsp-list">🦠 Proteins</label>
            <select id="nsp-list" multiple size="5">
        {nsp_options}
      </select>
    </div>
          <div class="control-group">
            <button id="back" style="display:none">⬅ Back</button>
            <button id="toggle">🔄 Protein View</button>
          </div>
          <div class="control-group">
            <label for="min-patients">Min patients per site</label>
            <div style="display: flex; align-items: center; gap: 5px;">
              <input id="min-patients" type="number" value="1" min="1" step="1" style="width:80px;padding:6px 8px;border:1px solid #e0e6ed;border-radius:6px;">
              <button id="min-patients-info" style="background: #3498db; color: white; border: none; border-radius: 50%; width: 16px; height: 16px; font-size: 10px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;" title="Click for information">ℹ</button>
            </div>
          </div>
          <div style="display:flex;gap:10px;align-items:center;margin-top:10px;">
            <button id="apply-filters" style="background: linear-gradient(135deg, #27ae60, #2ecc71); color: white; box-shadow: 0 4px 15px rgba(39, 174, 96, 0.3);">✅ Apply Filters</button>
            <button id="reset-filters" style="background: linear-gradient(135deg, #e67e22, #f39c12); color: white; box-shadow: 0 4px 15px rgba(230, 126, 34, 0.3);">🔄 Reset</button>
          </div>
        </div>
      </div>
      
      <div class="chart-section">
        <h2>📊 Selection Analysis</h2>
        <div class="heatmap-controls" style="text-align: center; margin-bottom: 15px;">
          <label style="font-weight: bold; margin-right: 10px;">Heatmap 1 View:</label>
          <select id="heatmap1-view" style="padding: 8px; border-radius: 4px; border: 1px solid #ccc; font-size: 14px;">
            <option value="mean">Mean Selection</option>
            <option value="selected">Selected/Total</option>
            <option value="non-variable">Non-Variable/Total</option>
          </select>
        </div>
        <div id="var-fig"></div>
        <p class="footnote">
          <b>Figure 1 – Selection analysis per patient versus gene/protein.</b> 
          Green indicates negative selection, white indicates neutral selection, orange indicates positive selection.
          Use the 'protein view / gene view' button (top-right) to swap the x-axis.
        </p>
      </div>
      
      <div class="chart-section">
        <h2>📊 β (slope) mean Analysis</h2>
        <div class="heatmap-controls" style="text-align: center; margin-bottom: 15px;">
          <label style="font-weight: bold; margin-right: 10px;">Heatmap 2 View:</label>
          <select id="heatmap2-view" style="padding: 8px; border-radius: 4px; border: 1px solid #ccc; font-size: 14px;">
            <option value="mean">Mean β</option>
            <option value="positive">Positive/Total</option>
            <option value="negative">Negative/Total</option>
          </select>
        </div>
        <div id="fig"></div>
        <p class="footnote">
          <b>Figure 2 – Mean β slope per patient versus gene/protein.</b> 
          Use the 'protein view / gene view' button (top-right) to swap the x-axis.
        </p>
      </div>
      
       <div class="chart-section">
        <h2>📈 Temporal Diversity</h2>
        <div id="var-plot" style="width:100%;height:350px;"></div>
        <p id="var-footnote" class="footnote">
          <b>Figure 3 – Temporal variability (diversity)</b> per genomic site, averaged across selected patients. 
          <b>Selected positions:</b> <span style="color:#00FF00">●</span> 1 patient, 
          <span style="color:#FFFF00">●</span> 2 patients, 
          <span style="color:#FFA500">●</span> 3 patients, 
          <span style="color:#FF8C00">●</span> 4-5 patients, 
          <span style="color:#FF4500">●</span> 6-8 patients, 
          <span style="color:#FF0000">●</span> 9+ patients.
        </p>
      </div>

       <div class="chart-section">
         <h2>🌳 UPGMA – Selection by genes</h2>
         <div id="upgma-genes-div" style="width:100%;height:380px;"></div>
         <p class="footnote">
           <b>Figure 4.</b> UPGMA dendrogram (binary Jaccard distance) by genes. Each leaf represents a patient, and the tree structure shows clustering based on temporal variability patterns across genes. 
           <b>How it's calculated:</b> Binary Jaccard distance is computed from presence/absence patterns of variants across genes for each patient. UPGMA clustering groups patients with similar variant profiles.
           <b>Branch heights:</b> Represent evolutionary distances between patient clusters. Shorter branches indicate more similar variant patterns.
           <b>Interpretation:</b> Patients clustered together share similar temporal variability patterns across genes. The tree structure reveals evolutionary relationships in viral variant selection.
         </p>
       </div>

       <div class="chart-section">
         <h2>🌳 UPGMA – Selection by proteins</h2>
         <div id="upgma-proteins-div" style="width:100%;height:380px;"></div>
         <p class="footnote">
           <b>Figure 4.</b> UPGMA dendrogram (binary Jaccard distance) by proteins. Each leaf represents a patient, and the tree structure shows clustering based on temporal variability patterns across proteins.
           <b>How it's calculated:</b> Binary Jaccard distance is computed from presence/absence patterns of variants across proteins for each patient. UPGMA clustering groups patients with similar variant profiles.
           <b>Branch heights:</b> Represent evolutionary distances between patient clusters. Shorter branches indicate more similar variant patterns.
           <b>Interpretation:</b> Patients clustered together share similar temporal variability patterns across proteins. The tree structure reveals evolutionary relationships in viral variant selection at the protein level.
         </p>
       </div>

       <div class="chart-section" style="text-align:center; margin-top: 10px;">
         <button id="upgma-filter-btn" style="background: #2ecc71; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; box-shadow: 0 4px 12px rgba(46, 204, 113, 0.3);">
           🌳 UPGMA: Filter controls
         </button>
       </div>
      
      <div class="chart-section" id="gene-section">
        <h2>🧬 Gene Analysis</h2>
        <div class="chart-container">
          <div id="gene-chart" style="width:100%;height:350px;"></div>
        </div>
        <p class="footnote">
          <b>Figure 5 – Per-site mutation frequency (number of patients).</b> 
          INTERACTIONS: • Hover to see β⁺/β⁻ counts • 
          Click any bar to open a detailed table with β, p-value and quick links to the per-site PDF plot and full-genome PDF.
        </p>
      </div>
      
      <div class="chart-section" id="protein-section" style="display:none;">
        <h2>🦠 Protein Analysis</h2>
        <div class="chart-container">
          <div id="nsp-chart" style="width:100%;height:350px;"></div>
        </div>
        <p class="footnote">
          <b>Figure 5 – Per-site mutation frequency (number of patients).</b> 
          INTERACTIONS: • Hover to see β⁺/β⁻ counts • 
          Click any bar to open a detailed table with β, p-value and quick links to the per-site PDF plot and full-genome PDF.
        </p>
      </div>
    </div>
    <!-- Modal for Min Patients Per Site Information -->
    <div id="min-patients-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 5% auto; padding: 20px; border-radius: 10px; width: 80%; max-width: 600px; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">ℹ️ Min Patients Per Site Filter</h3>
          <button id="close-modal" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div style="line-height: 1.6; color: #34495e;">
          <p><strong>What does this filter do?</strong></p>
          <p>This filter controls the minimum number of patients that must have data for a genomic position to be displayed in the visualizations. It helps focus on positions with sufficient data coverage across multiple patients.</p>
          
          <p><strong>How it works:</strong></p>
          <ul>
            <li><strong>Figure 3 (Temporal Variability):</strong> Filters which genomic positions appear as points in the scatter plot. Only positions with data in the minimum number of patients are shown as points. The trend line remains unchanged.</li>
            <li><strong>Figure 6 (Bar Charts):</strong> Filters which bars are displayed - bars with fewer patients than the threshold are hidden</li>
          </ul>
          
          <p><strong>Examples:</strong></p>
          <ul>
            <li><strong>Min = 1:</strong> Shows all positions that have data in at least 1 patient (default)</li>
            <li><strong>Min = 2:</strong> Only shows positions with data in 2 or more patients</li>
            <li><strong>Min = 3:</strong> Only shows positions with data in 3 or more patients</li>
          </ul>
          
          <p><strong>Note:</strong> This filter works in combination with the patient selection filter. If you select specific patients, the filter counts only those selected patients.</p>
        </div>
      </div>
    </div>
    
    <!-- Modal for Methodology -->
    <div id="methodology-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 2% auto; padding: 20px; border-radius: 10px; width: 90%; max-width: 800px; max-height: 90vh; overflow-y: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">📖 Methodology</h3>
          <button id="close-methodology-modal" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div style="line-height: 1.6; color: #34495e; white-space: pre-line;">{methodology_content}</div>
      </div>
    </div>

    <!-- Modal for UPGMA Visibility Filters -->
    <div id="upgma-filter-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 5% auto; padding: 20px; border-radius: 10px; width: 80%; max-width: 700px; max-height: 90vh; overflow-y: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">🌳 UPGMA Visibility Filters</h3>
          <button id="close-upgma-filter-modal" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div style="margin-bottom: 15px; color: #34495e;">
          <p style="margin: 0 0 10px 0;"><strong>Patient metadata categories required for inclusion (patients must have data):</strong></p>
          <div id="upgma-filter-categories" style="padding: 10px; border: 1px solid #dee2e6; border-radius: 6px; max-height: 300px; overflow-y: auto;"></div>
        </div>
        <div style="margin-bottom: 20px;">
          <label style="cursor: pointer;">
            <input type="checkbox" id="upgma-include-missing" style="margin-right: 8px;">
            Include patients not present in the filter TSV
          </label>
        </div>
        <div style="text-align: right;">
          <button id="apply-upgma-filter-btn" style="background: linear-gradient(135deg, #27ae60, #2ecc71); color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; font-weight: bold; box-shadow: 0 4px 15px rgba(39, 174, 96, 0.3);">Apply</button>
        </div>
      </div>
    </div>
    
    <!-- Modal for Advanced Filters Info -->
    <div id="advanced-filters-info-overlay" style="display: none; position: fixed; z-index: 100000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div id="advanced-filters-info-modal" style="background-color: white; margin: 5% auto; padding: 20px; border-radius: 10px; width: 80%; max-width: 700px; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 id="advanced-filters-info-title" style="margin: 0; color: #2c3e50;">ℹ️ Advanced Filters Info</h3>
          <button id="advanced-filters-info-close" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div id="advanced-filters-info-content" style="line-height: 1.6; color: #34495e; white-space: pre-line;"></div>
      </div>
    </div>
    
    <!-- Modal for More Filters -->
    <div id="more-filters-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 2% auto; padding: 20px; border-radius: 10px; width: 90%; max-width: 1800px; max-height: 90vh; overflow-y: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">🔍 More Filters</h3>
          <button id="close-more-filters-modal" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div style="margin-bottom: 20px; display: flex; gap: 10px; flex-wrap: wrap;">
          <button id="apply-filters-btn" style="background: linear-gradient(135deg, #27ae60, #2ecc71); color: white; border: none; border-radius: 6px; padding: 10px 20px; cursor: pointer; font-size: 16px; font-weight: bold; box-shadow: 0 4px 15px rgba(39, 174, 96, 0.3);">Apply Filters</button>
          <button id="reset-advanced-filters-btn" style="background: #95a5a6; color: white; border: none; border-radius: 6px; padding: 10px 20px; cursor: pointer; font-size: 16px; font-weight: bold;">Reset Advanced Filters</button>
        </div>
        <div id="filters-table-container" style="line-height: 1.6; color: #34495e;">
          <!-- La tabla se generará dinámicamente con JavaScript -->
        </div>
      </div>
    </div>
    
    <!-- Modal for Covarying Information -->
    <div id="covarying-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 5% auto; padding: 20px; border-radius: 10px; width: 80%; max-width: 600px; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">ℹ️ Covarying Filter Information</h3>
          <button id="close-covarying-modal" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div style="line-height: 1.6; color: #34495e;">
          <p><strong>What does this filter do?</strong></p>
          <p>This filter allows you to filter genomic positions based on their covariation status and the number of covarying positions.</p>
          
          <p><strong>Filter options:</strong></p>
          <ul>
            <li><strong>All positions:</strong> Shows all positions regardless of covariation status</li>
            <li><strong>Only covarying:</strong> Shows only positions that are covarying with other positions</li>
            <li><strong>Only non-covarying:</strong> Shows only positions that are not covarying</li>
          </ul>
          
          <p><strong>Min covarying count:</strong> When "Only covarying" is selected, this sets the minimum number of covarying positions required.</p>
          
          <p><strong>Date filter:</strong> Filter variants by their appearance date. Only variants that appeared on or before the selected date will be shown in the variability plot and bar charts.</p>
          
          <p><strong>Note:</strong> Click on "✓ Yes" in the Covarying column of tables to see detailed covariation information.</p>
        </div>
      </div>
    </div>
    
    <!-- Modal for Events Details -->
    <div id="events-modal" style="display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5);">
      <div style="background-color: white; margin: 2% auto; padding: 20px; border-radius: 10px; width: 90%; max-width: 900px; max-height: 90vh; overflow-y: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
          <h3 style="margin: 0; color: #2c3e50;">📊 Events Details</h3>
          <button onclick="document.getElementById('events-modal').style.display='none'" style="background: #e74c3c; color: white; border: none; border-radius: 50%; width: 30px; height: 30px; font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center; font-weight: bold;">×</button>
        </div>
        <div id="events-content" style="line-height: 1.6; color: #34495e;"></div>
      </div>
    </div>
    
    {BLOBS}
    {var_blob}
    {mapping_blob}
    {upgma_blob}
    {segments_blob}
    {selected_positions_blob}
    {selected_positions_with_flags_blob}
    {filter_data_blob}
    <script id="data-directory" type="application/json">"{args.dir}"</script>
    <script id="initial-mode" type="application/json">"gene"</script>
    <script id="toggle-target" type="application/json">""</script>
    <script>{JS}</script>
    </body></html>"""

    gene_view_blob = blob("gene", _public_view_payload(VIEWS["gene"]))
    protein_view_blob = blob("protein", _public_view_payload(VIEWS["protein"]))
    gene_upgma_blob = f'<script id="upgma-genes" type="application/json">{_json_dumps_compact(upgma_genes, cls=PlotlyJSONEncoder)}</script>'
    protein_upgma_blob = f'<script id="upgma-proteins" type="application/json">{_json_dumps_compact(upgma_proteins, cls=PlotlyJSONEncoder)}</script>'
    gene_segments_blob = f'<script id="segments" type="application/json">{_json_dumps_compact({"genes": segments["genes"], "proteins": []}, cls=PlotlyJSONEncoder)}</script>'
    protein_segments_blob = f'<script id="segments" type="application/json">{_json_dumps_compact({"genes": [], "proteins": segments["proteins"]}, cls=PlotlyJSONEncoder)}</script>'

    out_path = Path(args.out)
    stem = out_path.stem
    if stem.endswith("_genes"):
        gene_out = out_path
        protein_out = out_path.with_name(f"{stem[:-6]}_proteins{out_path.suffix or '.html'}")
    elif stem.endswith("_proteins"):
        protein_out = out_path
        gene_out = out_path.with_name(f"{stem[:-9]}_genes{out_path.suffix or '.html'}")
    else:
        suffix = out_path.suffix or ".html"
        gene_out = out_path.with_name(f"{stem}_genes{suffix}")
        protein_out = out_path.with_name(f"{stem}_proteins{suffix}")

    gene_html = (
        html.replace("__VIEW_BLOBS__", gene_view_blob)
            .replace("__UPGMA_BLOB__", gene_upgma_blob)
            .replace("__SEGMENTS_BLOB__", gene_segments_blob)
            .replace('<script id="initial-mode" type="application/json">"gene"</script>', '<script id="initial-mode" type="application/json">"gene"</script>')
            .replace('<script id="toggle-target" type="application/json">""</script>', f'<script id="toggle-target" type="application/json">{_json_dumps_compact(protein_out.name)}</script>')
    )
    protein_html = (
        html.replace("__VIEW_BLOBS__", protein_view_blob)
            .replace("__UPGMA_BLOB__", protein_upgma_blob)
            .replace("__SEGMENTS_BLOB__", protein_segments_blob)
            .replace('<script id="initial-mode" type="application/json">"gene"</script>', '<script id="initial-mode" type="application/json">"protein"</script>')
            .replace('<script id="toggle-target" type="application/json">""</script>', f'<script id="toggle-target" type="application/json">{_json_dumps_compact(gene_out.name)}</script>')
    )

    if args.minify_html_inline:
        gene_html = _minify_html_inline(gene_html)
        protein_html = _minify_html_inline(protein_html)

    gene_out.write_text(gene_html, encoding="utf-8")
    protein_out.write_text(protein_html, encoding="utf-8")
    print("[OK] HTML saved →", str(gene_out))
    print("[OK] HTML saved →", str(protein_out))

if __name__ == "__main__":
    try:
        if get_start_method(allow_none=True) != "fork":
            import multiprocessing as mp
            mp.set_start_method("fork")
    except RuntimeError:
        pass
    main()
