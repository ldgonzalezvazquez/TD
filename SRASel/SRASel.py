#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rank_selection_variants_fast.py   (v3  –  keep-ref + tidy plots)
────────────────────────────────────────────────────────────────
• Calcula β ponderado por profundidad + p-val permutación.
• No descarta la base mayoritaria (la “referencia local”).
• Genera:
    <id>_long.tsv
    <id>_variant_stats.tsv
    <id>_selected_ranked.tsv
    <id>_genome_selection.pdf
    <id>_plots/posXXXXX_ALL.pdf  (una por posición con rectas)
"""

import argparse, csv, re, pathlib, warnings, os, math, json, sys
import datetime as dt
from collections import defaultdict
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm
from scipy.ndimage import gaussian_filter1d
warnings.filterwarnings("ignore")
os.environ.setdefault("MKL_VERBOSE", "NO")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from typing import Tuple, List, Dict, Optional


TOKEN_RE = re.compile(r"^(\d+)")
AA_RE = re.compile(r"^([A-Z\*\?]+)(\d+)(?:\([^)]*\))?$|^(\d+)([A-Z\*\?]+)(?:\([^)]*\))?$|^(\d+)(amb)(?:\([^)]*\))?$|^(amb)(\d+)(?:\([^)]*\))?$")

# tablas de traducción sencillas ───────────────
CODON2AA = {
    #  ⇢ tabla estándar (solo tripletes que pueden aparecer en SARS-CoV-2)
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGG":"W","TGA":"*",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# ───────────────────── single final report helper ─────────────────────

def build_final_report_text(patient: str,
                            sigma: int = 350,
                            nperm: int = 2000,
                            freq_delta: float = 0.03,
                            thr_clear: float = 0.95,
                            eps_tie: float = 0.05,
                            depth_any: int = 100) -> str:
    """
    Returns a single multi-line string that documents, in English (Nature-like),
    all figures and tables produced for a given patient/run, with precise definitions.
    """
    return f"""
[REPORT — {patient}]

FIGURES
1) Genome-wide selection landscape (file: {patient}_genome_selection.pdf).
   This figure summarizes selection across the genome. Points mark labeled trajectories by genomic
   coordinate (x-axis) and weighted slope β (y-axis; time in days). Here β quantifies the per-day
   change in variant frequency f = k/n estimated by weighted least squares (WLS) regression of f on day t
   with weights w = n (read depth). Specifically:
       β = Sxy / Sxx,  where
       Sxy = Σ w_i (t_i − t̄_w)(f_i − f̄_w),   Sxx = Σ w_i (t_i − t̄_w)^2,
       t̄_w = (Σ w_i t_i) / (Σ w_i),          f̄_w = (Σ w_i f_i) / (Σ w_i).
   Variants are labeled ‘pos’ if p < 0.05 and β > 0, or ‘neg’ if p < 0.05 and β < 0. The p-values are
   obtained by a permutation test that shuffles the sampling days within each trajectory (n = {nperm} permutations).
   The right axis shows kernel-smoothed densities (Gaussian, σ = {sigma}) of selected positions, plotted upward for β > 0 (red)
   and downward for β < 0 (blue). Dashed horizontal lines denote neutral baselines: β̄_all (grey; all variants)
   and β̄_sel (black; only labeled variants). A colored bar above the x-axis delineates coding-sequence spans by gene.

2) Per-position temporal trajectories (files: {patient}_plots/posXXXXX_ALL.pdf).
   For each genomic position with at least one labeled trajectory, the panel displays the observed frequencies (points)
   across days for every (pos, var). Colored solid lines are WLS fits for labeled trajectories (β > 0, red; β < 0, blue);
   unlabeled trajectories are shown in light grey dashed. Each fit is drawn in its WLS form, i.e., a straight line
   with slope β that passes through the weighted means (t̄_w, f̄_w), which exactly represents the WLS solution without
   explicitly estimating an intercept. For context, genome-wide neutral reference lines are shown: β̄_all (grey) and β̄_sel (black).

TABLES
1) {patient}_long.tsv
   Long-format counts per (position, variant, day). Columns:
     pos, var, day, k (reads supporting the variant), n (total depth).
   Variants include nucleotides (A/C/G/T/N) and aggregated indels (INS, DEL).
   Rows are retained only for (pos, var) pairs that reach total depth ≥ {depth_any} on at least one day.

2) {patient}_variant_stats.tsv
   Per-(position, variant) statistics after a minimum temporal-change filter Δf ≥ {freq_delta:.3f},
   where Δf is computed from the fitted slope as Δf = β · (max(day) − min(day)).
   β is the WLS slope as defined above (weights w = n). Permutation p-values (n = {nperm}) are obtained by shuffling days
   within each trajectory; a trajectory is labeled ‘pos’ if p < 0.05 and β > 0, or ‘neg’ if p < 0.05 and β < 0.
   The score is |β| · (−log10 p) with the sign of β. Amino-acid annotation builds a codon by forcing the focal base at the
   position and choosing “clear” neighbors (frequency ≥ {thr_clear:.2f}) that differ from the focal position’s top frequency
   by more than {eps_tie:.2f}; ambiguity is reported as ‘ambiguous’ or ‘ambiguous (indel)’. Columns:
     gene, aa_idx, pos, var, beta, pval, score, tag, delta_f, codon, aa_change, ambiguous, annot_detail,
     depth (mean n), n_days, and JSON arrays for day, k, n.

3) {patient}_selected_ranked.tsv
   Ranked list of labeled variants (tag ∈ {{pos, neg}}), sorted by descending absolute score (|β| · −log10 p, signed by β).
   Provides a concise prioritization of trajectories for follow-up.

""".strip()

# ───────── Covariación (solo variantes seleccionadas, al final) ─────────
def _compute_covariation_annotations_selected(
    ranked_pre: pd.DataFrame,
    method: str = "beta_diff",
    threshold: float = 0.01,
    min_days: int = 3,
    min_range: float = 0.05,
    max_day_diff: int = 1
) -> pd.DataFrame:
    """
    Calcula co-variación entre variantes SELECCIONADAS (tag in {'pos','neg'}) dentro de un paciente.
    MÉTODO: Compara valores de beta directamente (no correlación de frecuencias).
    Dos variantes covarían si:
      1) |beta₁ - beta₂| ≤ threshold
      2) |día_inicio₁ - día_inicio₂| ≤ max_day_diff (aparecen en momentos similares)
    
    Devuelve un DataFrame con columnas:
        - pos, var
        - covarying (bool)
        - covarying_count (int)
        - covarying_with (str): lista formateada de pares con los que co-varía, incluyendo diff_beta y diff_day.
          Formato por pareja: nt=<pos>(<var>); aa=<gene>:<aa_change> [Δβ=0.xx, Δday=X]
    """
    if ranked_pre is None or ranked_pre.empty:
        return pd.DataFrame(columns=["pos","var","covarying","covarying_count","covarying_with"])

    # Trabajamos sólo con variantes seleccionadas
    sel = ranked_pre[ ranked_pre["tag"].isin(["pos","neg"]) ].copy()
    print(f"🔍 Variantes seleccionadas encontradas: {len(sel)}")
    if sel.empty:
        print("⚠️  No hay variantes seleccionadas para análisis de covariación")
        return pd.DataFrame(columns=["pos","var","covarying","covarying_count","covarying_with"])

    # Índice de metadatos por token
    # token = (pos:int, var:str)
    def _first_nonempty(row: pd.Series, keys: List[str]) -> Optional[str]:
        for k in keys:
            if k in row and pd.notna(row[k]) and str(row[k]).strip() != "":
                return str(row[k])
        return None

    meta_index: Dict[Tuple[int,str], Dict[str, Optional[str]]] = {}
    beta_dict: Dict[Tuple[int,str], float] = {}
    first_day_dict: Dict[Tuple[int,str], int] = {}
    
    for _, r in sel.iterrows():
        tok = (int(r["pos"]), r["var"])
        if tok in meta_index:
            continue
        meta_index[tok] = {
            "gene": _first_nonempty(r, ["gene","region","protein","prot"]),
            "aa_change": _first_nonempty(r, ["aa_change","aa","aa_change_str"]),
            "aa_pos": _first_nonempty(r, ["aa_pos","aa_position","aapos"]),
            "var_nt": _first_nonempty(r, ["var","nt_var","nt_change"]),
        }
        beta_dict[tok] = float(r["beta"])
        
        # Extraer día de primera aparición (mínimo de la lista 'day')
        if "day" in r and isinstance(r["day"], (list, tuple)) and len(r["day"]) > 0:
            first_day_dict[tok] = int(min(r["day"]))
        else:
            first_day_dict[tok] = 0  # Si no hay datos, asumir día 0

    def _fmt_partner(tok: Tuple[int,str], diff_beta: float, diff_day: int) -> str:
        pos, var = tok
        m = meta_index.get(tok, {})
        nt_var = m.get("var_nt") or str(var) or ""
        nt_part = f"nt={pos}" + (f"({nt_var})" if nt_var and nt_var != str(pos) else "")
        gene = m.get("gene") or "NA"
        aa_change = (m.get("aa_change") or "").strip()
        aa_pos = (m.get("aa_pos") or "").strip()
        if aa_change and aa_change.lower() != "nan":
            aa_part = f"aa={gene}:{aa_change}"
        elif aa_pos and aa_pos.lower() != "nan":
            aa_part = f"aa={gene}:{aa_pos}"
        else:
            aa_part = f"aa={gene}"
        return f"{nt_part}; {aa_part} [Δβ={diff_beta:.4f}, Δday={diff_day}]"

    # Comparación de betas Y temporalidad entre todas las variantes
    print(f"🔍 Analizando covariación por similitud de β (threshold: Δβ ≤ {threshold}, Δday ≤ {max_day_diff})...")
    print(f"📊 Variantes a analizar: {len(beta_dict)}")
    
    out_rows = []
    covarying_found = 0
    
    for tok1 in beta_dict.keys():
        beta1 = beta_dict[tok1]
        day1 = first_day_dict[tok1]
        partners = []
        
        # Comparar con todas las demás variantes
        for tok2 in beta_dict.keys():
            if tok1 == tok2:  # No comparar consigo misma
                continue
            
            beta2 = beta_dict[tok2]
            day2 = first_day_dict[tok2]
            
            diff_beta = abs(beta1 - beta2)
            diff_day = abs(day1 - day2)
            
            # Covariación si AMBOS criterios se cumplen:
            # 1) Beta similar
            # 2) Inicio temporal similar
            if diff_beta <= threshold and diff_day <= max_day_diff:
                partners.append((tok2, diff_beta, diff_day))
        
        # Formatear parejas
        partners_fmt = [_fmt_partner(tok, diff_b, diff_d) for tok, diff_b, diff_d in partners]
        
        is_covarying = bool(partners)
        if is_covarying:
            covarying_found += 1
            if covarying_found <= 3:  # Debug para las primeras 3
                print(f"  Covariación {covarying_found}: {tok1} con {len(partners)} parejas (β={beta1:.4f}, día_inicio={day1})")
        
        out_rows.append({
            "pos": tok1[0],
            "var": tok1[1],
            "covarying": is_covarying,
            "covarying_count": int(len(partners)),
            "covarying_with": " | ".join(partners_fmt) if partners_fmt else ""
        })

    print(f"📊 Covariaciones encontradas: {covarying_found} de {len(beta_dict)} variantes")
    return pd.DataFrame(out_rows)

# ───────────── helpers ─────────────────────────────────────────
def split_cell(cell: str):
    if cell == "NA":
        return None
    days, depth, bases, indels = cell.split("/", 3)
    days, depth = int(days), int(depth)
    counts = {m.group(1): int(m.group(2))
              for m in re.finditer(r"([ACTGN])(\d+)", bases)}
    # Extraer inserciones como diccionario {secuencia: count}
    ins_dict = {m.group(1): int(m.group(2)) 
                for m in re.finditer(r"\+([A-Z]+)(\d+)", "+"+indels)}
    # Extraer deleciones como diccionario {secuencia: count}
    dels_dict = {m.group(1): int(m.group(2)) 
                 for m in re.finditer(r"-([A-Z]+)(\d+)", "-"+indels)}
    return days, depth, counts, ins_dict, dels_dict

def process_line(line_data):
    """Procesa una línea del TSV y devuelve las filas generadas"""
    line, tokens, tp_num = line_data
    rows = []
    
    pos, gene = int(line[0]), line[1]
    # puede venir vacío ⇒ ValueError
    try:
        codon_pos = int(line[2])        # «1|2|3»
    except ValueError:
        codon_pos = None                # sin dato

    # aa_idx: si viene vacío/no numérico ⇒ marcar "ND" (no data)
    aa_raw = (line[3] or "").strip()
    if re.fullmatch(r"\d+", aa_raw):
        aa_idx_val = aa_raw
    else:
        # Cambia "ND" por "" si prefieres columna vacía
        aa_idx_val = "ND"

    by_tp = defaultdict(list)
    for tok, cell in zip(tokens, line[4:]):
        by_tp[tp_num[tok]].append(cell)

    for tp, clist in by_tp.items():
        n_tot = days = 0
        base_cnt = defaultdict(int)
        ins_cnt = defaultdict(int)
        dels_cnt = defaultdict(int)
        for c in clist:
            p = split_cell(c)
            if not p:
                continue
            d, n_i, bases, ins_dict, dels_dict = p
            days = d
            n_tot += n_i
            for b, x in bases.items():
                base_cnt[b] += x
            for seq, x in ins_dict.items():
                ins_cnt[seq] += x
            for seq, x in dels_dict.items():
                dels_cnt[seq] += x
        if n_tot == 0:
            continue

        # ► SNVs (A/C/T/G/N): igual que antes
        for b, c in base_cnt.items():
            if c > 0:
                rows.append((pos, b, days, c, n_tot))
        # ► INDELs como variantes separadas por secuencia
        #     Cada inserción única se trata como variante "+SEQ"
        for seq, c in ins_cnt.items():
            if c > 0:
                rows.append((pos, f"+{seq}", days, c, n_tot))
        # Cada deleción única se trata como variante "-SEQ"
        for seq, c in dels_cnt.items():
            if c > 0:
                rows.append((pos, f"-{seq}", days, c, n_tot))

    return rows, pos, gene, aa_idx_val, codon_pos

def process_chunk(chunk_data):
    """Procesa un chunk de filas ya procesadas del TSV en paralelo"""
    chunk_rows, tokens, tp_num, chunk_genes = chunk_data
    rows = []
    gene_of = {}
    aaidx_of = {}
    codonpos_of = {}
    
    # Las filas ya están procesadas, usar los mapeos de genes proporcionados
    for row in chunk_rows:
        pos, var, day, k, n = row
        rows.append(row)
        
        # Usar los mapeos de genes proporcionados
        if pos not in gene_of:
            gene_of[pos] = chunk_genes.get(pos, "unknown")
            aaidx_of[pos] = chunk_genes.get(f"{pos}_aa", "ND")
            codonpos_of[pos] = chunk_genes.get(f"{pos}_codon", None)
    
    return rows, gene_of, aaidx_of, codonpos_of

def make_long(tsv_path: pathlib.Path, depth_any=10, n_jobs=2, nperm=2000, freq_delta=0.01, thr_clear=0.95, eps_tie=0.05):
    """Nueva implementación: crea N archivos físicos _long separados y los procesa en paralelo"""
    pat = tsv_path.stem.split("_genes")[0]
    outdir = tsv_path.parent
    
    if n_jobs == 1:
        # Procesamiento secuencial (sin cambios)
        rows=[]; gene_of={}; aaidx_of={}; codonpos_of={}
        
        with tsv_path.open() as fh:
            rdr = csv.reader(fh, delimiter="\t")
            hdr = next(rdr); tokens = hdr[4:]
            tp_num = {tok: int(TOKEN_RE.match(tok).group(1)) for tok in tokens}

            # Leer todas las líneas en memoria
            fh.seek(0)
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr)  # saltar header
            all_lines = list(rdr)
            
            for line in tqdm(all_lines, desc=f"Leyendo {tsv_path.stem}", unit="pos"):
                result = process_line((line, tokens, tp_num))
                line_rows, pos, gene, aa_idx_val, codon_pos = result
                rows.extend(line_rows)
                gene_of[pos] = gene
                aaidx_of[pos] = aa_idx_val
                codonpos_of[pos] = codon_pos

        long = pd.DataFrame(rows, columns=["pos","var","day","k","n"])
        # mantener variantes que alcancen k ≥ depth_any en algún día
        keep = (long.groupby(["pos","var"])["k"].max() >= depth_any)
        long = long.merge(keep.reset_index().rename(columns={"k":"keep"}))
        return (long[long["keep"]].drop(columns="keep"),
                gene_of, aaidx_of, codonpos_of)
    
    else:
        # NUEVA IMPLEMENTACIÓN: Crear N archivos físicos separados
        print(f"🔄 Procesando {pat} con {n_jobs} cores...")
        
        # 1. Leer todas las líneas y aplicar filtro global de profundidad
        print("📖 Leyendo archivo TSV...")
        with tsv_path.open() as fh:
            rdr = csv.reader(fh, delimiter="\t")
            hdr = next(rdr); tokens = hdr[4:]
            tp_num = {tok: int(TOKEN_RE.match(tok).group(1)) for tok in tokens}

            fh.seek(0)
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr)  # saltar header
            all_lines = list(rdr)
        
        # Aplicar filtro global de profundidad ANTES de dividir en chunks
        print(f"🔍 Aplicando filtro global de profundidad (≥{depth_any} lecturas por variante)...")
        
        # Procesar todas las líneas para obtener datos completos
        all_rows = []
        genes_info = {}  # Diccionario para almacenar info de genes por posición
        for line in tqdm(all_lines, desc="Procesando líneas", unit="pos"):
            rows, pos, gene, aa_idx_val, codon_pos = process_line((line, tokens, tp_num))
            if rows:
                all_rows.extend(rows)
                # Guardar información de genes para esta posición
                if pos not in genes_info:
                    genes_info[pos] = gene
                    genes_info[f"{pos}_aa"] = aa_idx_val
                    genes_info[f"{pos}_codon"] = codon_pos
        
        # Crear DataFrame completo y aplicar filtro global
        print("📊 Aplicando filtros de calidad...")
        full_df = pd.DataFrame(all_rows, columns=["pos","var","day","k","n"])
        # CORREGIDO: Filtrar por lecturas de la variante específica (k), no cobertura total (n)
        keep_global = (full_df.groupby(["pos","var"])["k"].max() >= depth_any)
        valid_variants = set(keep_global[keep_global].index)
        
        print(f"✅ Filtro global: {len(valid_variants)} variantes válidas de {len(keep_global)} totales")
        
        # Filtrar filas procesadas manteniendo solo variantes válidas
        filtered_rows = []
        for row in all_rows:
            if (row[0], row[1]) in valid_variants:  # (pos, var)
                filtered_rows.append(row)
        
        # Agrupar filas por posición para mantener nucleótidos/aminoácidos juntos
        pos_groups = defaultdict(list)
        for row in filtered_rows:
            pos = row[0]
            pos_groups[pos].append(row)
        
        # 2. Crear chunks balanceados por posición
        print("📦 Dividiendo datos en chunks para procesamiento paralelo...")
        positions = sorted(pos_groups.keys())
        chunk_size = max(1, len(positions) // n_jobs)
        chunks = []
        
        for i in range(0, len(positions), chunk_size):
            chunk_positions = positions[i:i+chunk_size]
            chunk_rows = []
            chunk_genes = {}
            for pos in chunk_positions:
                chunk_rows.extend(pos_groups[pos])
                # Incluir información de genes para este chunk
                if pos in genes_info:
                    chunk_genes[pos] = genes_info[pos]
                    chunk_genes[f"{pos}_aa"] = genes_info[f"{pos}_aa"]
                    chunk_genes[f"{pos}_codon"] = genes_info[f"{pos}_codon"]
            # Crear un chunk que contenga las filas procesadas y info de genes
            chunks.append((chunk_rows, tokens, tp_num, chunk_genes))
        
        print(f"📁 Creando {len(chunks)} archivos temporales...")
        
        # 3. Crear archivos físicos separados
        chunk_files = []
        for i, chunk in enumerate(tqdm(chunks, desc="Creando chunks", unit="chunk")):
            chunk_file = outdir / f"{pat}_long_chunk_{i:02d}.tsv"
            chunk_files.append(chunk_file)
            
            # Procesar chunk y guardar como archivo físico
            chunk_rows, chunk_gene_of, chunk_aaidx_of, chunk_codonpos_of = process_chunk(chunk)
            
            if chunk_rows:
                chunk_df = pd.DataFrame(chunk_rows, columns=["pos","var","day","k","n"])
                # El filtro de profundidad ya se aplicó globalmente
                
                # Guardar archivo físico
                chunk_df.to_csv(chunk_file, sep="\t", index=False)
            else:
                # Crear archivo vacío para mantener consistencia
                pd.DataFrame(columns=["pos","var","day","k","n"]).to_csv(chunk_file, sep="\t", index=False)
        
        # 4. Procesar archivos en paralelo (análisis de permutaciones)
        print(f"🔬 Analizando variantes con {n_jobs} cores...")
        
        def process_chunk_file(chunk_file_path):
            """Procesa un archivo _long_chunk y devuelve estadísticas"""
            if not chunk_file_path.exists():
                return None, {}, {}, {}
            
            chunk_df = pd.read_csv(chunk_file_path, sep="\t")
            if chunk_df.empty:
                return None, {}, {}, {}
            
            # Reconstruir mapeos desde el chunk original
            chunk_idx = int(chunk_file_path.stem.split("_chunk_")[1])
            chunk_rows, tokens, tp_num, chunk_genes = chunks[chunk_idx]
            
            gene_of = {}; aaidx_of = {}; codonpos_of = {}
            # Usar los mapeos de genes del chunk
            for row in chunk_rows:
                pos = row[0]
                if pos not in gene_of:
                    gene_of[pos] = chunk_genes.get(pos, "unknown")
                    aaidx_of[pos] = chunk_genes.get(f"{pos}_aa", "ND")
                    codonpos_of[pos] = chunk_genes.get(f"{pos}_codon", None)
            
            return chunk_df, gene_of, aaidx_of, codonpos_of
        
        def analyse_chunk_variants(chunk_file_path, nperm, freq_delta, thr_clear, eps_tie, chunk_id):
            """Analiza todas las variantes de un chunk secuencialmente"""
            chunk_df, gene_of, aaidx_of, codonpos_of = process_chunk_file(chunk_file_path)
            if chunk_df is None or chunk_df.empty:
                return []
            
            # Construir resúmenes de frecuencias para este chunk
            _, fmean_pos = build_freq_summaries(chunk_df)
            
            # Agrupar variantes del chunk
            groups = list(chunk_df.groupby(["pos","var"]))
            
            # Analizar variantes del chunk secuencialmente (sin paralelización interna)
            stats = []
            for g in tqdm(groups, desc=f"Core {chunk_id}", unit="var", position=chunk_id, leave=False):
                result = analyse_variant(g, nperm, gene_of, aaidx_of, codonpos_of,
                                       chunk_df, freq_delta, fmean_pos, thr_clear, eps_tie)
                if result:
                    stats.append(result)
            
            return stats
        
        # Procesar archivos en paralelo (análisis de variantes)
        print(f"🔬 Procesando {len(chunk_files)} archivos en paralelo con {n_jobs} cores...")
        all_stats = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(analyse_chunk_variants)(chunk_file, nperm, freq_delta, thr_clear, eps_tie, i) 
            for i, chunk_file in enumerate(chunk_files)
        )
        
        # Aplanar lista de estadísticas
        stats = []
        for chunk_stats in all_stats:
            stats.extend(chunk_stats)
        
        # También necesitamos los mapeos para el reporte final
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(process_chunk_file)(chunk_file) for chunk_file in chunk_files
        )
        
        # 5. Unir resultados de todos los archivos
        print("🔗 Uniendo resultados de todos los cores...")
        all_dfs = []
        gene_of = {}; aaidx_of = {}; codonpos_of = {}
        
        for chunk_df, chunk_gene_of, chunk_aaidx_of, chunk_codonpos_of in results:
            if chunk_df is not None and not chunk_df.empty:
                all_dfs.append(chunk_df)
                gene_of.update(chunk_gene_of)
                aaidx_of.update(chunk_aaidx_of)
                codonpos_of.update(chunk_codonpos_of)
        
        if all_dfs:
            long = pd.concat(all_dfs, ignore_index=True)
        else:
            long = pd.DataFrame(columns=["pos","var","day","k","n"])
        
        # 6. Limpiar archivos temporales
        print("🧹 Limpiando archivos temporales...")
        for chunk_file in chunk_files:
            if chunk_file.exists():
                chunk_file.unlink()
        
        print(f"✅ Análisis completado: {len(long)} filas finales")
        return long, gene_of, aaidx_of, codonpos_of, stats

# — pendiente ponderada ————————————————————————————————
def beta_only(day,k,n):
    w = n
    W  = np.sum(w)
    x̄ = np.sum(w*day)   / W
    ȳ = np.sum(w*(k/n)) / W
    Sxx = np.sum(w*(day-x̄)**2)
    Sxy = np.sum(w*(day-x̄)*((k/n)-ȳ))
    beta = Sxy / Sxx if Sxx>0 else 0.
    return beta

def pval_perm(day,k,n,beta_real,nperm):
    """Versión secuencial simple para permutaciones"""
    greater = 0
    for _ in range(nperm):
        beta_star = beta_only(np.random.permutation(day),k,n)
        if abs(beta_star) >= abs(beta_real):
            greater += 1
    return (greater+1)/(nperm+1)

# ─────────────────── Frecuencias medias y anotación de codones/AA ──────────────
def build_freq_summaries(long_df: pd.DataFrame):
    """
    Devuelve:
      fmean_posvar: DataFrame con columnas [pos, var, freq_mean]
      fmean_pos:    dict {pos -> list[(var, freq_mean_desc)]}
    """
    tmp = long_df.assign(freq = long_df["k"] / long_df["n"])
    fmean_posvar = (tmp.groupby(["pos","var"])["freq"]
                       .mean()
                       .reset_index(name="freq_mean"))
    # dict rápido: pos -> [(var, freq_mean) ...] ordenado desc
    fmean_pos = {}
    for p, sub in fmean_posvar.groupby("pos"):
        lst = list(zip(sub["var"].values, sub["freq_mean"].values))
        lst.sort(key=lambda t: t[1], reverse=True)
        fmean_pos[p] = lst
    return fmean_posvar, fmean_pos

def infer_codon_aa_for_recta(pos: int,
                             var: str,
                             fmean_pos: dict,
                             codonpos_of: dict,
                             aaidx_of: dict,
                             thr_clear: float,
                             eps_tie: float,
                             selected_positions=None):
    """
    Construye el codón (y AA) de *esta recta (pos,var)* usando TODAS las rectas:
      - En 'pos' se fuerza la base = var.
      - En las otras posiciones del triplete se elige la base "clara"
        (freq >= thr_clear) y, además, suficientemente distinta de la frecuencia
        de la *base principal de la POSICIÓN* 'pos' (top de la posición):
        |f1 - f_pos_top| > eps_tie.
      - Si alguna de esas dos posiciones del triplete tiene rectas seleccionadas
        (tag ∈ {pos,neg}), se marca como ambiguo (vía `selected_positions`).
    Devuelve (codon, aa_label, ambiguous, detail_str).
    """
    # Manejar None como set vacío para evitar TypeError en "p in selected_positions"
    if selected_positions is None:
        selected_positions = set()
    
    if len(var) != 1:   # indel
        # Si la propia variante es un indel, marcamos explícitamente
        # la ambigüedad por indel.
        return None, "amb (indel)", True, "indel"
    if var == "N":
        # La base forzada en esta recta es desconocida → no podemos
        # determinar el codón/AA de manera fiable.
        # Marcamos ambigüedad explícita (no indel).
        return None, "amb", True, "forced_N"

    codon_pos = codonpos_of.get(pos)
    if codon_pos is None:
        return None, "amb", True, "no_codon_pos"

    codon_start = pos - (codon_pos - 1)
    triple = [codon_start, codon_start+1, codon_start+2]

    # base principal (top) de la POSICIÓN actual (no de la variante concreta)
    lst_here = fmean_pos.get(pos, [])
    if not lst_here:
        return None, "amb", True, "no_fmean_pos"
    v_pos_top, f_pos_top = lst_here[0]
    # normalizar parámetro opcional de posiciones seleccionadas
    if selected_positions is None:
        selected_positions = set()

    bases = []
    ambiguous_due_to_indel = False
    ambiguous = False
    detail_lines = []

    for p in triple:
        if p == pos:
            bases.append(var)
            detail_lines.append(f"{p}:FORCED={var};pos_top=({v_pos_top},{f_pos_top:.3f})")
            continue

        lst = fmean_pos.get(p)
        if not lst:
            bases.append("N")
            ambiguous = True
            p_selected = (p in selected_positions)
            detail_lines.append(f"{p}:NOINFO;pos_selected={p_selected}")
            continue

        # ── ignorar indels como “base clara” en posiciones vecinas ─────────
        lst_nt = [(v, f) for (v, f) in lst if len(v) == 1]   # A/C/G/T/N
        if not lst_nt:
            bases.append("N")
            ambiguous = True
            ambiguous_due_to_indel = True
            p_selected = (p in selected_positions)
            detail_lines.append(
                f"{p}:NO_NT_CANDIDATE(indel_only);pos_selected={p_selected}"
            )
            continue
        v1, f1 = lst_nt[0]
        f2 = lst_nt[1][1] if len(lst_nt) > 1 else 0.0

        # comparación frente a la *base principal de la POSICIÓN 'pos'*
        diff_vs_pos_top = abs(f1 - f_pos_top)
        ok_vs_pos_main = (diff_vs_pos_top > eps_tie) #Se calcula pero no se usa

        # ¿alguna recta seleccionada en esta OTRA posición?
        p_selected = (p in selected_positions)

        detail_lines.append(
            f"{p}:top=({v1},{f1:.3f}) second={f2:.3f} "
            f"pos_top=({v_pos_top},{f_pos_top:.3f}) |f1-f_pos_top|={diff_vs_pos_top:.3f} "
            f"pos_selected={p_selected}"
        )

        if (f1 < thr_clear) or p_selected:
            bases.append("N")
            ambiguous = True
        else:
            bases.append(v1)

    codon = "".join(bases)
    if ambiguous:
        # Etiqueta específica si la ambigüedad proviene de indels
        aa_label = "amb (indel)" if ambiguous_due_to_indel else "amb"
        return codon, aa_label, True, ";".join(detail_lines)
    aa = CODON2AA.get(codon, "?") if not ambiguous else "amb"
    # Devolver solo el aminoácido, sin el índice
    # El formateo completo se hará en _format_aa_change
    aa_label = aa

    return codon, aa_label, ambiguous, ";".join(detail_lines)

def analyse_variant(group, nperm, gene_of, aaidx_of, codonpos_of,
                    long_df, freq_delta, fmean_pos, thr_clear, eps_tie):
    (pos,var), grp = group
    if grp["day"].nunique() < 3:
        return None
    day, k, n = grp["day"].values, grp["k"].values, grp["n"].values
    freq = k/n
    # primero estimamos la pendiente
    beta = beta_only(day, k, n)
    # Δf por «fit»: pendiente × (t_max - t_min)
    delta_f = beta * (day.max() - day.min())
    if abs(delta_f) < freq_delta:
        return None

    p = pval_perm(day,k,n,beta,nperm)
    tag   = ("pos" if (p<0.05 and beta>0)
             else "neg" if (p<0.05 and beta<0) else "")
    score = math.copysign(abs(beta)*(-math.log10(p)), beta)

    # ─── anotación nt/aa por frecuencias globales ─────
    codon, aa_label, ambiguous, detail = infer_codon_aa_for_recta(
        pos, var, fmean_pos, codonpos_of, aaidx_of,
        thr_clear=thr_clear, eps_tie=eps_tie, selected_positions=None)

    # Formatear aa_change según el signo del beta
    aa_idx = aaidx_of[pos]
    is_indel = (len(var) != 1) or (aa_label == "amb (indel)")
    
    # Extraer solo el aminoácido del aa_label (sin el índice)
    if aa_label == "amb (indel)":
        aa = "amb"
    elif aa_label == "amb":
        aa = "amb"
    else:
        # Extraer solo la parte del aminoácido (sin el índice)
        aa = aa_label.rstrip('0123456789') if aa_label else "amb"
    
    aa_change = _format_aa_change(aa, aa_idx, beta, is_indel)

    return dict(gene=gene_of[pos], aa_idx=aaidx_of[pos], pos=pos, var=var,
                beta=beta, pval=p, score=score, tag=tag,
                delta_f=float(delta_f),
                codon=codon,
                aa_change=aa_change, ambiguous=bool(ambiguous),
                annot_detail=detail,
                depth=float(np.mean(n)), n_days=int(grp["day"].nunique()),
                day=day.tolist(),            # ← listas "puras"
                k=k.tolist(),
                n=n.tolist())

# ───────────────────── resumen por posición ─────────────────────────────
def _format_aa_change(aa: str, aa_idx: str, beta: float, is_indel: bool = False) -> str:
    """
    Formatea aa_change según el signo del beta:
    - Beta negativo: cambio+pos (ej: K242)
    - Beta positivo: pos+cambio (ej: 242K)
    - Ambiguo: amb+pos o pos+amb
    - Indel con beta < 0: (indel)amb242 o (indel)K242
    - Indel con beta > 0: 242amb(indel) o 242K(indel)
    """
    if aa is None or aa == "":
        aa = "amb"
    
    # Determinar si es indel
    if is_indel or aa == "amb (indel)":
        aa = "amb"
        indel_suffix = "(indel)"
    else:
        indel_suffix = ""
    
    # Si no hay aa_idx o no es numérico, usar "ND"
    if not aa_idx or not str(aa_idx).isdigit():
        aa_idx = "ND"
    
    # Formatear según el signo del beta
    if beta < 0:  # Beta negativo: cambio+pos
        if aa == "amb":
            return f"{indel_suffix}amb{aa_idx}" if indel_suffix else f"amb{aa_idx}"
        else:
            return f"{indel_suffix}{aa}{aa_idx}" if indel_suffix else f"{aa}{aa_idx}"
    else:  # Beta positivo: pos+cambio
        if aa == "amb":
            return f"{aa_idx}amb{indel_suffix}"
        else:
            return f"{aa_idx}{aa}{indel_suffix}"

def _split_aa_label(lbl: str) -> Tuple[str, str]:
    """'K100' -> ('K','100'); '100K' -> ('K','100'); 'ambiguous' -> ('amb','')"""
    # Considera "ambiguous" y variantes como "ambiguous (indel)" -> convertir a "amb"
    if lbl is None or (isinstance(lbl, str) and lbl.startswith("ambiguous")):
        return "amb", ""
    
    # Casos especiales para ambND, NDamb, y amb seguido de números
    if lbl == "ambND":
        return "amb", "ND"
    elif lbl == "NDamb":
        return "amb", "ND"
    elif lbl.startswith("amb") and lbl[3:].isdigit():
        return "amb", lbl[3:]
    
    m = AA_RE.match(lbl)
    if not m:
        return lbl, ""
    # Verificar si es formato aminoácido+número o número+aminoácido
    if m.group(1) and m.group(2):  # aminoácido+número (ej: D5152)
        return m.group(1), m.group(2)
    elif m.group(3) and m.group(4):  # número+aminoácido (ej: 5152D)
        return m.group(4), m.group(3)
    elif m.group(5) and m.group(6):  # número+amb (ej: 5624amb)
        return m.group(6), m.group(5)
    elif m.group(7) and m.group(8):  # amb+número (ej: amb614)
        return m.group(7), m.group(8)
    else:
        return lbl, ""

def summarise_positions(df_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Resumen por posición, guardando TODAS las rectas seleccionadas (pos/neg) como JSON:
      - vars, tags, betas, pvals, codons, aa_labels, ambiguous, annot_detail
    Además:
      - label:  D100K (↓ primero, ↑ después) si hay una subida y una bajada; o K100 si solo hay una.
      - synonymous: True/False si se puede decidir; False por defecto si no es determinable.
      - why_ambiguous_pos: concatenación de los detalles de anotación de las rectas ambiguas.
      - beta_up/p_up y beta_dn/p_dn: de las mejores rectas (mayor |score|) por signo.
      - score_up/score_dn: el score de esas “mejores” rectas.
    """
    sel = df_stats[df_stats["tag"].isin(["pos","neg"])]
    if sel.empty:
        return pd.DataFrame(columns=[
            "pos","aa_idx","label","synonymous","n_sel",
            "vars","tags","betas","pvals","codons","aa_labels","ambiguous","annot_detail",
            "why_ambiguous_pos","beta_up","p_up","beta_dn","p_dn","score_up","score_dn"
        ])

    rows = []

    for pos, sub in sel.groupby("pos"):
        # Listas JSON con TODAS las rectas seleccionadas
        vars_l        = sub["var"].tolist()
        tags_l        = sub["tag"].tolist()
        betas_l       = sub["beta"].tolist()
        pvals_l       = sub["pval"].tolist()
        codons_l      = sub["codon"].tolist() if "codon" in sub else [None]*len(sub)
        aa_labels_l   = sub["aa_change"].tolist() if "aa_change" in sub else [None]*len(sub)
        ambiguous_l   = sub["ambiguous"].tolist() if "ambiguous" in sub else [False]*len(sub)
        details_l     = sub["annot_detail"].tolist() if "annot_detail" in sub else [""]*len(sub)

        why_amb = [d for a, d in zip(ambiguous_l, details_l) if a]

        # top por score en cada signo
        pos_sel = sub[sub["tag"]=="pos"].sort_values("score", ascending=False)
        neg_sel = sub[sub["tag"]=="neg"].sort_values("score", ascending=False)

        beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        label = "amb"
        synonymous = None
        
        # Obtener aa_idx directamente desde la columna aa_idx del DataFrame
        aa_idx = sub["aa_idx"].iloc[0] if not sub.empty else ""

        # Verificar si hay múltiples rectas del mismo signo (debería ser ambiguo)
        multiple_pos = len(pos_sel) > 1
        multiple_neg = len(neg_sel) > 1
        
        if multiple_pos and not neg_sel.empty:
            # Múltiples subidas + una bajada → aminoácido de bajada{aa_idx}amb
            dn = neg_sel.iloc[0]
            aa_d, idx_d = _split_aa_label(dn["aa_change"])
            if aa_d == "amb":
                label = f"amb{aa_idx}amb" if aa_idx else "ambamb"
            else:
                # Usar el aminoácido extraído directamente (ya sin índice)
                label = f"{aa_d}{aa_idx}amb" if aa_idx else f"{aa_d}amb"
            synonymous = None
            beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        elif multiple_neg and not pos_sel.empty:
            # Múltiples bajadas + una subida → amb{aa_idx}aminoácido de subida
            up = pos_sel.iloc[0]
            aa_u, idx_u = _split_aa_label(up["aa_change"])
            if aa_u == "amb":
                label = f"amb{aa_idx}amb" if aa_idx else "ambamb"
            else:
                # Usar el aminoácido extraído directamente (ya sin índice)
                label = f"amb{aa_idx}{aa_u}" if aa_idx else f"amb{aa_u}"
            synonymous = None
            beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        elif multiple_pos and neg_sel.empty:
            # Solo múltiples subidas → {aa_idx}amb
            label = f"{aa_idx}amb" if aa_idx else "amb"
            synonymous = None
            beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        elif multiple_neg and pos_sel.empty:
            # Solo múltiples bajadas → amb{aa_idx}
            label = f"amb{aa_idx}" if aa_idx else "amb"
            synonymous = None
            beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        elif multiple_pos and multiple_neg:
            # Múltiples de ambos signos → amb{aa_idx}amb
            label = f"amb{aa_idx}amb" if aa_idx else "ambamb"
            synonymous = None
            beta_up = p_up = beta_dn = p_dn = score_up = score_dn = None
        elif not pos_sel.empty and not neg_sel.empty and not multiple_pos and not multiple_neg:
            up  = pos_sel.iloc[0]
            dn  = neg_sel.iloc[0]
            aa_u, idx_u = _split_aa_label(up["aa_change"])
            aa_d, idx_d = _split_aa_label(dn["aa_change"])
            # Usar aa_idx del DataFrame si no se puede extraer desde aa_change
            if not aa_idx and (idx_d or idx_u):
                aa_idx = idx_d or idx_u
            # Construir etiqueta: ↓ primero, ↑ después
            if aa_d == "amb" and aa_u == "amb":
                label = f"amb{aa_idx}amb" if aa_idx else "ambamb"
            else:
                # Extraer solo los aminoácidos sin los índices
                aa_d_only = aa_d.rstrip('0123456789') if aa_d else "amb"
                aa_u_only = aa_u.rstrip('0123456789') if aa_u else "amb"
                
                # Mostrar ambas direcciones: ↓ primero, ↑ después
                label = f"{aa_d_only}{aa_idx}{aa_u_only}"  # ↓ primero, ↑ después
            # Etiquetado sinónimo como booleano; None si no se puede decidir
            if aa_d not in ("amb","?") and aa_u not in ("amb","?"):
                synonymous = (aa_d == aa_u)   # True si AA iguales, False si distintas
            else:
                synonymous = None
            beta_up, p_up, score_up = up["beta"], up["pval"], up["score"]
            beta_dn, p_dn, score_dn = dn["beta"], dn["pval"], dn["score"]
        else:
            # sólo una recta seleccionada
            r = pos_sel.iloc[0] if not pos_sel.empty else neg_sel.iloc[0]
            aa, idx = _split_aa_label(r["aa_change"])
            # Usar aa_idx del DataFrame si no se puede extraer desde aa_change
            if not aa_idx and idx:
                aa_idx = idx
            if aa != "amb":
                # Usar el aminoácido extraído directamente (ya sin índice)
                # Si es subida ⇒ idx primero (100K); si es bajada ⇒ AA primero (K100)
                label = f"{aa_idx}{aa}" if r["tag"] == "pos" else f"{aa}{aa_idx}"
            else:
                # Si es ambiguo pero tenemos aa_idx, incluirlo
                if aa_idx:
                    label = f"{aa_idx}amb" if r["tag"] == "pos" else f"amb{aa_idx}"
                else:
                    # Si no hay aa_idx, usar ND y respetar la regla de negativo/positivo
                    label = "NDamb" if r["tag"] == "pos" else "ambND"
            if r["tag"] == "pos":
                beta_up, p_up, score_up = r["beta"], r["pval"], r["score"]
            else:
                beta_dn, p_dn, score_dn = r["beta"], r["pval"], r["score"]

        rows.append(dict(
            pos=pos,
            aa_idx=aa_idx,
            label=label,
            synonymous=synonymous,
            n_sel=len(sub),
            vars=json.dumps(vars_l),
            tags=json.dumps(tags_l),
            betas=json.dumps(betas_l),
            pvals=json.dumps(pvals_l),
            codons=json.dumps(codons_l),
            aa_labels=json.dumps(aa_labels_l),
            ambiguous=json.dumps(ambiguous_l),
            annot_detail=json.dumps(details_l),
            why_ambiguous_pos=" | ".join(why_amb),
            beta_up=beta_up,  p_up=p_up,
            beta_dn=beta_dn,  p_dn=p_dn,
            score_up=score_up, score_dn=score_dn,
        ))

    return pd.DataFrame(rows)

# ───────────────────── figura única por posición ────────────────────────
def plot_position(patient: str,
                  pos: int,
                  df_stats: pd.DataFrame,
                  long_df: pd.DataFrame,
                  beta_mean_sel: float,
                  beta_mean_all: float,
                  pos_row: pd.Series,
                  plot_dir: pathlib.Path,
                  gene_of: dict = None,
                  aaidx_of: dict = None):
    """
    Una figura por posición con alguna recta seleccionada:
      - rectas seleccionadas en rojo (β>0) / azul (β<0)
      - resto en gris
    """
    mates = long_df[ long_df["pos"] == pos ]
    if mates.empty:
        return

    # mapas rápidos
    tag_map  = {(r["pos"], r["var"]): r["tag"]  for _, r in df_stats.iterrows()}

    fig, ax = plt.subplots(figsize=(5,3.2))

    # Calcular escalas para k (punto central) y n (círculo exterior)
    k_global_max = float(mates["k"].max()) if not mates.empty else 0.0
    n_global_max = float(mates["n"].max()) if not mates.empty else 0.0

    for var, sub in mates.groupby("var"):
        d  = sub["day"].values
        fr = sub["k"].values / sub["n"].values
        b  = beta_only(d, sub["k"].values, sub["n"].values)

        tag = tag_map.get((pos,var), "")
        if tag == "pos":
            color, ec = "red", "k"
            z = 3
        elif tag == "neg":
            color, ec = "blue", "k"
            z = 3
        else:
            color, ec = "lightgrey", "lightgrey"
            z = 1

        # Calcular tamaños para k (punto central) y n (círculo exterior)
        if k_global_max > 0:
            s_k = 20.0 + 180.0 * (sub["k"].values.astype(float) / k_global_max)
        else:
            s_k = np.full_like(fr, 100.0, dtype=float)
        
        if n_global_max > 0:
            s_n = 30.0 + 200.0 * (sub["n"].values.astype(float) / n_global_max)
        else:
            s_n = np.full_like(fr, 130.0, dtype=float)

        # Dibujar círculos exteriores (n) primero - mismo color que k pero más claro
        ax.scatter(d, fr*100, s=s_n, c=color, alpha=0.2, edgecolors=color, linewidth=0.5, zorder=z-1)
        
        # Dibujar puntos centrales (k) encima
        ax.scatter(d, fr*100, s=s_k, c=color, edgecolors=ec,
                   linewidths=0.25, alpha=.7 if tag else .3, zorder=z)
        # recta
        xs = np.array(sorted(d))
        w   = sub["n"].values
        xw  = np.average(d,  weights=w)
        yw  = np.average(fr, weights=w)
        yfit = yw + b*(xs - xw)
        ax.plot(xs, yfit*100, color=color if tag else "grey",
                lw=1.0 if tag else .8, ls="-" if tag else "--", alpha=.8 if tag else .3, zorder=z)

    # líneas de neutralidad
    xspan = np.array([mates["day"].min(), mates["day"].max()])
    freq_all = (mates["k"].sum() / mates["n"].sum())
    # all
    ax.plot(xspan,
            (freq_all + beta_mean_all*(xspan-xspan.mean()))*100,
            color="grey", ls="--", lw=.8, alpha=.6, label=r"$\bar{\beta}_{all}$")
    # sel
    if beta_mean_sel is not None:
        ax.plot(xspan,
                (freq_all + beta_mean_sel*(xspan-xspan.mean()))*100,
                color="k", ls="--", lw=.8, alpha=.6, label=r"$\bar{\beta}_{sel}$")

    ax.set_xlabel("Day"); ax.set_ylabel("Frequency (%)")
    ax.grid(True, ls=":", alpha=0.3)

    # título con el resumen pos_row (si lo tenemos)
    if pos_row is not None:
        def _get(d, keys):
            for k in keys:
                if k in d and d[k] is not None:
                    s = str(d[k]).strip()
                    if s != "" and s.lower() not in {"nan", "none"}:
                        return d[k]
            return None

        def _safe_aa(x):
            if x is None:
                return "amb"
            s = str(x).strip()
            # Unificar todas las etiquetas de ambigüedad a "amb"
            if s == "" or s.lower() in {"nan", "na", "none", "nd", "ambiguous"} or s.startswith("ambiguous"):
                return "amb"
            return s

        def _safe_idx(x):
            try:
                return str(int(float(x)))
            except Exception:
                # Si faltase el índice (no debería), mostramos 'ND'
                return "ND"

        # Obtener aa_idx desde los datos originales si está disponible
        if aaidx_of is not None and pos in aaidx_of:
            aa_idx = aaidx_of[pos]
        else:
            aa_idx = _get(pos_row, ("aa_idx", "aa_pos", "aa_position", "aa_index"))
        
        # Obtener información de aminoácidos desde df_stats para esta posición
        pos_stats = df_stats[df_stats["pos"] == pos]
        if not pos_stats.empty:
            # Verificar si hay indels seleccionados en esta posición
            # Indels son variantes que empiezan con + (inserción) o - (deleción)
            indel_stats = pos_stats[pos_stats["var"].astype(str).str.startswith(('+', '-'))]
            has_selected_indel = not indel_stats.empty and any(tag in ["pos", "neg"] for tag in indel_stats["tag"].tolist())
            
            if has_selected_indel:
                # Si hay indels seleccionados, usar formato amb con (indel)
                indel_beta = indel_stats.iloc[0]["beta"]  # Usar el primer indel seleccionado
                if indel_beta < 0:  # Beta negativo: (indel)amb + pos
                    change_lbl = f"(indel)amb{aa_idx}"
                else:  # Beta positivo: pos + amb(indel)
                    change_lbl = f"{aa_idx}amb(indel)"
            else:
                # Usar la etiqueta de aminoácido del resumen si está disponible
                if pos_row is not None and "label" in pos_row and not pos_row["label"].startswith("ambiguous"):
                    change_lbl = pos_row["label"]
                else:
                    # Fallback: usar aa_change de las estadísticas
                    aa_changes = pos_stats["aa_change"].tolist()
                    if aa_changes and not aa_changes[0].startswith("ambiguous"):
                        change_lbl = aa_changes[0]
                    else:
                        change_lbl = f"amb{aa_idx}amb"
        else:
            change_lbl = f"amb{aa_idx}amb"

        syn_flag = pos_row.get("synonymous", None)
        syn_txt = ("syn" if syn_flag is True else
                   "nonsyn" if syn_flag is False and syn_flag is not None else
                   "amb")
        n_sel = pos_row.get('n_sel', 1)
        ttl = f"{patient} · pos {pos} · {change_lbl} [{syn_txt}]  (n_sel={n_sel})"
    else:
        ttl = f"{patient} · pos {pos}"
    
    # Agregar información del gen al principio del título
    if gene_of is not None and pos in gene_of:
        gene_info = gene_of[pos]
    else:
        # Fallback: intentar obtener desde df_stats
        gene_info = df_stats[df_stats["pos"] == pos]["gene"].iloc[0] if not df_stats[df_stats["pos"] == pos].empty else "unknown"
    ttl = f"{gene_info} · {ttl}"
    ax.set_title(ttl, fontsize=9)
    fig.tight_layout()

    out = plot_dir / f"pos{pos:05d}_ALL.pdf"
    fig.savefig(out)
    plt.close(fig)

def plot_genome(patient: str,
                df_stats: pd.DataFrame,
                long_df: pd.DataFrame,
                outdir: pathlib.Path,
                beta_mean_sel: float,
                beta_mean_all: float,
                tsv_path: pathlib.Path) -> None:
    """
    Resumen genómico:
      • Puntos seleccionados (pos/neg).
      • Curvas de densidad independientes    β>0 (↑ rojo)   β<0 (↓ azul)
        – se dibujan y rellenan en el eje derecho, ocupando todo su alto
      • Dos líneas horizontales de referencia: β̄_all y β̄_sel.
      • Barra de genes coloreada justo ENCIMA del eje-X (no tapa los datos)
      • Impresión en stdout de los rangos génicos
    """

    # ───────────── 1 · filtrar variantes seleccionadas  ─────────────────
    sel = df_stats[df_stats["tag"].isin(["pos", "neg"])].copy()
    if sel.empty:
        print(f"[WARN] {patient}: no hay variantes etiquetadas; omito figura")
        return

    sel["color"] = sel["tag"].map({"pos": "red", "neg": "blue"}).values
    # --- Usar el *_genes.tsv para fijar el rango genómico completo ---
    genes_tsv = tsv_path
    if not genes_tsv.is_file():
        genes_tsv = (outdir.parent / f"{patient}_genes.tsv")
        if not genes_tsv.is_file():
            genes_tsv = outdir / f"{patient}_genes.tsv"

    if genes_tsv.is_file():
        genes_raw = pd.read_csv(genes_tsv, sep="\t",
                                usecols=["position", "gene"])
    else:
        # Fallback si no existe el *_genes.tsv
        genes_raw = (df_stats[["pos", "gene"]]
                     .rename(columns={"pos": "position"})
                     .drop_duplicates())

    genome_min = int(genes_raw["position"].min())
    genome_max = int(genes_raw["position"].max())
    x_full     = np.arange(genome_min, genome_max + 1)


    # ───────────── 2 · figura base + eje de densidad  ───────────────────
    fig, ax  = plt.subplots(figsize=(13, 5))
    ax2      = ax.twinx()                     # ← densidad
    ax.set_xlabel("Genomic position\n")
    ax.set_ylabel(r"$\beta$")
    ax2.set_ylabel("Relative density")

    # ───────────── 3 · Sólo puntos (sin barras de IC)
    ax.scatter(sel["pos"], sel["beta"], s=46, c=sel["color"],
               edgecolors='k', linewidth=.25, zorder=3, alpha=.9)

    # ───────────── 4 · curvas de densidad (sólo en eje derecho) ─────────
    # ───—─ 3.5 · fija margen inferior y posición futura de barra genes ──
    # Asegurar que y=0 esté dentro del rango y añadir margen simétrico
    ymin_d, ymax_d = ax.get_ylim()            # límites tras trazar puntos
    span = (ymax_d - ymin_d)
    if not np.isfinite(span) or span <= 0:
        span = 1.0
    pad = 0.05 * span                         # 5 % del rango visible
    ymin_new = min(ymin_d, 0.0) - pad
    ymax_new = max(ymax_d, 0.0) + pad
    ax.set_ylim(ymin_new, ymax_new)
    y_bar = ymin_new + 0.40 * pad             # barra ≈ 40 % del padding inferior

    dens_max = 0.0
    for tag, sign, col in (("pos", +1, "red"),
                           ("neg", -1, "blue")):
        sub = sel[sel["tag"] == tag]
        if sub.empty:
            continue

        dens = (sub.groupby("pos").size()
                    .groupby(level=0).sum()          # evita duplicados
                    .reindex(x_full, fill_value=0)
                    .astype(float)
                    .sort_index())
        smooth = gaussian_filter1d(dens.values, sigma=350)
        dens_max = max(dens_max, smooth.max())    # ← guarda máximo global


        smooth_signed = sign * smooth               # arriba / abajo
        # relleno + línea en eje derecho
        ax2.fill_between(x_full, 0, smooth_signed,
                         color=col, alpha=.20, lw=0, zorder=1)
        ax2.plot(x_full, smooth_signed,
                 color=col, lw=1.1, alpha=.7, zorder=1)

    # ----- ajusta ax₂ para que y = 0 quede EXACTAMENTE alineado -----
    # fracción de la altura del eje izquierdo donde cae y=0
    frac0  = (0.0 - ymin_new) / (ymax_new - ymin_new)
    # Evitar extremos 0 o 1 por estabilidad numérica
    frac0 = min(max(frac0, 1e-6), 1 - 1e-6)
    upper2 =  1.05 * dens_max
    lower2 = -upper2 * frac0 / (1 - frac0)             # calcula para alinear y=0
    ax2.set_ylim(lower2, upper2)

    # Ticks cada 25 % del rango real del eje derecho
    yticks = np.linspace(lower2, upper2, 5)
    ax2.set_yticks(yticks)

    def fmt_tick(val):
        return "0" if abs(val) < 1e-100 else f"{val:.3g}"

    ax2.set_yticklabels([fmt_tick(t) for t in yticks])
    # Ocultar valores (numeración) del eje Y derecho pero mantener el eje
    ax2.tick_params(right=True, labelright=False)

    # asegurar coincidencia visual del 0 en ambos ejes (límites ya fijados)

    # ───────────── 5 · líneas de neutralidad β̄ ─────────────────────────
    #   · todas las variantes
    ax.axhline(beta_mean_all, color="grey", ls="--", lw=.9,
               label=r"$\bar{\beta}_{all}$")
    #   · solo seleccionadas (si existen)
    if beta_mean_sel is not None:
        ax.axhline(beta_mean_sel, color="k", ls="--", lw=.9,
                   label=r"$\bar{\beta}_{sel}$")

    # ───────────── 6 · ─────────────────────────────────────────────────
    genes_tbl = genes_raw.copy()
    genes_tbl["position"] = pd.to_numeric(genes_tbl["position"], errors="coerce")
    genes_tbl = genes_tbl.dropna(subset=["position"])
    genes_tbl["position"] = genes_tbl["position"].astype(int)
    genes_tbl = genes_tbl.sort_values("position")

    def _split_genes(val):
        if pd.isna(val):
            return []
        s = str(val).strip()
        if s == "":
            return []
        return [g.strip() for g in s.split(",") if g.strip() != ""]

    genes_tbl["gene_list"] = genes_tbl["gene"].map(_split_genes)
    genes_by_pos = (genes_tbl.groupby("position")["gene_list"]
                    .agg(lambda lists: sorted({g for lst in lists for g in lst})))

    pos_list = genes_by_pos.index.to_list()
    pos_genes = genes_by_pos.to_list()

    gene_order = []
    seen = set()
    for glist in pos_genes:
        for g in glist:
            if g not in seen:
                seen.add(g)
                gene_order.append(g)

    spans_rows = []
    active = {}
    prev_pos = None
    for pos, glist in zip(pos_list, pos_genes):
        current = set(glist)
        if prev_pos is not None:
            for g in list(active.keys()):
                if g not in current:
                    spans_rows.append((g, active[g], prev_pos))
                    del active[g]
        for g in current:
            if g not in active:
                active[g] = pos
        prev_pos = pos
    for g, start_pos in active.items():
        spans_rows.append((g, start_pos, prev_pos))

    spans = (pd.DataFrame(spans_rows, columns=["gene", "min", "max"])
             .sort_values("min"))

    # ── barra fina pegada justo ENCIMA del eje-X dentro del mismo ax ───────
    cmap = plt.cm.get_cmap("tab20")
    gene2color = {g: cmap(i % cmap.N) for i, g in enumerate(gene_order)}

    print("\n[GENE RANGES]")
    for g, mn, mx in spans[["gene", "min", "max"]].itertuples(index=False):
        ax.hlines(y_bar, mn, mx, lw=8,
                  color=gene2color[g], clip_on=False, zorder=4)
        print(f"  {g:10s} : {mn:>6d}–{mx:>6d}")

    # ───────────── 7 · leyenda (una línea, con marco) ──────────────────
    handles = [Patch(facecolor=gene2color[g], edgecolor='none')
               for g in gene_order]
    legend = ax.legend(handles,
                       gene_order,
                       title="Genes",
                       loc="upper center",
                       bbox_to_anchor=(0.5, -0.20),   # debajo del gráfico
                       ncol=len(gene_order),          # una única fila
                       fontsize=8,
                       title_fontsize=9,
                       frameon=True)
    legend.get_frame().set_linewidth(0.6)
    legend.get_frame().set_edgecolor('0.4')

    # ───────────── 9 · acabado y guardado  ─────────────────────────────
    ax.grid(True, ls=":", alpha=.3)
    ax.set_title(f"Genomic selection – {patient}")
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])   # deja hueco barra + leyenda
    fig.savefig(outdir / f"{patient}_genome_selection.pdf")
    plt.close(fig)

# ───────────── pipeline ────────────────────────────────────────
def process(tsv_path, outdir, nperm, n_jobs, resume, freq_delta,
            thr_clear, eps_tie, times_map=None, min_first_k=10,
            covar_enable=True, covar_method="beta_diff", covar_threshold=0.01,
            covar_min_days=3, covar_min_range=0.05, covar_max_day_diff=1):
    pat = tsv_path.stem.split("_genes")[0]
    stats_out = outdir/f"{pat}_variant_stats.tsv"
    if resume and stats_out.exists():
        print(f"⏭️  {pat}: ya procesado (skip)"); return
    print(f"\n🧬 Procesando paciente: {pat}")
    print("=" * 50)

    # Si n_jobs > 1, make_long ya hace el análisis de variantes en paralelo
    if n_jobs > 1:
        long, gene_of, aaidx_of, codonpos_of, stats = make_long(tsv_path, depth_any=100, n_jobs=n_jobs, 
                                                               nperm=nperm, freq_delta=freq_delta,
                                                               thr_clear=thr_clear, eps_tie=eps_tie)
        long.to_csv(outdir/f"{pat}_long.tsv",sep="\t",index=False)
        
        if not stats:
            print("   (nada pasa filtros)"); return
        df_full = pd.DataFrame(stats)
    else:
        # Procesamiento secuencial
        long,gene_of,aaidx_of,codonpos_of = make_long(tsv_path, depth_any=100, n_jobs=n_jobs,
                                                       nperm=nperm, freq_delta=freq_delta,
                                                       thr_clear=thr_clear, eps_tie=eps_tie)
        long.to_csv(outdir/f"{pat}_long.tsv",sep="\t",index=False)

        # --- sumarizamos frecuencias medias para TODAS las rectas (pos,var)
        _, fmean_pos = build_freq_summaries(long)

        groups = list(long.groupby(["pos","var"]))
        print(f"🔬 Analizando {len(groups)} variantes con {n_jobs} cores...")
        work = Parallel(n_jobs=n_jobs,backend="loky")(
            delayed(analyse_variant)(
                        g, nperm, gene_of, aaidx_of, codonpos_of,
                        long, freq_delta, fmean_pos, thr_clear, eps_tie)
                for g in tqdm(groups,desc=f"Analizando {pat}",unit="var"))
        stats=[d for d in work if d]
        
        if not stats:
            print("   (nada pasa filtros)"); return
        df_full = pd.DataFrame(stats)


    # ── RE-ANOTACIÓN POSTERIOR (solo en process):
    #     Si en un MISMO codón (triplete) hay ≥2 POSICIONES seleccionadas
    #     (alguna recta con tag ∈ {pos,neg}), marcamos ese codón como ambiguo
    #     para TODAS las rectas (pos,var) que caen en ese triplete.
    
    # Calcular fmean_pos si no está disponible (modo multi-core)
    if n_jobs > 1:
        _, fmean_pos = build_freq_summaries(long)
    
    selected_positions = set(df_full.loc[df_full["tag"].isin(["pos","neg"]), "pos"])

    def _codon_triple(p: int):
        cpos = codonpos_of.get(int(p))
        if cpos is None:
            return None
        start = int(p) - (cpos - 1)
        return (start, start+1, start+2)

    def _is_multisel_codon(p: int) -> bool:
        triple = _codon_triple(p)
        if triple is None:
            return False
        return sum((q in selected_positions) for q in triple) >= 2

    def _detail_reason(p: int) -> str:
        triple = _codon_triple(p)
        if triple is None:
            return "amb_due_to_multi_sel: no_codon_pos"
        sel = [q for q in triple if q in selected_positions]
        return (f"amb_due_to_multi_sel: codon={triple[0]},{triple[1]},{triple[2]} "
                f" selected={','.join(map(str, sel)) if sel else '[]'}")

    mask_multi = df_full["pos"].apply(_is_multisel_codon)
    if mask_multi.any():
        print(f"🔄 Re-calculando {mask_multi.sum()} variantes en codones con múltiples selecciones...")
        
        # Re-calcular codón y aa_change con selected_positions conocidas
        for idx in df_full[mask_multi].index:
            row = df_full.loc[idx]
            pos = int(row["pos"])
            var = row["var"]
            beta = row["beta"]
            aa_idx = row["aa_idx"]
            
            # Re-calcular con selected_positions
            codon_new, aa_label_new, ambiguous_new, detail_new = infer_codon_aa_for_recta(
                pos, var, fmean_pos, codonpos_of, aaidx_of,
                thr_clear=thr_clear, eps_tie=eps_tie, selected_positions=selected_positions
            )
            
            # Extraer aminoácido del aa_label
            if aa_label_new == "amb (indel)":
                aa = "amb"
                is_indel = True
            elif aa_label_new == "amb":
                aa = "amb"
                is_indel = len(var) != 1
            else:
                aa = aa_label_new.rstrip('0123456789') if aa_label_new else "amb"
                is_indel = len(var) != 1
            
            # Actualizar campos
            df_full.at[idx, "codon"] = codon_new
            df_full.at[idx, "ambiguous"] = True  # Forzar ambiguous=True por multi-selección
            df_full.at[idx, "aa_change"] = _format_aa_change(aa, aa_idx, beta, is_indel)
            
            # Añadir motivo al detalle (sin perder lo previo)
            prev_detail = str(row.get("annot_detail", ""))
            reason = _detail_reason(pos)
            new_detail = reason if (not prev_detail or prev_detail == "nan") else f"{prev_detail} | {reason}"
            df_full.at[idx, "annot_detail"] = new_detail

    # ───────── Construir ranked primero (con las listas day/k/n aún vivas) ─────────
    print("📋 Construyendo ranking de variantes...")
    ranked_pre = df_full[df_full["tag"]!=""].sort_values("score", ascending=False).copy()
    
    # ───────── Añadir columna con la parte después del guión de la variante ─────────
    ranked_pre["variant_suffix"] = ranked_pre["var"].str.split("-", n=1).str[1].fillna("")

    # ───────── Calcular fecha real para TODAS las variantes con n_days > 0 ─────────
    first_day = None
    if times_map and pat in times_map:
        print("📅 Calculando fechas de primera observación...")
        t0 = times_map.get(pat)  # datetime.date
        if t0 is not None:
            # Filtrar TODAS las filas con n_days > 0 (no solo las positivas)
            valid_idx = ranked_pre.index[ranked_pre["n_days"] > 0]
            print(f"📊 Total de variantes con n_days > 0: {len(valid_idx)}")
            
            if len(valid_idx) > 0:
                # Construir una tabla "long" con (pos,var,day,k) de todas las variantes válidas
                recs = []
                debug_samples = []
                for idx in tqdm(valid_idx, desc="Calculando fechas", unit="var"):
                    try:
                        days = ranked_pre.at[idx, "day"]
                        ks   = ranked_pre.at[idx, "k"]
                        p    = ranked_pre.at[idx, "pos"]
                        v    = ranked_pre.at[idx, "var"]
                        tag  = ranked_pre.at[idx, "tag"]
                    except Exception as e:
                        print(f"❌ Error extrayendo datos para idx {idx}: {e}")
                        continue
                    
                    if not isinstance(days, (list, tuple)) or not isinstance(ks, (list, tuple)):
                        print(f"⚠️  Datos no válidos para {p}-{v}: days={type(days)}, ks={type(ks)}")
                        continue
                    
                    # Encontrar el primer día donde k >= min_first_k
                    first_valid_day = None
                    
                    # Crear lista de (día_relativo, k) ordenada por día
                    # NOTA: Los días en las celdas YA son días relativos desde t0 (0, 5, 15, 20...)
                    # NO son índices de muestra
                    day_k_pairs = []
                    for d, k in zip(days, ks):
                        try:
                            if d is None: 
                                continue
                            kf = float(k) if k is not None else 0.0
                            # Los días ya vienen como días relativos desde t0
                            day_k_pairs.append((int(d), kf))
                        except Exception as e:
                            print(f"❌ Error procesando día {d}, k {k}: {e}")
                            continue
                    
                    # Ordenar por día y encontrar el primer día válido
                    day_k_pairs.sort(key=lambda x: x[0])  # Ordenar por día
                    for real_day, kf in day_k_pairs:
                        if kf >= float(min_first_k):
                            first_valid_day = real_day
                            
                            # Debug: imprimir primeras muestras
                            if len(debug_samples) < 5:
                                debug_samples.append({
                                    'pos': p, 'var': v, 'tag': tag,
                                    'day': real_day,
                                    'k': kf, 'min_first_k': min_first_k,
                                    'days_list': days, 'ks_list': ks,
                                    'day_k_pairs': day_k_pairs
                                })
                            break  # Tomar el primer día que cumple la condición
                    
                    if first_valid_day is not None:
                        recs.append((p, v, first_valid_day))
                    else:
                        print(f"⚠️  No se encontró día válido para {p}-{v} (tag: {tag})")
                
                print(f"🔍 Muestras de debug (primeras 5):")
                for sample in debug_samples:
                    print(f"  {sample}")
                
                if recs:
                    pos_long = pd.DataFrame(recs, columns=["pos","var","day"])
                    first_day = (pos_long
                                 .groupby(["pos","var"], as_index=False)["day"]
                                 .min()
                                 .rename(columns={"day":"first_seen_day"}))
                    # Convertir a fecha de calendario
                    first_day["date"] = first_day["first_seen_day"] \
                        .apply(lambda d: (t0 + dt.timedelta(days=int(d))).isoformat())
                    
                    print(f"📅 Fechas calculadas: {len(first_day)} variantes")
                    
                    # Unir solo la columna date a ranked_pre
                    ranked_pre = ranked_pre.merge(first_day[["pos","var","date"]],
                                                  on=["pos","var"], how="left")
                else:
                    print("❌ No se encontraron variantes válidas para calcular fechas")

    # ───────── Guardar variant_stats con day/k/n como JSON (como antes) ─────────
    df_full_out = df_full.copy()
    df_full_out["day"] = df_full_out["day"].apply(json.dumps)
    df_full_out["k"]   = df_full_out["k"].apply(json.dumps)
    df_full_out["n"]   = df_full_out["n"].apply(json.dumps)
    df_full_out.to_csv(outdir/f"{pat}_variant_stats.tsv", sep="\t", index=False)

    # ───────── Guardar selected_ranked con 'date' y sin listas ─────────
    ranked_out = ranked_pre.drop(columns=["day","k","n"], errors="ignore")

    # ───────── Covariación (solo variantes seleccionadas) DESPUÉS del threading ─────────
    if covar_enable:
        print(f"🔍 Calculando covariación por similitud de β (Δβ ≤ {covar_threshold}, Δday ≤ {covar_max_day_diff})...")
        cov_df = _compute_covariation_annotations_selected(
            ranked_pre,
            method=covar_method,
            threshold=covar_threshold,
            min_days=covar_min_days,
            min_range=covar_min_range,
            max_day_diff=covar_max_day_diff
        )
        print(f"📊 Resultados de covariación: {len(cov_df)} variantes analizadas")
        if not cov_df.empty:
            covarying_count = cov_df["covarying"].sum()
            print(f"✅ Encontradas {covarying_count} variantes covariantes")
            ranked_out = ranked_out.merge(cov_df, on=["pos","var"], how="left")
        else:
            print("⚠️  No se encontraron variantes para análisis de covariación")
            ranked_out["covarying"] = False
            ranked_out["covarying_count"] = 0
            ranked_out["covarying_with"] = ""
    else:
        ranked_out["covarying"] = False
        ranked_out["covarying_count"] = 0
        ranked_out["covarying_with"] = ""
    
    # ───────── Añadir columnas de fechas de primera aparición ─────────
    if first_day is not None:
        # Crear mapeo de (pos, var) a first_seen_day y date
        first_day_dict = first_day.set_index(['pos', 'var'])['first_seen_day'].to_dict()
        if 'date' in first_day.columns:
            date_dict = first_day.set_index(['pos', 'var'])['date'].to_dict()
        else:
            date_dict = {}
        
        ranked_out['days_from_t0'] = ranked_out.apply(
            lambda row: first_day_dict.get((row['pos'], row['var']), None), axis=1
        )
        ranked_out['date'] = ranked_out.apply(
            lambda row: date_dict.get((row['pos'], row['var']), None), axis=1
        )
    elif times_map and pat in times_map:
        # Si tenemos times_map pero no first_day, algo salió mal
        ranked_out['days_from_t0'] = None
        ranked_out['date'] = None
    else:
        ranked_out['days_from_t0'] = "No time file provided"
        ranked_out['date'] = "No time file provided"
    
    ranked_out.to_csv(outdir/f"{pat}_selected_ranked.tsv", sep="\t", index=False)

    df = df_full.drop(columns=["day","k","n"])

    beta_mean_all = df["beta"].mean()
    sel_mask = df["tag"].isin(["pos","neg"])
    beta_mean_sel = (df.loc[sel_mask, "beta"].mean()
                     if sel_mask.any() else None)
    print(f"   → {len(ranked_out)} sitios etiquetados")

    # ─ resumen por posición + gráficos ───────────────────────────────
    print("📊 Generando resúmenes por posición...")
    plot_dir = outdir/f"{pat}_plots"; plot_dir.mkdir(exist_ok=True)

    pos_summary = summarise_positions(df)
    pos_summary.to_csv(outdir/f"{pat}_per_position_summary.tsv",
                       sep="\t", index=False)

    # Mapa rápido pos -> fila resumen (para títulos)
    pos2row = {r.pos: r for _, r in pos_summary.iterrows()}

    # 1 figura por posición con alguna recta seleccionada
    print(f"🎨 Generando {len(pos_summary['pos'].unique())} gráficos de posición...")
    for pos in tqdm(pos_summary["pos"].unique(), desc="Gráficos posición", unit="pos"):
        plot_position(pat, pos, df, long,
                      beta_mean_sel, beta_mean_all,
                      pos2row.get(pos), plot_dir, gene_of, aaidx_of)

    # Figura genómica global
    print("🌍 Generando gráfico genómico global...")
    plot_genome(pat, df, long, outdir, beta_mean_sel, beta_mean_all, tsv_path)

    print("\n" + "=" * 50)
    print(f"✅ Análisis completado para {pat}")
    print(f"   📊 {len(ranked_out)} sitios etiquetados")
    print(f"   📁 Archivos guardados en: {outdir}")
    print("=" * 50)

    print(build_final_report_text(
        patient=pat,
        sigma=350,
        nperm=nperm,
        freq_delta=freq_delta,
        thr_clear=thr_clear,
        eps_tie=eps_tie,
        depth_any=100  # <-- ajusta si cambias el umbral en make_long()
    ))

# ───────────── main ────────────────────────────────────────────
if __name__ == "__main__":
    # ----------- modo normal (pipeline completo) -------------
    parser = argparse.ArgumentParser(
        description="Ranking de variantes (mantiene referencia)")
    parser.add_argument("tsv", nargs='+', type=pathlib.Path,
                        help="ficheros *_genes.tsv")
    parser.add_argument("-o","--outdir", default="ranking_out", type=pathlib.Path)
    parser.add_argument("-n","--nperm",  default=2000, type=int,
                        help="permutaciones por variante (def 2000)")
    parser.add_argument("-j","--jobs",   default=-1, type=int,
                        help="núcleos (-1 = todos)")
    parser.add_argument("--freq-delta",  default=0.03, type=float,
                        help="Δf mínimo (p.ej. 0.03 = 3 %)")
    parser.add_argument("--resume", action="store_true",
                        help="omitir pacientes ya procesados")
    parser.add_argument("--thr-clear", default=0.95, type=float,
                        help="Umbral de frecuencia para base 'clara'")
    parser.add_argument("--eps-tie", default=0.05, type=float,
                        help="Mínima diferencia entre la 1ª y 2ª base para no considerarlo empate")
    parser.add_argument("--times", type=pathlib.Path,
                        help="TSV con columnas: sample<TAB>date (formato d/m/yy o d/m/yyyy)")
    parser.add_argument("--first-k-min", default=10, type=int,
                        help="Lecturas mínimas para considerar la primera observación (def 10)")
    # ───────── Opciones de co-variación (activadas por defecto) ─────────
    parser.add_argument("--no-covar", dest="covar", action="store_false",
                        help="Desactiva la anotación de covariación intra-paciente (solo variantes seleccionadas)")
    parser.add_argument("--covar", dest="covar", action="store_true",
                        help="Activa la anotación de covariación intra-paciente (defecto)")
    parser.set_defaults(covar=True)
    parser.add_argument("--covar-method", choices=["beta_diff"], default="beta_diff",
                        help="Método de covariación: beta_diff (diferencia de β)")
    parser.add_argument("--covar-threshold", type=float, default=0.01,
                        help="Diferencia máxima de β para declarar covariación (defecto: 0.01)")
    parser.add_argument("--covar-min-days", type=int, default=3,
                        help="Mínimo nº de días observados por variante para entrar en el análisis")
    parser.add_argument("--covar-min-range", type=float, default=0.05,
                        help="Rango dinámico mínimo de frecuencia (max-min) por variante (defecto: 0.05)")
    parser.add_argument("--covar-max-day-diff", type=int, default=1,
                        help="Diferencia máxima de días de inicio para declarar covariación (defecto: 10)")
    args = parser.parse_args()

    args.outdir.mkdir(exist_ok=True)

    # Cargar mapa paciente -> fecha_real_min (t0) y mapeo de muestras
    times_map = None
    if args.times is not None and args.times.exists():
        tdf = pd.read_csv(args.times, sep="\t", header=None, names=["sample","date"], dtype=str)
        # parseo robusto de fechas (día/mes/año)
        tdf["date"] = pd.to_datetime(tdf["date"], dayfirst=True, errors="coerce")
        # paciente = todo antes del último '-'
        tdf["patient"] = tdf["sample"].str.rsplit("-", n=1).str[0]
        base = (tdf.dropna(subset=["date"])
                   .groupby("patient")["date"]
                   .min())
        # convertir a datetime.date
        times_map = {p: d.to_pydatetime().date() for p, d in base.items()}

    for tsv_file in args.tsv:
        process(tsv_file, args.outdir, args.nperm,
                args.jobs, args.resume, args.freq_delta,
                args.thr_clear, args.eps_tie,
                times_map=times_map, min_first_k=args.first_k_min,
                covar_enable=args.covar,
                covar_method=args.covar_method,
                covar_threshold=args.covar_threshold,
                covar_min_days=args.covar_min_days,
                covar_min_range=args.covar_min_range,
                covar_max_day_diff=args.covar_max_day_diff)
