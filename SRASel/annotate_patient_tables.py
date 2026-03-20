#!/usr/bin/env python3
"""
annotate_patient_tables.py
──────────────────────────
Add gene, codon-position (1|2|3) and amino-acid index to each genomic
position in every patient TSV produced by build_time_mutation_tsv.py.

Usage:
    ./annotate_patient_tables.py  --tsv-dir  tsv_conjuntos_tiempos \
                                  --gb       MN908947.3.gb
The script writes   <patient>_genes.tsv   files in the same directory.
"""
import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

def print_caveat_notice(ref_name: str) -> None:
    msg = f"""
[IMPORTANT CAVEAT] Annotation assumptions (reference: {ref_name})
This script is tuned for SARS-CoV-2 (NC_045512.2): coding sequences (CDS)
are unspliced (no introns), lie on the positive strand, and some CDS overlap. ORF1ab is
produced via a programmed −1 ribosomal frameshift and is often annotated in GenBank with
a CompoundLocation (join) plus a /ribosomal_slippage note. Under these conditions, assigning
codon_pos (1|2|3) and aa_idx by iterating genomic coordinates in reading order is appropriate.

If adapting to other genomes, please review:
• Strand: CDS on the negative (−) strand must be traversed in reading order (reverse
  genomic coordinates) before assigning codon_pos/aa_idx.
• Splicing: Eukaryotic genes frequently have introns (true splicing). Joined multi-exon CDS
  require explicit reading-order traversal across exons; frame continuity must be verified.
• Overlaps & alternative ORFs: Overlapping CDS or internal ORFs may map multiple (gene, frame)
  assignments to the same genomic position; downstream logic should decide how to collapse them.
• Coordinate convention: GenBank feature locations are 0-based, end-exclusive; this script
  converts them to 1-based, inclusive coordinates.

These assumptions are safe for canonical SARS-CoV-2, but reassess them if you change organism
or reference annotation.
""".strip()
    print(msg, file=sys.stderr)

# ─────────────────────────── helpers ──────────────────────────
def load_coding_table(gb_path: Path):
    """Return dict {genomic_pos (int) → (gene, codon_pos, aa_index)}."""
    record = SeqIO.read(gb_path, "genbank")

    # pos → OrderedDict{gene → list[(codon_pos, aa_idx)]}
    pos2info = defaultdict(lambda: OrderedDict())

    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene = feat.qualifiers.get("gene", ["CDS"])[0]

        # Build the list of genomic coordinates IN READING ORDER
        # ───────────── construir coordenadas en orden de lectura ─────────────
        coords = []
        parts = (feat.location.parts
                 if isinstance(feat.location, CompoundLocation)
                 else [feat.location])

        for i, part in enumerate(parts):
            start = int(part.start) + 1          # BioPython 0-based, end-exclusive
            end   = int(part.end)                # inclusive (1-based)

            coords.extend(range(start, end + 1))


        # Assign codon position & amino-acid index
        for idx, gpos in enumerate(coords):
            codon_pos = idx % 3 + 1              # 1,2,3
            aa_idx    = idx // 3 + 1             # 1-based

            gene_list = pos2info[gpos].setdefault(gene, [])
            # ── evita guardar dos veces el mismo codon_pos para el mismo gen
            if all(cp != codon_pos for cp, _ in gene_list):
                gene_list.append((codon_pos, aa_idx))

        # tras construir coords (en orden de lectura)
        last_cpos = (len(coords) - 1) % 3 + 1
        if last_cpos != 3:
            sys.stderr.write(
                f"⚠  CDS {gene} ends with codon_pos {last_cpos} "
                f"(length={len(coords)}; not a multiple of 3)\n"
            )

    return pos2info


def annotate_tsv(tsv_path: Path, pos2info: dict):
    """Write <patient>_genes.tsv with three extra columns."""
    out_path = tsv_path.with_name(f"{tsv_path.stem}_genes.tsv")

    with tsv_path.open() as fin, out_path.open("w", newline="") as fout:
        rd = csv.reader(fin, delimiter="\t")
        wr = csv.writer(fout, delimiter="\t")

        header = next(rd)
        wr.writerow([header[0], "gene", "codon_pos", "aa_idx", *header[1:]])

        for row in rd:
            try:
                pos = int(row[0])
            except ValueError:
                wr.writerow(row + ["", "", ""])  # keep malformed lines
                continue

            info = pos2info.get(pos)
            if info:
                genes  = []
                cpos   = []
                aindex = []

                for g, lst in info.items():          # Ordered: first seen first
                    genes.append(g)

                    # ---- condensar pares (codon_pos, aa_idx) sin perder orden ----
                    uniq_pairs = []
                    seen_cp = set()
                    for cp, ai in lst:
                        if cp not in seen_cp:
                            seen_cp.add(cp)
                            uniq_pairs.append((cp, ai))

                    # ── Mantener duplicado solo si llega exactamente en orden 3,1
                    #    (frameshift real); cualquier otro caso se colapsa.
                    if (
                        len(uniq_pairs) == 2
                        and uniq_pairs[0][0] == 3
                        and uniq_pairs[1][0] == 1
                    ):
                        keep_pairs = uniq_pairs         # (3,ai)(1,ai)
                    else:
                        keep_pairs = [uniq_pairs[0]]    # solo la primera pareja

                    cpos.extend(str(cp)  for cp, _ in keep_pairs)
                    aindex.extend(str(ai) for _,  ai in keep_pairs)

                wr.writerow([
                    row[0],
                    ",".join(genes),
                    ",".join(cpos),
                    ",".join(aindex),
                    *row[1:],
                ])
            else:
                wr.writerow([row[0], "", "", "", *row[1:]])
    return out_path


# ─────────────────────────── main ─────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Add gene & codon information to patient mutation tables")
    ap.add_argument("--tsv-dir", default="tsv_conjuntos_tiempos",
                    help="Directory containing the <patient>.tsv files")
    ap.add_argument("--gb", default="NC_045512.2/NC_045512.2.gb",
                    help="Path to MN908947.3 GenBank file (default: ./MN908947.3.gb)")
    args = ap.parse_args()

    tsv_dir = Path(args.tsv_dir).resolve()
    gb_path = Path(args.gb).resolve()

    if not tsv_dir.is_dir():
        sys.exit(f"❌  TSV directory not found: {tsv_dir}")
    if not gb_path.is_file():
        sys.exit(f"❌  GenBank file not found: {gb_path}")

    print_caveat_notice(gb_path.name)

    print(f"Reading coding information from {gb_path.name} …", file=sys.stderr)
    pos2info = load_coding_table(gb_path)

    tsv_files = sorted(p for p in tsv_dir.glob("*.tsv")
                       if not p.name.endswith("_genes.tsv"))

    if not tsv_files:
        sys.exit("❌  No *.tsv files found to annotate.")

    for tsv in tsv_files:
        out = annotate_tsv(tsv, pos2info)
        print(f"✔  {out.name} written ({tsv.name} + annotations)")

if __name__ == "__main__":
    main()
