#!/usr/bin/env python3
"""
build_table.py  –  Create a position×time matrix for one patient folder.

Usage:
    ./build_table.py  /path/to/patient_dir
"""
import argparse
import csv
import re
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# ───────────────────────── helpers ──────────────────────────────
DATE_FMT = "%d/%m/%y"
TOKEN_RE = re.compile(r"-(\d+)([a-z]?)$")        # →  ('1','a')  ('3','')
PLUS_RE  = re.compile(r"^\+([A-Z]+):(\d+)")
MINUS_RE = re.compile(r"^-([A-Z]+):(\d+)")
BASES    = "ACTGN"

def parse_time_file(tfile: Path):
    """
    Return:
        id2date   dict{sample_id→date}
        id2token  dict{sample_id→token}
        tokens_in_order  list of tokens as they appear (unique, keep first)
    """
    id2date, id2token = {}, {}
    tokens_order, seen = [], set()
    with tfile.open() as fh:
        for line in fh:                      # ─── clean NULL bytes ───
            if line.strip():
                sid_raw, d = line.rstrip().split("\t")
                sid   = sid_raw.replace('\x00', '')
                token = sid.rsplit("-", 1)[1]
                id2date[sid]  = datetime.strptime(d, DATE_FMT).date()
                id2token[sid] = token
                if token not in seen:
                    tokens_order.append(token)
                    seen.add(token)
    return id2date, id2token, tokens_order

def earliest_day0(id2date, id2token):
    t1_dates = []
    for sid, dt in id2date.items():
        m = TOKEN_RE.search(id2token[sid])
        if m and m.group(1) == "1":
            t1_dates.append(dt)
    # si no hay “1”, toma la fecha más antigua de todas
    return min(t1_dates) if t1_dates else min(id2date.values())

def parse_readcount_line(fields):
    """Return depth, dict{base:count}, list(+), list(-) for ONE position."""
    depth = int(fields[3])
    ref_base = fields[2]
    counts = {b: 0 for b in BASES}
    plus, minus = [], []

    for token in fields[4:]:
        if token.startswith("="):                # reference matches
            n = int(token.split(":")[1])
            counts[ref_base] += n
        elif token[0] in BASES:                  # A:…, C:…, …
            base = token[0]
            n = int(token.split(":")[1])
            counts[base] += n
        elif token.startswith("+"):
            m = PLUS_RE.match(token)
            if m: plus.append(f"+{m.group(1)}{m.group(2)}")
        elif token.startswith("-"):
            m = MINUS_RE.match(token)
            if m: minus.append(f"-{m.group(1)}{m.group(2)}")
    return depth, counts, plus, minus

# ───────────────────────── main ─────────────────────────────────
def main(pdir: Path):
    # 1· sanity & I/O paths
    if not pdir.is_dir():
        raise SystemExit(f"❌  {pdir} is not a directory")
    time_file = next(pdir.glob("*.time"), None)
    if not time_file:
        raise SystemExit("❌  No *.time file found in the directory")

    # 2· load dates & tokens
    id2date, id2token, tokens = parse_time_file(time_file)
    day0 = earliest_day0(id2date, id2token)

    # 3· find readcount files we actually have
    #     stem → "009A-A-1.readcounts",  necesitamos quitar también ".readcounts"
    rc_files = {}
    for f in pdir.glob("*.readcounts.tsv"):
        sid = re.sub(r"\.readcounts$", "", f.stem).replace('\x00', '')
        rc_files[sid] = f

    if not rc_files:
        raise SystemExit("❌  No *.readcounts.tsv files found")

    # 4· matrix: pos → token → code
    matrix = defaultdict(dict)

    for sid, fpath in rc_files.items():
        token = id2token.get(sid)
        if not token:                # sample without date/token – skip
            continue
        days = (id2date[sid] - day0).days

        with fpath.open(errors="ignore") as fh:
            tsv = csv.reader(fh, delimiter="\t")
            for row in tsv:
                if not row:
                    continue
                # limpia NULs en cada campo
                row = [c.replace('\x00', '') for c in row]
                try:
                    pos = int(row[1])
                except ValueError:
                    continue                 # línea malformada
                depth, counts, plus, minus = parse_readcount_line(row)

                code = (
                    f"{days}/{depth}/"
                    f"A{counts['A']}C{counts['C']}T{counts['T']}"
                    f"G{counts['G']}N{counts['N']}/"
                    f"{''.join(plus) if plus else '+0'}"
                    f"{''.join(minus) if minus else '-0'}"
                )
                matrix[pos][token] = code

    # 6· write output
    out = pdir / f"{pdir.name}.tsv"
    with out.open("w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["position"] + tokens)
        for pos in sorted(matrix):
            row = [matrix[pos].get(tok, "NA") for tok in tokens]
            wr.writerow([pos] + row)

    print(f"✔  Wrote {out}  ({len(matrix)} positions × {len(tokens)} time-points)")

# ───────────────────────── runner ───────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build position×time matrix for one patient")
    ap.add_argument("patient_dir", type=Path, help="Path to patient folder")
    args = ap.parse_args()
    main(args.patient_dir)
