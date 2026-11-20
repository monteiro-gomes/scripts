"""
13b_remap_enrichr_symbols.py
Fix for script 13: Enrichr symbols are upper case (possibly with ',weight'
suffixes); the DE database uses mouse-case symbols. This script
  1. diagnoses the symbol format of each library,
  2. normalises symbols (strip suffix, upper-case) and recomputes overlaps with
     the paper background (per cell type) and the frozen machinery list,
  3. lists terms for the regulators/marks of interest (primary, controls),
  4. looks for MECP2 / HDAC / KDM6B in other ChEA / consensus libraries.
INSPECTION ONLY: nothing is tested and nothing is frozen here.

Run (from repo root):
  python scripts\13b_remap_enrichr_symbols.py 2>&1 | Out-File -Encoding utf8 results\13b_output.txt

Outputs:
  results/13b_enrichr_terms_mapped.csv   (commit this)
"""
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"
CACHE = DATA / "enrichr_mouse"
CACHE.mkdir(parents=True, exist_ok=True)
RES.mkdir(exist_ok=True)

MAIN_LIBS = ["ChEA_2022", "ENCODE_Histone_Modifications_2015", "Epigenomics_Roadmap_HM_ChIP-seq"]
EXTRA_LIBS = ["ChEA_2016", "ChEA_2015", "ChEA_2013",
              "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]
REGULATORS = ["EZH2", "SUZ12", "JARID2", "MECP2", "REST", "TET1", "HDAC", "BRD4", "KDM6B"]
MISSING_REGS = ["MECP2", "HDAC", "KDM6B"]
MARKS = ["H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac", "H3K36me3"]
ROADMAP_KEEP = re.compile(r"^H3K\w+ (Brain|Neurosphere|Fetal Brain|H1$|H9$)", re.I)
MAX_EXTRA_PRINT = 25


def hr(t):
    print("\n" + "=" * 100)
    print(t)
    print("=" * 100)


def norm(g):
    return re.split(r"[,;\t]", str(g).strip())[0].strip().upper()


# ------------------------------------------------------------------ inputs
def load_backgrounds():
    p = DATA / "de_results_background.sqlite"
    con = sqlite3.connect(p)
    try:
        rows = con.execute("SELECT DISTINCT celltype, gene_name FROM de_results").fetchall()
    finally:
        con.close()
    bg = defaultdict(set)
    for ct, g in rows:
        bg[str(ct)].add(str(g))
    return dict(bg)


def load_machinery():
    p = RES / "epigenetic_machinery_geneset.csv"
    df = pd.read_csv(p)
    print("machinery CSV columns:", list(df.columns), "| rows:", len(df))
    col = next((c for c in df.columns if c.lower() in
                ("gene", "gene_name", "symbol", "gene_symbol", "genes")), df.columns[0])
    genes = set(df[col].astype(str))
    print(f"machinery gene column used: '{col}' -> {len(genes)} genes (expected 177)")
    return genes


def get_names():
    import gseapy as gp
    try:
        return set(gp.get_library_name(organism="Mouse"))
    except Exception as e:
        print(f"[!] get_library_name failed: {e}")
        return set()


def load_library(name):
    path = CACHE / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    import gseapy as gp
    try:
        lib = gp.get_library(name=name, organism="Mouse")
    except Exception as e:
        print(f"[!] get_library failed for {name}: {e}")
        return None
    path.write_text(json.dumps(lib), encoding="utf-8")
    return lib


# ------------------------------------------------------------------ matching
def token_regex(name):
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(name) + r"(?![0-9])", re.I)


def hdac_regex():
    return re.compile(r"(?<![A-Za-z0-9])HDAC\d*", re.I)


def term_matches(term, pat):
    rx = hdac_regex() if pat == "HDAC" else token_regex(pat)
    return rx.search(term) is not None


def species_of(lib_name, term):
    last = re.split(r"[\s_]", term.strip())[-1]
    if lib_name.startswith("ENCODE_Histone"):
        return last  # mm9 / hg19
    if lib_name.startswith("Epigenomics_Roadmap"):
        return "human(Roadmap)"
    return last  # ChEA: Mouse / Human / Rat / tissue word


# ------------------------------------------------------------------ diagnostics
def diagnose(lib_name, lib, union_bg, bg_up, mach_up):
    hr(f"DIAGNOSTIC {lib_name}")
    raw = [g for genes in lib.values() for g in genes]
    n = len(raw)
    has_comma = sum(("," in str(g)) for g in raw) / n
    all_upper = sum((str(g) == str(g).upper()) for g in raw) / n
    exact = sum((str(g) in union_bg) for g in raw) / n
    mapped = sum((norm(g) in bg_up["__union__"]) for g in raw) / n
    print("sample raw symbols (first term):", list(lib[next(iter(lib))])[:10])
    print(f"gene tokens: {n} | with ',' suffix: {has_comma:.1%} | all-upper: {all_upper:.1%}")
    print(f"matching background exactly: {exact:.1%} | after normalisation: {mapped:.1%}")
    universe = {norm(g) for g in raw}
    print(f"unique normalised genes in library: {len(universe)} | "
          f"in background (union): {len(universe & bg_up['__union__'])}")
    print(f"machinery genes present in library universe: {len(mach_up & universe)}/{len(mach_up)}")


# ------------------------------------------------------------------ rows
def make_row(lib_name, pat, term, genes, bg_up, mach_up):
    gn = {norm(g) for g in genes}
    row = {"library": lib_name, "pattern": pat, "term": term,
           "species_tag": species_of(lib_name, term), "n_genes": len(gn)}
    for ct, up in bg_up.items():
        if ct == "__union__":
            continue
        row[f"n_in_bg_{ct}"] = len(gn & up)
    m = sorted(gn & mach_up)
    row["n_machinery"] = len(m)
    row["machinery_genes"] = ";".join(m)
    return row


def fmt_row(r, cts):
    n = r["n_genes"]
    bgtxt = " ".join(f"{ct[:3]}={r[f'n_in_bg_{ct}']}({r[f'n_in_bg_{ct}'] / max(n, 1):.0%})"
                     for ct in cts)
    mg = r["machinery_genes"]
    mtxt = f"mach={r['n_machinery']}"
    if 0 < r["n_machinery"] <= 12:
        mtxt += f" [{mg}]"
    return f"   {r['term']} | n={n} | {bgtxt} | {mtxt}"


def main():
    bg = load_backgrounds()
    cts = list(bg.keys())
    bg_up = {ct: {g.upper() for g in genes} for ct, genes in bg.items()}
    bg_up["__union__"] = set().union(*bg_up.values())
    union_bg = set().union(*bg.values())
    # ambiguity check for upper-casing
    up2mouse = defaultdict(set)
    for g in union_bg:
        up2mouse[g.upper()].add(g)
    amb = {k: v for k, v in up2mouse.items() if len(v) > 1}
    print(f"background genes (union): {len(union_bg)} | upper-case collisions: {len(amb)}")
    if amb:
        print("  examples:", list(amb.items())[:5])
    mach_up = {g.upper() for g in load_machinery()}
    print("machinery in background (normalised):",
          {ct: len(mach_up & bg_up[ct]) for ct in cts}, "(notes say 171/177 and 169/177)")

    names = get_names()
    rows = []

    # ---------------- main libraries
    for lib_name in MAIN_LIBS:
        lib = load_library(lib_name)
        if not lib:
            continue
        diagnose(lib_name, lib, union_bg, bg_up, mach_up)

        if lib_name == "ChEA_2022":
            hr("ChEA_2022: regulator terms (normalised overlaps)")
            for pat in REGULATORS:
                hits = [t for t in lib if term_matches(t, pat)]
                print(f"\n--- {pat}: {len(hits)} terms")
                for t in hits:
                    r = make_row(lib_name, pat, t, lib[t], bg_up, mach_up)
                    rows.append(r)
                    print(fmt_row(r, cts) + f" | {r['species_tag']}")

        elif lib_name.startswith("ENCODE_Histone"):
            hr("ENCODE_Histone_Modifications_2015: mm9 terms only (hg19 excluded)")
            for pat in MARKS:
                hits = [t for t in lib if term_matches(t, pat) and t.strip().endswith("mm9")]
                print(f"\n--- {pat}: {len(hits)} mm9 terms")
                for t in sorted(hits):
                    r = make_row(lib_name, pat, t, lib[t], bg_up, mach_up)
                    rows.append(r)
                    print(fmt_row(r, cts))

        else:  # Roadmap
            hr("Epigenomics_Roadmap_HM_ChIP-seq: brain / neurosphere / ES-like terms")
            for pat in MARKS:
                hits = [t for t in lib if term_matches(t, pat)]
                if not hits:
                    print(f"\n--- {pat}: not in this library")
                    continue
                keep = [t for t in hits if ROADMAP_KEEP.search(t)]
                print(f"\n--- {pat}: {len(hits)} terms, {len(keep)} shown")
                for t in sorted(keep):
                    r = make_row(lib_name, pat, t, lib[t], bg_up, mach_up)
                    rows.append(r)
                    print(fmt_row(r, cts))

    # ---------------- other libraries for regulators missing in ChEA_2022
    hr("Other libraries: MECP2 / HDAC / KDM6B")
    for lib_name in EXTRA_LIBS:
        if names and lib_name not in names:
            print(f"\n[{lib_name}] not in mouse library list, skipped")
            continue
        lib = load_library(lib_name)
        if not lib:
            continue
        print(f"\n[{lib_name}] terms: {len(lib)}; sample term: {next(iter(lib))}")
        for pat in MISSING_REGS:
            hits = [t for t in lib if term_matches(t, pat)]
            print(f"  {pat}: {len(hits)} terms")
            for t in hits[:MAX_EXTRA_PRINT]:
                r = make_row(lib_name, pat, t, lib[t], bg_up, mach_up)
                rows.append(r)
                print(fmt_row(r, cts) + f" | {r['species_tag']}")
            if len(hits) > MAX_EXTRA_PRINT:
                print(f"   ... {len(hits) - MAX_EXTRA_PRINT} more (see CSV)")
                for t in hits[MAX_EXTRA_PRINT:]:
                    rows.append(make_row(lib_name, pat, t, lib[t], bg_up, mach_up))

    if rows:
        df = pd.DataFrame(rows).drop_duplicates(["library", "pattern", "term"])
        out = RES / "13b_enrichr_terms_mapped.csv"
        df.to_csv(out, index=False, encoding="utf-8")
        hr("SUMMARY")
        print(df.groupby(["library", "pattern"]).size().to_string())
        print(f"\nwritten: {out}")
    print("\nDONE. Nothing frozen: freezing happens in script 14 after reading this output.")


if __name__ == "__main__":
    main()