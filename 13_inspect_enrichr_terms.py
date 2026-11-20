"""
13_inspect_enrichr_terms.py
Inspect term names in ChEA_2022, ENCODE_Histone_Modifications_2015 and
Epigenomics_Roadmap_HM_ChIP-seq (Enrichr, mouse) for the regulators/marks of
interest. INSPECTION ONLY: nothing is tested and nothing is frozen here.

Run (from repo root):
  python scripts\13_inspect_enrichr_terms.py 2>&1 | Out-File -Encoding utf8 results\13_output.txt

Outputs:
  data/enrichr_mouse/<library>.json         cached libraries (not in Git)
  results/13_enrichr_term_candidates.csv    all matching terms (commit this)
"""
import json
import re
import sqlite3
import sys
from collections import Counter
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

LIBS = {
    "ChEA_2022": "regulator",
    "ENCODE_Histone_Modifications_2015": "mark",
    "Epigenomics_Roadmap_HM_ChIP-seq": "mark",
}
REGULATORS = ["EZH2", "SUZ12", "JARID2", "MECP2", "REST", "TET1", "HDAC", "BRD4", "KDM6B"]
MARKS = ["H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac", "H3K36me3"]
MAX_PRINT = 80  # terms printed per pattern per library (CSV always has all)


def hr(title):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


# ---------------------------------------------------------------- libraries
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


def list_available_names():
    import gseapy as gp
    try:
        return gp.get_library_name(organism="Mouse")
    except Exception as e:
        print(f"[!] get_library_name failed: {e}")
        return []


# ---------------------------------------------------------------- gene sets
def first_col_genes(df):
    for c in df.columns:
        if c.lower() in ("gene", "gene_name", "symbol", "gene_symbol", "genes"):
            return set(df[c].astype(str))
    return set(df.iloc[:, 0].astype(str))


def load_machinery():
    p = RES / "epigenetic_machinery_geneset.csv"
    if not p.exists():
        print(f"[!] {p} not found; machinery overlap skipped")
        return set()
    return first_col_genes(pd.read_csv(p))


def load_backgrounds():
    """Background genes per cell type from data/de_results_background.sqlite (made by 07b)."""
    p = DATA / "de_results_background.sqlite"
    out = {}
    if not p.exists():
        print(f"[!] {p} not found; background overlap skipped")
        return out
    con = sqlite3.connect(p)
    try:
        tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")]
        print("background DB tables:", tables)
        tbl = "de_results" if "de_results" in tables else tables[0]
        cols = [r[1] for r in con.execute(f"PRAGMA table_info({tbl})")]
        print(f"columns in {tbl}:", cols)
        gcol = next((c for c in cols if c.lower() in
                     ("gene", "gene_name", "symbol", "gene_symbol", "gene_id")), None)
        ccol = next((c for c in cols if "cell" in c.lower() or "cluster" in c.lower()
                     or "type" in c.lower()), None)
        if gcol is None:
            print("[!] could not identify a gene column; background overlap skipped")
            return out
        if ccol is None:
            out["all"] = {r[0] for r in con.execute(f"SELECT DISTINCT {gcol} FROM {tbl}")}
        else:
            for (ct,) in con.execute(f"SELECT DISTINCT {ccol} FROM {tbl}").fetchall():
                genes = {r[0] for r in con.execute(
                    f"SELECT DISTINCT {gcol} FROM {tbl} WHERE {ccol}=?", (ct,))}
                out[str(ct)] = genes
        for k, v in out.items():
            print(f"  background[{k}]: {len(v)} genes")
    except Exception as e:
        print(f"[!] background load failed: {e}")
    finally:
        con.close()
    return out


# ---------------------------------------------------------------- matching
def token_regex(name):
    # name must not be embedded in a longer alphanumeric token (H3K4me3 vs H3K4me3x)
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(name) + r"(?![0-9])", re.I)


def hdac_regex():
    return re.compile(r"(?<![A-Za-z0-9])HDAC\d*", re.I)


def match_terms(terms, pattern_name):
    rx = hdac_regex() if pattern_name == "HDAC" else token_regex(pattern_name)
    hits = []
    for t in terms:
        m = rx.search(t)
        if m:
            first_tok = re.split(r"[\s_]", t.strip())[0].upper()
            hits.append((t, first_tok == m.group(0).upper() or
                         (pattern_name == "HDAC" and first_tok.startswith("HDAC"))))
    return hits


def main():
    machinery = load_machinery()
    backgrounds = load_backgrounds()

    hr("Available Enrichr mouse library names matching our targets")
    names = list_available_names()
    print(f"{len(names)} mouse libraries listed")
    for n in names:
        if any(k in n for k in ("ChEA", "ENCODE_Histone", "Roadmap", "Epigenomics")):
            print("  ", n)

    rows = []
    for lib_name, kind in LIBS.items():
        hr(f"LIBRARY {lib_name}")
        lib = load_library(lib_name)
        if not lib:
            print("  -> not loaded, skipping")
            continue
        terms = list(lib.keys())
        sizes = pd.Series({t: len(set(g)) for t, g in lib.items()})
        print(f"terms: {len(terms)} | genes/term median {sizes.median():.0f}, "
              f"min {sizes.min()}, max {sizes.max()}")
        print("first 8 terms (format check):")
        for t in terms[:8]:
            print("   ", t, f"[{sizes[t]} genes]")
        last_tok = Counter(re.split(r"[\s_]", t.strip())[-1] for t in terms)
        print("most common last token (species/cell hint):", last_tok.most_common(10))

        patterns = REGULATORS if kind == "regulator" else MARKS
        # regulators also searched in histone libs? no: marks are not TF names
        for pat in patterns:
            hits = match_terms(terms, pat)
            print(f"\n--- {lib_name} :: {pat}: {len(hits)} matching terms")
            for t, is_first in hits[:MAX_PRINT]:
                genes = set(lib[t])
                n_mach = len(genes & machinery) if machinery else -1
                bg_txt = " ".join(f"bg[{k}]={len(genes & v)}" for k, v in backgrounds.items())
                print(f"   {'*' if is_first else ' '} {t} | n={len(genes)} "
                      f"| machinery={n_mach} {bg_txt}")
                row = {"library": lib_name, "pattern": pat, "kind": kind, "term": t,
                       "first_token_match": is_first, "n_genes": len(genes),
                       "n_machinery": n_mach}
                for k, v in backgrounds.items():
                    row[f"n_in_bg_{k}"] = len(genes & v)
                rows.append(row)
            if len(hits) > MAX_PRINT:
                print(f"   ... {len(hits) - MAX_PRINT} more (see CSV)")
                for t, is_first in hits[MAX_PRINT:]:
                    genes = set(lib[t])
                    rows.append({"library": lib_name, "pattern": pat, "kind": kind, "term": t,
                                 "first_token_match": is_first, "n_genes": len(genes),
                                 "n_machinery": len(genes & machinery) if machinery else -1})
            if not hits:
                print("   (none)")

        if kind == "mark":
            hr(f"{lib_name}: cell/tissue annotation of mark terms (to check brain relevance)")
            mark_terms = [t for t in terms if any(token_regex(m).search(t) for m in MARKS)]
            print(f"{len(mark_terms)} terms contain one of {MARKS}")
            for t in mark_terms[:300]:
                print("   ", t, f"[{sizes[t]}]")
            if len(mark_terms) > 300:
                print(f"   ... {len(mark_terms) - 300} more (see CSV)")
            brain_rx = re.compile(r"brain|cortex|hippocamp|neuron|neural|astro|NPC|ES[-_ ]|mESC|embryo",
                                  re.I)
            print("\nterms with brain/neural/ES-like words:")
            for t in mark_terms:
                if brain_rx.search(t):
                    print("   ", t, f"[{sizes[t]}]")

    if rows:
        df = pd.DataFrame(rows).drop_duplicates(["library", "pattern", "term"])
        out = RES / "13_enrichr_term_candidates.csv"
        df.to_csv(out, index=False, encoding="utf-8")
        hr("SUMMARY")
        print(df.groupby(["library", "pattern"]).size().to_string())
        print(f"\nwritten: {out}")

    print("\nDONE. Nothing frozen yet: next step is to choose terms (script 14) and tag `enrichr-sets-frozen`.")


if __name__ == "__main__":
    main()