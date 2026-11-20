"""
14_freeze_enrichr_sets.py
Freeze the regulator / histone-mark gene sets (Enrichr, cached mouse libraries)
for analyses 2b and 3, using a selection rule fixed BEFORE any test.
NOTHING IS TESTED HERE.

Rules (pre-specified):
  - Regulators (ChEA_2022, Mouse-labelled terms only):
      primary_neural  : EZH2 NPCs, SUZ12 NPCs, REST NEURONS
      secondary_proxy : per regulator, consensus of the mouse ESC/mESC terms
                        (genes present in >= max(2, ceil(k/2)) of the k terms)
      exploratory     : BRD4 (only mouse term is AML cells)
  - Marks (ENCODE_Histone_Modifications_2015, mm9 only):
      primary_neural  : brain (5 marks); cortical plate (H3K4me3, H3K27ac)
      secondary_proxy : cerebellum; Roadmap human hippocampus / mid-frontal
                        (H3K27me3, H3K36me3)
      negative_control: liver, kidney, heart, spleen, thymus, small intestine
  - A set is frozen only if >= MIN_BG of its genes are in the background of
    BOTH cell types (MIN_BG = 100).
  - Machinery genes stay in the sets but are flagged (is_machinery); the
    default for 2b/3 is to drop them from sets and DEG lists.

Run (from repo root):
  python scripts\14_freeze_enrichr_sets.py 2>&1 | Out-File -Encoding utf8 results\14_output.txt

Outputs (commit all, then tag `enrichr-sets-frozen`):
  results/enrichr_frozen_sets_manifest.csv   one row per candidate set
  results/enrichr_frozen_sets_genes.csv      long table: frozen sets x genes
  results/enrichr_frozen_sets_jaccard.csv    overlap between regulator/neural sets
  results/14_enrichr_provenance.json         library hashes, versions, date
"""
import gzip
import hashlib
import json
import math
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import date
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
RES.mkdir(exist_ok=True)

MIN_BG = 100
CHEA, ENC, ROAD = "ChEA_2022", "ENCODE_Histone_Modifications_2015", "Epigenomics_Roadmap_HM_ChIP-seq"
MARKS = ["H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac", "H3K36me3"]

# ---- explicit term lists (frozen by name, copied from the 13b output)
NEURAL_REG = [
    ("EZH2", "EZH2_NPC", "EZH2 27294783 Chip-Seq NPCs Mouse"),
    ("SUZ12", "SUZ12_NPC", "SUZ12 27294783 Chip-Seq NPCs Mouse"),
    ("REST", "REST_NEURONS", "REST 19997604 ChIP-ChIP NEURONS Mouse"),
]
ESC_CONSENSUS = {
    "EZH2": ["EZH2 27304074 Chip-Seq ESCs Mouse", "EZH2 18974828 ChIP-Seq MESCs Mouse",
             "EZH2 27294783 Chip-Seq ESCs Mouse"],
    "SUZ12": ["SUZ12 18555785 ChIP-Seq MESCs Mouse", "SUZ12 16625203 ChIP-ChIP MESCs Mouse",
              "SUZ12 18692474 ChIP-Seq MESCs Mouse", "SUZ12 18974828 ChIP-Seq MESCs Mouse",
              "SUZ12 18555785 Chip-Seq ESCs Mouse", "SUZ12 27294783 Chip-Seq ESCs Mouse"],
    "JARID2": ["JARID2 20064375 ChIP-Seq MESCs Mouse", "JARID2 20075857 ChIP-Seq MESCs Mouse"],
    "TET1": ["TET1 21451524 ChIP-Seq MESCs Mouse", "TET1 21490601 ChIP-Seq MESCs Mouse"],
    "REST": ["REST 21632747 ChIP-Seq MESCs Mouse", "REST 18959480 ChIP-ChIP MESCs Mouse"],
}
EXPLORATORY_REG = [("BRD4", "BRD4_AML", "BRD4 27068464 Chip-Seq AML-cells Mouse")]
EXCLUDED_NOTE = ("excluded by rule: Human/Rat/lymphoma terms, EZH2 myoblasts, SUZ12 MEFs, "
                 "hg19 ENCODE, HDAC2 ENCODE consensus; MECP2/HDAC/KDM6B absent from ChEA")

NEG_TISSUES = ["liver", "kidney", "heart", "spleen", "thymus", "small intestine"]
ROADMAP_TISSUES = ["Brain Hippocampus Middle", "Brain Mid Frontal Lobe"]


def hr(t):
    print("\n" + "=" * 100)
    print(t)
    print("=" * 100)


def norm(g):
    return re.split(r"[,;\t]", str(g).strip())[0].strip().upper()


def sha256(p):
    h = hashlib.sha256()
    h.update(Path(p).read_bytes())
    return h.hexdigest()


# ------------------------------------------------------------------ inputs
def load_backgrounds():
    con = sqlite3.connect(DATA / "de_results_background.sqlite")
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
    hr("MACHINERY CSV CHECK (notes say 177 genes, 17 groups)")
    print("path:", p, "| sha256:", sha256(p)[:16])
    print("columns:", list(df.columns))
    print(f"rows: {len(df)} | unique gene_name: {df['gene_name'].nunique()} | "
          f"unique after upper(): {df['gene_name'].str.upper().nunique()}")
    dups = df[df["gene_name"].duplicated(keep=False)]
    if len(dups):
        print("DUPLICATED gene_name rows:\n", dups.to_string())
    if "group" in df.columns:
        print(f"groups: {df['group'].nunique()}")
        print(df["group"].value_counts().to_string())
    if "in_db" in df.columns:
        print("in_db value counts:", df["in_db"].value_counts(dropna=False).to_dict())
        print("genes with in_db False/0:", df.loc[~df["in_db"].astype(bool), "gene_name"].tolist())
    print("-> compare with the frozen tag:  git show geneset-frozen:results/"
          "epigenetic_machinery_geneset.csv | Measure-Object -Line")
    return {g.upper() for g in df["gene_name"].astype(str)}


def load_lib(name):
    p = CACHE / f"{name}.json"
    if not p.exists():
        sys.exit(f"[!] cached library missing: {p} (run scripts 13/13b first)")
    lib = json.loads(p.read_text(encoding="utf-8"))
    return {t.strip(): {norm(g) for g in genes} for t, genes in lib.items()}, p


def gencode_symbols():
    p = DATA / "gencode.vM25.annotation.gtf.gz"
    if not p.exists():
        print(f"[!] {p} not found: GENCODE check skipped")
        return None
    rx = re.compile(r'gene_name "([^"]+)"')
    names = set()
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = rx.search(parts[8])
            if m:
                names.add(m.group(1).upper())
    return names


# ------------------------------------------------------------------ spec builder
def consensus(libn, terms):
    k = len(terms)
    need = max(2, math.ceil(k / 2))
    cnt = Counter(g for t in terms for g in libn[t])
    return {g for g, c in cnt.items() if c >= need}, need


def build_specs(libs):
    specs, missing = [], []

    def add(set_id, role, kind, name, lib, terms, mode, species, tissue, rule):
        absent = [t for t in terms if t not in libs[lib]]
        if absent:
            missing.append((set_id, absent))
            return
        if mode == "single":
            genes, need = set(libs[lib][terms[0]]), 1
        else:
            genes, need = consensus(libs[lib], terms)
        specs.append(dict(set_id=set_id, role=role, kind=kind, name=name, library=lib,
                          source_terms=";".join(terms), n_terms=len(terms), min_terms=need,
                          species=species, tissue=tissue, rule=rule, genes=genes))

    for reg, tag, term in NEURAL_REG:
        add(f"reg_{tag}", "primary_neural", "regulator", reg, CHEA, [term], "single",
            "mouse", tag.split("_")[1], "neural ChEA term")
    for reg, terms in ESC_CONSENSUS.items():
        add(f"reg_{reg}_mESC_cons", "secondary_proxy", "regulator", reg, CHEA, terms, "consensus",
            "mouse", "mESC/ESC", "consensus >= max(2, ceil(k/2)) terms")
    for reg, tag, term in EXPLORATORY_REG:
        add(f"reg_{tag}", "exploratory", "regulator", reg, CHEA, [term], "single",
            "mouse", "AML", "single non-neural term")

    for mark in MARKS:
        add(f"mark_{mark}_brain_mm9", "primary_neural", "mark", mark, ENC, [f"{mark} brain mm9"],
            "single", "mouse", "brain", "ENCODE mm9 brain")
    for mark in ("H3K4me3", "H3K27ac"):
        add(f"mark_{mark}_cortical_plate_mm9", "primary_neural", "mark", mark, ENC,
            [f"{mark} cortical plate mm9"], "single", "mouse", "cortical plate", "ENCODE mm9")
    for mark in MARKS:
        t = f"{mark} cerebellum mm9"
        if t in libs[ENC]:
            add(f"mark_{mark}_cerebellum_mm9", "secondary_proxy", "mark", mark, ENC, [t],
                "single", "mouse", "cerebellum", "ENCODE mm9")
    for mark in ("H3K27me3", "H3K36me3"):
        for tis in ROADMAP_TISSUES:
            add(f"mark_{mark}_{tis.replace(' ', '_')}_roadmap", "secondary_proxy", "mark", mark,
                ROAD, [f"{mark} {tis}"], "single", "human(Roadmap)", tis, "Roadmap brain")
    for mark in MARKS:
        for tis in NEG_TISSUES:
            t = f"{mark} {tis} mm9"
            if t in libs[ENC]:
                add(f"mark_{mark}_{tis.replace(' ', '_')}_mm9_NEG", "negative_control", "mark",
                    mark, ENC, [t], "single", "mouse", tis, "non-brain tissue control")
    return specs, missing


# ------------------------------------------------------------------ main
def main():
    bg = load_backgrounds()
    cts = list(bg)
    bg_up = {ct: {g.upper() for g in v} for ct, v in bg.items()}
    up2mouse = {g.upper(): g for v in bg.values() for g in v}
    mach = load_machinery()

    libs, paths = {}, {}
    for name in (CHEA, ENC, ROAD):
        libs[name], paths[name] = load_lib(name)

    # ---- are unmatched library genes non-expressed, or symbol failures?
    hr("SYMBOL DIAGNOSTIC: unmatched library genes vs GENCODE vM25")
    gsym = gencode_symbols()
    if gsym:
        for name, libn in libs.items():
            uni = set().union(*libn.values())
            in_bg = uni & set(up2mouse)
            rest = uni - in_bg
            in_gc = rest & gsym
            neither = rest - gsym
            print(f"{name}: unique genes {len(uni)} | in background {len(in_bg)} | "
                  f"not in bg but valid mouse symbol (not expressed/tested) {len(in_gc)} | "
                  f"not a GENCODE symbol (human-only/LOC/alias) {len(neither)}")
            print("   examples not in GENCODE:", sorted(neither)[:15])

    specs, missing = build_specs(libs)
    if missing:
        hr("MISSING TERMS (spec lists a term that is not in the library)")
        for sid, absent in missing:
            print(sid, "->", absent)

    # ---- manifest + gene table
    manifest, gene_rows = [], []
    for s in specs:
        genes = s["genes"]
        row = {k: v for k, v in s.items() if k != "genes"}
        row["n_raw"] = len(genes)
        for ct in cts:
            row[f"n_bg_{ct}"] = len(genes & bg_up[ct])
        row["n_machinery"] = len(genes & mach)
        row["frozen"] = all(row[f"n_bg_{ct}"] >= MIN_BG for ct in cts)
        manifest.append(row)
        if row["frozen"]:
            for g in sorted(genes):
                gene_rows.append({"set_id": s["set_id"], "role": s["role"], "gene_upper": g,
                                  "gene_mouse": up2mouse.get(g, ""),
                                  **{f"in_bg_{ct}": g in bg_up[ct] for ct in cts},
                                  "is_machinery": g in mach})
    man = pd.DataFrame(manifest)
    gen = pd.DataFrame(gene_rows)

    hr(f"CANDIDATE SETS (frozen requires >= {MIN_BG} background genes in BOTH cell types)")
    order = {"primary_neural": 0, "secondary_proxy": 1, "exploratory": 2, "negative_control": 3}
    man["_o"] = man["role"].map(order)
    man = man.sort_values(["_o", "kind", "set_id"]).drop(columns="_o")
    for _, r in man.iterrows():
        bgtxt = " ".join(f"{ct[:3]}={r[f'n_bg_{ct}']}" for ct in cts)
        print(f"{'FROZEN ' if r['frozen'] else 'DROPPED'} {r['role']:<16} {r['set_id']:<48} "
              f"raw={r['n_raw']:<5} {bgtxt} mach={r['n_machinery']} "
              f"(terms={r['n_terms']}, need>={r['min_terms']})")
    print("\n" + EXCLUDED_NOTE)
    print("\nframed counts by role:")
    print(man.groupby(["role", "frozen"]).size().to_string())

    # ---- overlap between regulator / primary / secondary sets (independence check)
    hr("JACCARD overlap (genes in background of both cell types) among non-control frozen sets")
    keep = man[(man["frozen"]) & (man["role"].isin(["primary_neural", "secondary_proxy"]))]
    inter_bg = set.intersection(*[bg_up[ct] for ct in cts])
    sets = {r["set_id"]: {g for g in next(s["genes"] for s in specs if s["set_id"] == r["set_id"])
                          if g in inter_bg} for _, r in keep.iterrows()}
    ids = list(sets)
    jac = pd.DataFrame(index=ids, columns=ids, dtype=float)
    for a in ids:
        for b in ids:
            u = len(sets[a] | sets[b])
            jac.loc[a, b] = len(sets[a] & sets[b]) / u if u else float("nan")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    print(jac.round(2).to_string())
    jac.to_csv(RES / "enrichr_frozen_sets_jaccard.csv", encoding="utf-8")

    # ---- write
    man.to_csv(RES / "enrichr_frozen_sets_manifest.csv", index=False, encoding="utf-8")
    gen.to_csv(RES / "enrichr_frozen_sets_genes.csv", index=False, encoding="utf-8")
    try:
        import gseapy
        gv = gseapy.__version__
    except Exception:
        gv = "unknown"
    prov = {"date": date.today().isoformat(), "gseapy": gv, "min_bg": MIN_BG,
            "libraries": {n: {"cache_file": str(p.relative_to(ROOT)), "sha256": sha256(p),
                              "n_terms": len(libs[n])} for n, p in paths.items()},
            "machinery_csv_sha256": sha256(RES / "epigenetic_machinery_geneset.csv"),
            "background_db_sha256": sha256(DATA / "de_results_background.sqlite"),
            "rule": "see docstring of scripts/14_freeze_enrichr_sets.py"}
    (RES / "14_enrichr_provenance.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")

    hr("WRITTEN")
    for f in ("enrichr_frozen_sets_manifest.csv", "enrichr_frozen_sets_genes.csv",
              "enrichr_frozen_sets_jaccard.csv", "14_enrichr_provenance.json"):
        print("  results\\" + f)
    print("\nNext: resolve the 177 vs 178 question, then commit these files and tag enrichr-sets-frozen.")
    print("Nothing has been tested.")


if __name__ == "__main__":
    main()