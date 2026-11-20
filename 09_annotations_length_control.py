# %%
# 09: mm10 annotations (GENCODE vM25, UCSC CpG islands); length- and sex-linked-gene controls
# Inputs : data/gencode.vM25.annotation.gtf.gz, data/cpgIslandExt.txt.gz,
#          data/de_results_background.sqlite, results/epigenetic_machinery_geneset.csv
# Outputs: results/gene_annotations.csv, results/annotation_summary.csv, results/length_matched_null.csv
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
DATA, RES = ROOT / "data", ROOT / "results"
GTF = DATA / "gencode.vM25.annotation.gtf.gz"
CGI = DATA / "cpgIslandExt.txt.gz"
DB_BG = DATA / "de_results_background.sqlite"
for f in (GTF, CGI, DB_BG):
    assert f.exists(), f"Missing file: {f}"

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
HOURS = [1, 2, 4, 24, 72]
N_EXPR_BINS, N_LEN_BINS, NDRAW = 10, 5, 5000
rng = np.random.default_rng(42)

# %% 1. Gene table from GENCODE vM25
gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                  names=["chrom", "feature", "start", "end", "strand", "attr"])
genes = gtf[gtf.feature == "gene"].copy()
genes["gene_name"] = genes.attr.str.extract(r'gene_name "([^"]+)"')[0]
genes["gene_type"] = genes.attr.str.extract(r'gene_type "([^"]+)"')[0]
genes = genes.dropna(subset=["gene_name"])
genes["gene_span"] = genes.end - genes.start + 1
genes = genes.sort_values("gene_span", ascending=False).drop_duplicates("gene_name")  # keep longest if duplicated symbol
genes["tss"] = np.where(genes.strand == "+", genes.start, genes.end)
genes["is_sex_linked"] = genes.chrom.isin(["chrX", "chrY"])
print("genes in annotation:", len(genes))

# %% 2. Promoter (TSS +/- 1 kb) overlap with CpG islands
cgi = pd.read_csv(CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])

def promoter_in_cgi(g_all, cgi_all, flank=1000):
    out = pd.Series(False, index=g_all.index)
    for chrom, g in g_all.groupby("chrom"):
        c = cgi_all[cgi_all.chrom == chrom].sort_values("start")
        if c.empty:
            continue
        cs = c.start.to_numpy()
        ce = np.maximum.accumulate(c.end.to_numpy())
        p_start = (g.tss - flank).to_numpy()
        p_end = (g.tss + flank).to_numpy()
        idx = np.searchsorted(cs, p_end, side="left") - 1
        hit = (idx >= 0) & (ce[np.clip(idx, 0, None)] > p_start)
        out.loc[g.index] = hit
    return out

genes["promoter_cgi"] = promoter_in_cgi(genes, cgi)
print("fraction of genes with promoter CGI:", round(genes.promoter_cgi.mean(), 3))
ann = genes.set_index("gene_name")[["chrom", "gene_type", "gene_span", "promoter_cgi", "is_sex_linked"]]
ann.reset_index().to_csv(RES / "gene_annotations.csv", index=False)

# %% 3. Load DE results (paper background) and machinery set
con = sqlite3.connect(DB_BG)
de = pd.read_sql("""SELECT gene_name, celltype, treatment_time, treatment_drug AS drug,
                           baseMean, log2FoldChange, stat, padj
                    FROM de_results WHERE control_time = '0h'""", con)
de["hours"] = de["treatment_time"].str.extract(r"(\d+)")[0].astype(int)
de["ct"] = de["celltype"].map(TARGETS)

geneset = pd.read_csv(RES / "epigenetic_machinery_geneset.csv")
machinery = set(geneset[geneset["in_db"]].gene_name)

for ct in TARGETS.values():
    bg_genes = set(de[de.ct == ct].gene_name)
    mapped = bg_genes & set(ann.index)
    print(f"{ct}: background {len(bg_genes)}, found in annotation {len(mapped)} "
          f"({100 * len(mapped) / len(bg_genes):.1f}%); machinery mapped: "
          f"{len(machinery & mapped)}/{len(machinery & bg_genes)}")
sex_m = sorted(g for g in machinery if g in ann.index and ann.loc[g, "is_sex_linked"])
print("machinery genes on X/Y:", sex_m)

# %% 4. Descriptive comparison: background vs machinery vs DEG classes (Psilo)
rows = []
for ct in TARGETS.values():
    sub = de[de.ct == ct]
    P1 = sub[(sub.drug == "Psilo") & (sub.hours == 1)].set_index("gene_name")
    P72 = sub[(sub.drug == "Psilo") & (sub.hours == 72)].set_index("gene_name")
    bgset = set(sub.gene_name)
    classes = {
        "background": bgset,
        "machinery": machinery & bgset,
        "Psilo up 1h": set(P1[(P1.padj < 0.05) & (P1.log2FoldChange > 0)].index),
        "Psilo up 72h": set(P72[(P72.padj < 0.05) & (P72.log2FoldChange > 0)].index),
        "Psilo down 72h": set(P72[(P72.padj < 0.05) & (P72.log2FoldChange < 0)].index),
    }
    for name, gs in classes.items():
        a = ann.reindex(sorted(gs)).dropna(subset=["gene_span"])
        rows.append(dict(ct=ct, group=name, n=len(gs), n_annotated=len(a),
                         median_span_kb=round(a.gene_span.median() / 1e3, 1),
                         pct_promoter_cgi=round(100 * a.promoter_cgi.astype(float).mean(), 1),
                         pct_X_or_Y=round(100 * a.is_sex_linked.astype(float).mean(), 1)))
summary = pd.DataFrame(rows)
summary.to_csv(RES / "annotation_summary.csv", index=False)
print("\n", summary.to_string(index=False))

# %% 5. Matched null for the whole machinery: expression only vs expression + length vs no X/Y
def make_strata(d, use_length):
    e = pd.qcut(np.log10(d.baseMean + 1), N_EXPR_BINS, labels=False, duplicates="drop").astype(str)
    if not use_length:
        return e
    l = pd.qcut(np.log10(d.gene_span), N_LEN_BINS, labels=False, duplicates="drop").astype(str)
    return e + "_" + l

def null_z(d, strata, genes):
    score = d["stat"]
    inset = score.index.intersection(genes)
    k = len(inset)
    obs = score[inset].mean()
    pool_mask = ~score.index.isin(inset)
    sv, xv = strata.to_numpy(), score.to_numpy()
    null_sum = np.zeros(NDRAW)
    for s, n in strata[inset].value_counts().items():
        pos = np.flatnonzero(pool_mask & (sv == s))
        assert len(pos) > 0, f"empty stratum {s}"
        idx = rng.integers(0, len(pos), size=(NDRAW, n))
        null_sum += xv[pos][idx].sum(axis=1)
    null = null_sum / k
    z = (obs - null.mean()) / null.std()
    p1 = (1 + (null >= obs).sum()) / (NDRAW + 1) if z > 0 else (1 + (null <= obs).sum()) / (NDRAW + 1)
    return dict(n_set=k, z=z, p_emp=min(1.0, 2 * p1))

rows = []
for ct in TARGETS.values():
    for drug in ["Psilo", "Ket"]:
        for h in HOURS:
            d = (de[(de.ct == ct) & (de.drug == drug) & (de.hours == h)]
                 .dropna(subset=["stat"]).drop_duplicates("gene_name").set_index("gene_name")
                 .join(ann[["gene_span", "promoter_cgi", "is_sex_linked"]], how="inner"))
            d = d[d.gene_span > 0]
            d_noxy = d[~d.is_sex_linked.astype(bool)]
            variants = {
                "expression only": (d, make_strata(d, False), machinery),
                "expression + length": (d, make_strata(d, True), machinery),
                "expression + length, no X/Y": (d_noxy, make_strata(d_noxy, True), machinery),
            }
            for vname, (dd, strata, gset) in variants.items():
                r = null_z(dd, strata, gset)
                rows.append(dict(ct=ct, drug=drug, hours=h, variant=vname, **r))
    print("done:", ct)

res = pd.DataFrame(rows)
res.to_csv(RES / "length_matched_null.csv", index=False)
print("\nz-scores (whole machinery vs matched random sets):")
print(res.pivot_table(index=["ct", "drug", "variant"], columns="hours", values="z").round(2).to_string())
print("\nempirical p (two-sided):")
print(res.pivot_table(index=["ct", "drug", "variant"], columns="hours", values="p_emp").round(4).to_string())