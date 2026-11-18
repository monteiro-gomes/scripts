# %%
# 04: Trajectories and set-level tests of the frozen epigenetic machinery gene set.
# Input : data/differential_expression_results.sqlite, results/epigenetic_machinery_geneset.csv
# Output: results/machinery_gene_stats_long.csv, results/machinery_set_tests.csv,
#         figures/machinery_heatmap_*.png, figures/machinery_set_NES.png
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
DB = ROOT / "data" / "differential_expression_results.sqlite"
RES, FIG = ROOT / "results", ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
HOURS = [1, 2, 4, 24, 72]
DRUGS = ["Psilo", "Ket"]
COLS = pd.MultiIndex.from_tuples([(d, h) for d in DRUGS for h in HOURS], names=["drug", "hours"])
MIN_SET = 5   # GSEA min_size: groups with fewer genes in the universe are not tested

# %% Frozen gene set
geneset = pd.read_csv(RES / "epigenetic_machinery_geneset.csv")
geneset = geneset[geneset["in_db"]].reset_index(drop=True)
groups = geneset.groupby("group", sort=False)["gene_name"].apply(list).to_dict()
machinery = set(geneset["gene_name"])
print("Genes:", len(machinery), "| groups:", len(groups))
print("Groups below MIN_SET (not tested):", [g for g, v in groups.items() if len(v) < MIN_SET])

# %% Load DE results for the two cell types (all genes: needed for the ranked lists)
con = sqlite3.connect(DB)
ph = ",".join("?" * len(TARGETS))
de = pd.read_sql(
    f"""SELECT gene_name, celltype, treatment_time, treatment_drug AS drug,
               baseMean, log2FoldChange, stat, padj
        FROM de_results
        WHERE celltype IN ({ph}) AND control_time = '0h'""",
    con, params=list(TARGETS))
de["hours"] = de["treatment_time"].str.extract(r"(\d+)")[0].astype(int)
de["ct"] = de["celltype"].map(TARGETS)
print(de.shape)

# %% Long table for machinery genes (for later scripts)
long = de[de.gene_name.isin(machinery)].merge(geneset[["gene_name", "group"]], on="gene_name")
long.drop(columns=["celltype", "treatment_time"]).to_csv(RES / "machinery_gene_stats_long.csv", index=False)

# %% A. Gene-level heatmaps (colour = Wald stat, dot = padj < 0.05)
def plot_heatmap(ct):
    sub = de[(de.ct == ct) & de.gene_name.isin(machinery)]
    Z = (sub.pivot_table(index="gene_name", columns=["drug", "hours"], values="stat")
            .reindex(index=geneset.gene_name.values, columns=COLS))
    P = (sub.pivot_table(index="gene_name", columns=["drug", "hours"], values="padj")
            .reindex(index=geneset.gene_name.values, columns=COLS))
    fig, ax = plt.subplots(figsize=(6, 0.16 * len(Z) + 2))
    sns.heatmap(Z, cmap="vlag", center=0, vmin=-6, vmax=6, ax=ax, yticklabels=True,
                cbar_kws={"label": "Wald stat vs 0h", "shrink": 0.4, "location": "top"})
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if P.iat[i, j] < 0.05:
                ax.text(j + 0.5, i + 0.5, "•", ha="center", va="center", fontsize=6)
    grp = geneset["group"].values
    starts = [0] + [i for i in range(1, len(grp)) if grp[i] != grp[i - 1]] + [len(grp)]
    for b in starts[1:-1]:
        ax.axhline(b, color="k", lw=0.6)
    for s, e in zip(starts[:-1], starts[1:]):
        ax.text(len(COLS) + 0.3, (s + e) / 2, grp[s], fontsize=6, va="center", clip_on=False)
    ax.axvline(len(HOURS), color="k", lw=1.5)
    ax.set_xticklabels([f"{d} {h}h" for d, h in COLS], rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"{ct}: epigenetic machinery", fontsize=9, pad=40)
    fig.savefig(FIG / f"machinery_heatmap_{ct.replace('/', '').replace(' ', '_')}.png",
                dpi=200, bbox_inches="tight")
    plt.show()

for ct in TARGETS.values():
    plot_heatmap(ct)

# %% B. Set-level tests: GSEA prerank (threshold-free) + Fisher up/down (padj < 0.05)
# One universe for both tests: genes with non-null stat AND padj (tested genes).
sets = {"ALL_MACHINERY": sorted(machinery)} | groups
rows = []
for ct in TARGETS.values():
    for drug in DRUGS:
        for h in HOURS:
            d = (de[(de.ct == ct) & (de.drug == drug) & (de.hours == h)]
                 .dropna(subset=["stat", "padj"]).drop_duplicates("gene_name"))
            uni_genes = set(d.gene_name)
            up = set(d[(d.padj < 0.05) & (d.log2FoldChange > 0)].gene_name)
            dn = set(d[(d.padj < 0.05) & (d.log2FoldChange < 0)].gene_name)

            # Ranked list: tiny reproducible jitter to break ties (low-count genes sit at stat ~ 0)
            rnk = d.set_index("gene_name")["stat"]
            jitter = np.random.default_rng(42).normal(0, 1e-6, len(rnk))
            rnk = (rnk + jitter).sort_values(ascending=False)

            pre = gp.prerank(rnk=rnk, gene_sets=sets, min_size=MIN_SET, max_size=1000,
                             permutation_num=1000, seed=42, threads=4,
                             no_plot=True, verbose=False, outdir=None)
            r = pre.res2d.set_index("Term")
            for name, genes in sets.items():
                inset = set(genes) & uni_genes
                n = len(inset)
                out = dict(ct=ct, drug=drug, hours=h, set=name, n_set=n,
                           n_up=len(inset & up), n_dn=len(inset & dn))
                for lab, sig in (("up", up), ("dn", dn)):
                    a = len(inset & sig)
                    b = n - a
                    c = len(sig) - a
                    dd = len(uni_genes) - n - c
                    out[f"p_fisher_{lab}"] = fisher_exact([[a, b], [c, dd]], alternative="greater")[1]
                if name in r.index:
                    out["NES"] = float(r.loc[name, "NES"])
                    out["p_gsea"] = float(r.loc[name, "NOM p-val"])
                else:
                    out["NES"] = out["p_gsea"] = np.nan
                rows.append(out)
    print("done:", ct)

res = pd.DataFrame(rows)

# Benjamini-Hochberg within two families: whole set vs functional groups
for pcol in ["p_gsea", "p_fisher_up", "p_fisher_dn"]:
    qcol = pcol.replace("p_", "q_")
    res[qcol] = np.nan
    for is_all in (True, False):
        m = ((res.set == "ALL_MACHINERY") == is_all) & res[pcol].notna()
        res.loc[m, qcol] = multipletests(res.loc[m, pcol], method="fdr_bh")[1]
res.to_csv(RES / "machinery_set_tests.csv", index=False)

# %% Summary: whole-set results
cols = ["ct", "drug", "hours", "n_set", "n_up", "n_dn", "NES", "p_gsea", "q_gsea"]
print(res[res.set == "ALL_MACHINERY"][cols].round(4).to_string(index=False))

# %% NES heatmap per cell type (star = q < 0.05)
order = ["ALL_MACHINERY"] + list(groups)
fig, axes = plt.subplots(1, 2, figsize=(11, 0.35 * len(order) + 2), sharey=True)
for ax, ct in zip(axes, TARGETS.values()):
    sub = res[res.ct == ct]
    NES = sub.pivot_table(index="set", columns=["drug", "hours"], values="NES").reindex(index=order, columns=COLS)
    Q = sub.pivot_table(index="set", columns=["drug", "hours"], values="q_gsea").reindex(index=order, columns=COLS)
    sns.heatmap(NES, cmap="vlag", center=0, vmin=-3, vmax=3, ax=ax, cbar=(ct == "L2/3 IT"),
                cbar_kws={"label": "NES"})
    for i in range(Q.shape[0]):
        for j in range(Q.shape[1]):
            if Q.iat[i, j] < 0.05:
                ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=10)
    ax.axvline(len(HOURS), color="k", lw=1.5)
    ax.set_xticklabels([f"{d} {h}h" for d, h in COLS], rotation=90)
    ax.set_title(ct)
    ax.set_xlabel("")
    ax.set_ylabel("")
fig.savefig(FIG / "machinery_set_NES.png", dpi=200, bbox_inches="tight")
plt.show()