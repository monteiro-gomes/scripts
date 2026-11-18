# %%
# 05: Specificity checks for the machinery signal
# Input : sqlite, results/epigenetic_machinery_geneset.csv
# Output: results/machinery_null_tests.csv, results/machinery_gene_patterns.csv,
#         figures/machinery_null_z.png
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
DB = ROOT / "data" / "differential_expression_results.sqlite"
RES, FIG = ROOT / "results", ROOT / "figures"

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
HOURS = [1, 2, 4, 24, 72]
NBINS, NDRAW = 20, 5000
rng = np.random.default_rng(42)

# %% Frozen gene set and data (universe = tested genes, i.e. non-null padj)
geneset = pd.read_csv(RES / "epigenetic_machinery_geneset.csv")
geneset = geneset[geneset["in_db"]]
machinery = set(geneset.gene_name)
groups = geneset.groupby("group", sort=False)["gene_name"].apply(list).to_dict()
sets = {"ALL_MACHINERY": sorted(machinery)} | groups

con = sqlite3.connect(DB)
ph = ",".join("?" * len(TARGETS))
de = pd.read_sql(
    f"""SELECT gene_name, celltype, treatment_time, treatment_drug AS drug,
               baseMean, log2FoldChange, lfcSE, stat, padj
        FROM de_results WHERE celltype IN ({ph}) AND control_time = '0h'""",
    con, params=list(TARGETS))
de["hours"] = de["treatment_time"].str.extract(r"(\d+)")[0].astype(int)
de["ct"] = de["celltype"].map(TARGETS)
de = de.dropna(subset=["stat", "padj"])

# %% Null model: random gene sets with the same baseMean-bin composition
def null_test(score, bins, genes, ndraw=NDRAW):
    inset = score.index.intersection(genes)
    k = len(inset)
    if k < 5:
        return None
    obs = score[inset].mean()
    pool, pool_bins = score.drop(inset), bins.drop(inset)
    null_sum = np.zeros(ndraw)
    for b, n in bins[inset].value_counts().items():
        vals = pool[pool_bins == b].to_numpy()
        assert len(vals) > 0, f"empty bin {b}"
        idx = rng.integers(0, len(vals), size=(ndraw, n))   # sampling with replacement (pools are large)
        null_sum += vals[idx].sum(axis=1)
    null = null_sum / k
    z = (obs - null.mean()) / null.std()
    p_one = (1 + (null >= obs).sum()) / (ndraw + 1) if obs > null.mean() \
        else (1 + (null <= obs).sum()) / (ndraw + 1)
    return dict(n_set=k, obs_mean=obs, null_mean=null.mean(), z=z, p_emp=min(1.0, 2 * p_one))

rows = []
for ct in TARGETS.values():
    for h in HOURS:
        sub = lambda drug: (de[(de.ct == ct) & (de.drug == drug) & (de.hours == h)]
                            .drop_duplicates("gene_name").set_index("gene_name"))
        P, K = sub("Psilo"), sub("Ket")
        common = P.index.intersection(K.index)
        contrast = (P.loc[common, "log2FoldChange"] - K.loc[common, "log2FoldChange"]) / \
                   np.sqrt(P.loc[common, "lfcSE"] ** 2 + K.loc[common, "lfcSE"] ** 2)
        contrast_bm = (P.loc[common, "baseMean"] + K.loc[common, "baseMean"]) / 2
        analyses = {"Psilo": (P["stat"], P["baseMean"]),
                    "Ket": (K["stat"], K["baseMean"]),
                    "Psilo_minus_Ket": (contrast, contrast_bm)}
        for an, (score, bm) in analyses.items():
            bins = pd.qcut(np.log10(bm + 1), NBINS, labels=False, duplicates="drop")
            for sname, genes in sets.items():
                r = null_test(score, bins, set(genes))
                if r:
                    rows.append(dict(ct=ct, hours=h, analysis=an, set=sname, **r))
    print("done:", ct)

res = pd.DataFrame(rows)
res["q"] = np.nan
for _, idx in res.groupby(["analysis", res["set"] == "ALL_MACHINERY"]).groups.items():
    res.loc[idx, "q"] = multipletests(res.loc[idx, "p_emp"], method="fdr_bh")[1]
res.to_csv(RES / "machinery_null_tests.csv", index=False)

allm = res[res["set"] == "ALL_MACHINERY"]
print("\nz-score vs matched random sets (whole machinery):")
print(allm.pivot_table(index=["ct", "analysis"], columns="hours", values="z").round(2).to_string())
print("\nq-values:")
print(allm.pivot_table(index=["ct", "analysis"], columns="hours", values="q").round(4).to_string())

# %% Figure
fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
for ax, ct in zip(axes, TARGETS.values()):
    for an, c in [("Psilo", "tab:red"), ("Ket", "tab:blue"), ("Psilo_minus_Ket", "k")]:
        s = allm[(allm.ct == ct) & (allm.analysis == an)].sort_values("hours")
        ax.plot(range(len(s)), s["z"], "-o", color=c, label=an)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xticks(range(len(HOURS)))
    ax.set_xticklabels([f"{h}h" for h in HOURS])
    ax.set_title(ct)
axes[0].set_ylabel("z vs matched random sets")
axes[0].legend(fontsize=7)
fig.savefig(FIG / "machinery_null_z.png", dpi=200, bbox_inches="tight")
plt.show()

# %% Which genes? Patterns at 1h vs 72h (psilocybin) + genome-wide reactivation test
def get(ct, h):
    return (de[(de.ct == ct) & (de.drug == "Psilo") & (de.hours == h)]
            .drop_duplicates("gene_name").set_index("gene_name")
            [["log2FoldChange", "stat", "padj"]].add_suffix(f"_{h}h"))

out = []
for ct in TARGETS.values():
    w = get(ct, 1).join(get(ct, 72), how="inner")
    rho, _ = spearmanr(w["stat_1h"], w["stat_72h"])
    mach = w.loc[w.index.intersection(machinery)]
    rho_m, _ = spearmanr(mach["stat_1h"], mach["stat_72h"])
    g1 = (w.padj_1h < 0.05) & (w.log2FoldChange_1h > 0)
    g72 = (w.padj_72h < 0.05) & (w.log2FoldChange_72h > 0)
    a, b, c, d = (g1 & g72).sum(), (g1 & ~g72).sum(), (~g1 & g72).sum(), (~g1 & ~g72).sum()
    odds, p = fisher_exact([[a, b], [c, d]])
    print(f"\n{ct}: Spearman(stat 1h, stat 72h) genome-wide = {rho:.3f}, machinery = {rho_m:.3f}")
    print(f"  up@1h = {g1.sum()}, up@72h = {g72.sum()}, both = {a}, odds ratio = {odds:.2f}, p = {p:.2e}")
    m = mach.copy()
    u1 = (m.padj_1h < 0.05) & (m.log2FoldChange_1h > 0)
    u72 = (m.padj_72h < 0.05) & (m.log2FoldChange_72h > 0)
    m["pattern"] = np.select([u1 & u72, u1, u72],
                             ["up_1h_and_72h", "up_1h_only", "up_72h_only"], "other")
    m["ct"] = ct
    out.append(m.join(geneset.set_index("gene_name")["group"]).reset_index())

pat = pd.concat(out)
pat.to_csv(RES / "machinery_gene_patterns.csv", index=False)
for ct, g in pat.groupby("ct"):
    print(f"\n=== {ct} ===")
    print(g["pattern"].value_counts().to_string())
    for p_ in ["up_1h_and_72h", "up_72h_only", "up_1h_only"]:
        names = g[g.pattern == p_].sort_values("group")
        print(f"{p_} ({len(names)}):", ", ".join(names.gene_name))