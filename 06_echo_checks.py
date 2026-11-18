# %%
# 06: Is the 72h wave an "echo" of 1h, and is the machinery special in it?
# Output: results/echo_matched_null.csv, figures/echo_corr_matrices.png
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata, spearmanr

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
COLS = pd.MultiIndex.from_product([["Psilo", "Ket"], HOURS], names=["drug", "hours"])
NBINS, NDRAW = 20, 5000
rng = np.random.default_rng(42)

geneset = pd.read_csv(RES / "epigenetic_machinery_geneset.csv")
machinery = set(geneset[geneset["in_db"]].gene_name)

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

def wide(ct, value):
    sub = de[de.ct == ct]
    return sub.pivot_table(index="gene_name", columns=["drug", "hours"], values=value).reindex(columns=COLS)

# %% Matched-null test for concordance between two score vectors
def matched_rho_null(x, y, bm, genes, ndraw=NDRAW):
    bins = pd.Series(pd.qcut(np.log10(bm + 1), NBINS, labels=False, duplicates="drop"), index=bm.index)
    inset = x.index.intersection(genes)
    obs = spearmanr(x[inset], y[inset])[0]
    in_mask = x.index.isin(inset)
    bvals = bins.to_numpy()
    xv, yv = x.to_numpy(), y.to_numpy()
    cols = []
    for b, n in bins[inset].value_counts().items():
        pos = np.flatnonzero(~in_mask & (bvals == b))
        cols.append(pos[rng.integers(0, len(pos), size=(ndraw, n))])
    sel = np.concatenate(cols, axis=1)
    rx, ry = rankdata(xv[sel], axis=1), rankdata(yv[sel], axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    null = (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
    return dict(n_set=len(inset), rho_set=obs, rho_null_mean=null.mean(), rho_null_sd=null.std(),
                z=(obs - null.mean()) / null.std(), p_one_sided=(1 + (null >= obs).sum()) / (ndraw + 1))

# %% Run per cell type
rows = []
fig, axes = plt.subplots(2, 2, figsize=(11, 10), gridspec_kw={"hspace": 0.55})
for r, ct in enumerate(TARGETS.values()):
    S, B = wide(ct, "stat"), wide(ct, "baseMean")
    L, E = wide(ct, "log2FoldChange"), wide(ct, "lfcSE")
    ok = S.dropna().index
    for t in (B, L, E):
        ok = ok.intersection(t.dropna().index)
    S, B, L, E = S.loc[ok], B.loc[ok], L.loc[ok], E.loc[ok]
    print(f"\n{ct}: {len(ok)} genes with all 10 conditions")

    # A. correlation matrix of Wald stats (all 10 conditions)
    C = S.corr(method="spearman")
    C.index = C.columns = [f"{d} {h}h" for d, h in C.columns]
    print("Spearman matrix, genome-wide (Psilo rows):")
    print(C.iloc[:5].round(2).to_string())
    sns.heatmap(C, vmin=-0.2, vmax=1, cmap="viridis", ax=axes[r, 0], annot=True, fmt=".2f",
                annot_kws={"size": 5}, cbar=False)
    axes[r, 0].set_title(f"{ct}: Wald stat vs 0h (shared control)")

    # B. control-free contrast: (Psilo - Ket) / se, the shared 0h group cancels
    D = (L["Psilo"] - L["Ket"]) / np.sqrt(E["Psilo"] ** 2 + E["Ket"] ** 2)
    CD = D.corr(method="spearman")
    CD.index = CD.columns = [f"{h}h" for h in CD.columns]
    print("Control-free (Psilo - Ket) correlation:")
    print(CD.round(2).to_string())
    sns.heatmap(CD, vmin=-0.2, vmax=1, cmap="viridis", ax=axes[r, 1], annot=True, fmt=".2f", cbar=False)
    axes[r, 1].set_title(f"{ct}: Psilo - Ket contrast (control cancels)")

    # Partial correlation of 1h and 72h given 24h (contrast version)
    rho3 = D[[1, 24, 72]].corr(method="spearman")
    r12, r13, r23 = rho3.loc[1, 72], rho3.loc[1, 24], rho3.loc[24, 72]
    pcor = (r12 - r13 * r23) / np.sqrt((1 - r13 ** 2) * (1 - r23 ** 2))
    print(f"partial rho(1h, 72h | 24h), contrast: {pcor:.2f}")

    # C. is the machinery more concordant than expression-matched random sets?
    bm = B.mean(axis=1)
    tests = {
        "stat 1h vs 72h": (S[("Psilo", 1)], S[("Psilo", 72)]),
        "stat 1h vs 24h (reference)": (S[("Psilo", 1)], S[("Psilo", 24)]),
        "contrast 1h vs 72h": (D[1], D[72]),
        "contrast 1h vs 24h (reference)": (D[1], D[24]),
        "contrast 24h vs 72h (reference)": (D[24], D[72]),
    }
    for name, (x, y) in tests.items():
        res = matched_rho_null(x, y, bm, machinery)
        rows.append(dict(ct=ct, test=name, **res))

fig.savefig(FIG / "echo_corr_matrices.png", dpi=200, bbox_inches="tight")
plt.show()

out = pd.DataFrame(rows)
out.to_csv(RES / "echo_matched_null.csv", index=False)
print("\nMachinery concordance vs expression-matched random sets:")
print(out.round(4).to_string(index=False))