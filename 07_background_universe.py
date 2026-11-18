# %%
# 07: Rebuild the paper's background universe (baseMean > 10 at any timepoint/drug)
# and check it against their reported sizes (L4/5 IT: 14951, L2/3 IT: 14557).
# Output: results/background_genes.csv
import sqlite3
from pathlib import Path
import pandas as pd

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
DB = ROOT / "data" / "differential_expression_results.sqlite"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
EXPECTED = {"L4/5 IT": 14951, "L2/3 IT": 14557}

geneset = pd.read_csv(RES / "epigenetic_machinery_geneset.csv")
machinery = set(geneset[geneset["in_db"]].gene_name)

# %% Load baseMean and padj for all comparisons
con = sqlite3.connect(DB)
ph = ",".join("?" * len(TARGETS))
df = pd.read_sql(
    f"SELECT gene_name, celltype, comparison_id, baseMean, padj FROM de_results WHERE celltype IN ({ph})",
    con, params=list(TARGETS))
df["ct"] = df["celltype"].map(TARGETS)
print("genes per cell type and comparison (should be 26100):")
print(df.groupby(["ct", "comparison_id"]).gene_name.nunique().groupby("ct").agg(["min", "max"]))

# %% Candidate definitions of "expressed"
g = df.groupby(["ct", "gene_name"])
summary = g.agg(max_bm=("baseMean", "max"), mean_bm=("baseMean", "mean"),
                n_nan_padj=("padj", lambda x: x.isna().sum())).reset_index()
rows = []
for ct, exp in EXPECTED.items():
    s = summary[summary.ct == ct]
    rows.append(dict(ct=ct, paper=exp,
                     max_gt10=(s.max_bm > 10).sum(),
                     max_ge10=(s.max_bm >= 10).sum(),
                     mean_gt10=(s.mean_bm > 10).sum()))
print("\nBackground size under each definition:")
print(pd.DataFrame(rows).to_string(index=False))

# %% Save the 'max baseMean > 10' background and check machinery coverage
bg = summary[summary.max_bm > 10].copy()
bg.to_csv(RES / "background_genes.csv", index=False)
for ct in EXPECTED:
    b = bg[bg.ct == ct]
    inb = machinery & set(b.gene_name)
    print(f"\n{ct}: background = {len(b)} | machinery in background = {len(inb)}/{len(machinery)}"
          f" | background genes with NaN padj in >=1 comparison = {(b.n_nan_padj > 0).sum()}")
    print("  machinery genes outside background:", sorted(machinery - inb))