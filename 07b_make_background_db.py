# %%
# 07b: Copy of the sqlite restricted to the paper's background genes (two target cell types)
# Input : data/differential_expression_results.sqlite, results/background_genes.csv
# Output: data/de_results_background.sqlite (table de_results, same columns)
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
OUT_DB = ROOT / "data" / "de_results_background.sqlite"
RES = ROOT / "results"

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
EXPECTED = {"L4/5 IT": 14951, "L2/3 IT": 14557}

bg = pd.read_csv(RES / "background_genes.csv")
bg["celltype"] = bg["ct"].map({v: k for k, v in TARGETS.items()})
bg = bg[["celltype", "gene_name"]]

con = sqlite3.connect(DB)
ph = ",".join("?" * len(TARGETS))
df = pd.read_sql(f"SELECT * FROM de_results WHERE celltype IN ({ph})", con, params=list(TARGETS))
print("rows before:", len(df))
df = df.merge(bg, on=["celltype", "gene_name"], how="inner")
print("rows after :", len(df))

n_genes = df.groupby(["celltype", "comparison_id"]).gene_name.nunique()
for full, short in TARGETS.items():
    assert set(n_genes.loc[full]) == {EXPECTED[short]}, f"unexpected gene count in {short}"
print("gene counts per comparison match the paper")

if OUT_DB.exists():
    OUT_DB.unlink()
dst = sqlite3.connect(OUT_DB)
df.to_sql("de_results", dst, index=False)
dst.close()
print("saved:", OUT_DB)