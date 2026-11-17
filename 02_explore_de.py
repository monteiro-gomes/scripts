# %%
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
con = sqlite3.connect(DB)

# %% 1. What comparisons and cell types exist?
comparisons = pd.read_sql("""
    SELECT comparison_id, control_time, control_drug, treatment_time, treatment_drug,
           COUNT(*) AS n_rows, COUNT(DISTINCT celltype) AS n_celltypes
    FROM de_results
    GROUP BY comparison_id, control_time, control_drug, treatment_time, treatment_drug
""", con)
print(comparisons.to_string())

celltypes = pd.read_sql(
    "SELECT celltype, COUNT(*) AS n_rows FROM de_results GROUP BY celltype ORDER BY celltype", con)
print(celltypes.to_string())

# %% 2. Pick our two target cell types
targets = celltypes[celltypes.celltype.str.contains(r"L2/3 IT|L4/5 IT", regex=True)]
print(targets.to_string())
TARGETS = targets.celltype.tolist()

# %% 3. Load only those cell types
placeholders = ",".join("?" * len(TARGETS))
de = pd.read_sql(
    f"SELECT * FROM de_results WHERE celltype IN ({placeholders})", con, params=TARGETS)
de["hours"] = de["treatment_time"].str.extract(r"(\d+)").astype(float)
print(de.shape)
print("fraction of NaN padj:", de["padj"].isna().mean().round(3))

# %% 4. Sanity check: number of DEGs per cell type, drug and time (padj < 0.05)
base = de[de["control_time"] == "0h"].copy()
sig = base[base["padj"] < 0.05]
counts = (sig.groupby(["celltype", "treatment_drug", "hours"])
             .agg(n_deg=("gene_name", "size"),
                  n_up=("log2FoldChange", lambda x: (x > 0).sum()),
                  n_down=("log2FoldChange", lambda x: (x < 0).sum()))
             .reset_index())
print(counts.pivot_table(index=["celltype", "treatment_drug"], columns="hours",
                         values="n_deg", fill_value=0).to_string())
tested = (base.dropna(subset=["padj"]).groupby(["celltype", "hours"])["gene_name"]
              .nunique().unstack())
print(tested.to_string())