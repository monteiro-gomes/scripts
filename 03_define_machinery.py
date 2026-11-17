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
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

TARGETS = ["006 L4/5 IT CTX Glut", "007 L2/3 IT CTX Glut"]

# %% A priori gene set (mouse symbols), grouped by function
MACHINERY = {
    "DNMT writers/maintenance": ["Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Uhrf1", "Uhrf2"],
    "DNA demethylation": ["Tet1", "Tet2", "Tet3", "Tdg", "Gadd45a", "Gadd45b", "Gadd45g"],
    "Methyl-DNA readers": ["Mecp2", "Mbd1", "Mbd2", "Mbd3", "Mbd4", "Zbtb33"],
    "H3K4 methylation": ["Kmt2a", "Kmt2b", "Kmt2c", "Kmt2d", "Kmt2e", "Setd1a", "Setd1b",
                         "Cxxc1", "Wdr5", "Ash2l", "Rbbp5", "Kdm5a", "Kdm5b", "Kdm5c",
                         "Kdm1a", "Kdm1b"],
    "H3K27 methylation (PRC2/PRC1/KDM6)": ["Ezh1", "Ezh2", "Suz12", "Eed", "Jarid2", "Kdm6a",
                                            "Kdm6b", "Rnf2", "Ring1", "Bmi1", "Cbx2", "Cbx4", "Cbx7"],
    "H3K9 methylation / HP1": ["Ehmt1", "Ehmt2", "Setdb1", "Atf7ip", "Suv39h1", "Suv39h2",
                               "Kdm4a", "Kdm4b", "Kdm4c", "Kdm3a", "Kdm3b", "Kdm7a", "Phf8",
                               "Phf2", "Cbx1", "Cbx3", "Cbx5", "Trim28"],
    "H3K36/H3K79 methylation": ["Nsd1", "Nsd2", "Nsd3", "Setd2", "Ash1l", "Kdm2a", "Kdm2b", "Dot1l"],
    "Arginine methylation": ["Prmt1", "Prmt5", "Carm1"],
    "Histone acetyltransferases": ["Ep300", "Crebbp", "Kat2a", "Kat2b", "Kat5", "Kat6a",
                                   "Kat6b", "Kat7", "Kat8"],
    "HDACs": [f"Hdac{i}" for i in range(1, 12)],
    "Sirtuins": [f"Sirt{i}" for i in range(1, 8)],
    "Acetyl-lysine readers": ["Brd1", "Brd2", "Brd3", "Brd4", "Brd8", "Brd9", "Brpf1", "Brpf3"],
    "Corepressor complexes": ["Rest", "Rcor1", "Rcor2", "Ncor1", "Ncor2", "Sin3a", "Sin3b",
                              "Chd3", "Chd4", "Mta1", "Mta2", "Gatad2a"],
    "BAF/SWI-SNF remodelers": ["Smarca2", "Smarca4", "Smarcb1", "Smarcc1", "Smarcc2", "Smarcd1",
                               "Arid1a", "Arid1b", "Arid2", "Pbrm1", "Actl6a", "Actl6b",
                               "Bcl7a", "Dpf1", "Dpf3", "Ss18", "Ss18l1"],
    "Other remodelers (CHD/ISWI/INO80)": ["Chd1", "Chd2", "Chd5", "Chd7", "Chd8", "Smarca5",
                                          "Smarca1", "Ino80", "Ep400", "Srcap"],
    "Histone variants/chaperones": ["H3f3a", "H3f3b", "H2afz", "H2az1", "H2afx", "Hira", "Atrx",
                                    "Daxx", "Asf1a", "Asf1b", "Chaf1a", "Chaf1b", "Nap1l1", "Nap1l2"],
    "Metabolic-epigenetic link": ["Acly", "Acss2", "Mat2a", "Mat2b", "Ahcy", "Nampt", "Nmnat1",
                                  "Idh1", "Idh2", "Ogdh", "Fh1", "Ogt", "Oga"],
}

geneset = pd.DataFrame(
    [(g, grp) for grp, genes in MACHINERY.items() for g in genes],
    columns=["gene_name", "group"])
dups = geneset[geneset.duplicated("gene_name", keep=False)]
print("Duplicated genes:\n", dups.to_string() if len(dups) else "none")
geneset = geneset.drop_duplicates("gene_name").reset_index(drop=True)
print("Genes in list:", len(geneset))

# %% Check symbols and detection (baseMean only, no effect sizes)
con = sqlite3.connect(DB)
ph_c = ",".join("?" * len(TARGETS))
ph_g = ",".join("?" * len(geneset))
q = f"""SELECT gene_name, celltype, comparison_id, baseMean, padj
        FROM de_results
        WHERE celltype IN ({ph_c}) AND gene_name IN ({ph_g})"""
df = pd.read_sql(q, con, params=TARGETS + geneset.gene_name.tolist())

det = (df.groupby(["gene_name", "celltype"])
         .agg(baseMean_mean=("baseMean", "mean"),
              n_tested=("padj", lambda x: x.notna().sum()))
         .unstack("celltype"))
det.columns = [f"{a}__{b}" for a, b in det.columns]
det = det.reset_index()

geneset = geneset.merge(det, on="gene_name", how="left")
geneset["in_db"] = geneset.iloc[:, 2:].notna().any(axis=1)

missing = geneset[~geneset.in_db]
print("\nNot found in the database (old/new symbol or not detected):")
print(missing[["gene_name", "group"]].to_string(index=False))

# %% Save the frozen list
geneset.to_csv(OUT / "epigenetic_machinery_geneset.csv", index=False)
print("\nSaved:", OUT / "epigenetic_machinery_geneset.csv")
print(geneset[geneset.in_db].groupby("group").size())