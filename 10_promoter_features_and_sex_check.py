# %%
# 10: Promoter features of DEG classes (adjusted) + sex-composition diagnostic
# Inputs : data/de_results_background.sqlite, results/gene_annotations.csv
# Outputs: results/promoter_feature_models.csv, results/sex_gene_check.csv
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
DATA, RES = ROOT / "data", ROOT / "results"
DB_BG = DATA / "de_results_background.sqlite"

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
HOURS = [1, 2, 4, 24, 72]
MIN_EVENTS = 30

ann = pd.read_csv(RES / "gene_annotations.csv").drop_duplicates("gene_name").set_index("gene_name")
con = sqlite3.connect(DB_BG)
de = pd.read_sql("""SELECT gene_name, celltype, treatment_time, treatment_drug AS drug,
                           baseMean, log2FoldChange, padj
                    FROM de_results WHERE control_time = '0h'""", con)
de["hours"] = de["treatment_time"].str.extract(r"(\d+)")[0].astype(int)
de["ct"] = de["celltype"].map(TARGETS)

z = lambda s: (s - s.mean()) / s.std()

# %% A. Logistic models: is a promoter CGI / gene length associated with being a DEG, given expression?
rows = []
for ct in TARGETS.values():
    for drug in ["Psilo", "Ket"]:
        for h in HOURS:
            d = (de[(de.ct == ct) & (de.drug == drug) & (de.hours == h)]
                 .dropna(subset=["padj"]).drop_duplicates("gene_name").set_index("gene_name")
                 .join(ann[["gene_span", "promoter_cgi", "is_sex_linked"]], how="inner"))
            d = d[d.gene_span > 0]
            expr = z(np.log10(d.baseMean + 1))
            X = sm.add_constant(pd.DataFrame({
                "CGI_promoter": d.promoter_cgi.astype(float),
                "length_z": z(np.log10(d.gene_span)),
                "expr_z": expr,
                "expr_z_sq": expr ** 2,
                "sex_linked": d.is_sex_linked.astype(float)}, index=d.index))
            outcomes = {"up": (d.padj < 0.05) & (d.log2FoldChange > 0),
                        "down": (d.padj < 0.05) & (d.log2FoldChange < 0)}
            for outcome, y in outcomes.items():
                y = y.astype(int)
                if y.sum() < MIN_EVENTS:
                    continue
                try:
                    fit = sm.Logit(y, X).fit(disp=0)
                except Exception as e:
                    print("fit failed:", ct, drug, h, outcome, e)
                    continue
                ci = fit.conf_int()
                for term in ["CGI_promoter", "length_z", "sex_linked"]:
                    rows.append(dict(ct=ct, drug=drug, hours=h, outcome=outcome, n_events=int(y.sum()),
                                     term=term, OR=np.exp(fit.params[term]),
                                     OR_lo=np.exp(ci.loc[term, 0]), OR_hi=np.exp(ci.loc[term, 1]),
                                     p=fit.pvalues[term]))
res = pd.DataFrame(rows)
res.to_csv(RES / "promoter_feature_models.csv", index=False)

res["OR_text"] = res.apply(lambda r: f"{r.OR:.2f} [{r.OR_lo:.2f}-{r.OR_hi:.2f}]", axis=1)
for term, label in [("CGI_promoter", "CGI promoter (OR vs no CGI)"),
                    ("length_z", "gene length (OR per 1 SD of log10 length)"),
                    ("sex_linked", "X/Y-linked (OR)")]:
    t = res[(res.term == term) & (res.drug == "Psilo")]
    print(f"\n{label}, Psilo; adjusted for expression, length, CGI, X/Y:")
    print(t.pivot_table(index=["ct", "outcome"], columns="hours", values="OR_text", aggfunc="first").to_string())
print("\nKetamine, CGI promoter (outcomes with >= 30 DEGs only):")
t = res[(res.term == "CGI_promoter") & (res.drug == "Ket")]
print(t.pivot_table(index=["ct", "outcome"], columns="hours", values="OR_text", aggfunc="first").to_string())

# %% B. Sex-composition diagnostic: Xist and Y-linked genes vs expected shift
F_FRAC = {("Psilo", 1): 2/3, ("Psilo", 2): 2/3, ("Psilo", 4): 2/4, ("Psilo", 24): 1/3, ("Psilo", 72): 2/4,
          ("Ket", 1): 2/3, ("Ket", 2): 2/4, ("Ket", 4): 2/4, ("Ket", 24): 2/3, ("Ket", 72): 2/3}
CONTROL_F = 6 / 8
SEX_GENES = ["Xist", "Ddx3y", "Eif2s3y", "Uty", "Kdm5d"]

sg = de[de.gene_name.isin(SEX_GENES)]
print("\nsex genes present in the background:", sorted(sg.gene_name.unique()))
obs = sg.pivot_table(index=["ct", "gene_name"], columns=["drug", "hours"], values="log2FoldChange")
print("\nobserved log2FC (shrunk) vs 0h:")
print(obs.round(2).to_string())

exp_x = {k: np.log2(f / CONTROL_F) for k, f in F_FRAC.items()}
exp_y = {k: np.log2((1 - f) / (1 - CONTROL_F)) for k, f in F_FRAC.items()}
exp = pd.DataFrame({"expected_Xist": exp_x, "expected_Y_genes": exp_y}).T
exp.columns = pd.MultiIndex.from_tuples(exp.columns, names=["drug", "hours"])
print("\nexpected log2FC if the shift reflects only the sex ratio (animal counts):")
print(exp.round(2).to_string())
obs.to_csv(RES / "sex_gene_check.csv")