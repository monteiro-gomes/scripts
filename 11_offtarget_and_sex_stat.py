# %%
# 11: Off-target (ambient-RNA proxy) fraction vs DEG classes; sex genes with unshrunk stat
# Inputs : data/differential_expression_results.sqlite (gene_expression table),
#          data/de_results_background.sqlite, results/gene_annotations.csv
# Outputs: results/offtarget_models.csv
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

TARGETS = {"006 L4/5 IT CTX Glut": "L4/5 IT", "007 L2/3 IT CTX Glut": "L2/3 IT"}
MIN_EVENTS = 100
z = lambda s: (s - s.mean()) / s.std()

# %% Baseline expression across all 16 cell types (untreated 0h)
con_full = sqlite3.connect(DATA / "differential_expression_results.sqlite")
expr = pd.read_sql("""SELECT gene_name, celltype, mean_expression FROM gene_expression
                      WHERE "time" = '0h' AND drug = 'none'""", con_full)
wide = expr.pivot_table(index="gene_name", columns="celltype", values="mean_expression")
print("baseline expression table:", wide.shape, "| min/max:", wide.min().min(), wide.max().max())
total = wide.sum(axis=1)

# %% DE results and annotations
ann = pd.read_csv(RES / "gene_annotations.csv").drop_duplicates("gene_name").set_index("gene_name")
con_bg = sqlite3.connect(DATA / "de_results_background.sqlite")
de = pd.read_sql("""SELECT gene_name, celltype, treatment_time, treatment_drug AS drug,
                           baseMean, log2FoldChange, stat, padj
                    FROM de_results WHERE control_time = '0h'""", con_bg)
de["hours"] = de["treatment_time"].str.extract(r"(\d+)")[0].astype(int)
de["ct"] = de["celltype"].map(TARGETS)

# %% A. Off-target fraction by DEG class, then models with and without it
CASES = [("Psilo", 1), ("Psilo", 72), ("Ket", 72)]
rows, desc = [], []
for full, ct in TARGETS.items():
    offtarget = (1 - wide[full] / (total + 1e-12)).rename("offtarget")
    for drug, h in CASES:
        d = (de[(de.ct == ct) & (de.drug == drug) & (de.hours == h)]
             .dropna(subset=["padj"]).drop_duplicates("gene_name").set_index("gene_name")
             .join(ann[["gene_span", "promoter_cgi", "is_sex_linked"]], how="inner")
             .join(offtarget, how="inner"))
        d = d[(d.gene_span > 0) & d.offtarget.notna()]
        up = (d.padj < 0.05) & (d.log2FoldChange > 0)
        dn = (d.padj < 0.05) & (d.log2FoldChange < 0)
        desc.append(dict(ct=ct, case=f"{drug} {h}h", background=round(d.offtarget.median(), 3),
                         up=round(d.offtarget[up].median(), 3) if up.sum() else np.nan,
                         down=round(d.offtarget[dn].median(), 3) if dn.sum() else np.nan,
                         n_up=int(up.sum()), n_down=int(dn.sum())))
        expr_z = z(np.log10(d.baseMean + 1))
        base_X = pd.DataFrame({"CGI_promoter": d.promoter_cgi.astype(float),
                               "length_z": z(np.log10(d.gene_span)),
                               "expr_z": expr_z, "expr_z_sq": expr_z ** 2,
                               "sex_linked": d.is_sex_linked.astype(float)}, index=d.index)
        for outcome, y in {"up": up, "down": dn}.items():
            if y.sum() < MIN_EVENTS:
                continue
            for model, X in {"without off-target": base_X,
                             "with off-target": base_X.assign(offtarget_z=z(d.offtarget))}.items():
                fit = sm.Logit(y.astype(int), sm.add_constant(X)).fit(disp=0)
                ci = fit.conf_int()
                for term in [t for t in ["CGI_promoter", "offtarget_z"] if t in fit.params.index]:
                    rows.append(dict(ct=ct, case=f"{drug} {h}h", outcome=outcome, model=model, term=term,
                                     OR=np.exp(fit.params[term]), lo=np.exp(ci.loc[term, 0]),
                                     hi=np.exp(ci.loc[term, 1])))

print("\nMedian off-target fraction (share of baseline expression outside the focal cell type):")
print(pd.DataFrame(desc).to_string(index=False))
res = pd.DataFrame(rows)
res["OR [95% CI]"] = res.apply(lambda r: f"{r.OR:.2f} [{r.lo:.2f}-{r.hi:.2f}]", axis=1)
res.to_csv(RES / "offtarget_models.csv", index=False)
print("\nCGI promoter OR, without vs with off-target adjustment:")
print(res[res.term == "CGI_promoter"].pivot_table(index=["ct", "case", "outcome"], columns="model",
      values="OR [95% CI]", aggfunc="first").to_string())
print("\nOff-target OR (per 1 SD):")
print(res[res.term == "offtarget_z"][["ct", "case", "outcome", "OR [95% CI]"]].to_string(index=False))

# %% B. Sex genes: unshrunk Wald stat (expected: Xist negative, Y genes positive if treated groups
#      have a lower female fraction than the 75%-female controls)
sg = de[de.gene_name.isin(["Xist", "Ddx3y", "Eif2s3y", "Uty", "Kdm5d"])]
print("\nUnshrunk Wald stat for sex genes:")
print(sg.pivot_table(index=["ct", "gene_name"], columns=["drug", "hours"], values="stat").round(2).to_string())
print("\npadj:")
print(sg.pivot_table(index=["ct", "gene_name"], columns=["drug", "hours"], values="padj").round(3).to_string())