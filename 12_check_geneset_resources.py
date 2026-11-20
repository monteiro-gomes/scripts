# %%
# 12: Check library availability and API versions (needs internet)
# Output: results/enrichr_mouse_libraries.txt (+ printed summary)
import sys
from pathlib import Path
import gseapy as gp
import decoupler as dc

try:
    ROOT = Path(__file__).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()
    if ROOT.name == "scripts":
        ROOT = ROOT.parent
RES = ROOT / "results"

print("python", sys.version.split()[0])
print("gseapy", gp.__version__)
print("decoupler", getattr(dc, "__version__", "?"))

# %% Enrichr mouse libraries mentioning TFs / chromatin / epigenetic terms
names = gp.get_library_name(organism="Mouse")
(RES / "enrichr_mouse_libraries.txt").write_text("\n".join(names), encoding="utf-8")
print("\nmouse libraries available:", len(names), "(full list saved to results/enrichr_mouse_libraries.txt)")
kw = ["chea", "encode", "histone", "chromatin", "polycomb", "epigen", "methyl", "roadmap",
      "h3k", "mecp", "jaspar", "trrust", "tf_", "tf-"]
hits = [n for n in names if any(k in n.lower() for k in kw)]
print("\nlibraries matching chromatin/TF keywords:")
for n in hits:
    print("  ", n)

# %% CollecTRI (TF-target network) for mouse
net = None
try:
    if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
        net = dc.op.collectri(organism="mouse")
    elif hasattr(dc, "get_collectri"):
        net = dc.get_collectri(organism="mouse", split_complexes=False)
    else:
        print("\nNo known CollecTRI function found in this decoupler version")
except Exception as e:
    print("\nCollecTRI download failed:", repr(e))

if net is not None:
    print("\nCollecTRI mouse:", net.shape, "| columns:", list(net.columns))
    print("TFs:", net["source"].nunique(), "| targets:", net["target"].nunique())
    for tf in ["Ezh2", "Suz12", "Mecp2", "Rest", "Tet1", "Brd4", "Hdac1", "Dnmt1",
               "Creb1", "Fos", "Jun", "Npas4", "Mef2c", "Egr1", "Smad3"]:
        print(f"  {tf}: {(net['source'] == tf).sum()} targets")