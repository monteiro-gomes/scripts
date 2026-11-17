import sqlite3
from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parents[1]   # ANALYSIS folder, when run as a script
except NameError:
    ROOT = Path.cwd()                            # fallback for interactive cells
    if ROOT.name == "scripts":
        ROOT = ROOT.parent

DB = ROOT / "data" / "differential_expression_results.sqlite"
print("Looking for:", DB, "| exists:", DB.exists())
assert DB.exists(), "File not found: check the name and that it is inside ANALYSIS\\data"

con = sqlite3.connect(DB)
cur = con.cursor()
for (name,) in cur.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall():
    print("==", name)
    for col in cur.execute(f"PRAGMA table_info({name})").fetchall():
        print("  ", col[1], col[2])
    print("  rows:", cur.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0])
    print("  sample:", cur.execute(f"SELECT * FROM {name} LIMIT 2").fetchall())