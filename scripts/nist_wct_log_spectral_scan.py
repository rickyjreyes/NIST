from pathlib import Path
import pandas as pd
import re

src = Path("/mnt/data/Pasted text.txt")
df = pd.read_csv(src, dtype=str)

def clean_excel_formula(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    # NIST CSV often arrives as ="123.45"
    if s.startswith('="') and s.endswith('"'):
        s = s[2:-1]
    elif s.startswith("="):
        s = s[1:]
    s = s.strip('"')
    return s

df = df.applymap(clean_excel_formula)

# Save full cleaned Fe table.
full_out = Path("/mnt/data/Fe_NIST_lines_cleaned.csv")
df.to_csv(full_out, index=False)

# Save Fe I only. In NIST, sp_num=1 is neutral atom, Fe I.
fe1 = df[df["sp_num"].astype(str).str.strip() == "1"].copy()
fe1_out = Path("/mnt/data/FeI_NIST_lines_cleaned.csv")
fe1.to_csv(fe1_out, index=False)

# Save Fe II and Fe III too, because they are useful follow-up controls.
fe2 = df[df["sp_num"].astype(str).str.strip() == "2"].copy()
fe2_out = Path("/mnt/data/FeII_NIST_lines_cleaned.csv")
fe2.to_csv(fe2_out, index=False)

fe3 = df[df["sp_num"].astype(str).str.strip() == "3"].copy()
fe3_out = Path("/mnt/data/FeIII_NIST_lines_cleaned.csv")
fe3.to_csv(fe3_out, index=False)

print({
    "input_rows": len(df),
    "FeI_rows": len(fe1),
    "FeII_rows": len(fe2),
    "FeIII_rows": len(fe3),
    "full": str(full_out),
    "FeI": str(fe1_out),
    "FeII": str(fe2_out),
    "FeIII": str(fe3_out),
})
