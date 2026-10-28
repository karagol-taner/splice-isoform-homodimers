import re, pandas as pd
s = pd.read_csv("data/human_sp_subunit.tsv.gz", sep="\t", dtype=str).fillna(""); s.columns = ["acc", "subunit"]
cx = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t").set_index("acc")
p = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
can = sorted(set(p.canonical_acc))
s = s[s.acc.isin(can)].copy()
s["general"] = s.subunit.apply(lambda t: re.sub(r"\[Isoform[^\]]*\]:[^\[]*", " ", t))
s["homodimer"] = s.general.str.contains(r"\bHomodimer|\bhomodimer", regex=True)
s["homo_any"] = s.general.str.contains(r"[Hh]omo(?:di|tri|tetra|penta|hexa|hepta|octa|multi|oligo)mer", regex=True)
s["monomer_only"] = s.general.str.contains(r"^\s*SUBUNIT:\s*Monomer\b|\bMonomer\.", regex=True) & ~s.homo_any
s["hc"] = s.acc.map(cx.hc_homodimer).fillna(False).astype(bool)
s.to_csv("results/canonical_subunit_calibration.tsv", sep="\t", index=False)
print("canonicals:", len(s))
