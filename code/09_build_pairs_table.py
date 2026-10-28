import pandas as pd
d = pd.read_csv("results/isoform_pairs_classified.tsv", sep="\t")
cx = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t")
hc = dict(zip(cx.acc, cx.hc_homodimer))
r = d[d.in_afdb_length_range == True].copy()
r["canon_hc"] = r.canonical_acc.map(hc).fillna(False).astype(bool)
r["iso_hc"] = r.isoform_acc.map(hc).fillna(False).astype(bool)
r.to_csv("results/pairs_with_hc.tsv", sep="\t", index=False)

p = pd.read_csv("results/pairs_with_hc.tsv", sep="\t")
mo = pd.read_csv("data/afdb_monomer_per_accession.tsv", sep="\t").set_index("acc")
mf = pd.read_csv("results/manifest_per_accession.tsv", sep="\t").set_index("acc")
for side, col in [("canon", "canonical_acc"), ("iso", "isoform_acc")]:
    p[f"{side}_mono_plddt"] = p[col].map(mo.mono_plddt)
    p[f"{side}_vlow"] = p[col].map(mo.mono_frac_vlow)
    p[f"{side}_pred"] = p[col].isin(mf.index)
    p[f"{side}_model_len"] = p[col].map(mf.chain_len)
    p[f"{side}_ipsae"] = p[col].map(mf.ipsae)
    p[f"{side}_pdockq2"] = p[col].map(mf.pdockq2)
    p[f"{side}_iptm"] = p[col].map(mf.ipTM)
p["model_len_ok"] = (p.canon_model_len == p.ref_len) & (p.iso_model_len == p.iso_len)
p["both_pred"] = p.canon_pred & p.iso_pred & p.model_len_ok
def group3(row):
    if row["class"] in ("truncation", "internal deletion"): return "subtraction"
    if row["class"] == "addition": return "addition"
    if row["class"] == "replacement": return "replacement"
    return "other"
p["grp"] = p.apply(group3, axis=1)
p = p[p.grp != "other"].copy()
p.to_csv("results/pairs_analysis.tsv", sep="\t", index=False)
print("pairs_with_hc:", len(r), "| pairs_analysis:", len(p))
