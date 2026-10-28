import json

import pandas as pd

M = pd.read_csv("data/manifest_extract/manifest_human.tsv.gz", sep="\t").rename(columns={"uniprotAccession": "acc"})
if M.acc.duplicated().any():
    raise SystemExit("an accession appears more than once in the human manifest rows")
M["ipsae"] = M[["ipSAE_AB", "ipSAE_BA"]].max(axis=1)
M["ipsae_min"] = M[["ipSAE_AB", "ipSAE_BA"]].min(axis=1)
M["pdockq2"] = M[["pDockQ2_AB", "pDockQ2_BA"]].max(axis=1)
M["lis"] = M[["LIS_AB", "LIS_BA"]].max(axis=1)
M["chain_len"] = (M.n0chn // 2).astype(int)
A = M[["acc", "modelEntityId", "chain_len", "ipTM", "ipsae", "ipsae_min", "pdockq2", "pDockQ", "lis",
       "N_clash_backbone", "N_clash_heavyAtom"]].copy()
A["rule_hc"] = (A.ipsae >= 0.6) & (A.pdockq2 >= 0.23)
A.to_csv("results/manifest_per_accession.tsv", sep="\t", index=False)

cx = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t")
q = cx.merge(A, on="acc", how="left", suffixes=("_api", ""))
pred = q.ipsae.notna()
served = q.hc_homodimer
chk = {
    "manifest_human_rows": int(len(A)),
    "queried_accessions": int(len(q)),
    "queried_predicted": int(pred.sum()),
    "queried_not_predicted": int((~pred).sum()),
    "served": int(served.sum()),
    "served_not_in_manifest": int((served & ~pred).sum()),
    "rule_vs_served_mismatches": int((q.loc[pred, "rule_hc"] != q.loc[pred, "hc_homodimer"]).sum()),
    "served_api_ipsae_within_001_of_larger_direction": int(((q.ipSAE - q.ipsae).abs() <= 0.0101)[served].sum()),
    "served_api_pdockq2_within_001_of_larger_direction": int(((q.pDockQ2 - q.pdockq2).abs() <= 0.0101)[served].sum()),
    "served_smaller_direction_ipsae_below_06": int((served & (q.ipsae_min < 0.6)).sum()),
    "served_whole_model_plddt_below_70": int((served & (q.pLDDT < 70)).sum()),
    "served_backbone_clashes_above_10": int((served & (q.N_clash_backbone > 10)).sum()),
    "predicted_ipsae_ge06": int((pred & (q.ipsae >= 0.6)).sum()),
    "predicted_ipsae_ge06_pdockq2_lt023": int((pred & (q.ipsae >= 0.6) & (q.pdockq2 < 0.23)).sum()),
    "predicted_ipsae_lt06_pdockq2_ge023": int((pred & (q.ipsae < 0.6) & (q.pdockq2 >= 0.23)).sum()),
}
json.dump(chk, open("results/manifest_rule_check.json", "w"), indent=1)
print(json.dumps(chk, indent=1))
if chk["served_not_in_manifest"] or chk["rule_vs_served_mismatches"]:
    raise SystemExit("deposition is not reproduced by ipSAE >= 0.6 and pDockQ2 >= 0.23; stop and inspect")
