import os
import sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoform_blocks import blocks
m = pd.read_csv("results/pair_mechanism.tsv", sep="\t")
up = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t"); useq = dict(zip(up.acc, up.seq))
rows = []
for r in m[m.novel_model > 0].itertuples():
    d = pd.read_csv(f"results/struct/{r.isoform_acc}.tsv.gz", sep="\t", keep_default_na=False, dtype={"partners_A": str, "rsa_A": str})
    d["rsa_A"] = pd.to_numeric(d.rsa_A, errors="coerce")
    seq = "".join(d.aa); bl = blocks(useq[r.canonical_acc], seq)
    sh = np.zeros(len(seq), bool)
    for _, y, s in bl: sh[y:y + s] = True
    d["novel"] = ~sh; d["iface"] = (d.iface_A == 1) | (d.iface_B == 1)
    d["plddt"] = np.minimum(d.plddt_A, d.plddt_B)
    d["acc"] = r.isoform_acc; d["gain"] = not r.canon_hc
    rows.append(d[["acc", "pos", "aa", "novel", "iface", "plddt", "rsa_A", "gain"]])
R = pd.concat(rows, ignore_index=True)
R.to_csv("results/residues_novel_models.tsv.gz", sep="\t", index=False)
print("residues:", len(R), "models:", R.acc.nunique())
