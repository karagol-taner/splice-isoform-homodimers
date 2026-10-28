import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoform_blocks import blocks

PLDDT_OK = 70.0
pairs = pd.read_csv("results/pairs_with_hc.tsv", sep="\t").rename(columns={"class": "cls"})
up = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
useq = dict(zip(up.acc, up.seq))

def load(acc):
    fn = f"results/struct/{acc}.tsv.gz"
    if not os.path.exists(fn):
        return None
    d = pd.read_csv(fn, sep="\t", keep_default_na=False, dtype={"partners_A": str, "rsa_A": str})
    d["rsa_A"] = pd.to_numeric(d.rsa_A, errors="coerce")
    d["partners"] = d.partners_A.apply(lambda s: [int(x) for x in str(s).split(";") if x not in ("", "nan")])
    return d

def iface_positions(d):
    return set(d.loc[(d.iface_A == 1) | (d.iface_B == 1), "pos"])

rows = []
for p in pairs[pairs.iso_hc].itertuples():
    iso = load(p.isoform_acc)
    if iso is None:
        continue
    iso_seq = "".join(iso.aa)
    can_seq = useq[p.canonical_acc]
    bl = blocks(can_seq, iso_seq)
    shared_iso = np.zeros(len(iso_seq), bool); iso2can = {}
    for a, b, s in bl:
        shared_iso[b:b + s] = True
        for k in range(s): iso2can[b + k + 1] = a + k + 1
    novel = ~shared_iso
    pos = iso.pos.values
    ifc = ((iso.iface_A == 1) | (iso.iface_B == 1)).values
    plddt = np.minimum(iso.plddt_A.values, iso.plddt_B.values)
    conf = plddt >= PLDDT_OK
    surf = (iso.rsa_A >= 0.25).values
    novel_pos = set(pos[novel])
    nn = 0
    for k in np.where(novel & (iso.iface_A == 1).values)[0]:
        nn += sum(1 for j in iso.partners.iloc[k] if j in novel_pos)
    rec = {"canonical_acc": p.canonical_acc, "isoform_acc": p.isoform_acc, "gene": p.gene}
    rec.update({
        "class": p.cls, "subclass": p.subclass,
        "canon_hc": bool(p.canon_hc), "iso_len_model": len(iso_seq), "iso_len_uniprot": p.iso_len,
        "novel_model": int(novel.sum()), "novel_conf": int((novel & conf).sum()),
        "n_iface": int(ifc.sum()), "n_iface_novel": int((ifc & novel).sum()),
        "n_iface_novel_conf": int((ifc & novel & conf).sum()),
        "n_surface": int(surf.sum()), "n_surface_novel": int((surf & novel).sum()),
        "novel_novel_contacts": int(nn),
        "mean_plddt_novel": float(plddt[novel].mean()) if novel.any() else np.nan,
        "mean_plddt_iface_novel": float(plddt[ifc & novel].mean()) if (ifc & novel).any() else np.nan,
    })
    if p.canon_hc:
        can = load(p.canonical_acc)
        rec["canon_model_matches_uniprot"] = bool(can is not None and "".join(can.aa) == can_seq)
        if can is not None and rec["canon_model_matches_uniprot"]:
            can_if = iface_positions(can)
            can2iso = {c: i for i, c in iso2can.items()}
            can_if_present = {can2iso[c] for c in can_if if c in can2iso}
            iso_if_shared = {q for q in pos[ifc] if q in iso2can}
            inter = len(can_if_present & iso_if_shared); union = len(can_if_present | iso_if_shared)
            rec.update({"canon_n_iface": len(can_if), "canon_iface_present_in_iso": len(can_if_present),
                        "iface_shared_overlap": inter, "iface_jaccard_shared": inter / union if union else np.nan})
    rows.append(rec)
out = pd.DataFrame(rows)
out.to_csv("results/pair_mechanism.tsv", sep="\t", index=False)
print("pairs with isoform model features:", len(out))
