import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoform_blocks import blocks

PLDDT_OK = 70.0
p = pd.read_csv("results/pairs_analysis.tsv", sep="\t").rename(columns={"class": "cls"})
up = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
useq = dict(zip(up.acc, up.seq))
cache = {}


def canonical_model(acc):
    if acc not in cache:
        d = pd.read_csv(f"results/struct/{acc}.tsv.gz", sep="\t", keep_default_na=False,
                        usecols=["pos", "aa", "plddt_A", "plddt_B", "iface_A", "iface_B"])
        cache[acc] = d
    return cache[acc]


rows, skipped = [], {"model_sequence_differs": 0, "no_model": 0}
for r in p[p.canon_hc].itertuples():
    try:
        d = canonical_model(r.canonical_acc)
    except FileNotFoundError:
        skipped["no_model"] += 1
        continue
    can_seq = useq[r.canonical_acc]
    if "".join(d.aa) != can_seq:
        skipped["model_sequence_differs"] += 1
        continue
    iso_seq = useq[r.isoform_acc]
    bl = blocks(can_seq, iso_seq)
    kept = np.zeros(len(can_seq), bool)
    for a, b, s in bl:
        kept[a:a + s] = True
    iface = ((d.iface_A == 1) | (d.iface_B == 1)).values
    plddt = np.minimum(d.plddt_A.values, d.plddt_B.values)
    conf = plddt >= PLDDT_OK
    sites, structured = 0, 0
    prev_a_end, prev_b_end = None, 0
    edges = [(a, b, s) for a, b, s in bl]
    for a, b, s in edges:
        if b > prev_b_end:
            sites += 1
            left = plddt[prev_a_end - 1] if prev_a_end is not None and prev_a_end > 0 else None
            right = plddt[a] if a < len(can_seq) else None
            fl = [x for x in (left, right) if x is not None]
            if fl and all(x >= PLDDT_OK for x in fl):
                structured += 1
        prev_a_end, prev_b_end = a + s, b + s
    if edges and len(iso_seq) > prev_b_end:
        sites += 1
        left = plddt[prev_a_end - 1]
        if left >= PLDDT_OK:
            structured += 1
    rows.append({
        "canonical_acc": r.canonical_acc, "isoform_acc": r.isoform_acc, "gene": r.gene, "grp": r.grp,
        "class": r.cls, "subclass": r.subclass,
        "iso_hc": bool(r.iso_hc), "loss": not bool(r.iso_hc),
        "can_len": len(can_seq), "iso_len": len(iso_seq), "novel": int(r.novel), "deleted": int(r.deleted),
        "can_iface": int(iface.sum()),
        "iface_retained": float(kept[iface].mean()) if iface.any() else np.nan,
        "can_conf": int(conf.sum()),
        "conf_removed": float((~kept[conf]).mean()) if conf.any() else np.nan,
        "removed_nonconf": int((~kept & ~conf).sum()),
        "ins_sites": sites, "ins_in_structured": structured,
    })

out = pd.DataFrame(rows)
out.to_csv("results/canonical_retention.tsv", sep="\t", index=False)
print("canonical high-confidence pairs:", int(p.canon_hc.sum()), "| analysed:", len(out), "| skipped:", skipped)
