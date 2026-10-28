import os, sys, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from isoform_blocks import classify
up = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
seq = dict(zip(up.acc, up.seq))
gene = dict(zip(up.acc, up.gene))
events = dict(zip(up.base, up.events))
iso = up[up.is_isoform_record]
rows = []
t0 = time.time()
for r in iso.itertuples():
    ref = seq.get(r.base)
    if ref is None:
        continue
    in_range = (15 <= len(ref) <= 1500) and (15 <= len(r.seq) <= 1500)
    if in_range:
        k = classify(ref, r.seq)
    else:
        k = {"ref_len": len(ref), "iso_len": len(r.seq)}
    k.update({"canonical_acc": r.base, "isoform_acc": r.acc, "gene": gene.get(r.base, ""),
              "uniprot_events": events.get(r.base, ""), "in_afdb_length_range": in_range})
    rows.append(k)
out = pd.DataFrame(rows)
cols = ["canonical_acc", "isoform_acc", "gene", "uniprot_events", "class", "subclass", "length_relation",
        "ref_len", "iso_len", "shared", "novel", "deleted", "novel_N", "novel_C", "novel_I",
        "deleted_N", "deleted_C", "deleted_I", "n_blocks", "longest_block", "novel_segments",
        "deleted_segments", "in_afdb_length_range"]
out = out.reindex(columns=cols)
out.to_csv("results/isoform_pairs_classified.tsv", sep="\t", index=False)
inr = out[out.in_afdb_length_range == True]
open("data/accessions_monomer.txt", "w").write("\n".join(sorted(set(inr.canonical_acc) | set(inr.isoform_acc))) + "\n")
print(f"classified {len(out)} isoform records in {time.time()-t0:.0f}s")
