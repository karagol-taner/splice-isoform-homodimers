import pandas as pd
p = pd.read_csv("results/pairs_with_hc.tsv", sep="\t")
accs = sorted(set(p.loc[p.iso_hc & (p.novel > 0), "isoform_acc"]))
open("data/sasa_needed.txt", "w").write("\n".join(accs) + "\n")
print("isoform models needing SASA:", len(accs))
