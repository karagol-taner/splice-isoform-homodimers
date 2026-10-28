import numpy as np, pandas as pd
p = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
rng = np.random.default_rng(7)
can = p.drop_duplicates("canonical_acc")
hc = can[can.canon_hc].canonical_acc.tolist(); non = can[~can.canon_hc].canonical_acc.tolist()
picked = list(rng.choice(hc, 150, replace=False)) + list(rng.choice(non, 150, replace=False))
open("data/copy_sample_accs.txt", "w").write("\n".join(picked) + "\n")
print(len(picked))
