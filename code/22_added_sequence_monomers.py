import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoform_blocks import blocks

OUT = "results/numbers_new.json"
NUM = json.load(open(OUT))
r4 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)

PL = {}
with gzip.open("data/afdb_monomer_plddt.jsonl.gz", "rt") as fh:
    for line in fh:
        d = json.loads(line)
        PL[d["acc"]] = np.array(d["plddt"])

pa = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
pa = pa[pa["class"] != "unrelated"].copy()


def positions(segs):
    out = []
    for s_ in str(segs).split(";"):
        if s_ and s_ != "nan":
            rng, site = s_.split(":")
            a, b = (int(v) for v in rng.split("-"))
            out.extend((p, site) for p in range(a, b + 1))
    return out


def frac(x, lo=None, hi=None):
    x = np.asarray(x)
    if lo is not None:
        return float((x >= lo).mean())
    return float((x < hi).mean())


mono = {}
rows = []
for r in pa[pa.novel > 0].itertuples():
    p = PL.get(r.isoform_acc)
    if p is None:
        continue
    add = positions(r.novel_segments)
    mask = np.zeros(len(p), bool)
    site = np.array([""] * len(p), dtype=object)
    for pos, s_ in add:
        mask[pos - 1] = True
        site[pos - 1] = s_
    rows.append(dict(isoform_acc=r.isoform_acc, canonical_acc=r.canonical_acc, cls=r.grp,
                     n_added=int(mask.sum()), added_lt50=float((p[mask] < 50).mean()), added_ge70=float((p[mask] >= 70).mean()),
                     added_mean=float(p[mask].mean()), n_shared=int((~mask).sum()),
                     shared_lt50=float((p[~mask] < 50).mean()) if (~mask).any() else np.nan,
                     shared_ge70=float((p[~mask] >= 70).mean()) if (~mask).any() else np.nan,
                     n_added_conf=int((p[mask] >= 70).sum()),
                     added_N=int((site == "N").sum()), added_C=int((site == "C").sum()), added_I=int((site == "I").sum())))
AD = pd.DataFrame(rows)
res_add, res_sh, res_site = [], [], {"N": [], "C": [], "I": []}
for r in pa[pa.novel > 0].itertuples():
    p = PL.get(r.isoform_acc)
    if p is None:
        continue
    mask = np.zeros(len(p), bool)
    for pos, s_ in positions(r.novel_segments):
        mask[pos - 1] = True
        res_site[s_].append(p[pos - 1])
    res_add.append(p[mask])
    res_sh.append(p[~mask])
ra, rs = np.concatenate(res_add), np.concatenate(res_sh)
mono["isoforms_with_added"] = int((pa.novel > 0).sum())
mono["isoforms_with_added_monomer"] = len(AD)
mono["added_residues"] = int(len(ra))
mono["shared_residues"] = int(len(rs))
mono["added_lt50"], mono["added_ge70"] = r4(100 * frac(ra, hi=50)), r4(100 * frac(ra, lo=70))
mono["shared_lt50"], mono["shared_ge70"] = r4(100 * frac(rs, hi=50)), r4(100 * frac(rs, lo=70))
mono["added_by_site"] = {s_: {"n": len(v), "lt50": r4(100 * frac(v, hi=50)), "ge70": r4(100 * frac(v, lo=70))}
                         for s_, v in res_site.items()}
by_cls = {}
for g_, g in AD.groupby("cls"):
    by_cls[g_] = {"isoforms": len(g), "added_residues": int(g.n_added.sum()),
                  "added_lt50": r4(100 * (g.added_lt50 * g.n_added).sum() / g.n_added.sum()),
                  "added_ge70": r4(100 * (g.added_ge70 * g.n_added).sum() / g.n_added.sum()),
                  "isoforms_majority_lt50": int((g.added_lt50 > 0.5).sum())}
mono["added_by_class"] = by_cls
mono["isoforms_majority_added_lt50"] = int((AD.added_lt50 > 0.5).sum())
AD.to_csv("results/added_residues_monomer.tsv.gz", sep="\t", index=False)

rem, kept = [], []
for r in pa[pa.deleted > 0].itertuples():
    p = PL.get(r.canonical_acc)
    if p is None:
        continue
    mask = np.zeros(len(p), bool)
    for pos, _ in positions(r.deleted_segments):
        mask[pos - 1] = True
    rem.append(p[mask])
    kept.append(p[~mask])
rem, kept = np.concatenate(rem), np.concatenate(kept)
edges = np.arange(0, 102, 2)
pd.DataFrame({"bin_start": edges[:-1], "added": np.histogram(ra, edges)[0], "shared": np.histogram(rs, edges)[0],
              "removed": np.histogram(rem, edges)[0], "kept": np.histogram(kept, edges)[0]}).to_csv(
    "results/table_monomer_plddt_hist.tsv", sep="\t", index=False)
mono["removed_residues"] = int(len(rem))
mono["removed_lt50"], mono["kept_lt50"] = r4(100 * frac(rem, hi=50)), r4(100 * frac(kept, hi=50))
mono["removed_ge70"], mono["kept_ge70"] = r4(100 * frac(rem, lo=70)), r4(100 * frac(kept, lo=70))

SW = pd.read_csv("results/swap_segments.tsv", sep="\t")
mrem, madd = [], []
for r in SW.itertuples():
    pc, pi = PL.get(r.canonical_acc), PL.get(r.isoform_acc)
    mrem.append(float(pc[r.can_start - 1:r.can_end].mean()) if pc is not None else np.nan)
    madd.append(float(pi[r.iso_start - 1:r.iso_end].mean()) if pi is not None else np.nan)
SW["removed_mean_plddt"], SW["added_mean_plddt"] = mrem, madd
SW.to_csv("results/swap_segments.tsv", sep="\t", index=False)
sw = {}
for lab, g in [("homologous", SW[SW.homologous]), ("other", SW[~SW.homologous])]:
    g = g.dropna(subset=["removed_mean_plddt", "added_mean_plddt"])
    sw[lab] = {"n": len(g), "removed_median": r4(g.removed_mean_plddt.median()), "added_median": r4(g.added_mean_plddt.median()),
               "added_ge70": int((g.added_mean_plddt >= 70).sum()), "added_lt50": int((g.added_mean_plddt < 50).sum())}
mono["swaps"] = sw

B = pa[pa.both_pred == True].copy()
nc = B[(B.canon_hc == False)].copy()
nc["n_added_conf"] = nc.isoform_acc.map(AD.set_index("isoform_acc").n_added_conf).fillna(0).astype(int)
nc["bin"] = np.select([nc.novel == 0, nc.n_added_conf == 0, nc.n_added_conf < 10], ["none added", "0 confident", "1-9 confident"], "10+ confident")
gb = {}
for b_ in ["none added", "0 confident", "1-9 confident", "10+ confident"]:
    g = nc[nc.bin == b_]
    gb[b_] = {"n": len(g), "gains": int((g.iso_hc == True).sum()), "rate": r4(100 * (g.iso_hc == True).mean()) if len(g) else None}
a, b = gb["10+ confident"], gb["0 confident"]
OR, p = fisher_exact([[a["gains"], a["n"] - a["gains"]], [b["gains"], b["n"] - b["gains"]]])
gb["10plus_vs_0"] = {"OR": r4(OR), "p": float(f"{p:.4g}")}
mono["gains_by_confident_added"] = gb

seqt = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
SEQ = dict(zip(seqt.acc, seqt.seq))
sh_rows = []
for r in B.itertuples():
    pc, pi = PL.get(r.canonical_acc), PL.get(r.isoform_acc)
    if pc is None or pi is None:
        continue
    ci, ii = [], []
    for a, b, s_ in blocks(SEQ[r.canonical_acc], SEQ[r.isoform_acc]):
        ci.extend(range(a, a + s_))
        ii.extend(range(b, b + s_))
    ci, ii = np.array(ci), np.array(ii)
    sh_rows.append(dict(canonical_acc=r.canonical_acc, isoform_acc=r.isoform_acc, grp=r.grp, canon_hc=r.canon_hc, iso_hc=r.iso_hc,
                        n_shared=len(ci), d_shared=float((pi[ii] - pc[ci]).mean())))
SH = pd.DataFrame(sh_rows)
SH.to_csv("results/shared_fold_confidence.tsv.gz", sep="\t", index=False)
fold = {"pairs": len(SH)}
ret = pd.read_csv("results/canonical_retention.tsv", sep="\t").set_index("isoform_acc").iface_retained
hcS = SH[SH.canon_hc == True].copy()
hcS["ret"] = hcS.isoform_acc.map(ret)
intact = hcS[hcS.ret >= 1.0]
lo, ke = intact[intact.iso_hc == False], intact[intact.iso_hc == True]
fold["intact"] = {"n": len(intact), "losses": len(lo), "d_median_loss": r4(lo.d_shared.median()), "d_median_kept": r4(ke.d_shared.median()),
                  "p_mwu": float(f"{mannwhitneyu(lo.d_shared, ke.d_shared).pvalue:.4g}")}
intact = intact.assign(dbin=np.select([intact.d_shared >= 0, intact.d_shared >= -2, intact.d_shared >= -5], ["0+", "-2to0", "-5to-2"], "lt-5"))
fold["intact_loss_by_d"] = {b_: {"n": int((intact.dbin == b_).sum()), "losses": int((intact[intact.dbin == b_].iso_hc == False).sum()),
                                 "rate": r4(100 * (intact[intact.dbin == b_].iso_hc == False).mean()) if (intact.dbin == b_).any() else None}
                            for b_ in ["0+", "-2to0", "-5to-2", "lt-5"]}
mono["shared_fold"] = fold

NUM["monomers"] = mono
json.dump(NUM, open(OUT, "w"), indent=1)
print(json.dumps(mono, indent=1)[:6000])
