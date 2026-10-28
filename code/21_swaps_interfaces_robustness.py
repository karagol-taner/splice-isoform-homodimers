import json
import os
import random
import sys

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from isoform_blocks import blocks

OUT = "results/numbers_new.json"
NUM = json.load(open(OUT)) if os.path.exists(OUT) else {}
r4 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)


def wilson(k, n):
    if n == 0:
        return [None, None]
    from statsmodels.stats.proportion import proportion_confint
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return [r4(100 * lo), r4(100 * hi)]


seq = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
SEQ = dict(zip(seq.acc, seq.seq))
pa = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
pa = pa[pa["class"] != "unrelated"].copy()
B = pa[pa.both_pred == True].copy()

al = Align.PairwiseAligner()
al.mode = "local"
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -11, -1
rng = random.Random(11)


_ALPHA = set(al.substitution_matrix.alphabet)


def clean(x):
    return "".join(c if c in _ALPHA else "X" for c in x)


def align_stats(a, b):
    a, b = clean(a), clean(b)
    if al.score(a, b) <= 0:
        return 0.0, 0, 0
    aln = al.align(a, b)[0]
    cols = ident = 0
    for (sa, ea), (sb, _) in zip(*aln.aligned):
        n = ea - sa
        cols += n
        ident += sum(1 for i in range(n) if a[sa + i] == b[sb + i])
    return aln.score, cols, ident


rows = []
for r in pa[pa["class"] == "replacement"].itertuples():
    can, iso = SEQ[r.canonical_acc], SEQ[r.isoform_acc]
    bl = blocks(can, iso)
    pa_, pb_ = 0, 0
    gaps = []
    for k, (a, b, s) in enumerate(bl):
        gaps.append(("N" if k == 0 else "I", pa_, a, pb_, b))
        pa_, pb_ = a + s, b + s
    gaps.append(("C", pa_, len(can), pb_, len(iso)))
    for site, ca, cb, ia, ib in gaps:
        if cb > ca and ib > ia:
            rem, add = can[ca:cb], iso[ia:ib]
            score, cols, ident = align_stats(add, rem)
            shuf = []
            for _ in range(100):
                s_ = list(clean(add))
                rng.shuffle(s_)
                shuf.append(al.score("".join(s_), clean(rem)))
            mu, sd = float(np.mean(shuf)), float(np.std(shuf))
            rows.append(dict(canonical_acc=r.canonical_acc, isoform_acc=r.isoform_acc, site=site,
                             can_start=ca + 1, can_end=cb, iso_start=ia + 1, iso_end=ib,
                             len_removed=cb - ca, len_added=ib - ia, score=score, aligned=cols,
                             identity=ident / cols if cols else 0.0, z=(score - mu) / sd if sd > 0 else 0.0))
SW = pd.DataFrame(rows)
SW["homologous"] = (SW.z >= 5) & (SW.aligned >= 20) & (SW.identity >= 0.25)
SW.to_csv("results/swap_segments.tsv", sep="\t", index=False)

SW["swap_len"] = SW[["len_removed", "len_added"]].min(axis=1)
best = SW.sort_values("swap_len", ascending=False).drop_duplicates("isoform_acc").set_index("isoform_acc")
pa["swap_homologous"] = pa.isoform_acc.map(best.homologous)
B["swap_homologous"] = B.isoform_acc.map(best.homologous)
rep = B[B["class"] == "replacement"]
hc = rep[rep.canon_hc == True]
sw = {"swaps": len(SW), "replacements_with_swap": int(SW.isoform_acc.nunique()),
      "replacements": int((pa["class"] == "replacement").sum()),
      "homologous_swaps": int(SW.homologous.sum()),
      "replacements_homologous": int(pa.swap_homologous.fillna(False).sum()),
      "z_threshold": 5, "min_aligned": 20, "min_identity": 0.25}
for lab, g in [("homologous", hc[hc.swap_homologous == True]), ("other", hc[hc.swap_homologous != True])]:
    k = int((g.iso_hc == False).sum())
    sw[f"loss_{lab}"] = {"n": len(g), "losses": k, "rate": r4(100 * k / len(g)) if len(g) else None, "ci": wilson(k, len(g))}
NUM["swaps"] = sw

ret = pd.read_csv("results/canonical_retention.tsv", sep="\t")
hcB = B[(B.canon_hc == True) & (B.deleted > 0)]
STRUCT = "results/struct"
cache = {}


def canon_features(acc):
    if acc not in cache:
        f = pd.read_csv(f"{STRUCT}/{acc}.tsv.gz", sep="\t")
        f["iface"] = ((f.iface_A > 0) | (f.iface_B > 0)).astype(int)
        f["plddt"] = f[["plddt_A", "plddt_B"]].min(axis=1)
        cache[acc] = f[["pos", "iface", "plddt"]]
    return cache[acc]


def removed_positions(segs):
    out = set()
    for s_ in str(segs).split(";"):
        if s_ and s_ != "nan":
            a, b = s_.split(":")[0].split("-")
            out.update(range(int(a), int(b) + 1))
    return out


res_rows = []
for r in hcB.itertuples():
    f = canon_features(r.canonical_acc)
    if len(f) != r.ref_len:
        continue
    rm = removed_positions(r.deleted_segments)
    f = f.assign(removed=f.pos.isin(rm).astype(int), canonical_acc=r.canonical_acc, isoform_acc=r.isoform_acc,
                 grp=r.grp)
    res_rows.append(f)
RR = pd.concat(res_rows, ignore_index=True)
RR["band"] = pd.cut(RR.plddt, [-1, 50, 70, 90, 101], right=False, labels=["<50", "50-70", "70-90", ">=90"])


def mh_or(df):
    num = den = 0.0
    for _, g in df.groupby("band", observed=True):
        a = ((g.iface == 1) & (g.removed == 1)).sum(); b = ((g.iface == 1) & (g.removed == 0)).sum()
        c = ((g.iface == 0) & (g.removed == 1)).sum(); d = ((g.iface == 0) & (g.removed == 0)).sum()
        n = a + b + c + d
        if n:
            num += a * d / n
            den += b * c / n
    return num / den if den else float("nan")


ir = {"pairs": int(RR.isoform_acc.nunique()), "canonicals": int(RR.canonical_acc.nunique()), "residues": len(RR)}
by = RR.groupby(["band", "iface"], observed=True).removed.agg(["sum", "size"]).reset_index()
ir["by_band"] = {f"{b}|{'iface' if i else 'other'}": {"removed": int(s), "n": int(n), "pct": r4(100 * s / n)}
                 for b, i, s, n in by.itertuples(index=False)}
ir["mh_or"] = r4(mh_or(RR))
canons = RR.canonical_acc.unique()
groups = {c: g for c, g in RR.groupby("canonical_acc")}
brng = np.random.default_rng(7)
boots = []
for _ in range(1000):
    pick = brng.choice(canons, size=len(canons), replace=True)
    boots.append(mh_or(pd.concat([groups[c] for c in pick], ignore_index=True)))
ir["mh_or_ci"] = [r4(np.nanpercentile(boots, 2.5)), r4(np.nanpercentile(boots, 97.5))]
NUM["interface_removal"] = ir
RR[["canonical_acc", "isoform_acc", "grp", "pos", "iface", "plddt", "band", "removed"]].to_csv(
    "results/interface_removal_residues.tsv.gz", sep="\t", index=False)

man = pd.read_csv("results/manifest_per_accession.tsv", sep="\t").set_index("acc")
calls = {
    "deposition": lambda s, p, t: (s >= 0.6) & (p >= 0.23),
    "ipsae_only": lambda s, p, t: s >= 0.6,
    "pdockq2_only": lambda s, p, t: p >= 0.23,
    "iptm_ge08": lambda s, p, t: t >= 0.8,
    "iptm_ge06": lambda s, p, t: t >= 0.6,
}
for si in [0.5, 0.55, 0.65, 0.7]:
    for pq in [0.15, 0.23, 0.3]:
        calls[f"grid_ipsae{si}_pdockq2{pq}"] = (lambda si_, pq_: (lambda s, p, t: (s >= si_) & (p >= pq_)))(si, pq)
for pq in [0.15, 0.3]:
    calls[f"grid_ipsae0.6_pdockq2{pq}"] = (lambda pq_: (lambda s, p, t: (s >= 0.6) & (p >= pq_)))(pq)

retmap = ret.set_index("isoform_acc").iface_retained
rob = {}
for name, f in calls.items():
    ch = f(B.canon_ipsae, B.canon_pdockq2, B.canon_iptm)
    ih = f(B.iso_ipsae, B.iso_pdockq2, B.iso_iptm)
    d = {}
    for g_ in ["addition", "subtraction", "replacement", "all"]:
        sel = (B.grp == g_) if g_ != "all" else pd.Series(True, index=B.index)
        n_c, n_n = int((sel & ch).sum()), int((sel & ~ch).sum())
        lo, ga = int((sel & ch & ~ih).sum()), int((sel & ~ch & ih).sum())
        d[g_] = {"canon_hc": n_c, "losses": lo, "loss_rate": r4(100 * lo / n_c) if n_c else None,
                 "canon_not": n_n, "gains": ga, "gain_rate": r4(100 * ga / n_n) if n_n else None}
    rr = B[ch & B.isoform_acc.isin(retmap.index)].copy()
    rr["ret"] = rr.isoform_acc.map(retmap)
    rr["bin"] = np.select([rr.ret >= 1.0, rr.ret >= 0.9, rr.ret >= 0.5], ["all", "90to100", "50to90"], "lt50")
    rr["loss"] = ~f(rr.iso_ipsae, rr.iso_pdockq2, rr.iso_iptm)
    d["loss_by_retention"] = {b_: {"n": int((rr.bin == b_).sum()), "losses": int(rr[rr.bin == b_].loss.sum()),
                                   "rate": r4(100 * rr[rr.bin == b_].loss.mean()) if (rr.bin == b_).any() else None}
                              for b_ in ["all", "90to100", "50to90", "lt50"]}
    rob[name] = d
accs = set(open("data/accessions_to_query.txt").read().split())
M = man[man.index.isin(accs)].copy()
M["band"] = pd.cut(M.ipTM, [-1, 0.5, 0.6, 0.8, 1.01], right=False, labels=["<0.5", "0.5-0.6", "0.6-0.8", ">=0.8"])
ag = M.groupby("band", observed=True).rule_hc.agg(["sum", "size"])
rob["iptm_vs_call"] = {"n": len(M), "deposited": int(M.rule_hc.sum()),
                       "by_band": {b_: {"deposited": int(s), "n": int(n), "pct": r4(100 * s / n)} for b_, (s, n) in ag.iterrows()},
                       "iptm_ge08_not_deposited": int(((M.ipTM >= 0.8) & ~M.rule_hc.astype(bool)).sum()),
                       "deposited_iptm_lt06": int(((M.ipTM < 0.6) & M.rule_hc.astype(bool)).sum())}
NUM["robustness"] = rob

def origin(ev):
    ev = str(ev)
    if "initiation" in ev and "splicing" not in ev and "promoter" not in ev:
        return "alternative initiation only"
    if "promoter" in ev and "splicing" not in ev:
        return "alternative promoter without splicing"
    if "splicing" in ev and ("promoter" in ev or "initiation" in ev):
        return "splicing with promoter or initiation"
    if "splicing" in ev:
        return "alternative splicing only"
    return "other"


pa["origin"] = pa.uniprot_events.map(origin)
NUM["origin"] = {"classes_by_origin": {f"{o}|{c}": int(n) for (o, c), n in pa.groupby(["origin", "class"]).size().items()}}

json.dump(NUM, open(OUT, "w"), indent=1)
print(json.dumps({k: NUM[k] for k in ["swaps"]}, indent=1)[:3000])
print(json.dumps(NUM["interface_removal"], indent=1)[:2500])
