import gzip
import importlib.metadata as md
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

OUT_RES = "results/aiupred_isoforms.jsonl.gz"
OUT = "results/added_sequence_disorder.json"
r4 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)

pa = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
pa = pa[pa["class"] != "unrelated"].copy()
seqt = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
SEQ = dict(zip(seqt.acc, seqt.seq))
R = pd.read_csv("results/residues_novel_models.tsv.gz", sep="\t")

targets = sorted(set(pa[pa.novel > 0].isoform_acc) | set(R.acc), key=lambda a: (len(SEQ[a]), a))

DIS = {}
if os.path.exists(OUT_RES):
    with gzip.open(OUT_RES, "rt") as fh:
        for line in fh:
            d = json.loads(line)
            DIS[d["acc"]] = np.array(d["disorder"])
todo = [a for a in targets if a not in DIS]
if todo:
    import torch
    from aiupred import AIUPred
    torch.set_num_threads(os.cpu_count() or 1)
    pred = AIUPred(force_cpu=True)
    t0 = time.time()
    with gzip.open(OUT_RES, "at") as fh:
        for k, a in enumerate(todo, 1):
            with torch.no_grad():
                d = np.round(pred.predict_disorder(SEQ[a]), 3)
            if len(d) != len(SEQ[a]):
                sys.exit(f"STOP: {a}: {len(d)} scores for {len(SEQ[a])} residues")
            fh.write(json.dumps({"acc": a, "disorder": d.tolist()}) + "\n")
            DIS[a] = d
            if k % 500 == 0:
                fh.flush()
                print(f"{k}/{len(todo)} {time.time() - t0:.0f}s", flush=True)
missing = [a for a in targets if a not in DIS]
if missing:
    sys.exit(f"STOP: no disorder for {len(missing)} isoforms")

PL = {}
with gzip.open("data/afdb_monomer_plddt.jsonl.gz", "rt") as fh:
    for line in fh:
        d = json.loads(line)
        PL[d["acc"]] = np.array(d["plddt"])


def positions(segs):
    out = []
    for s_ in str(segs).split(";"):
        if s_ and s_ != "nan":
            rng, site = s_.split(":")
            a, b = (int(v) for v in rng.split("-"))
            out.extend((p, site) for p in range(a, b + 1))
    return out


prev = json.load(open(OUT)) if os.path.exists(OUT) else {}


def version(pkg):
    try:
        return md.version(pkg).split("+")[0]
    except md.PackageNotFoundError:
        return prev.get(f"{pkg}_version")


res = {"predictor": "AIUPred", "aiupred_version": version("aiupred"), "torch_version": version("torch"), "cutoff": 0.5}

cols = {k: [] for k in ("dis", "plddt", "added", "site", "cls")}
n_iso, iso_frac = 0, []
for r in pa[pa.novel > 0].itertuples():
    p = PL.get(r.isoform_acc)
    if p is None:
        continue
    d = DIS[r.isoform_acc]
    if len(p) != len(d):
        sys.exit(f"STOP: {r.isoform_acc}: monomer model and sequence differ in length")
    n_iso += 1
    mask = np.zeros(len(d), bool)
    site = np.array([""] * len(d), dtype=object)
    for pos, s_ in positions(r.novel_segments):
        mask[pos - 1] = True
        site[pos - 1] = s_
    cols["dis"].append(d)
    cols["plddt"].append(p)
    cols["added"].append(mask)
    cols["site"].append(site)
    cols["cls"].append(np.full(len(d), r.grp, dtype=object))
    iso_frac.append(float((d[mask] >= 0.5).mean()))
D = pd.DataFrame({k: np.concatenate(v) for k, v in cols.items()})
D["disordered"] = D.dis >= 0.5
res["isoforms"] = n_iso
for lab, g in [("added", D[D.added]), ("shared", D[~D.added])]:
    lo, hi = g[g.plddt < 50], g[g.plddt >= 70]
    res[lab] = {"residues": len(g), "disordered": r4(100 * g.disordered.mean()),
                "plddt_lt50": {"residues": len(lo), "disordered": r4(100 * lo.disordered.mean())},
                "plddt_ge70": {"residues": len(hi), "disordered": r4(100 * hi.disordered.mean())}}
res["added_by_site"] = {s_: {"residues": int((D.added & (D.site == s_)).sum()),
                             "disordered": r4(100 * D[D.added & (D.site == s_)].disordered.mean())} for s_ in ["N", "C", "I"]}
res["added_by_class"] = {c_: {"residues": int((D.added & (D.cls == c_)).sum()),
                              "disordered": r4(100 * D[D.added & (D.cls == c_)].disordered.mean())}
                         for c_ in sorted(D[D.added].cls.unique())}
res["isoforms_majority_added_disordered"] = int((np.array(iso_frac) > 0.5).sum())

match = {a_ for a_, g in R.groupby("acc") if "".join(g.sort_values("pos").aa) == SEQ[a_]}
res["interface_models_sequence_matched"] = {"models": int(R.acc.nunique()), "matched": len(match)}
R = R[R.acc.isin(match)].copy()
R["dis"] = [DIS[a_][p_ - 1] for a_, p_ in zip(R.acc, R.pos)]
S = R[R.rsa_A >= 0.25].copy()
S["iface_i"] = S.iface.astype(int)
S["novel_i"] = S.novel.astype(int)
S["pbin"] = pd.cut(S.plddt, [0, 50, 70, 90, 101], right=False, labels=["lt50", "50to70", "70to90", "ge90"])
S["dbin"] = pd.cut(S.dis, [0, 0.25, 0.5, 0.75, 1.001], right=False, labels=["lt25", "25to50", "50to75", "ge75"])
res["surface_added_disordered"] = r4(100 * (S[S.novel].dis >= 0.5).mean())
res["surface_shared_disordered"] = r4(100 * (S[~S.novel].dis >= 0.5).mean())
models = {}
for name, f in [("unadjusted", "iface_i ~ novel_i"), ("disorder_adjusted", "iface_i ~ novel_i + C(dbin)"),
                ("plddt_adjusted", "iface_i ~ novel_i + C(pbin)"), ("plddt_and_disorder", "iface_i ~ novel_i + C(pbin) + C(dbin)")]:
    m = smf.logit(f, data=S).fit(disp=0)
    ci = m.conf_int().loc["novel_i"]
    models[name] = {"OR": r4(np.exp(m.params["novel_i"])), "lo": r4(np.exp(ci.iloc[0])), "hi": r4(np.exp(ci.iloc[1])),
                    "p": float(f"{m.pvalues['novel_i']:.3g}"), "n_residues": int(m.nobs)}
res["interface_models"] = models
band = (S.groupby(["dbin", "novel"], observed=True).agg(n=("iface", "size"), iface=("iface", "sum")).reset_index())
band["iface_rate"] = (100 * band.iface / band.n).round(4)
res["interface_by_disorder"] = band.to_dict(orient="records")

json.dump(res, open(OUT, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps({k: res[k] for k in ("isoforms", "added", "shared", "added_by_site", "added_by_class", "interface_models")}, indent=1))
