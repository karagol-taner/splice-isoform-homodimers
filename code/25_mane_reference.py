import json
import os
import re
import sys
import urllib.request

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

r4 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)
MANE = "data/uniprot_mane_select.tsv"
if not os.path.exists(MANE):
    url = ("https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606"
           "&fields=accession,xref_mane-select&format=tsv")
    resp = urllib.request.urlopen(url, timeout=600)
    rel = resp.headers.get("X-UniProt-Release")
    if rel != "2026_03":
        sys.exit(f"STOP: UniProt release {rel}, the study used 2026_03")
    open(MANE, "w").write(resp.read().decode())

m = pd.read_csv(MANE, sep="\t").fillna("")


def tags(s):
    out = []
    for part in [x.strip() for x in s.split(";") if x.strip()]:
        t = re.findall(r"\[([A-Z0-9-]+)\]", part)
        out.append(t[0] if t else "displayed")
    return out


M = dict(zip(m.Entry, m["MANE-Select"].map(tags)))
seqt = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
DISP = dict(zip(seqt.base, seqt.displayed_isoid))


def status(canonical_acc, isoform_acc):
    base = canonical_acc.split("-")[0]
    t = M.get(base)
    if t is None:
        return "entry absent"
    if not t:
        return "no MANE"
    if "displayed" in t or DISP.get(base) in t:
        return "canonical"
    if isoform_acc in t:
        return "isoform"
    return "another isoform"


pa = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
pa = pa[pa["class"] != "unrelated"].copy()
pa["mane"] = [status(c, i) for c, i in zip(pa.canonical_acc, pa.isoform_acc)]
B = pa[pa.both_pred == True].copy()
res = {"pairs": {k: int(v) for k, v in pa.mane.value_counts().items()},
       "analysed": {k: int(v) for k, v in B.mane.value_counts().items()}}


def transitions(d):
    out = {}
    for g_, s in d.groupby("grp"):
        hc, nc = s[s.canon_hc == True], s[s.canon_hc == False]
        out[g_] = {"pairs": len(s), "canon_hc": len(hc), "losses": int((hc.iso_hc == False).sum()),
                   "loss_rate": r4(100 * (hc.iso_hc == False).mean()), "gains": int((nc.iso_hc == True).sum()),
                   "gain_rate": r4(100 * (nc.iso_hc == True).mean())}
    return out


res["transitions_all"] = transitions(B)
res["transitions_canonical_mane"] = transitions(B[B.mane == "canonical"])

cr = pd.read_csv("results/canonical_retention.tsv", sep="\t")
cr = cr.merge(pa[["canonical_acc", "isoform_acc", "both_pred", "mane"]], on=["canonical_acc", "isoform_acc"], how="left")
cr = cr[cr.both_pred == True].copy()
cr["loss_i"] = cr.loss.astype(int)
cr["grp"] = pd.Categorical(cr.grp, categories=["addition", "subtraction", "replacement"])
cr["iface_lost10"] = 10 * (1 - cr.iface_retained)
cr["conf_removed10"] = 10 * cr.conf_removed


def logit(formula, data):
    mo = smf.logit(formula, data=data).fit(disp=0)
    ci = mo.conf_int()
    return {"n": int(mo.nobs), "terms": {t: {"OR": r4(np.exp(mo.params[t])), "lo": r4(np.exp(ci.loc[t, 0])),
                                             "hi": r4(np.exp(ci.loc[t, 1]))} for t in mo.params.index if t != "Intercept"}}


for lab, d in [("all", cr), ("canonical_mane", cr[cr.mane == "canonical"])]:
    res[f"retention_{lab}"] = {
        "pairs": len(d),
        "loss_all_kept": r4(100 * d[d.iface_retained == 1].loss.mean()),
        "loss_lt50_kept": r4(100 * d[d.iface_retained < 0.5].loss.mean()),
        "group_only": logit("loss_i ~ C(grp)", d),
        "with_interface": logit("loss_i ~ C(grp) + iface_lost10", d),
        "with_interface_and_fold": logit("loss_i ~ C(grp) + iface_lost10 + conf_removed10", d)}
json.dump(res, open("results/mane_reference.json", "w"), indent=1)
print(json.dumps(res, indent=1))
