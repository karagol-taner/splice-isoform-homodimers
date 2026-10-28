import json
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
N = {}
D4 = 4


def r4(x, d=D4):
    return round(float(x), d)


def wilson(k, n):
    lo, hi = proportion_confint(k, n, method="wilson")
    return [r4(100 * lo), r4(100 * hi)]


def group_of(cls):
    if cls in ("truncation", "internal deletion"):
        return "subtraction"
    if cls in ("addition", "replacement"):
        return cls
    return "other"


cl = pd.read_csv("results/isoform_pairs_classified.tsv", sep="\t")
N["isoform_records"] = int(len(cl))
r = cl[cl.in_afdb_length_range == True].copy()
N["pairs_in_range"] = int(len(r))
N["class_counts"] = r["class"].value_counts().to_dict()
N["class_by_length_relation"] = {f"{a}|{b}": int(v) for (a, b), v in
                                 r.groupby(["class", "length_relation"]).size().items()}
N["subclass_counts"] = {f"{a}|{b}": int(v) for (a, b), v in r.groupby(["class", "subclass"]).size().items()}
N["novel_ge10_all"] = int((r.novel >= 10).sum())
N["novel_gt0_all"] = int((r.novel > 0).sum())
N["novel_median_by_class"] = {k: float(v) for k, v in r.groupby("class").novel.median().items()}
tab1 = (r.groupby(["class", "length_relation"])
        .agg(pairs=("isoform_acc", "size"), novel_ge10=("novel", lambda x: int((x >= 10).sum())),
             median_novel=("novel", "median"), median_deleted=("deleted", "median"))
        .reset_index())
tab1.to_csv("results/table_partition.tsv", sep="\t", index=False)
N["class_table"] = {c: {"pairs": int(len(s_)), "novel_ge10": int((s_.novel >= 10).sum()),
                        "median_novel": float(s_.novel.median()), "median_deleted": float(s_.deleted.median())}
                    for c, s_ in r.groupby("class")}
def site_group(row):
    parts = (row.novel_N + row.deleted_N > 0, row.novel_C + row.deleted_C > 0, row.novel_I + row.deleted_I > 0)
    if sum(parts) > 1:
        return "more than one site"
    return "N-terminus" if parts[0] else ("C-terminus" if parts[1] else ("interior" if parts[2] else "none"))
r["site"] = r.apply(site_group, axis=1)
N["site_by_class"] = {f"{c}|{g}": int(v) for (c, g), v in r.groupby(["class", "site"]).size().items()}
N["deleted_median_by_class"] = {k: float(v) for k, v in r.groupby("class").deleted.median().items()}

cx = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t")
hcset = cx[cx.hc_homodimer]
N["afdb_accessions_queried"] = int(len(cx))
N["afdb_hc_homodimers"] = int(len(hcset))
N["manifest_rule_check"] = json.load(open("results/manifest_rule_check.json"))
mf = pd.read_csv("results/manifest_per_accession.tsv", sep="\t")
N["manifest_chain_len_range"] = [int(mf.chain_len.min()), int(mf.chain_len.max())]
up = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
up["len"] = up.seq.str.len()
up["predicted"] = up.acc.isin(set(mf.acc))
N["sp_coverage_by_length"] = {lab: {"sequences": int(len(s)), "predicted": int(s.predicted.sum()),
                                    "rate": r4(100 * s.predicted.mean())}
                              for lab, s in [("lt30", up[up.len < 30]), ("30to1000", up[(up.len >= 30) & (up.len <= 1000)]),
                                             ("1001to1500", up[(up.len > 1000) & (up.len <= 1500)]),
                                             ("1001to1250", up[(up.len > 1000) & (up.len <= 1250)]),
                                             ("1251to1500", up[(up.len > 1250) & (up.len <= 1500)]), ("gt1500", up[up.len > 1500])]}

p_all = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
N["pairs_classified"] = int(len(p_all))
N["pairs_iso_not_predicted"] = int((~p_all.iso_pred).sum())
ni = p_all[~p_all.iso_pred]
N["iso_not_predicted_len"] = {"lt30": int((ni.iso_len < 30).sum()), "gt1000": int((ni.iso_len > 1000).sum()),
                              "30to1000": int(((ni.iso_len >= 30) & (ni.iso_len <= 1000)).sum())}
N["api_losses_iso_not_predicted"] = int((p_all.canon_hc & ~p_all.iso_pred).sum())
mm = p_all.canon_pred & p_all.iso_pred & ~p_all.model_len_ok
N["set_aside"] = {"canon_not_predicted": int((~p_all.canon_pred).sum()),
                  "iso_not_predicted_canon_predicted": int((p_all.canon_pred & ~p_all.iso_pred).sum()),
                  "model_length_differs": int(mm.sum())}
N["model_length_differs_canonicals"] = int(p_all.loc[mm & (p_all.canon_model_len != p_all.ref_len), "canonical_acc"].nunique())
N["model_length_differs_isoforms"] = int((mm & (p_all.iso_model_len != p_all.iso_len)).sum())


def transitions(df):
    out = []
    for g, sub in list(df.groupby("grp")) + [("all", df)]:
        nb_g = int((~sub.canon_hc).sum()); k_g = int((~sub.canon_hc & sub.iso_hc).sum())
        nb_l = int(sub.canon_hc.sum()); k_l = int((sub.canon_hc & ~sub.iso_hc).sum())
        out.append({"group": g, "pairs": int(len(sub)), "canon_not_hc": nb_g, "gains": k_g,
                    "gain_rate": r4(100 * k_g / nb_g), "gain_ci": wilson(k_g, nb_g),
                    "canon_hc": nb_l, "losses": k_l, "loss_rate": r4(100 * k_l / nb_l), "loss_ci": wilson(k_l, nb_l),
                    "both_hc": int((sub.canon_hc & sub.iso_hc).sum())})
    return pd.DataFrame(out).set_index("group")


N["transitions_all_pairs"] = json.loads(transitions(p_all).to_json(orient="index"))
p = p_all[p_all.both_pred].copy()
N["pairs_analysed"] = int(len(p))
N["canonicals_analysed"] = int(p.canonical_acc.nunique())
N["canonical_hc"] = int(p.drop_duplicates("canonical_acc").canon_hc.sum())
N["isoform_hc"] = int(p.iso_hc.sum())
T = transitions(p)
T.to_csv("results/table_transitions.tsv", sep="\t")
N["transitions"] = json.loads(T.to_json(orient="index"))


def fisher_vs_addition(df, base, event, tag):
    a = df[base & (df.grp == "addition")]
    for g in ["subtraction", "replacement"]:
        b_ = df[base & (df.grp == g)]
        kb, ka = int(event[b_.index].sum()), int(event[a.index].sum())
        orr, pv = fisher_exact([[kb, len(b_) - kb], [ka, len(a) - ka]])
        N[f"fisher_{tag}_{g}_vs_addition"] = {"OR": r4(orr), "p": float(f"{pv:.4g}")}


fisher_vs_addition(p, ~p.canon_hc, p.iso_hc, "gain")
fisher_vs_addition(p, p.canon_hc, ~p.iso_hc, "loss")

ipo = []
for g, s in list(p.groupby("grp")) + [("all", p)]:
    below, above = s[s.canon_ipsae < 0.6], s[s.canon_ipsae >= 0.6]
    kg, kl = int((below.iso_ipsae >= 0.6).sum()), int((above.iso_ipsae < 0.6).sum())
    ipo.append({"group": g, "canon_below": int(len(below)), "gains": kg, "gain_rate": r4(100 * kg / len(below)),
                "gain_ci": wilson(kg, len(below)), "canon_above": int(len(above)), "losses": kl,
                "loss_rate": r4(100 * kl / len(above)), "loss_ci": wilson(kl, len(above))})
IPO = pd.DataFrame(ipo).set_index("group")
IPO.to_csv("results/table_transitions_ipsae_only.tsv", sep="\t")
N["transitions_ipsae_only"] = json.loads(IPO.to_json(orient="index"))
pi = p.copy()
pi["g_ip"] = (pi.canon_ipsae < 0.6) & (pi.iso_ipsae >= 0.6)
fisher_vs_addition(pi, pi.canon_ipsae < 0.6, pi.g_ip, "ipsae_gain")
gains = p[~p.canon_hc & p.iso_hc]
losses = p[p.canon_hc & ~p.iso_hc]
N["gain_types"] = {g: {"n": int(len(s)), "canon_ipsae_ge06": int((s.canon_ipsae >= 0.6).sum()),
                       "canon_ipsae_05to06": int(((s.canon_ipsae >= 0.5) & (s.canon_ipsae < 0.6)).sum()),
                       "canon_ipsae_lt05": int((s.canon_ipsae < 0.5).sum())}
                   for g, s in list(gains.groupby("grp")) + [("all", gains)]}
N["loss_types"] = {g: {"n": int(len(s)), "iso_ipsae_ge06": int((s.iso_ipsae >= 0.6).sum()),
                       "iso_ipsae_05to06": int(((s.iso_ipsae >= 0.5) & (s.iso_ipsae < 0.6)).sum()),
                       "iso_ipsae_lt05": int((s.iso_ipsae < 0.5).sum())}
                   for g, s in list(losses.groupby("grp")) + [("all", losses)]}
for tag, s in [("gains_pdockq2_type", gains[gains.canon_ipsae >= 0.6]), ("gains_ipsae_type", gains[gains.canon_ipsae < 0.6])]:
    N[tag] = {"n": int(len(s)), "canon_mono_plddt_median": r4(s.canon_mono_plddt.median()),
              "iso_mono_plddt_median": r4(s.iso_mono_plddt.median()),
              "d_mono_plddt_median": r4((s.iso_mono_plddt - s.canon_mono_plddt).median()),
              "canon_pdockq2_median": r4(s.canon_pdockq2.median()),
              "iso_pdockq2_median": r4(s.iso_pdockq2.median())}
can_hi, can_lo = p.canon_ipsae >= 0.8, p.canon_ipsae <= 0.4
rob_l = can_hi & (p.iso_ipsae <= 0.4)
rob_g = can_lo & (p.iso_ipsae >= 0.8)
N["robust_losses"] = {g: {"k": int(rob_l[s.index].sum()), "n": int(can_hi[s.index].sum()),
                          "call_losses": int((rob_l & p.canon_hc & ~p.iso_hc)[s.index].sum())}
                      for g, s in list(p.groupby("grp")) + [("all", p)]}
N["robust_gains"] = {g: {"k": int(rob_g[s.index].sum()), "n": int(can_lo[s.index].sum()),
                         "call_gains": int((rob_g & ~p.canon_hc & p.iso_hc)[s.index].sum())}
                     for g, s in list(p.groupby("grp")) + [("all", p)]}

q = p.dropna(subset=["canon_mono_plddt", "iso_mono_plddt", "canon_vlow", "iso_vlow"]).copy()
q["log_canon_len"] = np.log10(q.ref_len)
q["d_plddt"] = q.iso_mono_plddt - q.canon_mono_plddt
q["canon_vlow10"] = 10 * q.canon_vlow
q["gain"] = (~q.canon_hc & q.iso_hc).astype(int)
q["loss"] = (q.canon_hc & ~q.iso_hc).astype(int)
q["grp"] = pd.Categorical(q.grp, categories=["addition", "subtraction", "replacement"])


def logit(formula, data):
    m = smf.logit(formula, data=data).fit(disp=0)
    ci = m.conf_int()
    return {"n": int(m.nobs), "aic": round(float(m.aic), 1),
            "terms": {t: {"OR": r4(np.exp(m.params[t]), 6), "lo": r4(np.exp(ci.loc[t, 0]), 6),
                          "hi": r4(np.exp(ci.loc[t, 1]), 6), "p": float(f"{m.pvalues[t]:.4g}")} for t in m.params.index}}


reg = {"gain": logit("gain ~ C(grp) + log_canon_len + canon_vlow10 + d_plddt", q[~q.canon_hc]),
       "loss": logit("loss ~ C(grp) + log_canon_len + canon_vlow10 + d_plddt", q[q.canon_hc])}
qg = q[~q.canon_hc].copy()
qg["gain_ip"] = ((qg.canon_ipsae < 0.6) & (qg.iso_ipsae >= 0.6)).astype(int)
reg["gain_ipsae_only"] = logit("gain_ip ~ C(grp) + log_canon_len + canon_vlow10 + d_plddt", qg[qg.canon_ipsae < 0.6])
N["regression"] = reg
pd.DataFrame([{"model": k, "n": d["n"], "term": t, **v} for k, d in reg.items() for t, v in d["terms"].items()]) \
    .to_csv("results/table_regression.tsv", sep="\t", index=False)
gq = q[~q.canon_hc].copy()
gq["dbin"] = pd.cut(gq.d_plddt, [-100, -5, 0, 5, 10, 100], labels=["lt-5", "-5to0", "0to5", "5to10", "gt10"])
gb = []
for b_, s in gq.groupby("dbin", observed=True):
    gb.append({"grp": "all", "bin": b_, "n": int(len(s)), "gains": int(s.gain.sum()),
               "rate": r4(100 * s.gain.mean()), "ci": wilson(int(s.gain.sum()), len(s)),
               "gains_pdockq2_type": int((s.gain.astype(bool) & (s.canon_ipsae >= 0.6)).sum())})
for (g, b_), s in gq.groupby(["grp", "dbin"], observed=True):
    gb.append({"grp": g, "bin": b_, "n": int(len(s)), "gains": int(s.gain.sum()),
               "rate": r4(100 * s.gain.mean()), "ci": wilson(int(s.gain.sum()), len(s)),
               "gains_pdockq2_type": int((s.gain.astype(bool) & (s.canon_ipsae >= 0.6)).sum())})
GB = pd.DataFrame(gb)
GB.to_csv("results/table_gain_by_dplddt.tsv", sep="\t", index=False)
N["gain_by_dplddt"] = {f"{r.grp}|{r.bin}": {"n": r.n, "gains": r.gains, "rate": r.rate, "ci": r.ci,
                                            "gains_pdockq2_type": r.gains_pdockq2_type} for r in GB.itertuples()}
N["d_plddt_median_by_grp"] = {g: r4(s.d_plddt.median()) for g, s in gq.groupby("grp", observed=True)}

cr0 = pd.read_csv("results/canonical_retention.tsv", sep="\t")
cr0 = cr0.merge(p_all[["canonical_acc", "isoform_acc", "iso_pred", "both_pred", "canon_ipsae", "iso_ipsae", "canon_pdockq2",
                       "iso_pdockq2", "canon_mono_plddt", "iso_mono_plddt", "canon_vlow"]], on=["canonical_acc", "isoform_acc"], how="left")
N["retention_pairs_all"] = int(len(cr0))
N["retention_skipped_model_sequence"] = int(p_all.canon_hc.sum() - len(cr0))
cr = cr0[cr0.both_pred].copy()
N["retention_pairs"] = int(len(cr))
cr["ret_bin"] = pd.cut(cr.iface_retained, [0, 0.5, 0.9, 1.0, 1.0001], right=False, labels=["lt50", "50to90", "90to100", "all"])
rb = []
for (g, b_), s in list(cr.groupby(["grp", "ret_bin"], observed=True)) + [(("all", b_), s) for b_, s in cr.groupby("ret_bin", observed=True)]:
    rb.append({"grp": g, "bin": b_, "n": int(len(s)), "losses": int(s.loss.sum()),
               "rate": r4(100 * s.loss.mean()), "ci": wilson(int(s.loss.sum()), len(s)),
               "iso_ipsae_median": r4(s.iso_ipsae.median())})
RB = pd.DataFrame(rb)
RB.to_csv("results/table_loss_by_retention.tsv", sep="\t", index=False)
N["loss_by_retention"] = {f"{r.grp}|{r.bin}": {"n": r.n, "losses": r.losses, "rate": r.rate, "ci": r.ci,
                                               "iso_ipsae_median": r.iso_ipsae_median} for r in RB.itertuples()}
rho, pv = spearmanr(cr.iface_retained, cr.iso_ipsae)
N["iso_ipsae_vs_retention_spearman"] = {"rho": r4(rho), "p": float(f"{pv:.3g}"), "n": int(len(cr)),
                                        "x": "fraction of canonical interface residues kept"}
lo_ = cr[cr.loss]
N["losses_intact_iso_ipsae_ge05"] = {"k": int(((lo_.iface_retained == 1) & (lo_.iso_ipsae >= 0.5)).sum()),
                                     "n": int((lo_.iface_retained == 1).sum())}
N["losses_intact_iso_ipsae_ge06"] = {"k": int(((lo_.iface_retained == 1) & (lo_.iso_ipsae >= 0.6)).sum()),
                                     "n": int((lo_.iface_retained == 1).sum())}
N["losses_halfremoved_iso_ipsae_lt02"] = {"k": int(((lo_.iface_retained < 0.5) & (lo_.iso_ipsae < 0.2)).sum()),
                                          "n": int((lo_.iface_retained < 0.5).sum())}
cr["loss_i"] = cr.loss.astype(int)
cr["grp"] = pd.Categorical(cr.grp, categories=["addition", "subtraction", "replacement"])
cr["iface_lost10"] = 10 * (1 - cr.iface_retained)
cr["conf_removed10"] = 10 * cr.conf_removed
cr["ins_struct"] = (cr.ins_in_structured > 0).astype(int)
cr["d_plddt"] = cr.iso_mono_plddt - cr.canon_mono_plddt
cr["canon_vlow10"] = 10 * cr.canon_vlow
cr["log_canon_len"] = np.log10(cr.can_len)
crq = cr.dropna(subset=["d_plddt", "canon_vlow"])
ret_models = {"group_only": logit("loss_i ~ C(grp)", cr),
              "with_interface": logit("loss_i ~ C(grp) + iface_lost10", cr),
              "with_interface_and_fold": logit("loss_i ~ C(grp) + iface_lost10 + conf_removed10", cr),
              "full": logit("loss_i ~ C(grp) + iface_lost10 + conf_removed10 + d_plddt + canon_vlow10 + log_canon_len", crq),
              "intact_interface": logit("loss_i ~ C(grp) + d_plddt + canon_vlow10 + log_canon_len", crq[crq.iface_retained == 1])}
N["retention_models"] = ret_models
pd.DataFrame([{"model": k, "n": d["n"], "term": t, **v} for k, d in ret_models.items() for t, v in d["terms"].items()]) \
    .to_csv("results/table_retention_models.tsv", sep="\t", index=False)
cr["ipsae_bin"] = pd.cut(cr.canon_ipsae, [0.595, 0.65, 0.70, 0.80, 1.0], labels=["0.60to0.65", "0.65to0.70", "0.70to0.80", "0.80to1.00"])
rib = []
for (b_, ib), s in cr.groupby(["ret_bin", "ipsae_bin"], observed=True):
    rib.append({"bin": b_, "ipsae_bin": ib, "n": int(len(s)), "losses": int(s.loss.sum()),
                "rate": r4(100 * s.loss.mean()), "ci": wilson(int(s.loss.sum()), len(s))})
RIB = pd.DataFrame(rib)
RIB.to_csv("results/table_loss_by_retention_and_ipsae.tsv", sep="\t", index=False)
N["loss_by_retention_ipsae"] = {f"{r.bin}|{r.ipsae_bin}": {"n": r.n, "losses": r.losses, "rate": r.rate, "ci": r.ci}
                                for r in RIB.itertuples()}
ad = cr[cr.grp == "addition"]
N["addition_insert_structured"] = {"n": int((ad.ins_struct == 1).sum()), "losses": int(ad.loss[ad.ins_struct == 1].sum()),
                                   "n_other": int((ad.ins_struct == 0).sum()), "losses_other": int(ad.loss[ad.ins_struct == 0].sum())}
cr.to_csv("results/canonical_retention_analysed.tsv", sep="\t", index=False)

U = [json.loads(l) for l in open("data/uniparc_raw.jsonl")]
CP = [json.loads(l) for l in open("data/copies_complex_raw.jsonl")]
hcmap = dict(zip(cx.acc, cx.hc_homodimer))
mfi = mf.set_index("acc")
MC = pd.read_csv("data/manifest_extract/manifest_copies.tsv.gz", sep="\t")
MC["ipsae"] = MC[["ipSAE_AB", "ipSAE_BA"]].max(axis=1)
MC["pdockq2"] = MC[["pDockQ2_AB", "pDockQ2_BA"]].max(axis=1)
cip, cpq = dict(zip(MC.uniprotAccession, MC.ipsae)), dict(zip(MC.uniprotAccession, MC.pdockq2))
SE = {json.loads(l)["acc"]: json.loads(l) for l in open("data/copies_served_entries.jsonl")}
api_ip = {a: r["homo"][0]["complexPredictionAccuracy_ipSAE"] for a, r in SE.items() if r["homo"]}
api_pq = {a: r["homo"][0]["complexPredictionAccuracy_pDockQ2"] for a, r in SE.items() if r["homo"]}
rows = []
for c in CP:
    h = c["human_source"]
    rows.append({"copy": c["acc"], "human": h, "taxid": c["taxid"], "copy_deposited": bool(c["homo"]),
                 "human_deposited": bool(hcmap.get(h, False)), "copy_in_manifest": c["acc"] in cip,
                 "copy_predicted": (c["acc"] in cip) or bool(c["homo"]),
                 "copy_ipsae": cip.get(c["acc"], np.nan), "human_ipsae": mfi.ipsae.get(h, np.nan),
                 "copy_pdockq2": cpq.get(c["acc"], np.nan), "human_pdockq2": mfi.pdockq2.get(h, np.nan),
                 "copy_api_ipsae": api_ip.get(c["acc"], np.nan), "copy_api_pdockq2": api_pq.get(c["acc"], np.nan)})
D = pd.DataFrame(rows)
D.to_csv("results/identical_copies.tsv", sep="\t", index=False)
N["copies_sample_human"] = int(len(U))
N["copies_sample_human_hc"] = int(sum(bool(hcmap.get(u["acc"], False)) for u in U))
N["copies_sample_human_not_predicted"] = int(sum(u["acc"] not in mfi.index for u in U))
N["copies_human_with_copy"] = int(sum(1 for u in U if u["copies"]))
N["copies_queried"] = int(len(D))
N["copies_species"] = int(D.taxid.nunique())
N["copies_human_taxon"] = int((D.taxid == 9606).sum())
dd = D[D.copy_in_manifest & D.human_ipsae.notna()].copy()
dd["absdiff"] = (dd.copy_ipsae - dd.human_ipsae).abs()
N["copies_predicted"] = int(D.copy_predicted.sum())
N["copies_in_manifest"] = int(D.copy_in_manifest.sum())
N["copies_predicted_human_sources"] = int(D[D.copy_predicted].human.nunique())
nr = D[D.copy_deposited & ~D.copy_in_manifest].copy()
N["copies_served_without_manifest_row"] = {
    "copies": int(len(nr)), "human_sources": int(nr.human.nunique()),
    "human_sources_deposited": int(nr.drop_duplicates("human").human_deposited.sum()),
    "api_ipsae_within001_of_human_manifest": int(((nr.copy_api_ipsae - nr.human_ipsae).abs() <= 0.01 + 1e-9).sum()),
    "api_pdockq2_within001_of_human_manifest": int(((nr.copy_api_pdockq2 - nr.human_pdockq2).abs() <= 0.01 + 1e-9).sum())}
N["copies_ipsae_absdiff"] = {"median": r4(dd.absdiff.median()), "q90": r4(dd.absdiff.quantile(0.9)),
                             "max": r4(dd.absdiff.max()), "le002": r4((dd.absdiff <= 0.02).mean()),
                             "le005": r4((dd.absdiff <= 0.05).mean()), "n": int(len(dd))}
rho, pv = spearmanr(dd.human_ipsae, dd.copy_ipsae)
N["copies_ipsae_spearman"] = {"rho": r4(rho), "p": float(f"{pv:.3g}")}
rule_c = (dd.copy_ipsae >= 0.6) & (dd.copy_pdockq2 >= 0.23)
rule_h = (dd.human_ipsae >= 0.6) & (dd.human_pdockq2 >= 0.23)
N["copies_call_agreement"] = {"agree": int((rule_c == rule_h).sum()), "n": int(len(dd)),
                              "both_hc": int((rule_c & rule_h).sum()), "neither": int((~rule_c & ~rule_h).sum())}
near = dd[(dd.human_ipsae >= 0.5) & (dd.human_ipsae < 0.7)]
N["copies_near_threshold"] = {"n": int(len(near)), "absdiff_median": r4(near.absdiff.median()),
                              "ipsae_side_agree": int(((near.copy_ipsae >= 0.6) == (near.human_ipsae >= 0.6)).sum())}
dis = dd[rule_c != rule_h]
N["copies_call_disagreements"] = {"copies": int(len(dis)), "human_sources": int(dis.human.nunique()),
                                  "sources_total": int(dd.human.nunique()),
                                  "by_pdockq2_only": int(((dis.copy_ipsae >= 0.6) == (dis.human_ipsae >= 0.6)).sum()),
                                  "by_ipsae": int(((dis.copy_ipsae >= 0.6) != (dis.human_ipsae >= 0.6)).sum())}
pq_d = (dd.copy_pdockq2 - dd.human_pdockq2).abs()
N["copies_pdockq2_absdiff"] = {"median": r4(pq_d.median()), "q90": r4(pq_d.quantile(0.9)),
                               "max": r4(pq_d.max()), "le005": r4((pq_d <= 0.05).mean())}
N["copies_rule_reproduces_deposition"] = {"agree": int((dd.copy_deposited == rule_c).sum()), "n": int(len(dd))}

mech = pd.read_csv("results/pair_mechanism.tsv", sep="\t")
nov = mech[mech.novel_model > 0]
N["hc_isoform_models"] = int(len(mech))
N["hc_isoform_models_with_added"] = int(len(nov))
bp_pairs = set(zip(p.canonical_acc, p.isoform_acc))
in_set = np.array([k in bp_pairs for k in zip(nov.canonical_acc, nov.isoform_acc)], dtype=bool)
g = nov[(~nov.canon_hc).values & in_set].copy()
N["gains_with_added_models"] = int(len(g))
N["gains_added_at_iface"] = int((g.n_iface_novel > 0).sum())
N["gains_added_ge5_conf_iface"] = int((g.n_iface_novel_conf >= 5).sum())
N["gains_added_novel_novel_contacts"] = int((g.novel_novel_contacts > 0).sum())
oth = g[g.n_iface_novel_conf < 5]
N["gains_other_iface_residues"] = {"pairs": int(len(oth)), "iface": int(oth.n_iface.sum()), "added": int(oth.n_iface_novel.sum()),
                                   "added_conf": int(oth.n_iface_novel_conf.sum())}
top = g[g.n_iface_novel_conf >= 5].merge(p[["canonical_acc", "isoform_acc", "canon_ipsae", "iso_ipsae", "canon_pdockq2", "iso_pdockq2"]],
                                         on=["canonical_acc", "isoform_acc"], how="left")
N["gains_ge5_conf_list"] = top[["gene", "canonical_acc", "isoform_acc", "class", "n_iface", "n_iface_novel", "n_iface_novel_conf",
                                "mean_plddt_iface_novel", "canon_ipsae", "iso_ipsae", "canon_pdockq2", "iso_pdockq2"]] \
    .round(4).to_dict(orient="records")
R = pd.read_csv("results/residues_novel_models.tsv.gz", sep="\t")
nv, sh = R[R.novel], R[~R.novel]
N["added_residues"] = int(len(nv))
N["added_plddt_lt50"] = r4((nv.plddt < 50).mean(), 6)
N["added_plddt_ge70"] = r4((nv.plddt >= 70).mean(), 6)
N["shared_plddt_lt50"] = r4((sh.plddt < 50).mean(), 6)
N["shared_plddt_ge70"] = r4((sh.plddt >= 70).mean(), 6)
S = R[R.rsa_A >= 0.25].copy()
N["iface_share_added"] = r4(S[S.iface & S.novel].shape[0] / S[S.iface].shape[0], 6)
N["surface_share_added"] = r4(S[S.novel].shape[0] / S.shape[0], 6)
S["iface_i"] = S.iface.astype(int)
S["novel_i"] = S.novel.astype(int)
S["pbin"] = pd.cut(S.plddt, [0, 50, 70, 90, 101], right=False, labels=["lt50", "50to70", "70to90", "ge90"])
for name, f in [("unadjusted", "iface_i ~ novel_i"), ("plddt_adjusted", "iface_i ~ novel_i + C(pbin)")]:
    m = smf.logit(f, data=S).fit(disp=0)
    ci = m.conf_int().loc["novel_i"]
    N[f"residue_logit_{name}"] = {"OR": r4(np.exp(m.params["novel_i"]), 6),
                                  "lo": r4(np.exp(ci.iloc[0]), 6), "hi": r4(np.exp(ci.iloc[1]), 6),
                                  "p": float(f"{m.pvalues['novel_i']:.3g}"), "n_residues": int(m.nobs)}
band = (S.groupby(["pbin", "novel"], observed=True).agg(n=("iface", "size"), iface=("iface", "sum")).reset_index())
band["iface_rate"] = (100 * band.iface / band.n).round(4)
band.to_csv("results/table_interface_by_plddt.tsv", sep="\t", index=False)
in_set_m = np.array([k in bp_pairs for k in zip(mech.canonical_acc, mech.isoform_acc)], dtype=bool)
hb = mech[mech.canon_hc.values & in_set_m]
N["hchc_canon_model_differs"] = int((hb.canon_model_matches_uniprot == False).sum())
hh = hb[hb.iface_jaccard_shared.notna()].copy()
hh["grp"] = hh["class"].map(group_of)
N["hchc_pairs_with_models"] = int(len(hh))
N["hchc_jaccard"] = {gg: {"n": int(len(s)), "median": r4(s.iface_jaccard_shared.median()),
                          "below_0.5": int((s.iface_jaccard_shared < 0.5).sum())}
                     for gg, s in hh.groupby("grp")}

cal = pd.read_csv("results/canonical_subunit_calibration.tsv", sep="\t")
cur_len = dict(zip(up.acc, up.len))
cal = cal[cal.acc.isin(set(mf.acc)) & (cal.acc.map(dict(zip(mf.acc, mf.chain_len))) == cal.acc.map(cur_len))]
N["calibration_canonicals_predicted"] = int(len(cal))
has_text = cal.subunit.fillna("") != ""
for label, mask in [("annotated_homodimer", cal.homodimer), ("annotated_other_homooligomer", cal.homo_any & ~cal.homodimer),
                    ("annotated_monomer", cal.monomer_only),
                    ("text_other", has_text & ~cal.homo_any & ~cal.monomer_only),
                    ("no_subunit_text", ~has_text)]:
    k, n = int(cal.hc[mask].sum()), int(mask.sum())
    N[f"calibration_{label}"] = {"hc": k, "n": n, "rate": r4(100 * k / n), "ci": wilson(k, n)}

json.dump(N, open("results/numbers.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps({"pairs_in_range": N["pairs_in_range"], "pairs_analysed": N["pairs_analysed"],
                  "gains": N["transitions"]["all"]["gains"], "losses": N["transitions"]["all"]["losses"]}, indent=1))
