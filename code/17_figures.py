import json
import os

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
    "pdf.fonttype": 42, "figure.dpi": 100,
    "svg.hashsalt": "splice-isoform-homodimers",
})
N = json.load(open("results/numbers.json"))
NN = json.load(open("results/numbers_new.json"))
GROUPS = ["addition", "subtraction", "replacement"]
GLAB = {"addition": "Addition only", "subtraction": "Removal only", "replacement": "Replacement"}
GCOL = {"addition": "#9ecae1", "subtraction": "#4292c6", "replacement": "#08519c"}
GMARK = {"addition": "o", "subtraction": "s", "replacement": "^"}
ADDED, SHARED = "#d94801", "#636363"
THR = dict(color="black", lw=0.6, ls="--")


_LETTERS = []


def letter(ax, s, **_):
    _LETTERS.append((ax, s))


def _groups(vals, tol=0.03):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out, cur = [], [order[0]]
    for i in order[1:]:
        if vals[i] - vals[cur[-1]] <= tol:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    return out


def _place_letters(fig):
    items = [(ax, s) for ax, s in _LETTERS if ax.figure is fig]
    if not items:
        return
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    bbs = [ax.get_tightbbox(rend, bbox_extra_artists=[]).transformed(inv) for ax, _ in items]
    pos = [ax.get_position() for ax, _ in items]
    ys, xs = [b.y1 for b in bbs], [b.x0 for b in bbs]
    for grp in _groups([p.y1 for p in pos]):
        top = max(bbs[i].y1 for i in grp)
        for i in grp:
            ys[i] = top
    for grp in _groups([b.x0 for b in bbs], tol=0.02):
        left = min(bbs[i].x0 for i in grp)
        for i in grp:
            xs[i] = left
    for (ax, s), x, y in zip(items, xs, ys):
        fig.text(x - 0.012, y + 0.004, s, fontsize=9, fontweight="bold", va="top", ha="left")
    _LETTERS[:] = [it for it in _LETTERS if it[0].figure is not fig]


def save(fig, name):
    _place_letters(fig)
    fig.savefig(f"figures/{name}.png", dpi=400, bbox_inches="tight", metadata={"Software": None})
    fig.savefig(f"figures/{name}.svg", bbox_inches="tight", metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote", name)


cl = pd.read_csv("results/isoform_pairs_classified.tsv", sep="\t")
r = cl[cl.in_afdb_length_range == True].copy()
classes = ["truncation", "internal deletion", "replacement", "addition"]
clab = {"truncation": "Terminal truncation", "internal deletion": "Internal deletion",
        "replacement": "Replacement", "addition": "Addition"}

fig = plt.figure(figsize=(7.0, 5.4))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.62)

ax = fig.add_subplot(gs[0, 0])
rel_cols = {"shorter": "#525252", "equal": "#969696", "longer": "#d9d9d9"}
y = np.arange(len(classes))[::-1]
for yi, c in zip(y, classes):
    left = 0
    for rel in ["shorter", "equal", "longer"]:
        v = int(((r["class"] == c) & (r.length_relation == rel)).sum())
        if v:
            ax.barh(yi, v, left=left, color=rel_cols[rel], edgecolor="black", linewidth=0.4, height=0.62)
            left += v
    ax.text(left + 150, yi, f"{left:,}", va="center", fontsize=6.5)
ax.set_yticks(y)
ax.set_yticklabels([clab[c] for c in classes])
ax.set_xlabel("Isoform-canonical pairs")
ax.set_xlim(0, 13200)
ax.legend(handles=[Patch(facecolor=rel_cols[k], edgecolor="black", linewidth=0.4, label=f"Isoform {k}")
                   for k in ["shorter", "equal", "longer"]], frameon=False, loc="lower right")
ax.set_title(f"Isoform classes ({N['pairs_in_range']:,} pairs)", loc="left")
letter(ax, "a", dx=0.13)

ax = fig.add_subplot(gs[0, 1])
SB = N["site_by_class"]
sites = ["N-terminus", "C-terminus", "interior", "more than one site"]
site_cols = {"N-terminus": "#08519c", "C-terminus": "#6baed6", "interior": "#fd8d3c", "more than one site": "#d9d9d9"}
for yi, c in zip(y, classes):
    left = 0
    for g in sites:
        v = SB.get(f"{c}|{g}", 0)
        if v:
            ax.barh(yi, v, left=left, color=site_cols[g], edgecolor="black", linewidth=0.4, height=0.62)
            left += v
ax.set_yticks(y)
ax.set_yticklabels([clab[c] for c in classes])
ax.set_xlabel("Isoform-canonical pairs")
ax.set_xlim(0, 13200)
ax.legend(handles=[Patch(facecolor=site_cols[g], edgecolor="black", linewidth=0.4, label=g[0].upper() + g[1:])
                   for g in sites], frameon=False, loc="lower right", ncol=2, columnspacing=0.8, handlelength=1.4,
          borderaxespad=0.2, fontsize=6.0)
ax.set_title("Where the isoform differs", loc="left")
letter(ax, "b", dx=0.13)

ax = fig.add_subplot(gs[1, 0])
nov = r[r.novel > 0]
bins = np.unique(np.round(np.logspace(0, np.log10(nov.novel.max() + 1), 32)).astype(int))
for c, col in [("replacement", "#08519c"), ("addition", "#9ecae1")]:
    v = nov.loc[nov["class"] == c, "novel"]
    ax.hist(v, bins=bins, histtype="stepfilled" if c == "addition" else "step",
            color=col, edgecolor=col, linewidth=1.0, alpha=0.9 if c == "addition" else 1,
            label=f"{clab[c]} (median {v.median():.0f})")
top = 1.5 * ax.get_ylim()[1]
ax.set_ylim(0, top)
hmax = max(np.histogram(nov.loc[nov["class"] == c, "novel"], bins=bins)[0].max() for c in ["replacement", "addition"])
ax.set_yticks(range(0, int(hmax) + 1, 200))
ax.vlines(10, 0, 0.60 * top, color="black", lw=0.7, ls="--")
n_ge10_plotted = int((nov["class"].isin(["replacement", "addition"]) & (nov.novel >= 10)).sum())
ax.text(10.8, 0.985 * top, f"{n_ge10_plotted:,} replacements and additions\nadd 10 or more residues", fontsize=6.3,
        va="top")
ax.set_xscale("log")
ax.set_xlabel("Residues absent from the canonical sequence")
ax.set_ylabel("Pairs")
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.80))
ax.set_title("Added residues per isoform", loc="left")
letter(ax, "c", dx=0.08)

ax = fig.add_subplot(gs[1, 1])
rem = r[r.deleted > 0]
bins_d = np.unique(np.round(np.logspace(0, np.log10(rem.deleted.max() + 1), 32)).astype(int))
for c, col, ht in [("replacement", "#08519c", "step"), ("truncation", "#525252", "step"),
                   ("internal deletion", "#9ecae1", "stepfilled")]:
    v = rem.loc[rem["class"] == c, "deleted"]
    ax.hist(v, bins=bins_d, histtype=ht, color=col, edgecolor=col, linewidth=1.0, alpha=0.9 if ht == "stepfilled" else 1,
            label=f"{clab[c]} (median {v.median():.0f})")
top_d = 1.35 * ax.get_ylim()[1]
ax.set_ylim(0, top_d)
hmax_d = max(np.histogram(rem.loc[rem["class"] == c, "deleted"], bins=bins_d)[0].max()
             for c in ["replacement", "truncation", "internal deletion"])
ax.set_yticks(range(0, int(hmax_d) + 1, 200))
ax.set_xscale("log")
ax.set_xlabel("Canonical residues absent from the isoform")
ax.set_ylabel("Pairs")
ax.legend(frameon=False, loc="upper left", fontsize=6.0)
ax.set_title("Removed residues per isoform", loc="left")
letter(ax, "d", dx=0.08)
save(fig, "Figure1")


MH = pd.read_csv("results/table_monomer_plddt_hist.tsv", sep="\t")
MO = NN["monomers"]
fig = plt.figure(figsize=(7.0, 5.4))
gs = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.60)

ax = fig.add_subplot(gs[0, 0])
wdt = int(MH.bin_start.iloc[1] - MH.bin_start.iloc[0])
for col, color, lab, fill in [("shared", SHARED, "Shared", False), ("added", ADDED, "Added", True)]:
    dens = MH[col] / (MH[col].sum() * wdt)
    if fill:
        ax.bar(MH.bin_start, dens, width=wdt, align="edge", color=color, alpha=0.55, linewidth=0, label=lab)
    else:
        ax.step(np.append(MH.bin_start.values, 100), np.append(dens.values, dens.values[-1]), where="post",
                color=color, lw=1.1, label=lab)
ax.vlines([50, 70], 0, ax.get_ylim()[1], color="black", lw=0.6, ls=":")
ax.set_xlim(0, 100)
ax.set_xlabel("pLDDT, isoform monomer model")
ax.set_ylabel("Density")
ax.legend(frameon=False, loc="upper left", handlelength=1.2)
ax.set_title(f"{MO['added_residues']:,} added residues", loc="left")
letter(ax, "a", dx=0.085)

ax = fig.add_subplot(gs[0, 1])
st = MO["added_by_site"]
cats = [("Added\nN-term", st["N"]["lt50"]), ("Added\nC-term", st["C"]["lt50"]), ("Added\ninterior", st["I"]["lt50"]),
        ("Added\nall", MO["added_lt50"]), ("Shared", MO["shared_lt50"])]
cols2 = [ADDED, ADDED, ADDED, ADDED, SHARED]
ax.bar(range(len(cats)), [v for _, v in cats], color=cols2, edgecolor="black", linewidth=0.4, width=0.66)
ax.set_xticks(range(len(cats)))
ax.set_xticklabels([l for l, _ in cats], fontsize=5.6)
ax.set_ylabel("Residues with pLDDT <50 (%)")
ax.set_ylim(0, 85)
ax.set_title("Low confidence by location", loc="left")
letter(ax, "b", dx=0.09)

ax = fig.add_subplot(gs[1, 0])
SW = pd.read_csv("results/swap_segments.tsv", sep="\t").dropna(subset=["removed_mean_plddt", "added_mean_plddt"])
oth, hom = SW[~SW.homologous], SW[SW.homologous]
ax.scatter(oth.removed_mean_plddt, oth.added_mean_plddt, s=2, color="#bdbdbd", linewidth=0, rasterized=True,
           label=f"Other ({len(oth):,})")
ax.scatter(hom.removed_mean_plddt, hom.added_mean_plddt, s=8, facecolor=ADDED, edgecolor="black", linewidth=0.2,
           label=f"Homologous ({len(hom)})")
ax.plot([0, 100], [0, 100], color="#969696", lw=0.5, ls=":")
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_xlabel("Removed segment pLDDT (canonical)")
ax.set_ylabel("Added segment pLDDT (isoform)")
ax.legend(frameon=False, loc="upper left", fontsize=5.8, handletextpad=0.2, markerscale=1.4)
ax.set_title("Replacement swaps", loc="left")
letter(ax, "c", dx=0.085)

ax = fig.add_subplot(gs[1, 1])
SWL = NN["swaps"]
pr = [("Homologous\nswap", SWL["loss_homologous"], ADDED), ("Other\nswap", SWL["loss_other"], "#4292c6")]
for i, (lab, dd, col) in enumerate(pr):
    v, (lo, hi) = dd["rate"], dd["ci"]
    ax.bar(i, v, color=col, edgecolor="black", linewidth=0.4, width=0.6)
    ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="black", lw=0.7, capsize=3)
    ax.text(i, hi + 2, f"{dd['losses']}/{dd['n']}", ha="center", fontsize=6.0)
ax.set_xticks(range(2))
ax.set_xticklabels([l for l, _, _ in pr])
ax.set_xlim(-0.6, 1.6)
ax.set_ylim(0, 75)
ax.set_ylabel("Isoform not high-confidence (%)")
ax.set_title("Loss after a swap", loc="left")
letter(ax, "d", dx=0.09)
save(fig, "Figure2")


MF = pd.read_csv("results/manifest_per_accession.tsv", sep="\t")
CXa = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t")[["acc", "hc_homodimer"]]
Q5 = CXa.merge(MF, on="acc")
fig = plt.figure(figsize=(7.0, 2.7))
gs = fig.add_gridspec(1, 3, wspace=0.5)

ax = fig.add_subplot(gs[0, 0])
dep = Q5.hc_homodimer
ax.scatter(Q5.loc[~dep, "ipsae"], Q5.loc[~dep, "pdockq2"], s=0.8, color="#bdbdbd", linewidth=0, rasterized=True, label="Not deposited")
ax.scatter(Q5.loc[dep, "ipsae"], Q5.loc[dep, "pdockq2"], s=0.8, color="#08519c", linewidth=0, rasterized=True, label="Deposited")
ax.axvline(0.6, **THR)
ax.axhline(0.23, **THR)
ax.set_xlim(-0.02, 1.0)
ax.set_ylim(-0.02, 1.0)
ax.set_xlabel("ipSAE")
ax.set_ylabel("pDockQ2")
ax.legend(frameon=False, loc="upper left", fontsize=5.8, markerscale=6, handletextpad=0.1)
ax.set_title("Predicted homodimers", loc="left")
letter(ax, "a", dx=0.085)

ax = fig.add_subplot(gs[0, 1])
UP = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
UP["len"] = UP.seq.str.len()
UP["pred"] = UP.acc.isin(set(MF.acc))
lbins = [(1, 29, "<30"), (30, 100, "30-100"), (101, 500, "101-500"), (501, 1000, "501-1,000"),
         (1001, 1250, "1,001-1,250"), (1251, 1500, "1,251-1,500"), (1501, 10 ** 6, ">1,500")]
vals = [100 * UP[(UP.len >= a) & (UP.len <= b)].pred.mean() for a, b, _ in lbins]
ns = [int(((UP.len >= a) & (UP.len <= b)).sum()) for a, b, _ in lbins]
ax.bar(range(len(lbins)), vals, color="#4292c6", edgecolor="black", linewidth=0.4, width=0.65)
for i, (v, n_) in enumerate(zip(vals, ns)):
    ax.text(i, v + 2, f"{n_:,}", ha="center", fontsize=5.0, rotation=90, va="bottom")
ax.set_xticks(range(len(lbins)))
ax.set_xticklabels([l for _, _, l in lbins], fontsize=5.6, rotation=35, ha="right")
ax.set_ylim(0, 125)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlabel("Sequence length (residues)")
ax.set_ylabel("Predicted homodimer (%)")
ax.set_title("Coverage by length", loc="left")
letter(ax, "b", dx=0.06)

ax = fig.add_subplot(gs[0, 2])
IV = NN["robustness"]["iptm_vs_call"]["by_band"]
ibands = ["<0.5", "0.5-0.6", "0.6-0.8", ">=0.8"]
ilab = ["<0.5", "0.5-0.6", "0.6-0.8", "≥0.8"]
vals = [IV[b]["pct"] for b in ibands]
ax.bar(range(4), vals, color="#08519c", edgecolor="black", linewidth=0.4, width=0.64)
for k, b in enumerate(ibands):
    ax.text(k, vals[k] + 2, f"{IV[b]['deposited']:,}/\n{IV[b]['n']:,}", ha="center", fontsize=4.8, va="bottom")
ax.set_xticks(range(4))
ax.set_xticklabels(ilab, fontsize=5.8)
ax.set_xlabel("Homodimer ipTM")
ax.set_ylabel("Deposited (%)")
ax.set_ylim(0, 114)
ax.set_title("Deposition by ipTM", loc="left")
letter(ax, "c", dx=0.06)
save(fig, "Figure3")


T = N["transitions"]
P = pd.read_csv("results/pairs_analysis.tsv", sep="\t")
P = P[P.both_pred].copy()
CR = pd.read_csv("results/canonical_retention_analysed.tsv", sep="\t")
RET_BINS = ["lt50", "50to90", "90to100", "all"]
RET_LAB = ["<50%", "50-90%", "90-99%", "100%"]
THR = dict(color="black", lw=0.6, ls="--")


def rate_points(ax, xs, keys, src, mk, fc, label=None, ms=5.5):
    vs = np.array([src[k]["rate"] for k in keys])
    lo = np.array([src[k]["ci"][0] for k in keys])
    hi = np.array([src[k]["ci"][1] for k in keys])
    ax.errorbar(xs, vs, yerr=[vs - lo, hi - vs], fmt=mk, color="black", mfc=fc, mec="black", ms=ms,
                capsize=2.2, lw=0.7, label=label, zorder=3)
    return vs, lo, hi


fig = plt.figure(figsize=(7.0, 7.0))
gs = fig.add_gridspec(5, 3, height_ratios=[1, 0.5, 1, 0.9, 1], hspace=0, wspace=0.72, top=0.95, bottom=0.07)

ax = fig.add_subplot(gs[0, 0])
for i, g in enumerate(GROUPS):
    v, (lo, hi) = T[g]["loss_rate"], T[g]["loss_ci"]
    ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], fmt=GMARK[g], color="black", mfc=GCOL[g], mec="black",
                ms=6, capsize=2.5, lw=0.8)
    ax.text(i, hi + 3, f"{T[g]['losses']}/{T[g]['canon_hc']}", ha="center", fontsize=6.0)
ax.set_xticks(range(3))
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylabel("Isoform not high-confidence (%)")
ax.set_ylim(0, 70)
ax.set_xlim(-0.5, 2.5)
ax.set_title("Loss by isoform class", loc="left")
letter(ax, "a", dx=0.075)

ax = fig.add_subplot(gs[0, 1])
LB = N["loss_by_retention"]
for off, g in zip([-0.14, 0.0, 0.14], GROUPS):
    keys = [f"{g}|{b}" for b in RET_BINS if f"{g}|{b}" in LB]
    xs = [RET_BINS.index(k.split("|")[1]) + off for k in keys]
    vs, _, _ = rate_points(ax, xs, keys, LB, GMARK[g], GCOL[g], label=GLAB[g], ms=4.6)
    ax.plot(xs, vs, color=GCOL[g], lw=0.9, zorder=1)
for k, b in enumerate(RET_BINS):
    ax.text(k, 104, f"{LB[f'all|{b}']['n']}", ha="center", fontsize=5.6)
ax.set_xticks(range(4))
ax.set_xticklabels(RET_LAB, fontsize=6.0)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 110)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("Canonical interface kept")
ax.set_ylabel("Isoform not high-confidence (%)")
ax.legend(frameon=False, loc="lower left", fontsize=5.6, handletextpad=0.2)
ax.set_title("Loss by interface kept", loc="left", pad=10)
letter(ax, "b", dx=0.075)

ax = fig.add_subplot(gs[0, 2])
CR["ret_bin"] = pd.cut(CR.iface_retained, [0, 0.5, 0.9, 1.0, 1.0001], right=False, labels=RET_BINS)
data = [CR.loc[CR.ret_bin == b, "iso_ipsae"].values for b in RET_BINS]
bpl = ax.boxplot(data, positions=range(4), widths=0.55, showfliers=False, patch_artist=True,
                 medianprops=dict(color="black", lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for patch in bpl["boxes"]:
    patch.set_facecolor("#c6dbef")
rng = np.random.default_rng(5)
for i, d in enumerate(data):
    ax.scatter(i + rng.uniform(-0.2, 0.2, len(d)), d, s=1.6, color="black", alpha=0.3, linewidth=0, rasterized=True)
ax.axhline(0.6, **THR)
ax.set_xticks(range(4))
ax.set_xticklabels(RET_LAB, fontsize=6.0)
ax.set_xlabel("Canonical interface kept")
ax.set_ylabel("Isoform homodimer ipSAE")
ax.set_ylim(-0.03, 1.0)
ax.set_title("Isoform interface score", loc="left")
letter(ax, "c", dx=0.075)

ax = fig.add_subplot(gs[2, 0])
LI = N["loss_by_retention_ipsae"]
ibins = ["0.60to0.65", "0.65to0.70", "0.70to0.80", "0.80to1.00"]
ilab = ["0.60\nto 0.65", ">0.65\nto 0.70", ">0.70\nto 0.80", ">0.80"]
rcol = {"lt50": "#08306b", "50to90": "#2171b5", "90to100": "#6baed6", "all": "#c6dbef"}
for off, rb_ in zip([-0.15, -0.05, 0.05, 0.15], RET_BINS):
    xs = np.arange(4) + off
    vs, _, _ = rate_points(ax, xs, [f"{rb_}|{ib}" for ib in ibins], LI, "o", rcol[rb_], label=RET_LAB[RET_BINS.index(rb_)], ms=3.8)
    ax.plot(xs, vs, color=rcol[rb_], lw=0.9, zorder=1)
ax.set_xticks(range(4))
ax.set_xticklabels(ilab, fontsize=5.8)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 105)
ax.set_xlabel("Canonical homodimer ipSAE")
ax.set_ylabel("Isoform not high-confidence (%)")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.34), fontsize=5.6, ncol=4, columnspacing=0.6,
          handletextpad=0.1, title="Canonical interface kept", title_fontproperties={"size": 5.6})
ax.set_title("Loss by canonical ipSAE", loc="left")
letter(ax, "d", dx=0.075)

ax = fig.add_subplot(gs[2, 1])
RM = N["retention_models"]
rows = [("C(grp)[T.subtraction]", "Removal vs\naddition"), ("C(grp)[T.replacement]", "Replacement\nvs addition"),
        ("iface_lost10", "Interface\nremoved, /10%"), ("conf_removed10", "Confident\nresidues\nremoved, /10%")]
yy = np.arange(len(rows))[::-1]
for off, model, mk, fc in [(0.14, "group_only", "o", "white"), (-0.14, "with_interface_and_fold", "s", "#08519c")]:
    for yi, (t, _) in zip(yy, rows):
        d = RM[model]["terms"].get(t)
        if d is None:
            continue
        ax.errorbar(d["OR"], yi + off, xerr=[[d["OR"] - d["lo"]], [d["hi"] - d["OR"]]], fmt=mk, color="black",
                    mfc=fc, mec="black", ms=4.2, capsize=2, lw=0.7)
ax.axvline(1, **THR)
ax.set_xscale("log")
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
ax.set_xlim(0.45, 10)
ax.set_yticks(yy)
ax.set_yticklabels([lab for _, lab in rows], fontsize=5.8)
ax.set_ylim(-0.6, len(rows) - 0.3)
ax.set_xlabel("Odds ratio for loss")
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc="white", mec="black", ms=4.2, label="Class only"),
                   plt.Line2D([], [], marker="s", ls="", mfc="#08519c", mec="black", ms=4.2, label="+ residues removed")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.3), fontsize=5.8, ncol=2, columnspacing=0.6, handletextpad=0.1)
ax.set_title(f"Logistic models, n = {RM['group_only']['n']:,}", loc="left")
letter(ax, "e", dx=0.135)

ax = fig.add_subplot(gs[2, 2])
LT = N["loss_types"]
cats = [("iso_ipsae_ge06", "ipSAE ≥0.6\n(pDockQ2 fails)", "#c6dbef"), ("iso_ipsae_05to06", "ipSAE 0.5-0.6", "#6baed6"),
        ("iso_ipsae_lt05", "ipSAE <0.5", "#08519c")]
for i, g in enumerate(GROUPS):
    bottom = 0
    for k, lab, col in cats:
        v = 100 * LT[g][k] / LT[g]["n"]
        ax.bar(i, v, bottom=bottom, color=col, edgecolor="black", linewidth=0.4, width=0.6)
        bottom += v
    ax.text(i, 102, f"n = {LT[g]['n']}", ha="center", fontsize=5.8)
ax.set_xticks(range(3))
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylim(0, 110)
ax.set_ylabel("Losses (%)")
ax.legend(handles=[Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=l) for _, l, c in cats],
          frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=5.6, ncol=2, columnspacing=0.6,
          title="Isoform homodimer", title_fontproperties={"size": 5.6})
ax.set_title("What fails in a loss", loc="left")
letter(ax, "f", dx=0.075)

ax = fig.add_subplot(gs[4, 0])
IRB = NN["interface_removal"]["by_band"]
bands = ["<50", "50-70", "70-90", ">=90"]
blab = ["<50", "50-70", "70-90", "≥90"]
xb = np.arange(4)
for off, kind, col, lab in [(-0.19, "other", "#969696", "Elsewhere"), (0.19, "iface", "#08519c", "Interface")]:
    v = [IRB[f"{b}|{kind}"]["pct"] for b in bands]
    ax.bar(xb + off, v, width=0.36, color=col, edgecolor="black", linewidth=0.4, label=lab)
ax.set_xticks(xb)
ax.set_xticklabels(blab, fontsize=6.0)
ax.set_xlabel("Canonical residue pLDDT")
ax.set_ylabel("Removed by the isoform (%)")
ax.set_ylim(0, 44)
ax.legend(frameon=False, loc="upper right", fontsize=6.0, handletextpad=0.4)
mh = NN["interface_removal"]["mh_or"]
mlo, mhi = NN["interface_removal"]["mh_or_ci"]
ax.text(0.03, 0.96, f"MH OR {mh:.2f}\n({mlo:.2f} to {mhi:.2f})", transform=ax.transAxes, va="top", fontsize=5.8)
ax.set_title("Interface residues removed", loc="left")
letter(ax, "g", dx=0.075)

ax = fig.add_subplot(gs[4, 1])
IL = NN["monomers"]["shared_fold"]["intact_loss_by_d"]
dbins = ["0+", "-2to0", "-5to-2", "lt-5"]
dlab = ["≥0", "−2 to 0", "−5 to −2", "<−5"]
v = [IL[b]["rate"] for b in dbins]
ax.bar(range(4), v, color="#6baed6", edgecolor="black", linewidth=0.4, width=0.62)
for k, b in enumerate(dbins):
    ax.text(k, v[k] + 1, f"{IL[b]['losses']}/{IL[b]['n']}", ha="center", fontsize=5.6)
ax.set_xticks(range(4))
ax.set_xticklabels(dlab, fontsize=5.6)
ax.set_xlabel("Isoform minus canonical\nshared-fold pLDDT")
ax.set_ylabel("Loss with interface intact (%)")
ax.set_ylim(0, 42)
ax.set_title("Loss when the shared fold weakens", loc="left")
letter(ax, "h", dx=0.075)
save(fig, "Figure4")

TI = N["transitions_ipsae_only"]
fig = plt.figure(figsize=(7.0, 5.2))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.55)

ax = fig.add_subplot(gs[0, 0])
for i, g in enumerate(GROUPS):
    for off, src, kn, kd, fc in [(-0.12, T, "gains", "canon_not_hc", GCOL[g]), (0.12, TI, "gains", "canon_below", "white")]:
        v, (lo, hi) = src[g]["gain_rate"], src[g]["gain_ci"]
        ax.errorbar(i + off, v, yerr=[[v - lo], [hi - v]], fmt=GMARK[g], color="black", mfc=fc, mec="black",
                    ms=5.5, capsize=2.2, lw=0.8)
        ax.text(i + off, hi + 0.12, f"{src[g][kn]}/{src[g][kd]:,}", ha="center", fontsize=5.2, rotation=90, va="bottom")
ax.set_xticks(range(3))
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylabel("Isoform gains the call (%)")
ax.set_ylim(0, 4.6)
ax.set_xlim(-0.5, 2.5)
ax.legend(handles=[plt.Line2D([], [], marker="s", ls="", mfc="#4292c6", mec="black", ms=5, label="Deposition call"),
                   plt.Line2D([], [], marker="s", ls="", mfc="white", mec="black", ms=5, label="ipSAE ≥0.6 alone")],
          frameon=False, loc="upper left", fontsize=6.0)
ax.set_title("Gain by isoform class", loc="left")
letter(ax, "a", dx=0.07)

ax = fig.add_subplot(gs[0, 1:])
G = P[~P.canon_hc & P.iso_hc]
pd2 = G.canon_ipsae >= 0.6
ax.scatter(G.loc[~pd2, "canon_ipsae"], G.loc[~pd2, "iso_ipsae"], s=7, facecolor=ADDED, edgecolor="black", linewidth=0.2,
           label=f"Canonical below ipSAE 0.6 ({int((~pd2).sum())})")
ax.scatter(G.loc[pd2, "canon_ipsae"], G.loc[pd2, "iso_ipsae"], s=7, facecolor="#9ecae1", edgecolor="black", linewidth=0.2,
           label=f"Canonical fails pDockQ2 only ({int(pd2.sum())})")
ax.vlines(0.6, 0.55, 1.0, **THR)
ax.axhline(0.6, **THR)
ax.plot([0, 1], [0, 1], color="#969696", lw=0.5, ls=":")
ax.set_xlim(-0.02, 1.0)
ax.set_ylim(0.55, 1.1)
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("Canonical homodimer ipSAE")
ax.set_ylabel("Isoform homodimer ipSAE")
ax.legend(frameon=False, loc="upper left", fontsize=5.8, handletextpad=0.2, borderaxespad=0.2)
ax.set_title(f"The {len(G)} gains", loc="left")
letter(ax, "b", dx=0.07)

ax = fig.add_subplot(gs[1, 0])
GBd = N["gain_by_dplddt"]
dbins = ["lt-5", "-5to0", "0to5", "5to10", "gt10"]
dlab = ["≤−5", ">−5\nto 0", ">0\nto 5", ">5\nto 10", ">10"]
tot = np.array([GBd[f"all|{b}"]["rate"] for b in dbins])
nn = np.array([GBd[f"all|{b}"]["n"] for b in dbins])
pdq = np.array([100 * GBd[f"all|{b}"]["gains_pdockq2_type"] / GBd[f"all|{b}"]["n"] for b in dbins])
ax.bar(range(5), tot - pdq, color=ADDED, edgecolor="black", linewidth=0.4, width=0.6, label="Canonical below ipSAE 0.6")
ax.bar(range(5), pdq, bottom=tot - pdq, color="#9ecae1", edgecolor="black", linewidth=0.4, width=0.6, label="Canonical fails pDockQ2 only")
for k, b in enumerate(dbins):
    ax.text(k, tot[k] + 0.15, f"{GBd[f'all|{b}']['gains']}/{nn[k]:,}", ha="center", fontsize=5.2, rotation=90, va="bottom")
ax.set_xticks(range(5))
ax.set_xticklabels(dlab, fontsize=5.8)
ax.set_ylim(0, 8.6)
ax.set_xlabel("Isoform minus canonical monomer pLDDT")
ax.set_ylabel("Isoform gains the call (%)")
ax.legend(frameon=False, loc="upper left", fontsize=5.8)
ax.set_title("Gain by monomer confidence", loc="left")
letter(ax, "c", dx=0.07)

ax = fig.add_subplot(gs[1, 1])
RG, RL = N["robust_gains"], N["robust_losses"]
x = np.arange(3)
for off, src, col, lab in [(-0.19, RL, "#08519c", "Fall: canonical ≥0.8 to isoform ≤0.4"),
                           (0.19, RG, ADDED, "Rise: canonical ≤0.4 to isoform ≥0.8")]:
    vals = [100 * src[g]["k"] / src[g]["n"] for g in GROUPS]
    ax.bar(x + off, vals, width=0.36, color=col, edgecolor="black", linewidth=0.4, label=lab)
    for xi, g, v in zip(x + off, GROUPS, vals):
        ax.text(xi, max(v, 0) + 0.8, f"{src[g]['k']}/{src[g]['n']:,}", ha="center", fontsize=5.2, rotation=90, va="bottom")
ax.set_xticks(x)
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylabel("Pairs (%)")
ax.set_ylim(0, 70)
ax.set_yticks([0, 10, 20, 30, 40, 50, 60])
ax.legend(frameon=False, loc="upper left", fontsize=5.8)
ax.set_title("Large changes in ipSAE", loc="left")
letter(ax, "d", dx=0.07)

ax = fig.add_subplot(gs[1, 2])
GCA = NN["monomers"]["gains_by_confident_added"]
cbins = ["none added", "0 confident", "1-9 confident", "10+ confident"]
clab = ["None\nadded", "0", "1-9", "≥10"]
v = [GCA[b]["rate"] for b in cbins]
bcols = ["#d9d9d9", "#fdd0a2", "#fd8d3c", ADDED]
ax.bar(range(4), v, color=bcols, edgecolor="black", linewidth=0.4, width=0.64)
for k, b in enumerate(cbins):
    ax.text(k, v[k] + 0.06, f"{GCA[b]['gains']}/{GCA[b]['n']:,}", ha="center", fontsize=5.0, rotation=90, va="bottom")
ax.set_xticks(range(4))
ax.set_xticklabels(clab, fontsize=5.6)
ax.set_xlabel("Confident added residues\nin the isoform monomer")
ax.set_ylabel("Isoform gains the call (%)")
ax.set_ylim(0, 3.3)
ax.set_title("Gain by confident added sequence", loc="left")
letter(ax, "e", dx=0.07)
save(fig, "Figure5")

R = pd.read_csv("results/residues_novel_models.tsv.gz", sep="\t")
cl_idx = cl.set_index("isoform_acc")
fig = plt.figure(figsize=(7.0, 5.0))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.95], hspace=0.62)
top = outer[0].subgridspec(1, 3, wspace=0.75)
bot = outer[1].subgridspec(1, 4, wspace=0.55)

ax = fig.add_subplot(top[0, 0])
hb = np.linspace(0, 100, 41)
ax.hist(R.loc[~R.novel, "plddt"], bins=hb, density=True, histtype="step", color=SHARED, lw=1.1, label="Shared")
ax.hist(R.loc[R.novel, "plddt"], bins=hb, density=True, histtype="stepfilled", color=ADDED, alpha=0.55, lw=0, label="Added")
ax.set_ylim(0, 0.092)
ax.vlines([50, 70], 0, 0.07, color="black", lw=0.6, ls=":")
ax.set_xlabel("pLDDT, isoform homodimer")
ax.set_ylabel("Density")
ax.legend(frameon=False, loc="upper left", ncol=2, columnspacing=0.8, handlelength=1.2)
ax.set_title(f"{N['added_residues']:,} added residues", loc="left")
letter(ax, "a", dx=0.085)

ax = fig.add_subplot(top[0, 1])
band = pd.read_csv("results/table_interface_by_plddt.tsv", sep="\t")
order = ["lt50", "50to70", "70to90", "ge90"]
xb = np.arange(4)
for off, nv, col, hatch, lab in [(-0.19, False, SHARED, "", "Shared"), (0.19, True, ADDED, "////", "Added")]:
    vals = [float(band[(band.pbin == b) & (band.novel == nv)].iface_rate.iloc[0]) for b in order]
    ax.bar(xb + off, vals, width=0.36, color=col if not nv else "white", edgecolor=col, hatch=hatch, linewidth=0.8, label=lab)
ax.set_xticks(xb)
ax.set_xticklabels(["<50", "50-70", "70-90", "≥90"])
ax.set_xlabel("pLDDT band")
ax.set_ylabel("Exposed residues at\nthe interface (%)")
ax.set_ylim(0, 45)
ax.legend(frameon=False, loc="upper left", fontsize=6.0)
o = N["residue_logit_plddt_adjusted"]
ax.text(0.03, 0.70, f"Adjusted OR {o['OR']:.2f}\n({o['lo']:.2f} to {o['hi']:.2f})", transform=ax.transAxes, ha="left", va="top", fontsize=5.8)
ax.set_title("Matched confidence", loc="left")
letter(ax, "b", dx=0.085)

ax = fig.add_subplot(top[0, 2])
steps = [("Gains, isoform\ncarries added\nresidues", N["gains_with_added_models"]),
         ("Added residue\nat the interface", N["gains_added_at_iface"]),
         ("Contacts between\nadded residues", N["gains_added_novel_novel_contacts"]),
         ("≥5 confident\nadded residues\nat interface", N["gains_added_ge5_conf_iface"])]
yy = np.arange(len(steps))[::-1]
ax.barh(yy, [v for _, v in steps], color=["#d9d9d9", "#fdae6b", "#fd8d3c", ADDED], edgecolor="black", linewidth=0.4, height=0.6)
for yi, (_, v) in zip(yy, steps):
    ax.text(v + 4, yi, str(v), va="center", fontsize=6.0)
ax.set_yticks(yy)
ax.set_yticklabels([s_ for s_, _ in steps], fontsize=5.6)
ax.set_xlabel("Pairs")
ax.set_xlim(0, 280)
ax.set_title("Added residues in gains", loc="left")
letter(ax, "c", dx=0.16)

cands = [("Q9H579-2", "MROH8 isoform 2"), ("O75638-2", "CTAG2 isoform 2"), ("Q9Y4E5-4", "ZNF451 isoform 4"), ("Q96RF0-2", "SNX18 isoform 2")]
for k, (acc, name) in enumerate(cands):
    ax = fig.add_subplot(bot[0, k])
    d = pd.read_csv(f"results/struct/{acc}.tsv.gz", sep="\t", keep_default_na=False, dtype={"partners_A": str})
    segs = [tuple(int(v) for v in s_.split(":")[0].split("-")) for s_ in str(cl_idx.loc[acc, "novel_segments"]).split(";") if s_]
    pts = np.array([(i + 1, int(j)) for i, s_ in enumerate(d.partners_A) for j in str(s_).split(";") if j])
    isadd = np.zeros(len(d) + 1, bool)
    for a0, a1 in segs:
        isadd[a0:a1 + 1] = True
    added = isadd[pts[:, 0]] | isadd[pts[:, 1]]
    ax.scatter(pts[~added, 0], pts[~added, 1], s=1.5, color=SHARED, linewidth=0, rasterized=True)
    ax.scatter(pts[added, 0], pts[added, 1], s=1.5, color=ADDED, linewidth=0, rasterized=True)
    for a0, a1 in segs:
        ax.axvspan(a0, a1, color=ADDED, alpha=0.12, lw=0)
        ax.axhspan(a0, a1, color=ADDED, alpha=0.12, lw=0)
    lo_, hi_ = max(1, pts.min() - 15), min(len(d), pts.max() + 10)
    ax.set_xlim(lo_, hi_)
    ax.set_ylim(lo_, hi_)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=5.6)
    ax.set_xlabel("Chain A residue", fontsize=6.0)
    if k == 0:
        ax.set_ylabel("Chain B residue", fontsize=6.0)
    row = [r_ for r_ in N["gains_ge5_conf_list"] if r_["isoform_acc"] == acc][0]
    ax.set_title(f"{name}\nipSAE {row['canon_ipsae']:.2f} to {row['iso_ipsae']:.2f}", loc="left", fontsize=6.3)
    letter(ax, "defg"[k], dx=0.045)
save(fig, "Figure6")


fig = plt.figure(figsize=(7.0, 5.2))
gs = fig.add_gridspec(2, 3, hspace=0.72, wspace=0.60)
D = pd.read_csv("results/identical_copies.tsv", sep="\t")
dd = D[D.copy_in_manifest & D.human_ipsae.notna()]
hd = dd.human_deposited

ax = fig.add_subplot(gs[0, 0])
ax.scatter(dd.loc[~hd, "human_ipsae"], dd.loc[~hd, "copy_ipsae"], s=5, facecolor="#d9d9d9", edgecolor="black",
           linewidth=0.2, label="Human not deposited")
ax.scatter(dd.loc[hd, "human_ipsae"], dd.loc[hd, "copy_ipsae"], s=5, facecolor="#4292c6", edgecolor="black",
           linewidth=0.2, label="Human deposited")
ax.vlines(0.6, -0.02, 0.78, **THR)
ax.axhline(0.6, **THR)
ax.set_xlim(-0.02, 1.0)
ax.set_ylim(-0.02, 1.0)
ax.set_xlabel("ipSAE, human entry")
ax.set_ylabel("ipSAE, identical sequence")
ax.legend(frameon=False, loc="upper left", fontsize=5.6, handletextpad=0.1)
ax.set_title(f"{len(dd)} identical-sequence pairs", loc="left")
letter(ax, "a", dx=0.1)

ax = fig.add_subplot(gs[0, 1])
ax.scatter(dd.loc[~hd, "human_pdockq2"], dd.loc[~hd, "copy_pdockq2"], s=5, facecolor="#d9d9d9", edgecolor="black", linewidth=0.2)
ax.scatter(dd.loc[hd, "human_pdockq2"], dd.loc[hd, "copy_pdockq2"], s=5, facecolor="#4292c6", edgecolor="black", linewidth=0.2)
ax.axvline(0.23, **THR)
ax.axhline(0.23, **THR)
ax.plot([0, 1], [0, 1], color="#969696", lw=0.5, ls=":")
ax.set_xlim(-0.02, 1.0)
ax.set_ylim(-0.02, 1.0)
ax.set_xlabel("pDockQ2, human entry")
ax.set_ylabel("pDockQ2, identical sequence")
ax.set_title("pDockQ2 of the same pairs", loc="left")
letter(ax, "b", dx=0.1)

ax = fig.add_subplot(gs[0, 2])
cal = [("Homodimer", N["calibration_annotated_homodimer"]), ("Other homo-oligomer", N["calibration_annotated_other_homooligomer"]),
       ("Monomer", N["calibration_annotated_monomer"]), ("Other subunit text", N["calibration_text_other"]),
       ("No subunit text", N["calibration_no_subunit_text"])]
for i, (lab, d_) in enumerate(cal):
    v, (lo, hi) = d_["rate"], d_["ci"]
    ax.bar(i, v, color="#4292c6" if i < 2 else "#d9d9d9", edgecolor="black", linewidth=0.4, width=0.62)
    ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="black", lw=0.7, capsize=2)
    ax.text(i, hi + 1.2, f"{d_['hc']}/{d_['n']:,}", ha="center", fontsize=5.0, rotation=90, va="bottom")
ax.set_xticks(range(len(cal)))
ax.set_xticklabels([l for l, _ in cal], fontsize=5.6, rotation=35, ha="right")
ax.set_ylabel("Canonical high-confidence (%)")
ax.set_ylim(0, 50)
ax.set_title("UniProt SUBUNIT text", loc="left")
letter(ax, "c", dx=0.1)

ax = fig.add_subplot(gs[1, 0])
mech = pd.read_csv("results/pair_mechanism.tsv", sep="\t")
in_p = np.array([k in set(zip(P.canonical_acc, P.isoform_acc)) for k in zip(mech.canonical_acc, mech.isoform_acc)], dtype=bool)
hh = mech[mech.canon_hc.values & in_p & mech.iface_jaccard_shared.notna().values].copy()
hh["grp"] = hh["class"].map(lambda c: "subtraction" if c in ("truncation", "internal deletion") else c)
data = [hh.loc[hh.grp == g, "iface_jaccard_shared"].values for g in GROUPS]
bpj = ax.boxplot(data, positions=range(3), widths=0.5, showfliers=False, patch_artist=True,
                 medianprops=dict(color="black", lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
for patch, g in zip(bpj["boxes"], GROUPS):
    patch.set_facecolor(GCOL[g])
rng = np.random.default_rng(3)
for i, d_ in enumerate(data):
    ax.scatter(i + rng.uniform(-0.18, 0.18, len(d_)), d_, s=1.8, color="black", alpha=0.35, linewidth=0, rasterized=True)
    ax.text(i, -0.1, f"{len(d_)}", ha="center", fontsize=5.8)
ax.set_xticks(range(3))
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylim(-0.16, 1.05)
ax.set_ylabel("Interface overlap (Jaccard)")
ax.set_title("Both calls high-confidence", loc="left")
letter(ax, "d", dx=0.1)

callseq = [("deposition", "Deposition", "#08519c", "o"), ("ipsae_only", "ipSAE ≥0.6", "#4292c6", "s"),
           ("pdockq2_only", "pDockQ2 ≥0.23", "#6baed6", "^"), ("iptm_ge08", "ipTM ≥0.8", "#fd8d3c", "D"),
           ("iptm_ge06", "ipTM ≥0.6", "#fdae6b", "v")]
RB = NN["robustness"]

ax = fig.add_subplot(gs[1, 1])
rbins = ["all", "90to100", "50to90", "lt50"]
rlab = ["100%", "90-99%", "50-90%", "<50%"]
for name, lab, col, mk in callseq:
    lb = RB[name]["loss_by_retention"]
    v = [lb[b]["rate"] for b in rbins]
    ax.plot(range(4), v, marker=mk, color=col, ms=3.4, lw=0.9, label=lab)
ax.set_xticks(range(4))
ax.set_xticklabels(rlab, fontsize=5.8)
ax.set_xlabel("Canonical interface kept")
ax.set_ylabel("Isoform fails the call (%)")
ax.set_ylim(0, 105)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.32), fontsize=5.2, ncol=3, columnspacing=0.6, handletextpad=0.2)
ax.set_title("Loss by interface kept, by call", loc="left")
letter(ax, "e", dx=0.1)

ax = fig.add_subplot(gs[1, 2])
for name, lab, col, mk in callseq:
    v = [RB[name][g]["loss_rate"] for g in GROUPS]
    ax.plot(range(3), v, marker=mk, color=col, ms=3.4, lw=0.9, label=lab)
ax.set_xticks(range(3))
ax.set_xticklabels(["Addition\nonly", "Removal\nonly", "Replace-\nment"])
ax.set_ylabel("Loss rate (%)")
ax.set_ylim(0, 70)
ax.set_title("Loss by class, by call", loc="left")
letter(ax, "f", dx=0.1)
save(fig, "FigureS1")

