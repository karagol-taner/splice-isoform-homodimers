import json
import os
import sys
import time
import urllib.request
from difflib import SequenceMatcher

import pandas as pd

RAW = "data/candidate_evidence"
OUT = "results/candidate_isoform_evidence.json"
os.makedirs(RAW, exist_ok=True)
N = json.load(open("results/numbers.json"))
CAND = [r["isoform_acc"] for r in N["gains_ge5_conf_list"]]
pa = pd.read_csv("results/pairs_analysis.tsv", sep="\t").set_index("isoform_acc")
seqt = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
SEQ = dict(zip(seqt.acc, seqt.seq))
ALL = list(zip(seqt.acc, seqt.seq))


def fetch(url, path, accept=None, tries=6):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return open(path).read()
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept} if accept else {})
            txt = urllib.request.urlopen(req, timeout=120).read().decode()
            open(path, "w").write(txt)
            return txt
        except Exception as e:
            time.sleep(5 * (k + 1))
            err = e
    sys.exit(f"STOP: {url}: {err}")


def positions(segs):
    out = set()
    for s_ in str(segs).split(";"):
        if s_ and s_ != "nan":
            rng, _ = s_.split(":")
            a, b = (int(v) for v in rng.split("-"))
            out.update(range(a, b + 1))
    return out


def covered(query, target):
    sm = SequenceMatcher(None, query, target, autojunk=False)
    return [(b.a, b.size) for b in sm.get_matching_blocks() if b.size >= 5]


out = {}
for acc in CAND:
    base = acc.split("-")[0]
    r = pa.loc[acc]
    can_acc = r.canonical_acc
    iso, can = SEQ[acc], SEQ[can_acc]
    added = positions(r.novel_segments)
    rec = {"gene": None, "isoform_length": len(iso), "canonical_length": len(can), "added_residues": len(added)}

    txt = fetch(f"https://rest.uniprot.org/uniprotkb/{base}.txt", f"{RAW}/uniprot_{base}.txt")
    rec["gene"] = next(l.split("Name=")[1].split(";")[0].split(" ")[0] for l in txt.splitlines() if l.startswith("GN   Name="))
    rec["release"] = next(l[5:].strip() for l in txt.splitlines() if l.startswith("DT") and "entry version" in l)
    xr = {"Ensembl": [], "RefSeq": [], "CCDS": [], "MANE-Select": []}
    for l in txt.splitlines():
        if l.startswith("DR   "):
            db = l[5:].split(";")[0]
            if db in xr:
                tag = l.split("[")[1].split("]")[0] if "[" in l else None
                ids = [x.strip().split(" [")[0].rstrip(".") for x in l[5:].split(";")[1:]]
                xr[db].append({"ids": [i for i in ids if i and i != "-"], "isoform_tag": tag})
    rec["uniprot_xrefs"] = xr
    rec["isoform_xrefs"] = {db: [x["ids"] for x in v if x["isoform_tag"] == acc] for db, v in xr.items()}
    rec["mane_select"] = ("isoform" if any(x["isoform_tag"] == acc for x in xr["MANE-Select"]) else
                          "canonical" if any(x["isoform_tag"] in (None, can_acc + "-1") for x in xr["MANE-Select"]) else
                          "other isoform" if xr["MANE-Select"] else "none")

    pj = json.loads(fetch(f"https://www.ebi.ac.uk/proteins/api/proteomics?accession={acc}", f"{RAW}/proteomics_{acc}.json",
                          accept="application/json"))
    e = pj if isinstance(pj, dict) else pj[0]
    if e["sequence"] != iso:
        sys.exit(f"STOP: {acc}: proteomics mapping is to a different sequence")
    spec = []
    for f in e.get("features", []):
        b, en, pep = int(f["begin"]), int(f["end"]), f["peptide"]
        if iso[b - 1:en] != pep:
            sys.exit(f"STOP: {acc}: peptide {pep} not at {b}-{en}")
        if not (added & set(range(b, en + 1))) or pep in can:
            continue
        genes = sorted({a.split("-")[0] for a, s in ALL if pep in s})
        srcs = [ev["source"] for ev in f["evidences"] if "source" in ev]
        pa_ev = next((x for x in srcs if x.get("name") == "PeptideAtlas"), {})
        spec.append({"begin": b, "end": en, "peptide": pep, "added_overlap": len(added & set(range(b, en + 1))),
                     "sources": sorted({x["name"] for x in srcs}), "samples": pa_ev.get("properties", {}).get("Number of sample"),
                     "observations": pa_ev.get("properties", {}).get("Times of observations"),
                     "other_genes": [g for g in genes if g != base]})
    gene_spec = [s for s in spec if not s["other_genes"]]
    n_samp = [int(s["samples"]) for s in gene_spec if s["samples"]]
    rec["peptides_mapped"] = len(e.get("features", []))
    rec["isoform_specific_peptides"] = len(gene_spec)
    rec["isoform_specific_peptides_within_added"] = sum(1 for s in gene_spec if s["added_overlap"] == s["end"] - s["begin"] + 1)
    rec["isoform_specific_max_samples"] = max(n_samp) if n_samp else None
    rec["isoform_specific_peptides_also_in_other_genes"] = len(spec) - len(gene_spec)
    rec["peptide_list"] = spec

    prots = {}
    lk = json.loads(fetch(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{rec['gene']}?expand=1", f"{RAW}/ensembl_{rec['gene']}.json",
                          accept="application/json"))
    for t in lk.get("Transcript", []):
        tl = t.get("Translation")
        if tl:
            sq = json.loads(fetch(f"https://rest.ensembl.org/sequence/id/{tl['id']}?type=protein", f"{RAW}/ensembl_{tl['id']}.json",
                                  accept="application/json"))["seq"]
            prots[tl["id"]] = {"seq": sq, "transcript": t["id"], "ensembl_canonical": bool(t.get("is_canonical"))}
    refseq = sorted({i for x in xr["RefSeq"] for i in x["ids"] if i.startswith("NP_")})
    if refseq:
        fa = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=" + ",".join(refseq) + "&rettype=fasta&retmode=text",
                   f"{RAW}/refseq_{base}.fasta")
        for block in fa.split(">")[1:]:
            h, *s = block.split("\n")
            prots[h.split()[0]] = {"seq": "".join(s), "transcript": None, "ensembl_canonical": False}
    rows = []
    for pid, p in prots.items():
        bl = covered(iso, p["seq"])
        cov_iso = sum(n for _, n in bl)
        cov_add = sum(len(added & set(range(a + 1, a + n + 1))) for a, n in bl)
        cov_can = sum(n for _, n in covered(can, p["seq"]))
        rows.append({"protein": pid, "length": len(p["seq"]), "transcript": p["transcript"], "ensembl_canonical": p["ensembl_canonical"],
                     "identical_to_isoform": p["seq"] == iso, "identical_to_canonical": p["seq"] == can,
                     "isoform_residues_matched": cov_iso, "added_residues_matched": cov_add, "canonical_residues_matched": cov_can})
    rec["gene_proteins"] = rows
    out[acc] = rec


def checked(rec):
    peps = sorted((p["begin"], p["end"], p["peptide"]) for p in rec["peptide_list"]
                  if "PeptideAtlas" in p["sources"] and not p["other_genes"])
    prots = sorted((x["protein"], x["length"], x["identical_to_isoform"], x["identical_to_canonical"], x["added_residues_matched"])
                   for x in rec["gene_proteins"])
    return {"mane_select": rec["mane_select"], "peptides": peps, "proteins": prots}


if os.path.exists(OUT):
    dep = json.load(open(OUT))
    if set(dep) != set(out):
        sys.exit(f"STOP: candidates differ from {OUT}")
    diff = [f"{a} {k}" for a in out for k, v in checked(out[a]).items() if checked(dep[a])[k] != v]
    if diff:
        sys.exit(f"STOP: fetched evidence differs from {OUT}: " + ", ".join(diff))
    print("fetched evidence matches", OUT)
else:
    json.dump(out, open(OUT, "w"), indent=1)
for acc, rec in out.items():
    print(acc, rec["gene"], "MANE Select:", rec["mane_select"], "| specific peptides", rec["isoform_specific_peptides"],
          "(max samples", rec["isoform_specific_max_samples"], ") | proteins containing the added residues:",
          [(x["protein"], x["length"], x["added_residues_matched"]) for x in rec["gene_proteins"] if x["added_residues_matched"] >= 0.9 * rec["added_residues"]])
