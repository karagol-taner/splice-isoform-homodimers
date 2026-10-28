import gzip, re, collections
import pandas as pd

FASTA = "data/human_sp_isoforms.fasta.gz"
TSV = "data/human_sp_altprod.tsv.gz"

def read_fasta(fn):
    seqs, head, buf = {}, None, []
    with gzip.open(fn, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if head: seqs[head[0]] = (head[1], "".join(buf))
                acc = line.split("|")[1]
                seqs_hdr = line[1:]
                head, buf = (acc, seqs_hdr), []
            else:
                buf.append(line.strip())
        if head: seqs[head[0]] = (head[1], "".join(buf))
    return seqs

seqs = read_fasta(FASTA)
rows = []
for acc, (hdr, s) in seqs.items():
    base = re.sub(r"-\d+$", "", acc)
    gn = re.search(r" GN=(\S+)", hdr)
    rows.append({"acc": acc, "base": base, "is_isoform_record": acc != base, "gene": gn.group(1) if gn else "",
                 "length": len(s), "seq": s})
df = pd.DataFrame(rows)
tsv = pd.read_csv(TSV, sep="\t", dtype=str).fillna("")
tsv.columns = ["base", "entry_name", "gene_primary", "length_entry", "alt_products", "var_seq", "protein_existence"]

def parse_alt(s):
    events = re.findall(r"Event=([^;]+);", s)
    disp = None
    for m in re.finditer(r"IsoId=([^;]+); Sequence=([^;]+);", s):
        ids = [x.strip() for x in m.group(1).split(",")]
        if m.group(2).strip() == "Displayed":
            disp = ids[0]
    return events, disp
ev = tsv["alt_products"].apply(parse_alt)
tsv["events"] = ev.apply(lambda x: "|".join(x[0]))
tsv["displayed_isoid"] = ev.apply(lambda x: x[1] or "")
df = df.merge(tsv[["base", "entry_name", "protein_existence", "events", "displayed_isoid", "var_seq"]], on="base", how="left")
df.drop(columns=["var_seq"]).to_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t", index=False)

n_entries = df["base"].nunique()
with_iso = df.loc[df.is_isoform_record, "base"].unique()
print("sequences:", len(df), "| entries:", n_entries, "| entries with >=1 isoform record:", len(with_iso),
      "| isoform records:", int(df.is_isoform_record.sum()))
print("displayed isoform not -1 (entries with isoforms):",
      int(((tsv.base.isin(with_iso)) & (tsv.displayed_isoid != "") & (~tsv.displayed_isoid.str.endswith("-1"))).sum()))
print("event types among entries with isoforms:")
c = collections.Counter()
for e in tsv.loc[tsv.base.isin(with_iso), "events"]:
    for x in set(e.split("|")):
        if x: c[x.strip()] += 1
for k, v in c.most_common(): print("   ", k, v)
query = df.loc[df.base.isin(with_iso), "acc"].tolist()
open("data/accessions_to_query.txt", "w").write("\n".join(query) + "\n")
print("accessions to query (bases with isoforms + their isoforms):", len(query))
