import gzip
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

OUT = "data/afdb_monomer_plddt.jsonl.gz"
LOG = "data/afdb_monomer_plddt_log.json"
URL = "https://alphafold.ebi.ac.uk/files/{entry}-confidence_v{ver}.json"

seq = pd.read_csv("data/uniprot_human_sp_sequences.tsv.gz", sep="\t")
SEQ = dict(zip(seq.acc, seq.seq))
mono = pd.read_csv("data/afdb_monomer_per_accession.tsv", sep="\t").set_index("acc")
cl = pd.read_csv("results/isoform_pairs_classified.tsv", sep="\t")
cl = cl[cl.in_afdb_length_range == True]
targets = sorted(set(cl.isoform_acc) | set(cl.canonical_acc))

todo, skipped = [], {"no_model": [], "sequence_differs": []}
for a in targets:
    if a not in mono.index or pd.isna(mono.loc[a, "mono_entry"]):
        skipped["no_model"].append(a)
        continue
    if hashlib.md5(SEQ[a].encode()).hexdigest() != mono.loc[a, "mono_seq_checksum"]:
        skipped["sequence_differs"].append(a)
        continue
    todo.append((a, mono.loc[a, "mono_entry"], int(mono.loc[a, "mono_version"])))

done = set()
if os.path.exists(OUT):
    with gzip.open(OUT, "rt") as fh:
        for line in fh:
            done.add(json.loads(line)["acc"])
todo = [t for t in todo if t[0] not in done]
print(f"targets {len(targets):,}; already done {len(done):,}; to fetch {len(todo):,}; "
      f"no model {len(skipped['no_model'])}; sequence differs {len(skipped['sequence_differs'])}", flush=True)

S = requests.Session()


def fetch(item):
    acc, entry, ver = item
    for attempt in range(4):
        try:
            r = S.get(URL.format(entry=entry, ver=ver), timeout=60)
            if r.status_code == 200:
                d = r.json()
                p = [round(float(x), 2) for x in d["confidenceScore"]]
                if len(p) != len(SEQ[acc]):
                    return acc, ver, None, f"length {len(p)} != {len(SEQ[acc])}"
                return acc, ver, p, None
            if r.status_code == 404:
                return acc, ver, None, "404"
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))
    return acc, ver, None, "failed after retries"


failed = {}
n_ok = 0
with gzip.open(OUT, "at") as fh, ThreadPoolExecutor(max_workers=8) as ex:
    futs = [ex.submit(fetch, t) for t in todo]
    for i, f in enumerate(as_completed(futs), 1):
        acc, ver, p, err = f.result()
        if p is None:
            failed[acc] = err
        else:
            fh.write(json.dumps({"acc": acc, "version": ver, "plddt": p}) + "\n")
            n_ok += 1
        if i % 1000 == 0:
            fh.flush()
            print(f"{i:,}/{len(todo):,} fetched, {len(failed)} failed", flush=True)

json.dump({"targets": len(targets), "fetched_now": n_ok, "already_done": len(done), "failed": failed,
           "skipped": {k: v for k, v in skipped.items()}}, open(LOG, "w"), indent=1)
if len(failed) > 0.02 * max(1, len(todo)):
    sys.exit(f"STOP: {len(failed)} of {len(todo)} requests failed; see {LOG}")
print(f"done: {n_ok:,} fetched now, {len(failed)} failed", flush=True)
