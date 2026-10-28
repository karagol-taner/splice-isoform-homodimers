import json
import os
import random
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

UA = {"User-Agent": "isoform-homodimer-study (karagolt@mit.edu)"}
OUT = "data/copies_served_entries.jsonl"
KEEP = ("modelEntityId", "latestVersion", "globalMetricValue", "complexPredictionAccuracy_ipTM",
        "complexPredictionAccuracy_ipSAE", "complexPredictionAccuracy_pDockQ", "complexPredictionAccuracy_pDockQ2",
        "complexPredictionAccuracy_LIS", "organismScientificName", "taxId", "sequenceStart", "sequenceEnd")


def get_json(url, tries=8):
    delay = 1.0
    for _ in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return 200, json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(delay + random.random()); delay = min(2 * delay, 60); continue
            return e.code, None
        except Exception:
            time.sleep(delay + random.random()); delay = min(2 * delay, 60)
    return -1, None


served = [json.loads(l) for l in open("data/copies_complex_raw.jsonl")]
served = [c for c in served if c["homo"]]
done = set()
if os.path.exists(OUT):
    done = {json.loads(l)["acc"] for l in open(OUT)}
todo = [c for c in served if c["acc"] not in done]
print("served copies:", len(served), "| to fetch:", len(todo), flush=True)


def fetch(c):
    acc = c["acc"]
    st, d = get_json(f"https://alphafold.ebi.ac.uk/api/complex/{acc}")
    homo = []
    if st == 200 and d:
        for e in d:
            comp = e.get("complexComposition", [])
            if (e.get("assemblyType") == "Homo" and e.get("oligomericState") == "dimer"
                    and [x.get("identifier") for x in comp] == [acc] and [x.get("stoichiometry") for x in comp] == [2]):
                homo.append({k: e.get(k) for k in KEEP})
    return {"acc": acc, "human_source": c["human_source"], "status": st, "homo": homo}


with open(OUT, "a") as fo, ThreadPoolExecutor(4) as ex:
    for rec in ex.map(fetch, todo):
        fo.write(json.dumps(rec) + "\n"); fo.flush()
recs = [json.loads(l) for l in open(OUT)]
bad = [r["acc"] for r in recs if r["status"] != 200 or not r["homo"]]
print("records:", len(recs), "| without a served homodimer now:", len(bad), bad[:10])
