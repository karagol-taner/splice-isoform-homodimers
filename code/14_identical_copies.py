import json, os, re, sys, time, random, urllib.request, urllib.error, urllib.parse
from concurrent.futures import ThreadPoolExecutor
LIST = sys.argv[1]; NW = int(sys.argv[2]) if len(sys.argv) > 2 else 4
accs = [a.strip() for a in open(LIST) if a.strip()]
UA = {"User-Agent": "isoform-homodimer-study (karagolt@mit.edu)"}

def get_json(url, tries=8):
    delay = 1.0
    for _ in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return 200, json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504): time.sleep(delay + random.random()); delay = min(2 * delay, 60); continue
            return e.code, None
        except Exception:
            time.sleep(delay + random.random()); delay = min(2 * delay, 60)
    return -1, None

def done_set(fn):
    s = set()
    if os.path.exists(fn):
        for line in open(fn):
            try: s.add(json.loads(line)["acc"])
            except Exception: pass
    return s

U = "data/uniparc_raw.jsonl"; done = done_set(U)
def uniparc(acc):
    q = urllib.parse.quote(f"isoform:{acc}")
    st, d = get_json(f"https://rest.uniprot.org/uniparc/search?query={q}&fields=upi,accession&format=json&size=5")
    copies = []
    base = acc.split("-")[0]
    if st == 200 and d:
        for res in d.get("results", []):
            for x in res.get("uniProtKBAccessions", []) or []:
                if re.search(r"\.\d+$", x):
                    continue
                if x.split("-")[0] == base:
                    continue
                copies.append({"id": x, "upi": res.get("uniParcId")})
    return {"acc": acc, "status": st, "copies": copies}
todo = [a for a in accs if a not in done]
print("uniparc to do:", len(todo), flush=True)
with open(U, "a") as fo, ThreadPoolExecutor(NW) as ex:
    for rec in ex.map(uniparc, todo):
        fo.write(json.dumps(rec) + "\n"); fo.flush()

copies = {}
for line in open(U):
    r = json.loads(line)
    if r["acc"] in set(accs):
        for c in r["copies"]:
            copies.setdefault(c["id"], r["acc"])
C = "data/copies_complex_raw.jsonl"; done = done_set(C)
todo = [c for c in copies if c not in done]
print("copies to query:", len(todo), flush=True)
def cx(acc):
    _, o = get_json(f"https://rest.uniprot.org/uniprotkb/{acc}?fields=organism_id,reviewed,sequence&format=json")
    taxid = ((o or {}).get("organism") or {}).get("taxonId") if o else None
    seqlen = ((o or {}).get("sequence") or {}).get("length") if o else None
    st, d = get_json(f"https://alphafold.ebi.ac.uk/api/complex/{acc}")
    homo = []
    if st == 200 and d:
        for e in d:
            if e.get("assemblyType") == "Homo" and [c.get("identifier") for c in e.get("complexComposition", [])] == [acc]:
                homo.append({k: e.get(k) for k in ("modelEntityId", "organismScientificName", "globalMetricValue",
                             "complexPredictionAccuracy_ipTM", "complexPredictionAccuracy_ipSAE", "sequenceEnd")})
    return {"acc": acc, "human_source": copies[acc], "taxid": taxid, "copy_len": seqlen, "status": st, "homo": homo}
with open(C, "a") as fo, ThreadPoolExecutor(NW) as ex:
    for rec in ex.map(cx, todo):
        fo.write(json.dumps(rec) + "\n"); fo.flush()
print("done", flush=True)
