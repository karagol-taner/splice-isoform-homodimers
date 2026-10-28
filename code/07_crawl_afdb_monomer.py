import json, os, sys, time, random, threading
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

OUT = "data/afdb_monomer_raw.jsonl"
API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else None

accs = [a.strip() for a in open("data/accessions_monomer.txt") if a.strip()]
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try: done.add(json.loads(line)["acc"])
        except Exception: pass
todo = [a for a in accs if a not in done]
if LIMIT: todo = todo[:LIMIT]
print(f"total {len(accs)}, done {len(done)}, to do {len(todo)}, workers {N_WORKERS}", flush=True)

lock = threading.Lock()
def fetch(acc):
    delay = 1.0
    for _ in range(8):
        try:
            req = urllib.request.Request(API.format(acc), headers={"User-Agent": "isoform-homodimer-study (karagolt@mit.edu)"})
            with urllib.request.urlopen(req, timeout=40) as r:
                body = r.read().decode("utf8")
            return {"acc": acc, "status": 200, "data": json.loads(body)}
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"acc": acc, "status": 404, "data": None}
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(delay + random.random()); delay = min(delay * 2, 60); continue
            return {"acc": acc, "status": e.code, "data": None}
        except Exception:
            time.sleep(delay + random.random()); delay = min(delay * 2, 60)
    return {"acc": acc, "status": -1, "data": None}

t0 = time.time(); n = 0
with open(OUT, "a") as fo, ThreadPoolExecutor(N_WORKERS) as ex:
    futs = {ex.submit(fetch, a): a for a in todo}
    for f in as_completed(futs):
        rec = f.result()
        with lock:
            fo.write(json.dumps(rec) + "\n"); fo.flush()
        n += 1
        if n % 500 == 0:
            el = time.time() - t0
            print(f"{n}/{len(todo)} in {el:.0f}s ({n/el:.1f}/s)", flush=True)
el = time.time() - t0
print(f"finished {n} in {el:.0f}s ({n/max(el,1e-9):.1f}/s)", flush=True)
