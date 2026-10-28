import gzip, os, sys, tarfile, time, random, urllib.request, urllib.error
import gemmi
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 4
cx = pd.read_csv("data/afdb_complex_per_accession.tsv", sep="\t")
pairs = pd.read_csv("results/pairs_with_hc.tsv", sep="\t")
need = set(pairs.loc[pairs.canon_hc, "canonical_acc"]) | set(pairs.loc[pairs.iso_hc, "isoform_acc"])
cx = cx[cx.hc_homodimer & cx.acc.isin(need)]
os.makedirs("models", exist_ok=True)
todo = [(r.acc, r.model_id, int(r.version)) for r in cx.itertuples() if not os.path.exists(f"models/{r.acc}.cif.gz")]
print("models needed:", len(cx), "to download:", len(todo), flush=True)
def get(item):
    acc, mid, v = item
    url = f"https://alphafold.ebi.ac.uk/files/{mid}-model_v{v}.cif"
    delay = 1.0
    for _ in range(8):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "isoform-homodimer-study (karagolt@mit.edu)"}), timeout=60) as r:
                data = r.read()
            with gzip.open(f"models/{acc}.cif.gz.part", "wb") as f: f.write(data)
            os.replace(f"models/{acc}.cif.gz.part", f"models/{acc}.cif.gz")
            return acc, "ok"
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504): time.sleep(delay + random.random()); delay = min(2 * delay, 60); continue
            return acc, f"http {e.code}"
        except Exception:
            time.sleep(delay + random.random()); delay = min(2 * delay, 60)
    return acc, "failed"
bad = []
t0 = time.time()
with ThreadPoolExecutor(N_WORKERS) as ex:
    for k, f in enumerate(as_completed([ex.submit(get, it) for it in todo]), 1):
        acc, st = f.result()
        if st != "ok": bad.append((acc, st))
        if k % 250 == 0: print(f"{k}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
print("done; failures:", len(bad), bad[:10], flush=True)
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
ref = {}
with tarfile.open("results/struct_features.tar.gz") as t:
    for m in t.getmembers():
        if m.name.endswith(".tsv.gz"):
            rows = gzip.decompress(t.extractfile(m).read()).decode().splitlines()[1:]
            ref[os.path.basename(m.name)[:-len(".tsv.gz")]] = [tuple(r.split("\t")[1:4]) for r in rows]
def check(acc):
    fn = f"models/{acc}.cif.gz"
    if not os.path.exists(fn):
        return acc, "missing"
    if acc not in ref:
        return acc, "no deposited features"
    st = gemmi.read_structure(fn)
    st.remove_hydrogens()
    chains = [ch for ch in st[0]]
    if len(chains) != 2 or len(chains[0]) != len(chains[1]):
        return acc, "chains differ"
    got = [(ONE.get(a.name, "X"), f"{a[0].b_iso:.2f}", f"{b[0].b_iso:.2f}") for a, b in zip(chains[0], chains[1])]
    return acc, "ok" if got == ref[acc] else "differs"
with ThreadPoolExecutor(N_WORKERS) as ex:
    res = list(ex.map(check, sorted(cx.acc)))
off = [r for r in res if r[1] != "ok"]
print(f"models checked: {len(res)}; matching the deposited features: {len(res) - len(off)}", flush=True)
if off:
    sys.exit(f"STOP: {len(off)} models are missing or differ from the models analysed: {off[:10]}")
