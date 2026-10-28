import glob, gzip, os, sys, tempfile
import numpy as np
from multiprocessing import Pool

CUTOFF = 5.0
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLU": 223, "GLN": 225, "GLY": 104,
          "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
          "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
SASA_NEEDED = set(open("data/sasa_needed.txt").read().split()) if os.path.exists("data/sasa_needed.txt") else set()


def rsa_chainA(fn):
    from Bio.PDB import MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley
    with gzip.open(fn, "rt") as f, tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as tmp:
        tmp.write(f.read())
        path = tmp.name
    s = MMCIFParser(QUIET=True).get_structure("m", path)
    os.unlink(path)
    A = list(s[0])[0]
    ShrakeRupley(n_points=100).compute(A, level="R")
    return [r.sasa / MAXASA.get(r.get_resname(), 200) for r in A if r.id[0] == " "]


def process(fn):
    import gemmi
    from scipy.spatial import cKDTree
    acc = os.path.basename(fn)[:-len(".cif.gz")]
    out = f"results/struct/{acc}.tsv.gz"
    if os.path.exists(out):
        return acc, "skip"
    try:
        st = gemmi.read_structure(fn)
        st.remove_hydrogens()
        chains = [ch for ch in st[0]]
        if len(chains) != 2:
            return acc, f"chains={len(chains)}"
        A, B = chains
        if len(A) != len(B):
            return acc, "chain length mismatch"

        def heavy(ch):
            xyz, owner = [], []
            for k, r in enumerate(ch):
                for a in r:
                    xyz.append((a.pos.x, a.pos.y, a.pos.z))
                    owner.append(k)
            return np.array(xyz), np.array(owner)

        xa, oa = heavy(A)
        xb, ob = heavy(B)
        tb = cKDTree(xb)
        partners = [set() for _ in range(len(A))]
        iface_B = np.zeros(len(B), bool)
        for i, nb in enumerate(tb.query_ball_point(xa, CUTOFF)):
            for j in nb:
                partners[oa[i]].add(int(ob[j]))
                iface_B[ob[j]] = True
        rsa = rsa_chainA(fn) if acc in SASA_NEEDED else None
        lines = ["pos\taa\tplddt_A\tplddt_B\trsa_A\tiface_A\tiface_B\tpartners_A"]
        for k in range(len(A)):
            ra, rb = A[k], B[k]
            rsa_s = "" if rsa is None else f"{rsa[k]:.3f}"
            lines.append(f"{k+1}\t{ONE.get(ra.name, 'X')}\t{ra[0].b_iso:.2f}\t{rb[0].b_iso:.2f}\t{rsa_s}\t"
                         f"{int(len(partners[k]) > 0)}\t{int(iface_B[k])}\t"
                         + ";".join(str(j + 1) for j in sorted(partners[k])))
        with gzip.open(out + ".part", "wt") as f:
            f.write("\n".join(lines) + "\n")
        os.replace(out + ".part", out)
        return acc, "ok"
    except Exception as e:
        return acc, f"error {type(e).__name__}: {e}"


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    os.makedirs("results/struct", exist_ok=True)
    files = sorted(glob.glob("models/*.cif.gz"))
    with Pool(n) as pool:
        res = pool.map(process, files, chunksize=8)
    bad = [r for r in res if r[1] not in ("ok", "skip")]
    print(f"models {len(files)}, processed {sum(r[1] == 'ok' for r in res)}, "
          f"skipped {sum(r[1] == 'skip' for r in res)}, problems {len(bad)}")
    for b in bad[:20]:
        print("  ", b)
