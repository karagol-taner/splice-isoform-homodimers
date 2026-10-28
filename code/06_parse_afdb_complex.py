import json
import pandas as pd
rows = []
for line in open("data/afdb_complex_raw.jsonl"):
    r = json.loads(line)
    acc, st, data = r["acc"], r["status"], r["data"]
    rec = {"acc": acc, "api_status": st, "hc_homodimer": False}
    if st == 200 and data:
        homo = [e for e in data if e.get("assemblyType") == "Homo" and e.get("oligomericState") == "dimer"
                and [c.get("identifier") for c in e.get("complexComposition", [])] == [acc]
                and [c.get("stoichiometry") for c in e.get("complexComposition", [])] == [2]]
        hetero = [e for e in data if e.get("assemblyType") == "Hetero"]
        rec["n_hetero_entries"] = len(hetero)
        if homo:
            e = sorted(homo, key=lambda x: -x.get("latestVersion", 0))[0]
            rec.update({"hc_homodimer": True, "model_id": e.get("modelEntityId"), "version": e.get("latestVersion"),
                        "pLDDT": e.get("globalMetricValue"), "ipTM": e.get("complexPredictionAccuracy_ipTM"),
                        "ipSAE": e.get("complexPredictionAccuracy_ipSAE"), "pDockQ": e.get("complexPredictionAccuracy_pDockQ"),
                        "pDockQ2": e.get("complexPredictionAccuracy_pDockQ2"), "LIS": e.get("complexPredictionAccuracy_LIS"),
                        "seq_end": e.get("sequenceEnd"), "provider": e.get("providerId"), "is_isoform_flag": e.get("isIsoform"),
                        "n_homo_entries": len(homo)})
    rows.append(rec)
df = pd.DataFrame(rows)
df.to_csv("data/afdb_complex_per_accession.tsv", sep="\t", index=False)
print("accessions:", len(df), "| status counts:", df.api_status.value_counts().to_dict())
hc = df[df.hc_homodimer]
print("served homodimers:", len(hc))
for c in ["ipTM", "ipSAE", "pLDDT", "pDockQ2", "LIS"]:
    print(f"  {c}: min {hc[c].min()}, 1st pct {hc[c].quantile(0.01):.3f}, median {hc[c].median():.3f}, max {hc[c].max()}")
print("served with ipSAE < 0.6:", int((hc.ipSAE < 0.6).sum()), "| pLDDT < 70:", int((hc.pLDDT < 70).sum()))
print("accessions with >1 homodimer entry:", int((hc.n_homo_entries > 1).sum()))
