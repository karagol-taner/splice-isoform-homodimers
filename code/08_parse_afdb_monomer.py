import json
import pandas as pd
rows = []
for line in open("data/afdb_monomer_raw.jsonl"):
    r = json.loads(line)
    acc, st, data = r["acc"], r["status"], r["data"]
    rec = {"acc": acc, "mono_status": st}
    if st == 200 and data:
        hit = [e for e in data if e.get("uniprotAccession") == acc and not e.get("isComplex")]
        if hit:
            e = hit[0]
            rec.update({"mono_entry": e.get("entryId"), "mono_version": e.get("latestVersion"),
                        "mono_plddt": e.get("globalMetricValue"),
                        "mono_frac_vlow": e.get("fractionPlddtVeryLow"), "mono_frac_low": e.get("fractionPlddtLow"),
                        "mono_frac_conf": e.get("fractionPlddtConfident"), "mono_frac_vhigh": e.get("fractionPlddtVeryHigh"),
                        "mono_len": (e.get("sequenceEnd") or 0) - (e.get("sequenceStart") or 1) + 1,
                        "mono_seq_checksum": e.get("sequenceChecksum")})
    rows.append(rec)
df = pd.DataFrame(rows)
df.to_csv("data/afdb_monomer_per_accession.tsv", sep="\t", index=False)
print("accessions:", len(df), "| with monomer summary:", int(df.mono_plddt.notna().sum()), "| status:", df.mono_status.value_counts().to_dict())
