import os, urllib.request
os.makedirs("data", exist_ok=True)
Q = "%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29"
FILES = {
    "data/human_sp_isoforms.fasta.gz": f"https://rest.uniprot.org/uniprotkb/stream?query={Q}&format=fasta&includeIsoform=true&compressed=true",
    "data/human_sp_altprod.tsv.gz": f"https://rest.uniprot.org/uniprotkb/stream?query={Q}&fields=accession%2Cid%2Cgene_primary%2Clength%2Ccc_alternative_products%2Cft_var_seq%2Cprotein_existence&format=tsv&compressed=true",
    "data/human_sp_subunit.tsv.gz": f"https://rest.uniprot.org/uniprotkb/stream?query={Q}&fields=accession%2Ccc_subunit&format=tsv&compressed=true",
}
for fn, url in FILES.items():
    with urllib.request.urlopen(url, timeout=600) as r:
        rel, date = r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")
        open(fn, "wb").write(r.read())
    print(fn, "UniProt release", rel, date)
