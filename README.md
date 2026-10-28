# Interface removal reshapes the predicted homodimers of human splice isoforms while added sequence seldom builds interfaces

**Taner Karagöl¹,\*,‡, Alper Karagöl²,³,\*,‡**

¹ McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA\
² Broad Institute of MIT and Harvard, Cambridge, MA, USA\
³ Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA\
\* These authors contributed equally\
‡ Corresponding authors

## Abstract

Alternative splicing removes, replaces and adds protein sequence, yet structural surveys of isoforms have modelled single chains, and isoform and canonical homodimers have been compared at proteome scale only for truncated isoforms. We classified 20,171 human Swiss-Prot isoforms by the sequence blocks they share with their canonical forms and compared canonical and isoform homodimers in the 2026 AlphaFold Database complex release, whose deposited human models had an ipSAE of at least 0.6 and a pDockQ2 of at least 0.23. Most isoforms (60.8%) carried residues absent from the canonical, and of 376,017 such residues in the isoform monomer models 66.8% had pLDDT below 50, against 25.4% of shared residues. Among 18,508 pairs with both homodimers predicted, isoforms lost the canonical's high-confidence call in 44.9% of removal-only and 54.3% of replacement pairs, against 19.4% of addition-only pairs; the class odds ratios of 3.39 and 4.94 fell to between 0.99 and 1.76 once interface removal was in the model, and the rates at which interface and non-interface residues were removed showed no detectable difference. Isoforms that kept every interface residue lost the call in 17.3 to 19.4% of pairs and those that removed more than half in 94.2%. Gains were rare (1.79 to 2.24%), with no detectable difference between addition and removal and no detectable rise with more confident added sequence; replacement swaps were seldom homologous to what they replaced. The pattern held under ipSAE alone, pDockQ2 alone and ipTM thresholds. In three isoforms, confidently placed added residues form a predicted homodimer interface the canonical lacks. Splicing reshapes predicted homodimers mainly by removing interfaces, and interfaces built from added sequence are rare.

**Keywords:** Alternative splicing, Protein isoforms, Homodimer interfaces, AlphaFold Database, ipSAE

## Repository layout

```
splice-isoform-homodimers/
├── README.md                                 this file
├── LICENSE                                   MIT licence for the code
├── code/                                     pipeline, numbered in run order
│   ├── 01_fetch_uniprot.py                   human Swiss-Prot entries, isoform sequences and annotation from the UniProt REST API
│   ├── 02_parse_uniprot.py                   sequence table and isoform-canonical pairs
│   ├── 03_classify_isoforms.py               shared-block classification of every pair
│   ├── 05_crawl_afdb_complex.py              AlphaFold Database complex endpoint for every accession
│   ├── 06_parse_afdb_complex.py              served homodimer of each accession
│   ├── 07_crawl_afdb_monomer.py              AlphaFold Database monomer endpoint
│   ├── 07b_fetch_monomer_plddt.py            per-residue pLDDT of the monomer models
│   ├── 08_parse_afdb_monomer.py              monomer model summary of each accession
│   ├── 08b_parse_manifest.py                 manifest scores per accession and the deposition rule
│   ├── 09_build_pairs_table.py               pairs joined with the high-confidence calls
│   ├── 10_download_models.py                 homodimer models, fetched and checked against the deposited features
│   ├── 10a_select_sasa.py                    models whose solvent accessibility step 11 computes
│   ├── 11_structure_features.py              per-residue interface flags, pLDDT and accessibility
│   ├── 12_mechanism.py                       interface composition of each served isoform homodimer
│   ├── 12b_canonical_interface_retention.py  fraction of each canonical interface an isoform keeps
│   ├── 13_residue_table.py                   residues of isoform homodimers that carry added residues
│   ├── 14_identical_copies.py                identical sequences in other UniProtKB entries
│   ├── 14a_select_copy_sample.py             canonical entries of the identical-sequence check
│   ├── 14b_copies_served_entries.py          API entries of the served identical copies
│   ├── 15_subunit_calibration.py             calls against UniProt SUBUNIT annotation
│   ├── 16_analysis_main.py                   main statistics, into results/numbers.json
│   ├── 17_figures.py                         Figures 1 to 6 and Supplementary Figure S1
│   ├── 19_supplementary_tables.py            Supplementary Tables S1 to S9
│   ├── 21_swaps_interfaces_robustness.py     swaps, interface removal, robustness and origin
│   ├── 22_added_sequence_monomers.py         added and removed sequence in the monomer models
│   ├── 22b_added_sequence_disorder.py        sequence-based disorder of added and shared residues (AIUPred)
│   ├── 24_candidate_isoform_evidence.py      evidence of the four candidates, fetched and checked
│   ├── 25_mane_reference.py                  transitions and loss models with the MANE Select product as reference
│   ├── isoform_blocks.py                     shared blocks between an isoform and its canonical
│   └── colab/manifest_extract.ipynb          human rows of the release manifest, on Google Colab
├── data/                                     UniProt and AlphaFold Database extracts
│   ├── accessions_monomer.txt                accessions of the pairs within the length range
│   ├── accessions_to_query.txt               accessions queried at the complex endpoint
│   ├── afdb_complex_per_accession.tsv        complex endpoint result per accession
│   ├── afdb_complex_raw.jsonl                complex endpoint responses
│   ├── afdb_monomer_per_accession.tsv        monomer model summary per accession
│   ├── afdb_monomer_plddt.jsonl.gz.part0-2   per-residue pLDDT of the monomer models, in three parts
│   ├── copies_complex_raw.jsonl              complex endpoint responses for identical copies
│   ├── copies_served_entries.jsonl           API entries of the served copies
│   ├── copy_sample_accs.txt                  canonical entries of the identical-sequence check
│   ├── human_sp_altprod.tsv.gz               UniProt alternative-products annotation
│   ├── human_sp_isoforms.fasta.gz            human Swiss-Prot sequences with isoforms
│   ├── human_sp_subunit.tsv.gz               UniProt SUBUNIT annotation
│   ├── manifest_extract/                     manifest rows of human accessions and identical copies
│   ├── sasa_needed.txt                       models whose solvent accessibility step 11 computes
│   ├── uniparc_raw.jsonl                     UniParc records of the sampled sequences
│   ├── uniprot_human_sp_sequences.tsv.gz     sequences of UniProtKB release 2026_03
│   ├── uniprot_mane_select.tsv               MANE Select cross-references of the release
│   └── uniprot_release_headers.txt           UniProt release of the download
├── results/                                  derived tables and statistics
│   ├── added_residues_monomer.tsv.gz         confidence of added sequence per isoform
│   ├── added_sequence_disorder.json          disorder of added and shared residues, interface models
│   ├── aiupred_isoforms.jsonl.gz             per-residue AIUPred disorder of the isoforms
│   ├── candidate_isoform_evidence.json       protein and transcript evidence of the candidates
│   ├── canonical_retention.tsv               canonical interface kept per high-confidence pair
│   ├── canonical_retention_analysed.tsv      the same for the analysed pairs, with scores
│   ├── canonical_subunit_calibration.tsv     calls against UniProt SUBUNIT annotation
│   ├── identical_copies.tsv                  identical copies and their scores
│   ├── interface_removal_residues.tsv.gz     canonical residues by interface, band and removal
│   ├── isoform_pairs_classified.tsv          every isoform-canonical pair and its class
│   ├── mane_reference.json                   the analysis with the MANE Select product as reference
│   ├── manifest_per_accession.tsv            manifest scores and deposition call per accession
│   ├── manifest_rule_check.json              deposition rule against the API
│   ├── numbers.json                          statistics of the main analysis
│   ├── numbers_new.json                      statistics of the monomer, swap and robustness analyses
│   ├── pair_mechanism.tsv                    interface composition of served isoform homodimers
│   ├── pairs_analysis.tsv                    classified pairs with scores, calls and confidence
│   ├── pairs_with_hc.tsv                     pairs within the length range with their calls
│   ├── residues_novel_models.tsv.gz          residues of isoform homodimers with added residues
│   ├── shared_fold_confidence.tsv.gz         confidence change of the shared residues
│   ├── struct_features.tar.gz                per-residue features of every homodimer model used
│   ├── swap_segments.tsv                     replacement swaps and their homology call
│   └── table_*.tsv                           tables behind the figures and supplementary tables (10 files)
├── figures/                                  Figures 1 to 6 and Supplementary Figure S1, PNG and SVG
└── supplementary/
    ├── Supplementary_Information.pdf         supplementary information
    └── Supplementary_Tables.xlsx             Supplementary Tables S1 to S9
```

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `results/isoform_pairs_classified.tsv` | 22,131 | Every human Swiss-Prot isoform paired with its canonical: class, added and removed residues by position, shared blocks, and whether both sequences are 15 to 1,500 residues long. |
| `results/pairs_with_hc.tsv` | 20,171 | The pairs within that length range, with the high-confidence call of each homodimer. |
| `results/pairs_analysis.tsv` | 20,163 | The classified pairs that share at least one block, with manifest scores, calls, model lengths and monomer confidence. |
| `results/manifest_per_accession.tsv` | 91,641 | One manifest row per human accession: chain length, ipTM, pDockQ, the larger-direction ipSAE, pDockQ2 and LIS, the smaller-direction ipSAE, clash counts and the deposition call. |
| `results/canonical_retention_analysed.tsv` | 1,662 | Fraction of each high-confidence canonical interface that the isoform keeps, for the analysed pairs. |
| `results/interface_removal_residues.tsv.gz` | 741,321 | Canonical residues of high-confidence pairs whose isoform removes residues, with interface label, pLDDT band and removal. |
| `results/swap_segments.tsv` | 11,336 | Every replacement swap, its local alignment to the sequence it replaces and the homology call. |
| `results/added_residues_monomer.tsv.gz` | 12,216 | Confidence of added and shared residues per isoform in the monomer models. |
| `results/residues_novel_models.tsv.gz` | 326,484 | Residues of the high-confidence isoform homodimers that carry added residues, with interface label, pLDDT and accessibility. |
| `results/aiupred_isoforms.jsonl.gz` | 12,263 | AIUPred disorder of every residue of the isoforms that carry added residues, one record per isoform. |
| `results/struct_features.tar.gz` | 2,180 files | Per-residue interface flags, pLDDT and relative accessibility of every homodimer model used. |
| `data/uniprot_human_sp_sequences.tsv.gz` | 42,562 | Sequences of the human Swiss-Prot entries and isoforms, UniProtKB release 2026_03. |
| `data/afdb_complex_per_accession.tsv` | 32,719 | AlphaFold Database complex endpoint result for every queried accession. |
| `data/manifest_extract/manifest_human.tsv.gz` | 91,641 | Human rows of the release manifest (model_entity_metadata_mapping.csv). |
| `data/afdb_monomer_plddt.jsonl.gz.part0-2` | 29,975 | Per-residue pLDDT of the monomer models, one record per accession. |
| `data/uniprot_mane_select.tsv` | 20,431 | MANE Select transcript of each human Swiss-Prot entry; an isoform identifier in brackets marks a MANE Select product other than the displayed isoform. |

`results/numbers.json` and `results/numbers_new.json` hold the computed statistics, `results/added_sequence_disorder.json`, `results/candidate_isoform_evidence.json` and `results/mane_reference.json` those of the disorder, candidate and MANE Select analyses, and `supplementary/Supplementary_Tables.xlsx` holds Supplementary Tables S1 to S9.

## Reproducing the analysis

All scripts use relative paths and run from the repository root. The steps below run offline on the deposited files and reproduce the tables and statistics in `results/`, the figures and the workbook.

```bash
cat data/afdb_monomer_plddt.jsonl.gz.part0 data/afdb_monomer_plddt.jsonl.gz.part1 \
    data/afdb_monomer_plddt.jsonl.gz.part2 > data/afdb_monomer_plddt.jsonl.gz   # MD5 72d488099daf296a5ef62a051daa3a19
tar xzf results/struct_features.tar.gz -C results     # per-residue features into results/struct/
python3 code/02_parse_uniprot.py
python3 code/03_classify_isoforms.py
python3 code/06_parse_afdb_complex.py
python3 code/08b_parse_manifest.py                    # stops unless the deposition rule reproduces the API
python3 code/09_build_pairs_table.py
python3 code/10a_select_sasa.py
python3 code/12_mechanism.py
python3 code/12b_canonical_interface_retention.py
python3 code/13_residue_table.py
python3 code/14a_select_copy_sample.py
python3 code/15_subunit_calibration.py
python3 code/16_analysis_main.py
python3 code/21_swaps_interfaces_robustness.py
python3 code/22_added_sequence_monomers.py
python3 code/22b_added_sequence_disorder.py           # reads the deposited AIUPred scores
python3 code/25_mane_reference.py
python3 code/17_figures.py
python3 code/19_supplementary_tables.py
```

Steps 01, 05, 07, 07b, 14 and 14b query UniProt or the AlphaFold Database and return the releases current when they run; their outputs used here are in `data/`. Step 08 parses the monomer summary crawl that step 07 writes into `data/afdb_monomer_per_accession.tsv`. `code/10_download_models.py` fetches each homodimer model by the model identifier and version in `data/afdb_complex_per_accession.tsv`, checks the sequence and per-residue pLDDT of both chains against `results/struct_features.tar.gz`, and stops if any model is missing or differs; step 11 computes the per-residue features from them. Step 14 takes the sample list as its argument: `python3 code/14_identical_copies.py data/copy_sample_accs.txt`. Step 22b runs AIUPred only for isoforms missing from `results/aiupred_isoforms.jsonl.gz`, and step 25 queries UniProt only if `data/uniprot_mane_select.tsv` is missing and stops unless the release is 2026_03. `code/24_candidate_isoform_evidence.py` fetches the UniProtKB entries, mapped peptides and Ensembl and RefSeq proteins of the four candidate isoforms, checks the MANE Select status, isoform-specific peptides and gene proteins derived from them against `results/candidate_isoform_evidence.json`, and stops if any of them differs.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Data handling | `pandas`, `numpy` |
| Statistics | `scipy`, `statsmodels` |
| Sequence and structure | `biopython`, `gemmi` |
| Database queries | `requests`, `urllib` (standard library) |
| Figures and workbook | `matplotlib`, `openpyxl` |
| Manifest extraction | `polars` (Google Colab) |
| Disorder prediction | `aiupred` (https://github.com/doszilab/AIUPred), `torch` |

Analysis tested with Python 3.13.16, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0, Biopython 1.88 and gemmi 0.7.5; figures and workbook with Python 3.10.12, pandas 2.3.3, matplotlib 3.10.9 and openpyxl 3.1.5; the manifest was read with Polars 1.35.2, and disorder was predicted with AIUPred 3.1.3 under PyTorch 2.14.1. AIUPred is available from its authors under an academic licence.

## Data availability

* UniProtKB release 2026_03, human Swiss-Prot entries with isoforms, retrieved through the UniProt REST API (`https://rest.uniprot.org/`); the extracts used are in `data/`.
* AlphaFold Protein Structure Database 2026 complex release: API entries (`https://alphafold.ebi.ac.uk/`) and the human and identical-copy rows of the release manifest are in `data/`; the homodimer models are fetched and checked by `code/10_download_models.py`.
* MANE Select cross-references of UniProtKB release 2026_03 are in `data/`. For the four candidate isoforms, the UniProtKB entries, the peptides of the EBI Proteins API proteomics mapping (`https://www.ebi.ac.uk/proteins/api/`) and the proteins of their genes from the Ensembl REST API (`https://rest.ensembl.org/`) and NCBI RefSeq (`https://eutils.ncbi.nlm.nih.gov/`) are fetched and checked by `code/24_candidate_isoform_evidence.py`; the evidence derived from them is in `results/candidate_isoform_evidence.json`.
* Disorder was predicted with AIUPred (Erdős and Dosztányi, 2024); the per-residue scores are in `results/`.
* The supplementary information and Supplementary Tables S1 to S9 are in `supplementary/`.

## Citation

Karagöl T., Karagöl A. *Interface removal reshapes the predicted homodimers of human splice isoforms while added sequence seldom builds interfaces.* (2026).

## Contact

Questions and issues are welcome through GitHub Issues or by email to the corresponding authors, Taner Karagöl (karagolt@mit.edu) and Alper Karagöl (akaragol@mit.edu).

## License

The code in this repository is released under the MIT License (`LICENSE`). The UniProt and AlphaFold Database extracts in `data/` and the tables derived from them are distributed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), the licence of those sources.
