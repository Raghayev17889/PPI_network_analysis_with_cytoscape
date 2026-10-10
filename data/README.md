# Data

Nothing in this folder is version-controlled. Everything is downloaded or regenerated:

```bash
python src/data_download.py   # data/raw/: Hetionet v1.0 nodes + edges (~15 MB, SHA-256 verified)
python src/hetionet.py        # data/processed/DOID_10652/: tables for the configured disease
```

| file (`data/processed/DOID_10652/`) | contents |
|---|---|
| `genes.tsv` | Entrez ID and symbol of all 20,945 Hetionet genes |
| `interactome.tsv` | 147,164 protein-protein interactions (`GiG`) |
| `disease_genes.tsv` | associated (`DaG`), up-regulated (`DuG`) and down-regulated (`DdG`) genes |
| `pathways.gmt`, `go_biological_process.gmt` | gene sets in GMT format (usable in Cytoscape EnrichmentMap) |
| `drugs.tsv`, `drug_targets.tsv` | approved drugs treating the disease and the proteins they bind |
| `druggable_genes.tsv` | proteins bound by any DrugBank compound |
| `metadata.json` | disease ID, name and counts |

Source: Himmelstein et al., *eLife* 2017, <https://github.com/hetio/hetionet>, commit
`8a6cc0c5604d0e2908786e4de38af67b8e46ee4a` (CC0).
