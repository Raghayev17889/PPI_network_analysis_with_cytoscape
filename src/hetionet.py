"""Extract a disease-centred dataset from Hetionet.

Hetionet node identifiers look like ``Gene::5663`` or ``Disease::DOID:10652``.
This module turns the node and edge files into plain tables keyed by Entrez
gene ID (as strings) and gene-set dictionaries that the rest of the project
uses:

=====================  ==========================================================
``interactions``       undirected protein-protein interactions (``GiG``)
``disease_genes``      genes associated with the disease (``DaG``) plus genes
                       up/down-regulated in disease tissue (``DuG``/``DdG``)
``pathways``           pathway gene sets (``GpPW``)
``go_processes``       GO biological-process gene sets (``GpBP``)
``drugs``              approved compounds that treat the disease (``CtD``)
``drug_targets``       proteins bound by those compounds (``CbG``)
``druggable``          proteins bound by any compound in Hetionet (``CbG``)
=====================  ==========================================================

Usage::

    python src/hetionet.py            # writes data/processed/
"""

from __future__ import annotations

import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from src.config import load_config

logger = logging.getLogger(__name__)

METAEDGES = ("GiG", "DaG", "DuG", "DdG", "GpPW", "GpBP", "CtD", "CbG")


def strip_kind(identifier: str) -> str:
    """``'Gene::5663'`` -> ``'5663'``; ``'Disease::DOID:10652'`` -> ``'DOID:10652'``."""
    return identifier.split("::", 1)[1] if "::" in identifier else identifier


def load_nodes(path: str | Path) -> pd.DataFrame:
    """Hetionet nodes with columns ``id``, ``name``, ``kind`` and the bare ``identifier``."""
    nodes = pd.read_csv(path, sep="\t", dtype=str)
    nodes["identifier"] = nodes["id"].map(strip_kind)
    return nodes


def load_edges(path: str | Path, metaedges: tuple[str, ...] | None = METAEDGES) -> pd.DataFrame:
    """Hetionet edges (``source``, ``metaedge``, ``target``), optionally filtered by type."""
    edges = pd.read_csv(path, sep="\t", dtype=str)
    if metaedges is not None:
        edges = edges[edges["metaedge"].isin(metaedges)]
    return edges.reset_index(drop=True)


@dataclass
class DiseaseData:
    """Everything the analysis needs about one disease, keyed by Entrez gene ID."""

    disease_id: str
    disease_name: str
    genes: pd.DataFrame                       # entrez, symbol
    interactions: pd.DataFrame                # gene_a, gene_b
    disease_genes: pd.DataFrame               # entrez, symbol, evidence
    pathways: dict[str, set[str]] = field(default_factory=dict)
    go_processes: dict[str, set[str]] = field(default_factory=dict)
    drugs: pd.DataFrame = field(default_factory=pd.DataFrame)          # drugbank_id, name
    drug_targets: pd.DataFrame = field(default_factory=pd.DataFrame)   # entrez, symbol, drug
    druggable: set[str] = field(default_factory=set)

    @property
    def symbols(self) -> dict[str, str]:
        return dict(zip(self.genes["entrez"], self.genes["symbol"]))

    def genes_with_evidence(self, evidence: str) -> list[str]:
        return sorted(self.disease_genes.loc[self.disease_genes["evidence"] == evidence, "entrez"])

    @property
    def seed_genes(self) -> list[str]:
        """Disease-associated genes (``DaG``): the seeds of the disease module."""
        return self.genes_with_evidence("associated")

    def summary(self) -> dict[str, int | str]:
        return {
            "disease": f"{self.disease_name} ({self.disease_id})",
            "genes": len(self.genes),
            "interactions": len(self.interactions),
            "associated_genes": len(self.seed_genes),
            "upregulated_genes": len(self.genes_with_evidence("upregulated")),
            "downregulated_genes": len(self.genes_with_evidence("downregulated")),
            "pathways": len(self.pathways),
            "go_processes": len(self.go_processes),
            "approved_drugs": len(self.drugs),
            "drug_target_genes": self.drug_targets["entrez"].nunique() if len(self.drug_targets) else 0,
            "druggable_genes": len(self.druggable),
        }


def _gene_sets(edges: pd.DataFrame, metaedge: str, nodes: pd.DataFrame) -> dict[str, set[str]]:
    """Gene sets named ``"<name> [<identifier>]"`` so that duplicate names stay distinct."""
    subset = edges[edges["metaedge"] == metaedge]
    names = dict(zip(nodes["id"], nodes["name"]))
    sets: dict[str, set[str]] = {}
    for target, group in subset.groupby("target"):
        label = f"{names.get(target, strip_kind(target))} [{strip_kind(target)}]"
        sets[label] = set(group["source"].map(strip_kind))
    return sets


def extract_disease_data(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    disease_id: str,
    disease_name: str | None = None,
) -> DiseaseData:
    """Build :class:`DiseaseData` for ``disease_id`` (e.g. ``"DOID:10652"``)."""
    disease_node = f"Disease::{disease_id}"
    if disease_node not in set(nodes["id"]):
        raise KeyError(f"Disease {disease_id} is not in Hetionet")
    names = dict(zip(nodes["id"], nodes["name"]))
    disease_name = disease_name or names[disease_node]

    gene_nodes = nodes[nodes["kind"] == "Gene"]
    genes = pd.DataFrame({"entrez": gene_nodes["identifier"], "symbol": gene_nodes["name"]}).reset_index(drop=True)
    symbol = dict(zip(genes["entrez"], genes["symbol"]))

    ppi = edges[edges["metaedge"] == "GiG"]
    interactions = pd.DataFrame(
        {"gene_a": ppi["source"].map(strip_kind), "gene_b": ppi["target"].map(strip_kind)}
    ).reset_index(drop=True)

    rows = []
    for metaedge, evidence in (("DaG", "associated"), ("DuG", "upregulated"), ("DdG", "downregulated")):
        targets = edges.loc[(edges["metaedge"] == metaedge) & (edges["source"] == disease_node), "target"]
        rows += [(strip_kind(t), evidence) for t in targets]
    disease_genes = pd.DataFrame(rows, columns=["entrez", "evidence"])
    disease_genes.insert(1, "symbol", disease_genes["entrez"].map(symbol))

    treats = edges[(edges["metaedge"] == "CtD") & (edges["target"] == disease_node)]
    drugs = pd.DataFrame({"drugbank_id": treats["source"].map(strip_kind), "name": treats["source"].map(names)})
    binds = edges[edges["metaedge"] == "CbG"]
    disease_binds = binds[binds["source"].isin(treats["source"])]
    drug_targets = pd.DataFrame(
        {
            "entrez": disease_binds["target"].map(strip_kind),
            "symbol": disease_binds["target"].map(strip_kind).map(symbol),
            "drug": disease_binds["source"].map(names),
        }
    ).sort_values(["entrez", "drug"]).reset_index(drop=True)

    data = DiseaseData(
        disease_id=disease_id,
        disease_name=disease_name,
        genes=genes,
        interactions=interactions,
        disease_genes=disease_genes.sort_values(["evidence", "entrez"]).reset_index(drop=True),
        pathways=_gene_sets(edges, "GpPW", nodes),
        go_processes=_gene_sets(edges, "GpBP", nodes),
        drugs=drugs.sort_values("name").reset_index(drop=True),
        drug_targets=drug_targets,
        druggable=set(binds["target"].map(strip_kind)),
    )
    logger.info("Extracted %s", data.summary())
    return data


# ---------------------------------------------------------------------------
# Persistence (data/processed/)
# ---------------------------------------------------------------------------


def write_gmt(gene_sets: dict[str, set[str]], path: str | Path) -> Path:
    """Write gene sets in GMT format (name, description, genes...)."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as handle:
        for name in sorted(gene_sets):
            label, _, identifier = name.rpartition(" [")
            members = "\t".join(sorted(gene_sets[name], key=lambda g: (len(g), g)))
            handle.write(f"{name}\t{identifier.rstrip(']') or label}\t{members}\n")
    return path


def read_gmt(path: str | Path) -> dict[str, set[str]]:
    sets = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = set(fields[2:])
    return sets


def save_disease_data(data: DiseaseData, directory: str | Path) -> Path:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    data.genes.to_csv(directory / "genes.tsv", sep="\t", index=False)
    data.interactions.to_csv(directory / "interactome.tsv", sep="\t", index=False)
    data.disease_genes.to_csv(directory / "disease_genes.tsv", sep="\t", index=False)
    data.drugs.to_csv(directory / "drugs.tsv", sep="\t", index=False)
    data.drug_targets.to_csv(directory / "drug_targets.tsv", sep="\t", index=False)
    pd.Series(sorted(data.druggable), name="entrez").to_csv(directory / "druggable_genes.tsv", sep="\t", index=False)
    write_gmt(data.pathways, directory / "pathways.gmt")
    write_gmt(data.go_processes, directory / "go_biological_process.gmt")
    meta = {"disease_id": data.disease_id, "disease_name": data.disease_name, "summary": data.summary()}
    (directory / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return directory


def load_disease_data(directory: str | Path) -> DiseaseData:
    directory = Path(directory)
    meta = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    read = lambda name: pd.read_csv(directory / name, sep="\t", dtype=str)  # noqa: E731
    return DiseaseData(
        disease_id=meta["disease_id"],
        disease_name=meta["disease_name"],
        genes=read("genes.tsv"),
        interactions=read("interactome.tsv"),
        disease_genes=read("disease_genes.tsv"),
        pathways=read_gmt(directory / "pathways.gmt"),
        go_processes=read_gmt(directory / "go_biological_process.gmt"),
        drugs=read("drugs.tsv"),
        drug_targets=read("drug_targets.tsv"),
        druggable=set(read("druggable_genes.tsv")["entrez"]),
    )


def prepare(config: dict | None = None, force: bool = False) -> DiseaseData:
    """Download Hetionet if needed, extract the configured disease and cache it."""
    from src.data_download import download_hetionet

    config = config or load_config()
    processed = Path(config["data"]["processed_dir"]) / config["disease"]["id"].replace(":", "_")
    if (processed / "metadata.json").exists() and not force:
        return load_disease_data(processed)
    files = download_hetionet(config["data"]["raw_dir"])
    data = extract_disease_data(
        load_nodes(files["nodes"]), load_edges(files["edges"]), config["disease"]["id"], config["disease"]["name"]
    )
    save_disease_data(data, processed)
    return data


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    result = prepare(force="--force" in sys.argv)
    for key, value in result.summary().items():
        print(f"{key:22s} {value}")
