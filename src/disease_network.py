"""Construct and filter the Alzheimer's disease interaction network.

Filtering steps (each recorded in :attr:`DiseaseNetwork.filter_log`):

1. **seed genes** - disease-associated genes from Hetionet;
2. **in interactome** - seeds present in the filtered human interactome;
3. **module expansion** - add the top DIAMOnD proteins;
4. **induced subgraph** - keep only interactions among the selected proteins;
5. **largest component** - drop proteins that are not connected to the main
   network (they cannot have meaningful centralities or clusters).

Every node is annotated with what Cytoscape needs for styling: gene symbol,
role (``seed`` or ``diamond``), DIAMOnD rank and p-value, interactome degree,
differential expression in disease tissue, and drug-target status.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx
import pandas as pd

from src.hetionet import DiseaseData

logger = logging.getLogger(__name__)

NODE_COLUMNS = [
    "id", "symbol", "role", "diamond_rank", "diamond_p", "interactome_degree", "expression",
    "ad_drug_target", "ad_drugs", "druggable",
]


@dataclass
class DiseaseNetwork:
    graph: nx.Graph
    filter_log: list[dict] = field(default_factory=list)

    def log_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.filter_log).set_index("step")

    def node_table(self) -> pd.DataFrame:
        """One row per protein with every node attribute (Cytoscape node table)."""
        rows = [{"id": n, **attrs} for n, attrs in self.graph.nodes(data=True)]
        table = pd.DataFrame(rows)
        ordered = [c for c in NODE_COLUMNS if c in table.columns]
        rest = [c for c in table.columns if c not in ordered]
        return table[ordered + rest].sort_values(["role", "symbol"], ascending=[False, True]).reset_index(drop=True)

    def edge_table(self) -> pd.DataFrame:
        """One row per interaction (Cytoscape edge table: source, target, interaction, ...)."""
        rows = [
            {"source": u, "target": v, "interaction": "pp", **attrs} for u, v, attrs in self.graph.edges(data=True)
        ]
        return pd.DataFrame(rows).sort_values(["source", "target"]).reset_index(drop=True)


def _log(log: list, step: str, description: str, proteins: int, interactions: int | None = None) -> None:
    log.append({"step": step, "description": description, "proteins": proteins, "interactions": interactions})


def build_disease_network(
    interactome: nx.Graph,
    data: DiseaseData,
    diamond_genes: pd.DataFrame,
    n_diamond: int = 100,
) -> DiseaseNetwork:
    """Seeds + top ``n_diamond`` DIAMOnD proteins, filtered and annotated.

    Parameters
    ----------
    interactome:
        Filtered human interactome (see :func:`src.interactome.build_interactome`).
    data:
        Disease data extracted from Hetionet.
    diamond_genes:
        Output of :func:`src.diamond.diamond`, ranked.
    n_diamond:
        How many DIAMOnD proteins to add.
    """
    log: list[dict] = []
    seeds_all = data.seed_genes
    _log(log, "seed genes", "disease-associated genes (Hetionet DaG)", len(seeds_all))
    seeds = [s for s in seeds_all if s in interactome]
    _log(log, "in interactome", "seeds present in the filtered interactome", len(seeds))

    added = diamond_genes.head(n_diamond)
    selected = list(dict.fromkeys(seeds + added["gene"].tolist()))
    _log(log, "module expansion", f"+ top {len(added)} DIAMOnD proteins", len(selected))

    graph = interactome.subgraph(selected).copy()
    _log(log, "induced subgraph", "interactions among selected proteins", graph.number_of_nodes(), graph.number_of_edges())

    giant = max(nx.connected_components(graph), key=len)
    dropped = sorted(set(graph) - giant)
    graph = graph.subgraph(giant).copy()
    _log(log, "largest component", "drop proteins disconnected from the main network",
         graph.number_of_nodes(), graph.number_of_edges())

    # ---- annotation -----------------------------------------------------
    symbols = data.symbols
    rank = dict(zip(added["gene"], added["rank"]))
    p_value = dict(zip(added["gene"], added["p_value"]))
    up = set(data.genes_with_evidence("upregulated"))
    down = set(data.genes_with_evidence("downregulated"))
    drugs_by_gene = data.drug_targets.groupby("entrez")["drug"].apply(lambda s: "; ".join(sorted(set(s)))).to_dict()
    seed_set = set(seeds)

    for node in graph:
        graph.nodes[node].update(
            {
                "symbol": symbols.get(node, node),
                "role": "seed" if node in seed_set else "diamond",
                "diamond_rank": int(rank[node]) if node in rank else 0,
                "diamond_p": float(p_value[node]) if node in p_value else 1.0,
                "interactome_degree": int(interactome.degree(node)),
                "expression": "up" if node in up else "down" if node in down else "none",
                "ad_drug_target": node in drugs_by_gene,
                "ad_drugs": drugs_by_gene.get(node, ""),
                "druggable": node in data.druggable,
            }
        )
    for u, v in graph.edges():
        roles = sorted((graph.nodes[u]["role"], graph.nodes[v]["role"]))
        graph.edges[u, v]["edge_type"] = "-".join(roles)

    network = DiseaseNetwork(graph, log)
    if dropped:
        logger.info("Dropped %d disconnected proteins: %s", len(dropped), [symbols.get(d, d) for d in dropped])
    logger.info("Disease network: %d proteins, %d interactions", graph.number_of_nodes(), graph.number_of_edges())
    return network


def export_tables(network: DiseaseNetwork, directory: str | Path, prefix: str = "ad_network") -> dict[str, Path]:
    """Write node, edge and filter-log tables (importable into Cytoscape)."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {
        "nodes": directory / f"{prefix}_nodes.tsv",
        "edges": directory / f"{prefix}_edges.tsv",
        "filter_log": directory / f"{prefix}_filtering.tsv",
    }
    network.node_table().to_csv(paths["nodes"], sep="\t", index=False)
    network.edge_table().to_csv(paths["edges"], sep="\t", index=False)
    network.log_frame().to_csv(paths["filter_log"], sep="\t")
    return paths
