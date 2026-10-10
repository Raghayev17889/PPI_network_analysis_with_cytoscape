"""Build the human interactome graph and map disease genes onto it."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)


def induced_subgraph(graph: nx.Graph, nodes) -> nx.Graph:
    """Copy of the subgraph induced by ``nodes``, in ``graph``'s own node and edge order.

    ``graph.subgraph(nodes)`` iterates over ``nodes`` as a set, and the order of a
    set of strings changes with Python's per-process hash seed. That order then
    leaks into edge orientation, layouts and floating-point sums (betweenness),
    so results would differ slightly between runs. Building the copy in the
    parent graph's order keeps every output byte-for-byte reproducible.
    """
    keep = set(nodes)
    sub = nx.Graph()
    sub.add_nodes_from((n, graph.nodes[n]) for n in graph if n in keep)
    sub.add_edges_from((u, v, d) for u, v, d in graph.edges(data=True) if u in keep and v in keep)
    return sub


def build_interactome(
    interactions: pd.DataFrame,
    symbols: dict[str, str] | None = None,
    largest_component: bool = True,
    exclude: list[str] | set[str] = (),
) -> nx.Graph:
    """Undirected protein-protein interaction graph keyed by Entrez ID.

    Self-interactions and duplicate pairs are dropped. Nodes carry a
    ``symbol`` attribute. ``exclude`` removes proteins before the component
    is computed (used for ubiquitin: UBC alone has 8,611 partners in
    Hetionet because proteomics screens record every ubiquitinated protein as
    an interaction, which would make it the bottleneck of any network). With
    ``largest_component`` only the largest connected component is kept, on
    which network distances are defined.
    """
    graph = nx.Graph()
    pairs = interactions[["gene_a", "gene_b"]].astype(str)
    pairs = pairs[pairs["gene_a"] != pairs["gene_b"]]
    graph.add_edges_from(pairs.itertuples(index=False, name=None))
    excluded = [g for g in map(str, exclude) if g in graph]
    graph.remove_nodes_from(excluded)
    if excluded:
        logger.info("Excluded %d proteins: %s", len(excluded), excluded)
    if largest_component and graph.number_of_nodes():
        giant = max(nx.connected_components(graph), key=len)
        removed = graph.number_of_nodes() - len(giant)
        graph = induced_subgraph(graph, giant)
        logger.info("Interactome LCC: %d proteins (%d outside removed)", len(giant), removed)
    if symbols:
        nx.set_node_attributes(graph, {n: symbols.get(n, n) for n in graph}, "symbol")
    return graph


@dataclass
class SeedMapping:
    """Which disease genes are present in the interactome."""

    mapped: list[str]
    missing: list[str]

    @property
    def coverage(self) -> float:
        total = len(self.mapped) + len(self.missing)
        return len(self.mapped) / total if total else 0.0


def symbols_to_entrez(symbols: list[str], mapping: dict[str, str]) -> list[str]:
    """Translate gene symbols to Entrez IDs using an ``entrez -> symbol`` mapping."""
    reverse = {v: k for k, v in mapping.items()}
    return [reverse[s] for s in symbols if s in reverse]


def map_genes(graph: nx.Graph, genes: list[str] | set[str]) -> SeedMapping:
    """Split ``genes`` into those present in ``graph`` and those absent."""
    genes = sorted(set(map(str, genes)))
    present = [g for g in genes if g in graph]
    absent = [g for g in genes if g not in graph]
    return SeedMapping(present, absent)


def interactome_summary(graph: nx.Graph) -> dict[str, float]:
    degrees = [d for _, d in graph.degree()]
    n = graph.number_of_nodes()
    return {
        "proteins": n,
        "interactions": graph.number_of_edges(),
        "mean_degree": 2 * graph.number_of_edges() / n if n else 0.0,
        "max_degree": max(degrees, default=0),
        "density": nx.density(graph) if n > 1 else 0.0,
    }
