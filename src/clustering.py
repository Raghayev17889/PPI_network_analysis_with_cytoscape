"""Functional clusters: community detection in the disease network.

Louvain modularity optimisation (Blondel et al., 2008) partitions the
network into densely connected groups of proteins. Clusters are relabelled
1, 2, ... by decreasing size; clusters smaller than ``min_size`` are merged
into cluster 0 ("other") so the Cytoscape colour legend stays readable.
Because Louvain is randomised, :func:`partition_stability` reruns it with
different seeds and reports the adjusted Rand index against the reference.
"""

from __future__ import annotations

import logging
from math import comb

import networkx as nx
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def detect_communities(
    graph: nx.Graph,
    method: str = "louvain",
    resolution: float = 1.0,
    seed: int = 42,
) -> list[set[str]]:
    """Communities sorted by decreasing size."""
    if method == "louvain":
        communities = nx.community.louvain_communities(graph, resolution=resolution, seed=seed)
    elif method == "greedy":
        communities = nx.community.greedy_modularity_communities(graph, resolution=resolution)
    else:
        raise ValueError(f"Unknown method {method!r}; use 'louvain' or 'greedy'")
    return sorted((set(c) for c in communities), key=lambda c: (-len(c), min(c)))


def membership(communities: list[set[str]], min_size: int = 1) -> dict[str, int]:
    """``protein -> cluster`` with clusters numbered from 1; small ones become 0."""
    labels: dict[str, int] = {}
    next_label = 1
    for community in communities:
        if len(community) >= min_size:
            label, next_label = next_label, next_label + 1
        else:
            label = 0
        for node in community:
            labels[node] = label
    return labels


def adjusted_rand_index(a: dict[str, int], b: dict[str, int]) -> float:
    """Adjusted Rand index between two partitions of the same proteins."""
    nodes = sorted(set(a) & set(b))
    table = pd.crosstab(pd.Series([a[n] for n in nodes]), pd.Series([b[n] for n in nodes])).to_numpy()
    sum_cells = sum(comb(int(x), 2) for x in table.ravel())
    sum_rows = sum(comb(int(x), 2) for x in table.sum(axis=1))
    sum_cols = sum(comb(int(x), 2) for x in table.sum(axis=0))
    total = comb(len(nodes), 2)
    expected = sum_rows * sum_cols / total if total else 0.0
    maximum = (sum_rows + sum_cols) / 2
    return float((sum_cells - expected) / (maximum - expected)) if maximum != expected else 1.0


def partition_stability(
    graph: nx.Graph, reference: dict[str, int], runs: int = 20, resolution: float = 1.0, seed: int = 0
) -> pd.DataFrame:
    """Adjusted Rand index of ``runs`` Louvain partitions against ``reference``."""
    rows = []
    for run in range(runs):
        other = membership(detect_communities(graph, "louvain", resolution, seed=seed + run + 1))
        rows.append({"run": run + 1, "adjusted_rand_index": adjusted_rand_index(reference, other)})
    return pd.DataFrame(rows)


def cluster_summary(graph: nx.Graph, labels: dict[str, int], top: int = 5) -> pd.DataFrame:
    """Size, composition, density and most connected proteins of every cluster."""
    rows = []
    for cluster in sorted(set(labels.values())):
        members = [n for n, c in labels.items() if c == cluster]
        sub = graph.subgraph(members)
        roles = [graph.nodes[n].get("role", "") for n in members]
        hubs = sorted(members, key=lambda n: (-graph.degree(n), n))[:top]
        rows.append(
            {
                "cluster": cluster,
                "size": len(members),
                "seed_genes": roles.count("seed"),
                "diamond_genes": roles.count("diamond"),
                "internal_edges": sub.number_of_edges(),
                "density": nx.density(sub) if len(members) > 1 else 0.0,
                "top_proteins": ", ".join(graph.nodes[n].get("symbol", n) for n in hubs),
            }
        )
    return pd.DataFrame(rows)


def cluster_network(
    graph: nx.Graph,
    method: str = "louvain",
    resolution: float = 1.0,
    seed: int = 42,
    min_size: int = 5,
) -> tuple[dict[str, int], float]:
    """Cluster ``graph``; returns ``(labels, modularity)`` and sets node attribute ``cluster``."""
    communities = detect_communities(graph, method, resolution, seed)
    labels = membership(communities, min_size)
    modularity = nx.community.modularity(graph, communities, resolution=resolution)
    nx.set_node_attributes(graph, labels, "cluster")
    logger.info(
        "%d clusters (>= %d proteins), modularity %.3f",
        len({c for c in labels.values() if c > 0}), min_size, modularity,
    )
    return labels, modularity
