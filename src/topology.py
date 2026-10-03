"""Centrality analysis: hubs, bottlenecks and network robustness.

Following Yu et al. (*PLoS Comput Biol* 3, e59, 2007), proteins are
classified by two independent properties:

* **hub** - many interaction partners (high degree);
* **bottleneck** - many shortest paths pass through it (high betweenness).

Hub-bottlenecks are the most central proteins: removing them fragments the
network. Non-hub bottlenecks connect modules without many partners of their
own and are easy to miss when ranking by degree alone.
"""

from __future__ import annotations

import logging

import networkx as nx
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

CLASS_ORDER = ["hub-bottleneck", "non-hub bottleneck", "hub non-bottleneck", "peripheral"]


def centrality_table(graph: nx.Graph) -> pd.DataFrame:
    """Degree, betweenness, closeness, eigenvector, clustering and k-core per protein."""
    nodes = list(graph.nodes())
    table = pd.DataFrame(
        {
            "id": nodes,
            "symbol": [graph.nodes[n].get("symbol", n) for n in nodes],
            "degree": [graph.degree(n) for n in nodes],
            "betweenness": pd.Series(nx.betweenness_centrality(graph, normalized=True)).reindex(nodes).to_numpy(),
            "closeness": pd.Series(nx.closeness_centrality(graph)).reindex(nodes).to_numpy(),
            "eigenvector": pd.Series(nx.eigenvector_centrality_numpy(graph)).reindex(nodes).abs().to_numpy(),
            "clustering": pd.Series(nx.clustering(graph)).reindex(nodes).to_numpy(),
            "k_core": pd.Series(nx.core_number(graph)).reindex(nodes).to_numpy(),
        }
    )
    return table


def classify_hubs_bottlenecks(
    table: pd.DataFrame, hub_quantile: float = 0.9, bottleneck_quantile: float = 0.9
) -> pd.DataFrame:
    """Add ``is_hub``, ``is_bottleneck`` and ``centrality_class`` columns.

    Thresholds are quantiles of the network's own degree and betweenness
    distributions (``0.9`` = top 10%).
    """
    if not (0 < hub_quantile < 1 and 0 < bottleneck_quantile < 1):
        raise ValueError("quantiles must be in (0, 1)")
    table = table.copy()
    hub_cut = table["degree"].quantile(hub_quantile)
    bottleneck_cut = table["betweenness"].quantile(bottleneck_quantile)
    table["is_hub"] = table["degree"] >= hub_cut
    table["is_bottleneck"] = table["betweenness"] >= bottleneck_cut
    table["centrality_class"] = np.select(
        [
            table["is_hub"] & table["is_bottleneck"],
            ~table["is_hub"] & table["is_bottleneck"],
            table["is_hub"] & ~table["is_bottleneck"],
        ],
        CLASS_ORDER[:3],
        default=CLASS_ORDER[3],
    )
    table.attrs["hub_degree_threshold"] = float(hub_cut)
    table.attrs["bottleneck_betweenness_threshold"] = float(bottleneck_cut)
    return table


def key_proteins(table: pd.DataFrame, n: int = 20) -> pd.DataFrame:
    """Top proteins by betweenness (ties broken by degree)."""
    columns = [c for c in ["symbol", "role", "centrality_class", "degree", "betweenness", "closeness",
                           "cluster", "ad_drug_target", "druggable"] if c in table.columns]
    return table.sort_values(["betweenness", "degree"], ascending=False).head(n)[columns].reset_index(drop=True)


def attack_curve(graph: nx.Graph, order: list[str], fractions: np.ndarray | None = None) -> pd.DataFrame:
    """Relative size of the largest component as proteins are removed in ``order``."""
    fractions = np.linspace(0, 0.3, 31) if fractions is None else np.asarray(fractions)
    n = graph.number_of_nodes()
    rows = []
    for fraction in fractions:
        removed = order[: int(round(fraction * n))]
        remaining = graph.copy()
        remaining.remove_nodes_from(removed)
        giant = max((len(c) for c in nx.connected_components(remaining)), default=0)
        rows.append({"fraction_removed": float(fraction), "lcc_fraction": giant / n})
    return pd.DataFrame(rows)


def robustness_analysis(graph: nx.Graph, table: pd.DataFrame, n_random: int = 20, seed: int = 42) -> pd.DataFrame:
    """Compare targeted removal (by betweenness, by degree) with random removal."""
    by_betweenness = table.sort_values("betweenness", ascending=False)["id"].tolist()
    by_degree = table.sort_values("degree", ascending=False)["id"].tolist()
    curves = [
        attack_curve(graph, by_betweenness).assign(strategy="betweenness"),
        attack_curve(graph, by_degree).assign(strategy="degree"),
    ]
    rng = np.random.default_rng(seed)
    nodes = list(graph.nodes())
    random_curves = [attack_curve(graph, list(rng.permutation(nodes))) for _ in range(n_random)]
    random_mean = random_curves[0][["fraction_removed"]].copy()
    random_mean["lcc_fraction"] = np.mean([c["lcc_fraction"] for c in random_curves], axis=0)
    curves.append(random_mean.assign(strategy="random"))
    return pd.concat(curves, ignore_index=True)
