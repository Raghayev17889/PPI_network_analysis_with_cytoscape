"""Is there an Alzheimer's disease module in the interactome?

Network medicine (Menche et al., Science 2015) expects the products of
disease genes to cluster in one neighbourhood of the interactome, the
*disease module*. Here that is tested by comparing the largest connected
component (LCC) formed by the seed genes with the LCCs of random gene sets
of the same size and degree distribution (degree-preserving sampling, so hubs
are not over-represented by chance).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def lcc_size(graph: nx.Graph, genes: list[str] | set[str]) -> int:
    """Size of the largest connected component induced by ``genes``."""
    subgraph = graph.subgraph([g for g in genes if g in graph])
    if subgraph.number_of_nodes() == 0:
        return 0
    return max(len(c) for c in nx.connected_components(subgraph))


def degree_bins(graph: nx.Graph, min_bin_size: int = 100) -> dict[str, list[str]]:
    """Group proteins with similar degree so each group has ``min_bin_size`` members.

    Returns ``node -> list of nodes in the same degree bin``. Degrees are
    visited in increasing order and a bin is closed once it holds enough
    proteins; the highest-degree remainder joins the last bin.
    """
    by_degree: dict[int, list[str]] = {}
    for node, degree in graph.degree():
        by_degree.setdefault(degree, []).append(node)
    bins: list[list[str]] = []
    current: list[str] = []
    for degree in sorted(by_degree):
        current.extend(sorted(by_degree[degree]))
        if len(current) >= min_bin_size:
            bins.append(current)
            current = []
    if current:
        if bins:
            bins[-1].extend(current)
        else:
            bins.append(current)
    return {node: members for members in bins for node in members}


def degree_matched_sample(
    genes: list[str], bins: dict[str, list[str]], rng: np.random.Generator
) -> list[str]:
    """Random genes with the same degree profile as ``genes`` (no repeats)."""
    chosen: set[str] = set()
    for gene in genes:
        pool = bins[gene]
        for _ in range(50):
            candidate = pool[int(rng.integers(len(pool)))]
            if candidate not in chosen:
                break
        chosen.add(candidate)
    return list(chosen)


@dataclass
class ModuleSignificance:
    """Observed seed LCC compared with degree-matched random gene sets."""

    n_genes: int
    observed_lcc: int
    random_lcc: np.ndarray

    @property
    def mean(self) -> float:
        return float(self.random_lcc.mean())

    @property
    def std(self) -> float:
        return float(self.random_lcc.std(ddof=1))

    @property
    def z_score(self) -> float:
        return (self.observed_lcc - self.mean) / self.std if self.std > 0 else float("inf")

    @property
    def empirical_p(self) -> float:
        """One-sided p-value with the +1 correction (never exactly zero)."""
        return float((np.sum(self.random_lcc >= self.observed_lcc) + 1) / (len(self.random_lcc) + 1))

    def summary(self) -> dict[str, float]:
        return {
            "genes_in_interactome": self.n_genes,
            "observed_lcc": self.observed_lcc,
            "observed_lcc_fraction": self.observed_lcc / self.n_genes if self.n_genes else 0.0,
            "random_lcc_mean": self.mean,
            "random_lcc_std": self.std,
            "z_score": self.z_score,
            "empirical_p": self.empirical_p,
            "n_random": len(self.random_lcc),
        }


def module_significance(
    graph: nx.Graph,
    genes: list[str],
    n_random: int = 1000,
    seed: int = 42,
    min_bin_size: int = 100,
) -> ModuleSignificance:
    """Test whether ``genes`` form a larger connected module than chance."""
    genes = [g for g in dict.fromkeys(genes) if g in graph]
    rng = np.random.default_rng(seed)
    bins = degree_bins(graph, min_bin_size)
    observed = lcc_size(graph, genes)
    random_sizes = np.array(
        [lcc_size(graph, degree_matched_sample(genes, bins, rng)) for _ in range(n_random)], dtype=float
    )
    result = ModuleSignificance(len(genes), observed, random_sizes)
    logger.info("Disease module: %s", result.summary())
    return result
