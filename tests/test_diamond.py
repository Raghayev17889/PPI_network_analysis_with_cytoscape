from __future__ import annotations

import networkx as nx
import numpy as np
import pytest
from scipy.stats import hypergeom

from src.diamond import diamond, diamond_cross_validation, hypergeom_log_sf


def test_hypergeom_matches_scipy():
    rng = np.random.default_rng(0)
    draws = rng.integers(1, 60, 200)
    x = np.array([rng.integers(0, d + 1) for d in draws])
    ours = np.exp(hypergeom_log_sf(x, 5000, 80, draws))
    np.testing.assert_allclose(ours, hypergeom.sf(x - 1, 5000, 80, draws), rtol=1e-9, atol=1e-300)
    # impossible counts have probability zero
    assert hypergeom_log_sf(np.array([10]), 100, 5, np.array([3]))[0] == -np.inf


def _toy():
    graph = nx.path_graph(30)
    graph = nx.relabel_nodes(graph, str)
    seeds = ["s1", "s2", "s3"]
    graph.add_edges_from([("hub", s) for s in seeds] + [("weak", "s1"), ("weak", "0")])
    graph.add_edges_from([("s1", "1"), ("s2", "5"), ("s3", "9")])
    return graph, seeds


def test_diamond_picks_most_connected_first():
    graph, seeds = _toy()
    result = diamond(graph, seeds, n_genes=3)
    assert result.loc[0, "gene"] == "hub"
    assert result.loc[0, "links_to_module"] == 3
    assert result["rank"].tolist() == [1, 2, 3]
    assert result["p_value"].between(0, 1).all()
    assert not set(result["gene"]) & set(seeds)


def test_diamond_stops_when_component_is_exhausted():
    graph = nx.Graph([("a", "b"), ("b", "c")])
    graph.add_node("lonely")
    result = diamond(graph, ["a"], n_genes=10)
    assert result["gene"].tolist() == ["b", "c"]


def test_diamond_errors():
    graph, seeds = _toy()
    with pytest.raises(ValueError):
        diamond(graph, ["not-there"])
    with pytest.raises(ValueError):
        diamond(graph, seeds, alpha=0)
    assert diamond(graph, seeds, n_genes=2, alpha=3).loc[0, "gene"] == "hub"


def test_cross_validation_recovers_planted_module():
    graph = nx.relabel_nodes(nx.gnm_random_graph(300, 600, seed=3), str)
    module = [str(i) for i in range(25)]
    graph.add_edges_from((a, b) for i, a in enumerate(module) for b in module[i + 1 :] if (i + int(b)) % 2 == 0)
    cv = diamond_cross_validation(graph, module, n_genes=20, folds=5)
    assert len(cv) == 5
    assert cv["recovered"].sum() > 3 * cv["expected"].sum()
