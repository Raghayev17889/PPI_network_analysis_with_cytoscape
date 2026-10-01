from __future__ import annotations

import networkx as nx
import numpy as np

from src.disease_module import degree_bins, degree_matched_sample, lcc_size, module_significance


def _graph_with_clique(seed=0):
    graph = nx.relabel_nodes(nx.gnm_random_graph(400, 800, seed=seed), str)
    clique = [str(i) for i in range(20)]
    graph.add_edges_from((a, b) for i, a in enumerate(clique) for b in clique[i + 1 :])
    return graph, clique


def test_lcc_size():
    graph = nx.Graph([("a", "b"), ("b", "c"), ("x", "y")])
    assert lcc_size(graph, ["a", "b", "c", "x", "y", "missing"]) == 3
    assert lcc_size(graph, ["a", "x"]) == 1
    assert lcc_size(graph, []) == 0


def test_degree_bins_cover_every_node():
    graph, _ = _graph_with_clique()
    bins = degree_bins(graph, min_bin_size=50)
    assert set(bins) == set(graph)
    assert all(len(members) >= 50 for members in bins.values())
    for node, members in bins.items():
        assert node in members


def test_degree_matched_sample():
    graph, clique = _graph_with_clique()
    bins = degree_bins(graph, 50)
    sample = degree_matched_sample(clique, bins, np.random.default_rng(1))
    assert len(sample) == len(set(sample)) == len(clique)
    allowed = {node for c in clique for node in bins[c]}
    assert set(sample) <= allowed  # every draw comes from a seed's degree bin


def test_clique_is_a_significant_module():
    graph, clique = _graph_with_clique()
    result = module_significance(graph, clique, n_random=100, min_bin_size=30)
    assert result.observed_lcc == 20
    assert result.z_score > 3
    assert result.empirical_p < 0.05
    summary = result.summary()
    assert summary["n_random"] == 100 and summary["observed_lcc_fraction"] == 1.0
