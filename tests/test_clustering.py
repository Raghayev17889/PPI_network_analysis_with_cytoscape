from __future__ import annotations

import networkx as nx
import pytest

from src.clustering import (
    adjusted_rand_index,
    cluster_network,
    cluster_summary,
    detect_communities,
    membership,
    partition_stability,
)


@pytest.fixture
def two_cliques():
    graph = nx.disjoint_union(nx.complete_graph(8), nx.complete_graph(6))
    graph.add_edge(0, 8)
    graph.add_edge(20, 21)  # a separate tiny pair
    graph = nx.relabel_nodes(graph, str)
    nx.set_node_attributes(graph, {n: ("seed" if int(n) % 2 else "diamond") for n in graph}, "role")
    return graph


def test_detects_the_two_cliques(two_cliques):
    communities = detect_communities(two_cliques)
    assert communities[0] == {str(i) for i in range(8)}
    assert communities[1] == {str(i) for i in range(8, 14)}
    assert detect_communities(two_cliques, "greedy")[0] == communities[0]
    with pytest.raises(ValueError):
        detect_communities(two_cliques, "spectral")


def test_membership_merges_small_clusters(two_cliques):
    labels = membership(detect_communities(two_cliques), min_size=3)
    assert labels["0"] == 1 and labels["8"] == 2 and labels["20"] == 0 == labels["21"]


def test_adjusted_rand_index():
    a = {"a": 1, "b": 1, "c": 2, "d": 2}
    assert adjusted_rand_index(a, a) == pytest.approx(1.0)
    assert adjusted_rand_index(a, {"a": 7, "b": 7, "c": 3, "d": 3}) == pytest.approx(1.0)  # label names don't matter
    assert adjusted_rand_index(a, {"a": 1, "b": 2, "c": 1, "d": 2}) < 0.1


def test_cluster_network_and_summary(two_cliques):
    labels, modularity = cluster_network(two_cliques, min_size=3)
    assert modularity > 0.3
    assert two_cliques.nodes["3"]["cluster"] == 1
    summary = cluster_summary(two_cliques, labels).set_index("cluster")
    assert summary.loc[1, "size"] == 8 and summary.loc[1, "density"] == 1.0
    assert summary.loc[1, "seed_genes"] + summary.loc[1, "diamond_genes"] == 8
    stability = partition_stability(two_cliques, labels, runs=3)
    assert (stability["adjusted_rand_index"] > 0.9).all()


def test_partition_is_identical_across_processes(tmp_path):
    """String hashing differs per process; the partition must not."""
    import subprocess
    import sys

    script = (
        "import networkx as nx; from src.clustering import detect_communities;"
        "g = nx.relabel_nodes(nx.les_miserables_graph(), lambda n: 'P_' + n);"
        "print(sorted(sorted(c) for c in detect_communities(g)))"
    )
    outputs = {
        subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True,
                       env={"PYTHONHASHSEED": str(seed), "PYTHONPATH": "."}).stdout
        for seed in (1, 2, 3)
    }
    assert len(outputs) == 1
