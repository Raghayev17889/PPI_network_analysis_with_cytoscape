from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from src.topology import attack_curve, centrality_table, classify_hubs_bottlenecks, key_proteins, robustness_analysis


@pytest.fixture
def barbell():
    """Two 6-cliques joined through a single bridge protein 'B'."""
    graph = nx.barbell_graph(6, 1)
    graph = nx.relabel_nodes(graph, {6: "B", **{i: f"n{i}" for i in range(13) if i != 6}})
    nx.set_node_attributes(graph, {n: f"SYM_{n}" for n in graph}, "symbol")
    return graph


def test_centrality_table_matches_networkx(barbell):
    table = centrality_table(barbell).set_index("id")
    assert table.loc["B", "betweenness"] == pytest.approx(nx.betweenness_centrality(barbell)["B"])
    assert table.loc["n0", "clustering"] == 1.0
    assert table.loc["B", "symbol"] == "SYM_B"
    assert (table["eigenvector"] >= 0).all()


def test_bridge_is_a_non_hub_bottleneck(barbell):
    table = classify_hubs_bottlenecks(centrality_table(barbell), 0.9, 0.9).set_index("id")
    assert table.loc["B", "centrality_class"] == "non-hub bottleneck"
    assert table["is_bottleneck"].sum() >= 1
    assert set(table["centrality_class"]) <= {"hub-bottleneck", "non-hub bottleneck", "hub non-bottleneck", "peripheral"}
    assert table.attrs["hub_degree_threshold"] > 0
    with pytest.raises(ValueError):
        classify_hubs_bottlenecks(centrality_table(barbell), 1.0, 0.5)


def test_key_proteins_sorted(barbell):
    top = key_proteins(classify_hubs_bottlenecks(centrality_table(barbell)), n=3)
    assert top.loc[0, "symbol"] == "SYM_B"
    assert top["betweenness"].is_monotonic_decreasing


def test_attack_curve_and_robustness(barbell):
    curve = attack_curve(barbell, ["B"], fractions=np.array([0.0, 1 / 13]))
    assert curve["lcc_fraction"].tolist() == [1.0, pytest.approx(6 / 13)]
    table = centrality_table(barbell)
    curves = robustness_analysis(barbell, table, n_random=3)
    assert set(curves["strategy"]) == {"betweenness", "degree", "random"}
    targeted = curves[curves.strategy == "betweenness"].set_index("fraction_removed")["lcc_fraction"]
    random = curves[curves.strategy == "random"].set_index("fraction_removed")["lcc_fraction"]
    assert (targeted <= random + 1e-9).all()
