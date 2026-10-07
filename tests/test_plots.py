from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd

from src import plots
from src.topology import centrality_table, classify_hubs_bottlenecks, robustness_analysis


def test_all_figures_render(tmp_path):
    graph = nx.relabel_nodes(nx.barbell_graph(6, 2), str)
    nx.set_node_attributes(graph, {n: f"G{n}" for n in graph}, "symbol")
    table = classify_hubs_bottlenecks(centrality_table(graph))
    cv = pd.DataFrame(
        {
            "held_out_seeds": [10, 10],
            "background": [1000, 1000],
            "recovered": [3, 2],
            "expected": [0.4, 0.4],
            "recovered_ranks": [[1, 5, 30], [2, 40]],
        }
    )
    enrichment = pd.DataFrame(
        {
            "cluster": [1, 1, 2],
            "source": ["pathway"] * 3,
            "term": ["Apoptosis [PC7_1]", "A very long pathway name " * 4 + "[X]", "Synapse [PC7_2]"],
            "p_value": [1e-8, 1e-4, 1e-6],
            "fdr": [1e-6, 1e-3, 1e-5],
        }
    )
    figures = [
        plots.plot_module_significance(np.random.default_rng(0).normal(80, 9, 200).round(), 135, 6.0, 183),
        plots.plot_diamond_recovery(cv, 50),
        plots.plot_centrality_scatter(table, n_labels=3),
        plots.plot_robustness(robustness_analysis(graph, table, n_random=2)),
        plots.plot_cluster_enrichment(enrichment, {1: "Apoptosis", 2: "Synapse"}),
    ]
    for i, fig in enumerate(figures):
        path = plots.save(fig, tmp_path / f"figure_{i}.png")
        assert path.stat().st_size > 1000
