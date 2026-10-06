from __future__ import annotations

import json

import networkx as nx
import pandas as pd
import pytest

from src.cytoscape_viz import (
    CLUSTER_COLORS,
    OTHER_COLOR,
    build_style,
    cluster_color,
    cytoscape_tables,
    to_cytoscape_js,
    write_offline_files,
    write_style_json,
)
from src.topology import centrality_table, classify_hubs_bottlenecks


@pytest.fixture
def annotated_graph():
    graph = nx.barbell_graph(5, 1)
    graph = nx.relabel_nodes(graph, str)
    for n in graph:
        i = int(n)
        graph.nodes[n].update(
            symbol=f"G{n}",
            role="seed" if i % 2 else "diamond",
            cluster=1 if i < 5 else 2 if i > 5 else 0,
            ad_drug_target=(n == "0"),
            druggable=(n in {"0", "1"}),
        )
    for u, v in graph.edges():
        graph.edges[u, v]["edge_type"] = "-".join(sorted((graph.nodes[u]["role"], graph.nodes[v]["role"])))
    return graph


def test_cluster_colors():
    assert cluster_color(1) == CLUSTER_COLORS[0]
    assert cluster_color(0) == OTHER_COLOR
    assert cluster_color(9) == CLUSTER_COLORS[0]


def test_cytoscape_tables(annotated_graph):
    table = classify_hubs_bottlenecks(centrality_table(annotated_graph))
    nodes, edges = cytoscape_tables(annotated_graph, table, {1: "Alpha", 2: "Beta"})
    nodes = nodes.set_index("id")
    assert nodes.loc["0", "drug_status"] == "AD drug target"
    assert nodes.loc["1", "drug_status"] == "druggable"
    assert nodes.loc["2", "drug_status"] == "none"
    assert nodes.loc["0", "cluster_name"] == "1. Alpha"
    assert nodes.loc["5", "cluster_name"] == "0. other"
    assert nodes.loc["5", "label"] == "G5"     # the bridge is a bottleneck
    assert nodes.loc["0", "label"] == "G0"     # AD drug targets are always labelled
    assert "betweenness" in nodes.columns and "G5" in nodes.loc["5", "tooltip"]
    assert {"source", "target", "interaction", "edge_type"} <= set(edges.columns)
    assert len(edges) == annotated_graph.number_of_edges()


def test_build_style_structure(tmp_path):
    defaults, mappings = build_style([0, 1, 2], max_betweenness=0.4)
    props = {m["visualProperty"]: m for m in mappings}
    assert props["NODE_FILL_COLOR"]["map"] == [
        {"key": "0", "value": OTHER_COLOR}, {"key": "1", "value": CLUSTER_COLORS[0]},
        {"key": "2", "value": CLUSTER_COLORS[1]},
    ]
    assert props["NODE_SHAPE"]["mappingColumn"] == "role"
    sizes = props["NODE_SIZE"]["points"]
    assert sizes[-1]["value"] == pytest.approx(0.4) and float(sizes[0]["equal"]) < float(sizes[-1]["equal"])
    assert props["NODE_LABEL"]["mappingType"] == "passthrough"
    assert any(d["visualProperty"] == "NETWORK_BACKGROUND_PAINT" for d in defaults)
    saved = json.loads(write_style_json(tmp_path / "style.json", [1], 0.1).read_text())
    assert saved["title"] == "Alzheimer PPI" and saved["mappings"]


def test_offline_exports(tmp_path, annotated_graph):
    table = classify_hubs_bottlenecks(centrality_table(annotated_graph))
    nodes, edges = cytoscape_tables(annotated_graph, table)
    positions = {n: (float(i), 2.0 * i) for i, n in enumerate(annotated_graph)}
    cyjs = to_cytoscape_js(nodes, edges, positions)
    assert len(cyjs["elements"]["nodes"]) == annotated_graph.number_of_nodes()
    assert cyjs["elements"]["nodes"][1]["position"] == {"x": 1.0, "y": 2.0}
    assert all("id" in e["data"] for e in cyjs["elements"]["edges"])

    files = write_offline_files(nodes, edges, tmp_path)
    reloaded = nx.read_graphml(files["graphml"])
    assert reloaded.number_of_edges() == annotated_graph.number_of_edges()
    assert reloaded.nodes["0"]["drug_status"] == "AD drug target"
    assert json.loads(files["cyjs"].read_text())["elements"]["nodes"]
    assert isinstance(pd.read_json(files["cyjs"], typ="series"), pd.Series)


def test_clustered_layout_is_deterministic_and_separates_clusters(annotated_graph):
    from src.cytoscape_viz import clustered_layout

    labels = nx.get_node_attributes(annotated_graph, "cluster")
    first = clustered_layout(annotated_graph, labels)
    assert first == clustered_layout(annotated_graph, labels)
    assert set(first) == set(annotated_graph)

    def centroid(cluster):
        pts = [first[n] for n, c in labels.items() if c == cluster]
        return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)

    (x1, y1), (x2, y2) = centroid(1), centroid(2)
    assert ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5 > 100
