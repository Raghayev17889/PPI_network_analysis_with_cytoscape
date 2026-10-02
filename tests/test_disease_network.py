from __future__ import annotations

import pandas as pd
import pytest

from src.diamond import diamond
from src.disease_network import build_disease_network, export_tables
from src.hetionet import extract_disease_data, load_edges, load_nodes
from src.interactome import build_interactome


@pytest.fixture
def toy_network(hetionet_files):
    nodes, edges = hetionet_files
    data = extract_disease_data(load_nodes(nodes), load_edges(edges), "DOID:0001")
    interactome = build_interactome(data.interactions, data.symbols)  # LCC = genes 1-5
    ranked = diamond(interactome, data.seed_genes, n_genes=3)
    return build_disease_network(interactome, data, ranked, n_diamond=2), ranked


def test_filter_log(toy_network):
    network, _ = toy_network
    log = network.log_frame()
    assert log.index.tolist() == [
        "seed genes", "in interactome", "module expansion", "induced subgraph", "largest component",
    ]
    assert log.loc["seed genes", "proteins"] == 3       # genes 1, 2, 6
    assert log.loc["in interactome", "proteins"] == 2   # gene 6 is outside the LCC
    assert log.loc["module expansion", "proteins"] == 4
    assert log.loc["largest component", "proteins"] == network.graph.number_of_nodes()


def test_node_annotations(toy_network):
    network, ranked = toy_network
    nodes = network.node_table().set_index("id")
    assert set(nodes.index) == {"1", "2", "3", "4"}
    assert nodes.loc["1", "role"] == "seed" and nodes.loc["1", "diamond_rank"] == 0
    assert nodes.loc["3", "role"] == "diamond" and nodes.loc["3", "diamond_rank"] == 1
    assert nodes.loc["3", "symbol"] == "G3"
    assert nodes.loc["4", "expression"] == "up"
    assert bool(nodes.loc["2", "ad_drug_target"]) and nodes.loc["2", "ad_drugs"] == "Drugamab"
    assert bool(nodes.loc["2", "druggable"]) and not bool(nodes.loc["1", "druggable"])
    assert nodes.loc["3", "interactome_degree"] == 3


def test_edges_and_export(tmp_path, toy_network):
    network, _ = toy_network
    edges = network.edge_table()
    assert {"source", "target", "interaction", "edge_type"} <= set(edges.columns)
    assert set(edges["edge_type"]) <= {"seed-seed", "diamond-seed", "diamond-diamond"}
    paths = export_tables(network, tmp_path)
    assert pd.read_csv(paths["nodes"], sep="\t").shape[0] == 4
    assert pd.read_csv(paths["edges"], sep="\t").shape[0] == network.graph.number_of_edges()
