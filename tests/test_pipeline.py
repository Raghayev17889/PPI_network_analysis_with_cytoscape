"""End-to-end run of the pipeline on the miniature Hetionet (no downloads, no Cytoscape)."""

from __future__ import annotations

import json

import pandas as pd

from src.hetionet import extract_disease_data, load_edges, load_nodes, save_disease_data
from src.pipeline import run


def _config(tmp_path):
    return {
        "disease": {"id": "DOID:0001", "name": "toy disease"},
        "data": {"raw_dir": tmp_path / "raw", "processed_dir": tmp_path / "processed"},
        "interactome": {"exclude_symbols": [], "largest_component": True},
        "disease_module": {"random_sets": 20, "min_degree_bin": 2},
        "diamond": {"n_genes": 2, "alpha": 1, "network_genes": 2, "cv_folds": 2},
        "topology": {"hub_quantile": 0.8, "bottleneck_quantile": 0.8, "random_attacks": 2},
        "clustering": {"method": "louvain", "resolution": 1.0, "min_cluster_size": 1, "stability_runs": 2},
        "enrichment": {"fdr": 1.0, "min_set_size": 1, "max_set_size": 500, "label_source": "pathway"},
        "cluster_names": {1: "First"},
        "output": {"results_dir": tmp_path / "results", "cytoscape_dir": tmp_path / "cytoscape"},
        "random_seed": 1,
    }


def test_pipeline_end_to_end(tmp_path, hetionet_files):
    nodes, edges = hetionet_files
    data = extract_disease_data(load_nodes(nodes), load_edges(edges), "DOID:0001")
    # Pre-populate the processed cache so prepare() does not download anything.
    save_disease_data(data, tmp_path / "processed" / "DOID_0001")

    summary = run(_config(tmp_path), use_cytoscape=False, cross_validate=True)

    results = tmp_path / "results"
    for table in ("diamond_ranking", "node_metrics", "key_proteins", "clusters", "cluster_enrichment",
                  "ad_network_nodes", "ad_network_edges", "ad_network_filtering", "robustness"):
        assert (results / "tables" / f"{table}.tsv").exists(), table
    for figure in ("module_significance", "diamond_recovery", "robustness", "hubs_bottlenecks", "cluster_enrichment"):
        assert (results / "figures" / f"{figure}.png").exists(), figure
    for name in ("network/alzheimer_network.cyjs", "network/alzheimer_network.graphml", "styles/alzheimer_style.json"):
        assert (tmp_path / "cytoscape" / name).exists(), name

    saved = json.loads((results / "summary.json").read_text())
    assert saved["network"]["proteins"] == summary["network"]["proteins"] == 4
    assert saved["seed_mapping"]["in_interactome"] == 2
    clusters = pd.read_csv(results / "tables" / "clusters.tsv", sep="\t")
    assert clusters.loc[clusters["cluster"] == 1, "name"].item() == "First"
