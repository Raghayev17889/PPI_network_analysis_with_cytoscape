"""End-to-end pipeline: from Hetionet to a styled Cytoscape session.

Steps
-----
1. Download Hetionet and extract the disease dataset (``src.hetionet``).
2. Build the human interactome, excluding ubiquitin (``src.interactome``).
3. Test whether the disease genes form a connected module (``src.disease_module``).
4. Rank candidate module genes with DIAMOnD and cross-validate (``src.diamond``).
5. Build and filter the disease network (``src.disease_network``).
6. Hubs, bottlenecks and robustness (``src.topology``).
7. Louvain clusters and their stability (``src.clustering``).
8. Pathway / GO enrichment of every cluster (``src.enrichment``).
9. Tables and figures in ``results/``.
10. Cytoscape: network files, visual style and - if Cytoscape is running -
    styled views, cluster subnetworks, exported images and a ``.cys`` session
    in ``cytoscape/`` (``src.cytoscape_viz``).

Usage::

    python src/pipeline.py                 # full run (start Cytoscape first)
    python src/pipeline.py --no-cytoscape  # analysis + offline network files only
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

import numpy as np
import pandas as pd

from src import cytoscape_viz as cy
from src import plots
from src.clustering import cluster_network, cluster_summary, partition_stability
from src.config import ROOT, load_config
from src.diamond import diamond, diamond_cross_validation
from src.disease_module import module_significance
from src.disease_network import build_disease_network, export_tables
from src.enrichment import cluster_labels, enrich, enrich_clusters
from src.hetionet import prepare
from src.interactome import build_interactome, interactome_summary, map_genes, symbols_to_entrez
from src.topology import centrality_table, classify_hubs_bottlenecks, key_proteins, robustness_analysis

logger = logging.getLogger(__name__)


def _relative(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def run(config: dict, use_cytoscape: bool = True, cross_validate: bool = True) -> dict:
    """Run every step; returns a JSON-serialisable summary."""
    t0 = time.perf_counter()
    seed = int(config.get("random_seed", 42))
    results_dir = Path(config["output"]["results_dir"])
    tables, figures = results_dir / "tables", results_dir / "figures"
    cyto_dir = Path(config["output"]["cytoscape_dir"])
    for directory in (tables, figures, cyto_dir):
        directory.mkdir(parents=True, exist_ok=True)
    summary: dict = {"disease": f"{config['disease']['name']} ({config['disease']['id']})"}

    # 1-2. Data and interactome --------------------------------------------
    data = prepare(config)
    symbols = data.symbols
    summary["hetionet"] = data.summary()
    excluded = symbols_to_entrez(config["interactome"].get("exclude_symbols", []), symbols)
    interactome = build_interactome(
        data.interactions, symbols, config["interactome"].get("largest_component", True), excluded
    )
    summary["interactome"] = interactome_summary(interactome)

    # 3. Disease module -----------------------------------------------------
    mapping = map_genes(interactome, data.seed_genes)
    summary["seed_mapping"] = {
        "seed_genes": len(mapping.mapped) + len(mapping.missing),
        "in_interactome": len(mapping.mapped),
        "missing": [symbols.get(g, g) for g in mapping.missing],
    }
    module_cfg = config["disease_module"]
    significance = module_significance(
        interactome, mapping.mapped, module_cfg["random_sets"], seed, module_cfg["min_degree_bin"]
    )
    summary["module_significance"] = significance.summary()
    plots.save(
        plots.plot_module_significance(
            significance.random_lcc, significance.observed_lcc, significance.z_score, significance.n_genes
        ),
        figures / "module_significance.png",
    )

    # 4. DIAMOnD ------------------------------------------------------------
    dcfg = config["diamond"]
    ranking = diamond(interactome, mapping.mapped, dcfg["n_genes"], dcfg["alpha"])
    ranking.insert(2, "symbol", ranking["gene"].map(symbols))
    ranking.to_csv(tables / "diamond_ranking.tsv", sep="\t", index=False)
    if cross_validate:
        cv = diamond_cross_validation(interactome, mapping.mapped, dcfg["n_genes"], dcfg["cv_folds"], seed)
        cv.drop(columns="recovered_ranks").to_csv(tables / "diamond_cross_validation.tsv", sep="\t", index=False)
        summary["diamond_cross_validation"] = {
            "recovered": int(cv["recovered"].sum()),
            "expected": float(cv["expected"].sum()),
            "fold_enrichment": float(cv["recovered"].sum() / cv["expected"].sum()),
        }
        plots.save(plots.plot_diamond_recovery(cv, dcfg["n_genes"]), figures / "diamond_recovery.png")

    # 5. Disease network ----------------------------------------------------
    network = build_disease_network(interactome, data, ranking, dcfg["network_genes"])
    graph = network.graph
    summary["filtering"] = network.log_frame().reset_index().to_dict(orient="records")

    # 6. Topology -----------------------------------------------------------
    tcfg = config["topology"]
    centrality = classify_hubs_bottlenecks(centrality_table(graph), tcfg["hub_quantile"], tcfg["bottleneck_quantile"])
    curves = robustness_analysis(graph, centrality, tcfg["random_attacks"], seed)
    curves.to_csv(tables / "robustness.tsv", sep="\t", index=False)
    plots.save(plots.plot_robustness(curves), figures / "robustness.png")

    # 7. Clusters -----------------------------------------------------------
    ccfg = config["clustering"]
    labels, modularity = cluster_network(graph, ccfg["method"], ccfg["resolution"], seed, ccfg["min_cluster_size"])
    stability = partition_stability(graph, labels, ccfg["stability_runs"], ccfg["resolution"], seed)
    clusters = cluster_summary(graph, labels)

    # 8. Enrichment ---------------------------------------------------------
    ecfg = config["enrichment"]
    annotated = set().union(*data.pathways.values()) | set().union(*data.go_processes.values())
    background = set(interactome) & annotated
    collections = {"pathway": data.pathways, "go_bp": data.go_processes}
    enrichment = enrich_clusters(
        labels, collections, background, symbols, ecfg["fdr"], ecfg["min_set_size"], ecfg["max_set_size"]
    )
    enrichment.to_csv(tables / "cluster_enrichment.tsv", sep="\t", index=False)
    whole = enrich(list(graph), data.pathways, background, ecfg["min_set_size"], ecfg["max_set_size"], symbols)
    whole[whole["fdr"] <= ecfg["fdr"]].to_csv(tables / "network_pathway_enrichment.tsv", sep="\t", index=False)
    auto = cluster_labels(enrichment, ecfg.get("label_source", "pathway"))
    curated = {int(k): v for k, v in (config.get("cluster_names") or {}).items()}
    names = {c: curated.get(c, auto.get(c, f"cluster {c}")) for c in sorted(set(labels.values())) if c > 0}
    clusters["name"] = clusters["cluster"].map(lambda c: names.get(c, "other"))
    clusters["top_pathway"] = clusters["cluster"].map(lambda c: auto.get(c, ""))
    clusters.to_csv(tables / "clusters.tsv", sep="\t", index=False)
    plots.save(plots.plot_cluster_enrichment(enrichment, names), figures / "cluster_enrichment.png")

    # 9. Node table, key proteins, figures ---------------------------------
    nodes_table = network.node_table().merge(
        centrality.drop(columns=["symbol"]), on="id"
    ).assign(cluster=lambda t: t["id"].map(labels))
    nodes_table["cluster_name"] = nodes_table["cluster"].map(lambda c: names.get(c, "other"))
    nodes_table.sort_values("betweenness", ascending=False).to_csv(tables / "node_metrics.tsv", sep="\t", index=False)
    key = key_proteins(nodes_table, 25)
    key.to_csv(tables / "key_proteins.tsv", sep="\t", index=False)
    centrality_plot = centrality.copy()
    plots.save(plots.plot_centrality_scatter(centrality_plot), figures / "hubs_bottlenecks.png")
    export_tables(network, tables)

    summary["network"] = {
        "proteins": graph.number_of_nodes(),
        "interactions": graph.number_of_edges(),
        "seed_genes": int((nodes_table["role"] == "seed").sum()),
        "diamond_genes": int((nodes_table["role"] == "diamond").sum()),
        "ad_drug_targets": int(nodes_table["ad_drug_target"].sum()),
        "druggable": int(nodes_table["druggable"].sum()),
    }
    summary["topology"] = {
        "hub_degree_threshold": centrality.attrs["hub_degree_threshold"],
        "bottleneck_betweenness_threshold": centrality.attrs["bottleneck_betweenness_threshold"],
        "classes": centrality["centrality_class"].value_counts().to_dict(),
        "top_bottlenecks": key["symbol"].head(10).tolist(),
    }
    summary["clustering"] = {
        "modularity": modularity,
        "clusters": {int(r.cluster): {"name": r.name, "size": int(r.size)} for r in clusters.itertuples()},
        "stability_ari_mean": float(stability["adjusted_rand_index"].mean()),
        "stability_ari_min": float(stability["adjusted_rand_index"].min()),
    }

    # 10. Cytoscape ---------------------------------------------------------
    cy_nodes, cy_edges = cy.cytoscape_tables(graph, centrality, names)
    positions = cy.clustered_layout(graph, labels)
    files = cy.write_offline_files(cy_nodes, cy_edges, cyto_dir / "network", positions)
    style_cap = float(centrality["betweenness"].quantile(0.98))
    cy.write_style_json(cyto_dir / "styles" / "alzheimer_style.json", sorted(set(labels.values())), style_cap)
    summary["cytoscape"] = {"offline_files": {k: _relative(v) for k, v in files.items()}}

    base_url = config.get("cytoscape", {}).get("base_url", cy.BASE_URL)
    if use_cytoscape and cy.cytoscape_available(base_url):
        summary["cytoscape"].update(
            build_cytoscape_session(graph, cy_nodes, cy_edges, labels, names, positions, style_cap, cyto_dir,
                                    base_url)
        )
    elif use_cytoscape:
        logger.warning("Cytoscape is not reachable at %s: skipped live session (offline files written)", base_url)

    summary["runtime_seconds"] = round(time.perf_counter() - t0, 1)
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=_json_default), "utf-8")
    logger.info("Done in %.0fs. Results in %s", summary["runtime_seconds"], results_dir)
    return summary


def build_cytoscape_session(graph, nodes, edges, labels, names, positions, style_cap, cyto_dir, base_url) -> dict:
    """Load, style, lay out and export every Cytoscape view; save the session."""
    import py4cytoscape as p4c

    p4c.close_session(save_before_closing=False, base_url=base_url)
    clusters = sorted(set(labels.values()))
    figures = cyto_dir / "figures"

    main = cy.load_network(nodes, edges, "Alzheimer's disease network", "Alzheimer's disease PPI", base_url)
    cy.apply_style(main, clusters, style_cap, base_url=base_url)
    cy.apply_positions(positions, main, base_url)
    cy.add_legend(main, positions, names, base_url=base_url)
    exported = cy.export_view(main, figures / "alzheimer_network", base_url=base_url)

    for cluster in (c for c in clusters if c > 0):
        members = [n for n, c in labels.items() if c == cluster]
        sub = cy.create_subnetwork(main, members, f"Cluster {cluster}: {names[cluster]}", graph, labels,
                                   base_url=base_url)
        cy.show_all_labels(sub, nodes[nodes["id"].isin(members)], base_url)
        exported += cy.export_view(sub, figures / f"cluster_{cluster}", formats=("PNG",), base_url=base_url)

    key = nodes[(nodes["centrality_class"] != "peripheral") | (nodes["drug_status"] == "AD drug target")]
    sub = cy.create_subnetwork(main, key["id"].tolist(), "Hubs, bottlenecks and drug targets", graph, labels,
                               base_url=base_url)
    cy.show_all_labels(sub, key, base_url)
    exported += cy.export_view(sub, figures / "key_proteins", formats=("PNG", "SVG"), base_url=base_url)

    p4c.set_current_network(main, base_url=base_url)
    style_xml = cy.export_style_xml(cyto_dir / "styles" / "alzheimer_style.xml", base_url=base_url)
    session = cy.save_session(cyto_dir / "alzheimer_network.cys", base_url)
    logger.info("Cytoscape session saved to %s", session)
    return {
        "session": _relative(session),
        "style_xml": _relative(style_xml),
        "figures": [_relative(p) for p in exported],
        "cytoscape_version": p4c.cytoscape_version_info(base_url=base_url).get("cytoscapeVersion"),
    }


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return _relative(value)
    raise TypeError(f"Not JSON serialisable: {type(value)}")


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description="Alzheimer's disease PPI network analysis with Cytoscape.")
    parser.add_argument("--config", type=Path, help="YAML configuration (default: config/config.yaml)")
    parser.add_argument("--no-cytoscape", action="store_true", help="skip the live Cytoscape session")
    parser.add_argument("--skip-cv", action="store_true", help="skip DIAMOnD cross-validation (faster)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    matplotlib.use("Agg")
    return run(load_config(args.config), use_cytoscape=not args.no_cytoscape, cross_validate=not args.skip_cv)


if __name__ == "__main__":
    main()
