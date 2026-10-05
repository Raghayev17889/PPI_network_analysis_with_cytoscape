"""Cytoscape automation with py4cytoscape.

Cytoscape desktop (3.10+) must be running; py4cytoscape talks to its REST API
(CyREST) on ``http://127.0.0.1:1234/v1``. The visual style is defined here as
plain CyREST JSON, so it is version-controlled, testable without Cytoscape and
identical on every run.

Visual encoding of the Alzheimer's disease network
---------------------------------------------------
=====================  =====================================================
node fill colour       functional cluster (Louvain), fixed palette order
node shape             ellipse = known disease gene, diamond = DIAMOnD gene
node size              betweenness centrality (bottlenecks are large)
node border            red = target of an approved AD drug, dark grey =
                       targeted by any DrugBank compound
node label             gene symbol of hubs, bottlenecks and AD drug targets
edge colour / opacity  seed-seed interactions darker than the rest
=====================  =====================================================

Without a running Cytoscape the network can still be exported as GraphML and
Cytoscape.js JSON (:func:`write_offline_files`) and opened in Cytoscape via
*File > Import > Network from File*.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)

STYLE_NAME = "Alzheimer PPI"
BASE_URL = "http://127.0.0.1:1234/v1"

#: Categorical palette (fixed order, colour-vision-deficiency checked) for clusters 1..8.
CLUSTER_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER_COLOR = "#c3c2b7"
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781", "surface": "#fcfcfb"}
DRUG_BORDER = {"AD drug target": ("#d03b3b", 7.0), "druggable": ("#52514e", 3.0), "none": ("#fcfcfb", 1.5)}


def cluster_color(cluster: int) -> str:
    return CLUSTER_COLORS[(cluster - 1) % len(CLUSTER_COLORS)] if cluster > 0 else OTHER_COLOR


# ---------------------------------------------------------------------------
# Node / edge tables for Cytoscape
# ---------------------------------------------------------------------------


def cytoscape_tables(
    graph: nx.Graph,
    centrality: pd.DataFrame,
    cluster_names: dict[int, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Node and edge DataFrames in the layout ``create_network_from_data_frames`` expects.

    Adds display columns: ``label`` (symbol for key proteins, empty otherwise),
    ``drug_status``, ``cluster_name`` and a ``tooltip``.
    """
    nodes = pd.DataFrame([{"id": n, **attrs} for n, attrs in graph.nodes(data=True)])
    metrics = [c for c in ("degree", "betweenness", "closeness", "eigenvector", "clustering", "k_core",
                           "is_hub", "is_bottleneck", "centrality_class") if c in centrality.columns]
    nodes = nodes.merge(centrality[["id"] + metrics], on="id", how="left")
    cluster_names = cluster_names or {}
    if "cluster" not in nodes:
        nodes["cluster"] = 0
    nodes["cluster_name"] = [
        f"{c}. {cluster_names.get(c, 'cluster ' + str(c))}" if c > 0 else "0. other" for c in nodes["cluster"]
    ]
    nodes["drug_status"] = [
        "AD drug target" if t else "druggable" if d else "none"
        for t, d in zip(nodes.get("ad_drug_target", False), nodes.get("druggable", False))
    ]
    key = (nodes.get("centrality_class", "peripheral") != "peripheral") | (nodes["drug_status"] == "AD drug target")
    nodes["label"] = [s if k else "" for s, k in zip(nodes["symbol"], key)]
    nodes["tooltip"] = [
        f"{r.symbol} ({r.role}) | cluster {r.cluster_name} | degree {r.degree} | betweenness {r.betweenness:.3f}"
        for r in nodes.itertuples()
    ]
    edges = pd.DataFrame(
        [{"source": u, "target": v, "interaction": "pp", **attrs} for u, v, attrs in graph.edges(data=True)]
    )
    return nodes, edges


# ---------------------------------------------------------------------------
# Visual style (CyREST JSON)
# ---------------------------------------------------------------------------


def _discrete(prop: str, column: str, column_type: str, mapping: dict) -> dict:
    return {
        "mappingType": "discrete",
        "mappingColumn": column,
        "mappingColumnType": column_type,
        "visualProperty": prop,
        "map": [{"key": str(k), "value": str(v)} for k, v in mapping.items()],
    }


def _continuous(prop: str, column: str, points: list[tuple[float, float]]) -> dict:
    return {
        "mappingType": "continuous",
        "mappingColumn": column,
        "mappingColumnType": "Double",
        "visualProperty": prop,
        "points": [{"value": x, "lesser": str(y), "equal": str(y), "greater": str(y)} for x, y in points],
    }


def _passthrough(prop: str, column: str) -> dict:
    return {"mappingType": "passthrough", "mappingColumn": column, "mappingColumnType": "String",
            "visualProperty": prop}


def build_style(clusters: list[int], max_betweenness: float) -> tuple[list[dict], list[dict]]:
    """Defaults and mappings of the network style.

    Node size grows with betweenness; the scale saturates at
    ``max_betweenness`` so a single dominant protein does not shrink the rest.
    """
    top = max(max_betweenness, 1e-6)
    defaults = {
        "NETWORK_BACKGROUND_PAINT": INK["surface"],
        "NODE_SHAPE": "ELLIPSE",
        "NODE_SIZE": 28,
        "NODE_FILL_COLOR": OTHER_COLOR,
        "NODE_BORDER_PAINT": INK["surface"],
        "NODE_BORDER_WIDTH": 1.5,
        "NODE_LABEL_COLOR": INK["primary"],
        "NODE_LABEL_FONT_SIZE": 13,
        "NODE_LABEL_FONT_FACE": "SansSerif.bold,Font.BOLD,13",
        "NODE_TRANSPARENCY": 235,
        "EDGE_STROKE_UNSELECTED_PAINT": INK["muted"],
        "EDGE_WIDTH": 0.8,
        "EDGE_TRANSPARENCY": 60,
    }
    mappings = [
        _discrete("NODE_FILL_COLOR", "cluster", "Integer", {c: cluster_color(c) for c in clusters}),
        _discrete("NODE_SHAPE", "role", "String", {"seed": "ELLIPSE", "diamond": "DIAMOND"}),
        _continuous("NODE_SIZE", "betweenness", [(0.0, 22), (top * 0.1, 40), (top * 0.5, 75), (top, 110)]),
        _discrete("NODE_BORDER_PAINT", "drug_status", "String", {k: v[0] for k, v in DRUG_BORDER.items()}),
        _discrete("NODE_BORDER_WIDTH", "drug_status", "String", {k: v[1] for k, v in DRUG_BORDER.items()}),
        _passthrough("NODE_LABEL", "label"),
        _passthrough("NODE_TOOLTIP", "tooltip"),
        _discrete("EDGE_TRANSPARENCY", "edge_type", "String",
                  {"seed-seed": 150, "diamond-seed": 70, "diamond-diamond": 45}),
        _discrete("EDGE_STROKE_UNSELECTED_PAINT", "edge_type", "String",
                  {"seed-seed": INK["secondary"], "diamond-seed": INK["muted"], "diamond-diamond": INK["muted"]}),
    ]
    return [{"visualProperty": k, "value": v} for k, v in defaults.items()], mappings


def write_style_json(path: str | Path, clusters: list[int], max_betweenness: float, name: str = STYLE_NAME) -> Path:
    """Save the style definition as JSON (version-controlled, reusable via CyREST)."""
    defaults, mappings = build_style(clusters, max_betweenness)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"title": name, "defaults": defaults, "mappings": mappings}, indent=2), "utf-8")
    return path


# ---------------------------------------------------------------------------
# Offline exports (no Cytoscape needed)
# ---------------------------------------------------------------------------


def to_cytoscape_js(nodes: pd.DataFrame, edges: pd.DataFrame, positions: dict | None = None) -> dict:
    """Cytoscape.js JSON (``.cyjs``), importable by Cytoscape desktop and cytoscape.js."""

    def clean(value):
        if isinstance(value, (bool, int, float, str)) or value is None:
            return value
        return value.item() if hasattr(value, "item") else str(value)

    elements = {"nodes": [], "edges": []}
    for row in nodes.to_dict(orient="records"):
        node = {"data": {"name": row["id"], **{k: clean(v) for k, v in row.items()}}}
        if positions and row["id"] in positions:
            x, y = positions[row["id"]]
            node["position"] = {"x": float(x), "y": float(y)}
        elements["nodes"].append(node)
    for i, row in enumerate(edges.to_dict(orient="records")):
        data = {k: clean(v) for k, v in row.items()}
        data["id"] = f"e{i}"
        elements["edges"].append({"data": data})
    return {"format_version": "1.0", "generated_by": "PPI_network_analysis_with_cytoscape",
            "data": {"name": STYLE_NAME}, "elements": elements}


def write_offline_files(nodes: pd.DataFrame, edges: pd.DataFrame, directory: str | Path,
                        positions: dict | None = None, prefix: str = "alzheimer_network") -> dict[str, Path]:
    """Write ``.cyjs`` and ``.graphml`` files that open directly in Cytoscape."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    cyjs = directory / f"{prefix}.cyjs"
    cyjs.write_text(json.dumps(to_cytoscape_js(nodes, edges, positions), indent=1), encoding="utf-8")
    graph = nx.Graph()
    for row in nodes.to_dict(orient="records"):
        graph.add_node(row["id"], **{k: v for k, v in row.items() if k != "id" and pd.notna(v)})
    for row in edges.to_dict(orient="records"):
        graph.add_edge(row["source"], row["target"], **{k: v for k, v in row.items() if k not in ("source", "target")})
    graphml = directory / f"{prefix}.graphml"
    nx.write_graphml(graph, graphml)
    return {"cyjs": cyjs, "graphml": graphml}


# ---------------------------------------------------------------------------
# Live Cytoscape session
# ---------------------------------------------------------------------------


def cytoscape_available(base_url: str = BASE_URL) -> bool:
    """``True`` if Cytoscape answers on its REST port."""
    try:
        import py4cytoscape as p4c

        p4c.cytoscape_ping(base_url=base_url)
        return True
    except Exception:  # noqa: BLE001 - any failure means "not reachable"
        return False


def load_network(nodes: pd.DataFrame, edges: pd.DataFrame, title: str, collection: str,
                 base_url: str = BASE_URL) -> int:
    """Create the network in Cytoscape from the node/edge tables; returns its SUID."""
    import py4cytoscape as p4c

    suid = p4c.create_network_from_data_frames(nodes, edges, title=title, collection=collection, base_url=base_url)
    logger.info("Loaded %s into Cytoscape (SUID %s)", title, suid)
    return suid


def apply_style(network: int, clusters: list[int], max_betweenness: float, name: str = STYLE_NAME,
                base_url: str = BASE_URL) -> str:
    """Create (or replace) the visual style and apply it to ``network``."""
    import py4cytoscape as p4c

    defaults, mappings = build_style(clusters, max_betweenness)
    if name in p4c.get_visual_style_names(base_url=base_url):
        p4c.delete_visual_style(name, base_url=base_url)
    p4c.create_visual_style(name, defaults=defaults, mappings=mappings, base_url=base_url)
    p4c.set_visual_style(name, network=network, base_url=base_url)
    # Node width/height are locked so NODE_SIZE controls both dimensions.
    p4c.lock_node_dimensions(True, style_name=name, base_url=base_url)
    return name
